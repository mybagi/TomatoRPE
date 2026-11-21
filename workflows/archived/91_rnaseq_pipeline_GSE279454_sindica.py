"""
91_rnaseq_pipeline_GSE279454_sindica.py
GSE279454 (Serendipita indica × 干旱番茄) RNA-seq 上游 —— 双端最小验证集.
  fastp QC -> HISAT2 比对 SL4.0 -> samtools 排序 -> featureCounts 定量.
输出 Solyc_ID raw count 矩阵 + log2(CPM+1) 矩阵.

设计: Solanum lycopersicum, 最小验证集 6 样本 (2 组 × 3 重复), 双端:
  Ctrl40FC = 干旱 + 空白 (Mock)      (SRR31025502,501,500)  ← 核心对比对照
  Sind40FC = 干旱 + 印度梨形孢        (SRR31025484,483,482)  ← 核心对比处理
核心逆转对比: Sind40FC vs Ctrl40FC (干旱下 ± S. indica), 拿真 padj 验证逆转评分.
番茄本体 -> 无需 ortholog 映射.
非破坏性 (BAM 已存在则跳过).
"""

import os
import subprocess
import sys
import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed

PROJECT_DIR = Path(os.environ.get("TOMATO_PROJECT_ROOT", "."))
DATA_DIR = PROJECT_DIR / "data"
OUTPUT_DIR = DATA_DIR / "elicitors/Fungi_Sindica_drought_tomato_GSE279454"
FASTQ_DIR = OUTPUT_DIR / "fastq_raw"
QC_DIR = OUTPUT_DIR / "qc"
BAM_DIR = OUTPUT_DIR / "bam"
COUNTS_DIR = OUTPUT_DIR / "counts"
LOG_FILE = PROJECT_DIR / f"logs/91_rnaseq_GSE279454_{datetime.now():%Y%m%d_%H%M%S}.log"

GENOME_INDEX = DATA_DIR / "genome/tomato_SL4.0/hisat2_index/SL4.0"
GFF_FILE = DATA_DIR / "genome/tomato_SL4.0/ITAG4.0_gene_models.gff"

FASTP = "/usr/bin/fastp"
SAMTOOLS = "/usr/bin/samtools"
HISAT2 = "hisat2"
FEATURECOUNTS = "featureCounts"

THREADS_PER_SAMPLE = 8
MAX_PARALLEL = 3

# SRR -> (label, layout)
SAMPLE_META = {
    "SRR31025502": ("Ctrl40FC_1", "PE"), "SRR31025501": ("Ctrl40FC_2", "PE"),
    "SRR31025500": ("Ctrl40FC_3", "PE"), "SRR31025484": ("Sind40FC_1", "PE"),
    "SRR31025483": ("Sind40FC_2", "PE"), "SRR31025482": ("Sind40FC_3", "PE"),
}


class Tee:
    def __init__(self, *f): self.files = f
    def write(self, o):
        for x in self.files: x.write(o); x.flush()
    def flush(self):
        for x in self.files: x.flush()


def setup_dirs():
    for d in [QC_DIR, BAM_DIR, COUNTS_DIR]:
        d.mkdir(parents=True, exist_ok=True)


def run_sample(srr):
    """单样本: fastp -> HISAT2 -> samtools sort (PE)."""
    label, layout = SAMPLE_META.get(srr, (srr, "PE"))
    bam_file = BAM_DIR / f"{srr}.sorted.bam"
    sam_file = BAM_DIR / f"{srr}.sam"

    if bam_file.exists() and bam_file.stat().st_size > 1000:
        print(f"  [{srr}] ({label}) BAM 已存在，跳过")
        return srr, True, None

    try:
        fq1 = FASTQ_DIR / f"{srr}_1.fastq.gz"
        fq2 = FASTQ_DIR / f"{srr}_2.fastq.gz"
        if not (fq1.exists() and fq2.exists()):
            print(f"  [{srr}] ({label}) 缺 PE FASTQ，跳过"); return srr, False, None
        c1 = QC_DIR / f"{srr}_1.clean.fq.gz"
        c2 = QC_DIR / f"{srr}_2.clean.fq.gz"
        cmd_fastp = [FASTP, "-i", str(fq1), "-I", str(fq2),
                     "-o", str(c1), "-O", str(c2),
                     "-h", str(QC_DIR/f"{srr}_fastp.html"),
                     "-j", str(QC_DIR/f"{srr}_fastp.json"),
                     "--thread", "4",
                     "--qualified_quality_phred", "20",
                     "--length_required", "36"]
        r = subprocess.run(cmd_fastp, capture_output=True, text=True, timeout=2400)
        if r.returncode != 0:
            print(f"  [{srr}] fastp 失败: {r.stderr[:200]}"); return srr, False, None
        cmd_hisat2 = [HISAT2, "-x", str(GENOME_INDEX),
                      "-1", str(c1), "-2", str(c2),
                      "-S", str(sam_file), "-p", str(THREADS_PER_SAMPLE),
                      "--no-unal", "--dta"]

        r = subprocess.run(cmd_hisat2, capture_output=True, text=True, timeout=5400)
        if r.returncode != 0:
            print(f"  [{srr}] HISAT2 失败: {r.stderr[:200]}"); return srr, False, None
        aln = ""
        for line in r.stderr.split("\n"):
            if "overall alignment rate" in line:
                aln = line.strip(); print(f"  [{srr}] ({label}) {aln}"); break

        cmd_sort = [SAMTOOLS, "sort", "-@", str(THREADS_PER_SAMPLE),
                    "-o", str(bam_file), str(sam_file)]
        r = subprocess.run(cmd_sort, capture_output=True, text=True, timeout=2400)
        if r.returncode != 0:
            print(f"  [{srr}] samtools sort 失败"); return srr, False, None

        sam_file.unlink(missing_ok=True)
        c1.unlink(missing_ok=True)
        c2.unlink(missing_ok=True)
        return srr, True, aln

    except subprocess.TimeoutExpired:
        print(f"  [{srr}] 超时!"); return srr, False, None
    except Exception as e:
        print(f"  [{srr}] 异常: {e}"); return srr, False, None


def run_featurecounts():
    """对所有 BAM 运行 featureCounts (双端 -p --countReadPairs)."""
    print("\n" + "=" * 60)
    print("Step 4: featureCounts 基因定量 (双端)")
    print("=" * 60)
    bam_files = sorted(BAM_DIR.glob("*.sorted.bam"))
    if not bam_files:
        raise RuntimeError("无 BAM 文件!")
    print(f"输入: {len(bam_files)} 个 BAM")
    output_file = COUNTS_DIR / "featureCounts_output.txt"
    cmd = [
        FEATURECOUNTS,
        "-a", str(GFF_FILE),
        "-o", str(output_file),
        "-F", "GFF", "-t", "exon", "-g", "Parent",
        "-T", "16", "--minOverlap", "10",
        "-p", "--countReadPairs",
    ] + [str(b) for b in bam_files]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=5400)
    if r.returncode != 0:
        print(f"featureCounts 失败!\n{r.stderr[:500]}")
        raise RuntimeError("featureCounts failed")
    summary_file = Path(str(output_file) + ".summary")
    if summary_file.exists():
        print("\n定量统计摘要:")
        print(pd.read_csv(summary_file, sep="\t").to_string(index=False))
    print(f"\n[OK] featureCounts 完成: {output_file}")
    return output_file


def extract_solyc(gene_id):
    if ":" in gene_id:
        gene_id = gene_id.split(":")[-1]
    parts = gene_id.split(".")
    if len(parts) >= 2 and parts[0].startswith("Solyc"):
        return f"{parts[0]}.{parts[1]}"
    return gene_id


def build_matrices(fc_output):
    """featureCounts 输出 -> Solyc raw count + log2(CPM+1) 矩阵."""
    print("\n" + "=" * 60)
    print("Step 5: 构建表达矩阵")
    print("=" * 60)
    df = pd.read_csv(fc_output, sep="\t", comment="#")
    meta_cols = ["Geneid", "Chr", "Start", "End", "Strand", "Length"]
    bam_cols = [c for c in df.columns if c not in meta_cols]
    rename = {}
    for c in bam_cols:
        srr = Path(c).name.replace(".sorted.bam", "")
        label = SAMPLE_META.get(srr, (srr, ""))[0]
        rename[c] = f"{label}|{srr}"
    df = df.rename(columns=rename)
    df["Solyc"] = df["Geneid"].apply(extract_solyc)
    data_cols = list(rename.values())
    counts = df.groupby("Solyc")[data_cols].sum()
    counts = counts[counts.sum(axis=1) > 0]
    print(f"基因数 (count>0): {len(counts)}")

    raw_out = OUTPUT_DIR / "GSE279454_Solyc_expr_counts.csv"
    counts.to_csv(raw_out)
    print(f"[OK] raw count: {raw_out}")

    cpm = counts / counts.sum(axis=0) * 1e6
    log2cpm = np.log2(cpm + 1)
    log_out = OUTPUT_DIR / "GSE279454_Solyc_expr_log2.csv"
    log2cpm.to_csv(log_out)
    print(f"[OK] log2(CPM+1): {log_out}")
    return counts


def main():
    setup_dirs()
    log_f = open(LOG_FILE, "w")
    sys.stdout = Tee(sys.__stdout__, log_f)
    print("=" * 60)
    print(f"GSE279454 S. indica 干旱番茄 RNA-seq 上游 | {datetime.now():%Y-%m-%d %H:%M}")
    print("=" * 60)

    srrs = list(SAMPLE_META.keys())
    print(f"\nStep 1-3: fastp -> HISAT2 -> sort ({len(srrs)} 样本, 并行 {MAX_PARALLEL})")
    ok = 0
    with ThreadPoolExecutor(max_workers=MAX_PARALLEL) as ex:
        futs = {ex.submit(run_sample, s): s for s in srrs}
        for fu in as_completed(futs):
            srr, success, aln = fu.result()
            if success: ok += 1
    print(f"\n比对完成: {ok}/{len(srrs)}")
    if ok == 0:
        raise RuntimeError("无成功样本")

    fc = run_featurecounts()
    build_matrices(fc)
    print("\n[全部完成]")


if __name__ == "__main__":
    main()
