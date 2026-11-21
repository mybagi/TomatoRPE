"""
90_rnaseq_pipeline_PRJNA925706_aba.py
PRJNA925706 (ABA 干旱番茄, "RNA-seq of tomato") RNA-seq 上游 —— 双端版 (330 单端).
  fastp QC (PE/SE) -> HISAT2 比对 SL4.0 (PE/SE) -> samtools 排序 -> featureCounts 定量.
输出 Solyc_ID raw count 矩阵 + log2(CPM+1) 矩阵.

设计: Solanum lycopersicum, 4 组 × 4 重复, 双端 150bp (SRR23283330 单端):
  WT = 常水对照            (38,37,30*,29)   *330单端
  PE = 干旱(PEG) 无ABA     (28,27,26,25)   ← 核心对比对照
  PA = 干旱 + ABA          (24,23,36,35)   ← 核心对比处理
  PC = 干旱 + 其他处理     (34,33,32,31)   [具体待原文确认]
核心逆转对比: PA vs PE (干旱 ± ABA), 情境匹配主线番茄干旱靶.
番茄本体 -> 无需 ortholog 映射.
分组来源: sample_alias (ENA experiment_alias/run_alias), 见 sample_groups.tsv.
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
OUTPUT_DIR = DATA_DIR / "elicitors/ABA_drought_tomato_PRJNA925706"
FASTQ_DIR = OUTPUT_DIR / "fastq"
QC_DIR = OUTPUT_DIR / "qc"
BAM_DIR = OUTPUT_DIR / "bam"
COUNTS_DIR = OUTPUT_DIR / "counts"
LOG_FILE = PROJECT_DIR / f"logs/90_rnaseq_PRJNA925706_{datetime.now():%Y%m%d_%H%M%S}.log"

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
    "SRR23283338": ("WT_1", "PE"), "SRR23283337": ("WT_2", "PE"),
    "SRR23283330": ("WT_3", "SE"), "SRR23283329": ("WT_4", "PE"),
    "SRR23283328": ("PE_1", "PE"), "SRR23283327": ("PE_2", "PE"),
    "SRR23283326": ("PE_3", "PE"), "SRR23283325": ("PE_4", "PE"),
    "SRR23283324": ("PA_1", "PE"), "SRR23283323": ("PA_2", "PE"),
    "SRR23283336": ("PA_3", "PE"), "SRR23283335": ("PA_4", "PE"),
    "SRR23283334": ("PC_1", "PE"), "SRR23283333": ("PC_2", "PE"),
    "SRR23283332": ("PC_3", "PE"), "SRR23283331": ("PC_4", "PE"),
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
    """单样本: fastp -> HISAT2 -> samtools sort. 支持 PE/SE."""
    label, layout = SAMPLE_META.get(srr, (srr, "PE"))
    bam_file = BAM_DIR / f"{srr}.sorted.bam"
    sam_file = BAM_DIR / f"{srr}.sam"

    if bam_file.exists() and bam_file.stat().st_size > 1000:
        print(f"  [{srr}] ({label}) BAM 已存在，跳过")
        return srr, True, None

    try:
        if layout == "PE":
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
        else:  # SE
            fq = FASTQ_DIR / f"{srr}.fastq.gz"
            if not fq.exists():
                print(f"  [{srr}] ({label}) 缺 SE FASTQ，跳过"); return srr, False, None
            c1 = QC_DIR / f"{srr}.clean.fq.gz"
            cmd_fastp = [FASTP, "-i", str(fq), "-o", str(c1),
                         "-h", str(QC_DIR/f"{srr}_fastp.html"),
                         "-j", str(QC_DIR/f"{srr}_fastp.json"),
                         "--thread", "4",
                         "--qualified_quality_phred", "20",
                         "--length_required", "36"]
            r = subprocess.run(cmd_fastp, capture_output=True, text=True, timeout=2400)
            if r.returncode != 0:
                print(f"  [{srr}] fastp 失败: {r.stderr[:200]}"); return srr, False, None
            cmd_hisat2 = [HISAT2, "-x", str(GENOME_INDEX),
                          "-U", str(c1),
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
        if layout == "PE":
            c2.unlink(missing_ok=True)
        return srr, True, aln

    except subprocess.TimeoutExpired:
        print(f"  [{srr}] 超时!"); return srr, False, None
    except Exception as e:
        print(f"  [{srr}] 异常: {e}"); return srr, False, None


def run_featurecounts():
    """对所有 BAM 运行 featureCounts (双端 -p --countReadPairs; 330单端由featureCounts自动按frag计)."""
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

    raw_out = OUTPUT_DIR / "PRJNA925706_Solyc_expr_counts.csv"
    counts.to_csv(raw_out)
    print(f"[OK] raw count: {raw_out}")

    cpm = counts / counts.sum(axis=0) * 1e6
    log2cpm = np.log2(cpm + 1)
    log_out = OUTPUT_DIR / "PRJNA925706_Solyc_expr_log2.csv"
    log2cpm.to_csv(log_out)
    print(f"[OK] log2(CPM+1): {log_out}")
    return counts


def main():
    setup_dirs()
    log_f = open(LOG_FILE, "w")
    sys.stdout = Tee(sys.__stdout__, log_f)
    print("=" * 60)
    print(f"PRJNA925706 ABA 干旱番茄 RNA-seq 上游 | {datetime.now():%Y-%m-%d %H:%M}")
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
