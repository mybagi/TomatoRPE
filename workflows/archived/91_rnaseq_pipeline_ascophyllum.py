#!/usr/bin/env python3
"""
Ascophyllum nodosum 海藻提取物 priming 番茄抗旱 RNA-seq 流程
PRJNA592330: 干旱背景下 SF(海藻) vs water(对照), 各6样本(S3/S4时间点各3rep), 单端75bp.
fastp QC -> HISAT2 比对 SL4.0 -> samtools 排序 -> featureCounts 定量.
"""
import os
import subprocess, sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

PROJECT_DIR = Path(os.environ.get("TOMATO_PROJECT_ROOT", "."))
DATA_DIR = PROJECT_DIR / "data"
OUTPUT_DIR = DATA_DIR / "elicitors/Ascophyllum_seaweed_drought_tomato_PRJNA592330"
FASTQ_DIR = OUTPUT_DIR / "fastq"
QC_DIR = OUTPUT_DIR / "qc"
BAM_DIR = OUTPUT_DIR / "bam"
COUNTS_DIR = OUTPUT_DIR / "counts"

GENOME_INDEX = DATA_DIR / "genome/tomato_SL4.0/hisat2_index/SL4.0"
GFF_FILE = DATA_DIR / "genome/tomato_SL4.0/ITAG4.0_gene_models.gff"
FASTP = "/usr/bin/fastp"
SAMTOOLS = "/usr/bin/samtools"
HISAT2 = "hisat2"
FEATURECOUNTS = "featureCounts"
THREADS_PER_SAMPLE = 8
MAX_PARALLEL = 3

# 干旱背景核心对比: SF+drought(6) vs water+drought(6)
SAMPLES = {
    "SRR10550655": "SF_drought_S3_r1", "SRR10550654": "SF_drought_S3_r2", "SRR10550643": "SF_drought_S3_r3",
    "SRR10550632": "SF_drought_S4_r1", "SRR10550631": "SF_drought_S4_r2", "SRR10550630": "SF_drought_S4_r3",
    "SRR10550648": "water_drought_S3_r1", "SRR10550647": "water_drought_S3_r2", "SRR10550646": "water_drought_S3_r3",
    "SRR10550645": "water_drought_S4_r1", "SRR10550644": "water_drought_S4_r2", "SRR10550642": "water_drought_S4_r3",
}


def setup_dirs():
    for d in [QC_DIR, BAM_DIR, COUNTS_DIR]:
        d.mkdir(parents=True, exist_ok=True)


def run_sample(srr):
    label = SAMPLES[srr]
    bam_file = BAM_DIR / f"{srr}.sorted.bam"
    sam_file = BAM_DIR / f"{srr}.sam"
    if bam_file.exists() and bam_file.stat().st_size > 1000:
        print(f"  [{srr}] ({label}) BAM 已存在，跳过", flush=True)
        return srr, True, None
    try:
        fq = FASTQ_DIR / f"{srr}.fastq.gz"
        if not fq.exists():
            print(f"  [{srr}] 缺 FASTQ，跳过"); return srr, False, None
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
                aln = line.strip(); print(f"  [{srr}] ({label}) {aln}", flush=True); break
        cmd_sort = [SAMTOOLS, "sort", "-@", str(THREADS_PER_SAMPLE),
                    "-o", str(bam_file), str(sam_file)]
        r = subprocess.run(cmd_sort, capture_output=True, text=True, timeout=2400)
        if r.returncode != 0:
            print(f"  [{srr}] samtools sort 失败"); return srr, False, None
        sam_file.unlink(missing_ok=True)
        c1.unlink(missing_ok=True)
        return srr, True, aln
    except subprocess.TimeoutExpired:
        print(f"  [{srr}] 超时!"); return srr, False, None
    except Exception as e:
        print(f"  [{srr}] 异常: {e}"); return srr, False, None


def run_featurecounts():
    print("\nStep: featureCounts 基因定量 (单端)", flush=True)
    bam_files = sorted(BAM_DIR.glob("*.sorted.bam"))
    if not bam_files:
        print("  无 BAM"); return False
    output_file = COUNTS_DIR / "featureCounts_output.txt"
    cmd = [FEATURECOUNTS, "-a", str(GFF_FILE), "-o", str(output_file),
           "-F", "GFF", "-t", "exon", "-g", "Parent",
           "-T", "8"] + [str(b) for b in bam_files]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
    if r.returncode != 0:
        print(f"  featureCounts 失败: {r.stderr[:300]}"); return False
    print(f"  完成: {output_file}", flush=True)
    return True


if __name__ == "__main__":
    setup_dirs()
    print(f"Step: 比对 {len(SAMPLES)} 单端样本 (并行 {MAX_PARALLEL})", flush=True)
    ok = 0
    with ProcessPoolExecutor(max_workers=MAX_PARALLEL) as ex:
        futs = {ex.submit(run_sample, srr): srr for srr in SAMPLES}
        for f in as_completed(futs):
            srr, success, _ = f.result()
            if success: ok += 1
    print(f"\n比对完成: {ok}/{len(SAMPLES)}", flush=True)
    if ok >= 4:
        run_featurecounts()
