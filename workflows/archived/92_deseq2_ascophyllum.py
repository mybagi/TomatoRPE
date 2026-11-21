#!/usr/bin/env python3
"""
Ascophyllum 海藻 DESeq2 差异分析: SF(海藻)+drought vs water+drought.
输出 DEG.csv (gene,log2FC,pvalue,padj,direction,solyc) 供脚本76逆转评分.
时间点(S3/S4)作为批次协变量纳入模型.
"""
import os
import pandas as pd, numpy as np
from pathlib import Path
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

BASE = Path(os.environ.get("TOMATO_PROJECT_ROOT", ".")) / "data/elicitors/Ascophyllum_seaweed_drought_tomato_PRJNA592330"
CNT = BASE / "counts/featureCounts_output.txt"
OUT = BASE / "DEG.csv"

# SRR -> (处理组, 时间点)
GROUPS = {
    "SRR10550655": ("SF", "S3"), "SRR10550654": ("SF", "S3"), "SRR10550643": ("SF", "S3"),
    "SRR10550632": ("SF", "S4"), "SRR10550631": ("SF", "S4"), "SRR10550630": ("SF", "S4"),
    "SRR10550648": ("water", "S3"), "SRR10550647": ("water", "S3"), "SRR10550646": ("water", "S3"),
    "SRR10550645": ("water", "S4"), "SRR10550644": ("water", "S4"), "SRR10550642": ("water", "S4"),
}

# 读featureCounts (跳过首行注释)
df = pd.read_csv(CNT, sep="\t", comment="#")
gene_col = df.columns[0]
# BAM列 -> SRR
bam_cols = [c for c in df.columns if c.endswith(".sorted.bam")]
srr_of = {c: Path(c).name.replace(".sorted.bam", "") for c in bam_cols}
counts = df[[gene_col] + bam_cols].copy()
counts.columns = [gene_col] + [srr_of[c] for c in bam_cols]
counts = counts.set_index(gene_col)

# 转置: 样本 x 基因
srrs = list(GROUPS.keys())
counts = counts[srrs].T
counts = counts.astype(int)

# 元数据
meta = pd.DataFrame({
    "condition": [GROUPS[s][0] for s in srrs],
    "timepoint": [GROUPS[s][1] for s in srrs],
}, index=srrs)

print("样本分组:")
print(meta.groupby(["condition", "timepoint"]).size())

# 过滤低表达基因
keep = (counts.sum(axis=0) >= 10)
counts = counts.loc[:, keep]
print(f"\n保留基因: {counts.shape[1]} (原始过滤后)")

# DESeq2: ~timepoint + condition (时间点作协变量), 比较 SF vs water
dds = DeseqDataSet(counts=counts, metadata=meta,
                   design="~timepoint + condition", refit_cooks=True)
dds.deseq2()
stat = DeseqStats(dds, contrast=["condition", "SF", "water"])
stat.summary()
res = stat.results_df.copy()

res = res.reset_index()
res = res.rename(columns={res.columns[0]: "gene", "log2FoldChange": "log2FC"})
# 清理基因ID: "mRNA:Solyc00g500047.1.1" -> "Solyc00g500047"
res["solyc"] = res["gene"].str.extract(r"(Solyc\d+g\d+)", expand=False)
res = res.dropna(subset=["solyc"])
res["direction"] = "ns"
res.loc[(res["padj"] < 0.05) & (res["log2FC"] > 0), "direction"] = "up"
res.loc[(res["padj"] < 0.05) & (res["log2FC"] < 0), "direction"] = "down"
res["gene"] = res["solyc"]
res = res[["gene", "log2FC", "pvalue", "padj", "direction", "solyc"]]
res.to_csv(OUT, index=False)

n_sig = (res["padj"] < 0.05).sum()
print(f"\n显著DEG (padj<0.05): {n_sig}")
print(f"  up: {(res['direction']=='up').sum()}, down: {(res['direction']=='down').sum()}")
print(f"输出: {OUT}")
