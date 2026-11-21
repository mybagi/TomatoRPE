#!/usr/bin/env python3
"""
GSE279454 S. indica DESeq2 差异分析: Sind40FC (印度梨形孢+干旱) vs Ctrl40FC (空白+干旱).
最小验证集 3v3, 单一时间点 -> design = ~condition (无协变量).
输出 DEG.csv (gene,log2FC,pvalue,padj,direction,solyc) 供下游逆转评分.
"""
import os
import pandas as pd, numpy as np
from pathlib import Path
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

BASE = Path(os.environ.get("TOMATO_PROJECT_ROOT", ".")) / "data/elicitors/Fungi_Sindica_drought_tomato_GSE279454"
CNT = BASE / "counts/featureCounts_output.txt"
OUT = BASE / "DEG.csv"

# SRR -> 处理组
GROUPS = {
    "SRR31025484": "Sind40FC", "SRR31025483": "Sind40FC", "SRR31025482": "Sind40FC",
    "SRR31025502": "Ctrl40FC", "SRR31025501": "Ctrl40FC", "SRR31025500": "Ctrl40FC",
}

df = pd.read_csv(CNT, sep="\t", comment="#")
gene_col = df.columns[0]
bam_cols = [c for c in df.columns if c.endswith(".sorted.bam")]
srr_of = {c: Path(c).name.replace(".sorted.bam", "") for c in bam_cols}
counts = df[[gene_col] + bam_cols].copy()
counts.columns = [gene_col] + [srr_of[c] for c in bam_cols]
counts = counts.set_index(gene_col)

srrs = list(GROUPS.keys())
counts = counts[srrs].T.astype(int)

meta = pd.DataFrame({"condition": [GROUPS[s] for s in srrs]}, index=srrs)
print("样本分组:")
print(meta.groupby("condition").size())

keep = (counts.sum(axis=0) >= 10)
counts = counts.loc[:, keep]
print(f"\n保留基因: {counts.shape[1]}")

dds = DeseqDataSet(counts=counts, metadata=meta,
                   design="~condition", refit_cooks=True)
dds.deseq2()
stat = DeseqStats(dds, contrast=["condition", "Sind40FC", "Ctrl40FC"])
stat.summary()
res = stat.results_df.copy()

res = res.reset_index()
res = res.rename(columns={res.columns[0]: "gene", "log2FoldChange": "log2FC"})
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
