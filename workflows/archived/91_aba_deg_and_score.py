"""
91_aba_deg_and_score.py
E-② ABA 逆转分析 (PRJNA925706, 番茄本体·干旱背景, 双端·4重复)

pydeseq2 差异表达 (带重复, 有 FDR):
  核心: PA vs PE (干旱+ABA vs 干旱)  -> ABA 在干旱背景的效应
  健全性: PE vs WT (干旱 vs 常水)    -> 阳性对照, 应有强信号(抓到生物学)
  辅助: PC vs PE (其他处理 vs 干旱)  -> 参照

输出:
  DEG.csv (gene/log2FC/padj, PA vs PE, 喂 scripts/76 逆转评分)
  DEG_PEvsWT.csv (健全性)
非破坏性。
"""
import os
import sys
import numpy as np, pandas as pd
from pathlib import Path
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

ROOT = Path(os.environ.get("TOMATO_PROJECT_ROOT", "."))
DDIR = ROOT / "data/elicitors/ABA_drought_tomato_PRJNA925706"
COUNTS = DDIR / "PRJNA925706_Solyc_expr_counts.csv"

counts = pd.read_csv(COUNTS, index_col=0)
print(f"count矩阵: {counts.shape[0]} 基因 × {counts.shape[1]} 样本")

groups = {c: c.split("_")[0] for c in counts.columns}


def run_deseq(trt, ctrl, tag):
    cols = [c for c in counts.columns if groups[c] in (trt, ctrl)]
    sub = counts[cols].T
    meta = pd.DataFrame({"condition": [groups[c] for c in cols]}, index=cols)
    # 过滤低表达
    keep = sub.sum(axis=0) >= 10
    sub = sub.loc[:, keep]
    dds = DeseqDataSet(counts=sub, metadata=meta, design_factors="condition",
                       ref_level=["condition", ctrl], quiet=True)
    dds.deseq2()
    st = DeseqStats(dds, contrast=["condition", trt, ctrl], quiet=True)
    st.summary()
    res = st.results_df.copy()
    res.index.name = "gene"
    n_sig = ((res["padj"] < 0.05) & (res["log2FoldChange"].abs() >= 1)).sum()
    print(f"\n[{tag}] {trt} vs {ctrl}: {len(res)} 基因检验, "
          f"|log2FC|>=1 & padj<0.05 = {n_sig} 个 DEG")
    return res


print("\n=== 健全性: PE vs WT (干旱本身, 阳性对照) ===")
res_pewt = run_deseq("PE", "WT", "健全性")
deg_pewt = pd.DataFrame({"gene": res_pewt.index,
                         "log2FC": res_pewt["log2FoldChange"].values,
                         "padj": res_pewt["padj"].values})
deg_pewt.to_csv(DDIR / "DEG_PEvsWT.csv", index=False)

print("\n=== 核心: PA vs PE (干旱+ABA vs 干旱) ===")
res_pape = run_deseq("PA", "PE", "核心")
deg = pd.DataFrame({"gene": res_pape.index,
                    "log2FC": res_pape["log2FoldChange"].values,
                    "padj": res_pape["padj"].values})
deg.to_csv(DDIR / "DEG.csv", index=False)
print(f"[OK] PA vs PE → {DDIR/'DEG.csv'}")

# 健全性判定
n_drought = ((res_pewt["padj"] < 0.05) & (res_pewt["log2FoldChange"].abs() >= 1)).sum()
n_aba = ((res_pape["padj"] < 0.05) & (res_pape["log2FoldChange"].abs() >= 1)).sum()
print(f"\n{'='*50}")
print(f"健全性: 干旱本身(PE vs WT) {n_drought} 个 DEG",
      "→ 抓到生物学 ✓" if n_drought > 500 else "→ ⚠ 干旱信号偏弱")
print(f"ABA干旱效应(PA vs PE) {n_aba} 个 DEG",
      "→ ABA有转录扰动" if n_aba > 100 else "→ ⚠ ABA扰动弱")
