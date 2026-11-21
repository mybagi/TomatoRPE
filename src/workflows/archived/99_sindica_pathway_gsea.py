#!/usr/bin/env python3
"""
99_sindica_pathway_gsea.py  [优化方向2: 通路层富集定靶点]

问题: S.indica 单基因 log2FC≈0 都不显著, 但作用是"全局微调".
思路: 单基因抓不到的信号, 通路层(GSEA)可能抓到——GSEA 用整条基因排序,
      对协同小变化敏感. 若真菌在某 GO_BP 通路(如 ABA 信号)整体显著富集,
      则从"被真菌显著推动的通路"里选靶点, 真菌关联证据比单基因硬.

方法: 以 S.indica 全基因 log2FC_MLE 排序 → gseapy prerank → GO_BP 基因集.
      GO_BP 基因集从 GAF (gene_association.itag4.0.nr, namespace P) 构建.
输出:
  results/sindica_pathway_gsea.tsv        (全部通路 NES/FDR)
  results/sindica_pathway_core_targets.tsv (核心基因落在显著富集通路的映射)
"""
import os
import re
from pathlib import Path
from collections import defaultdict
import pandas as pd
import gseapy as gp

ROOT = Path(os.environ.get("TOMATO_PROJECT_ROOT", "."))
sv = lambda s: re.sub(r"\.\d+$", "", str(s))

# --- 1. 真菌全基因排序 ---
deg = pd.read_csv(ROOT / "results/sindica_strict_DEG.tsv", sep="\t")
deg["g"] = deg["Solyc_ID"].map(sv)
rnk = (deg.dropna(subset=["log2FC_MLE"])
          .groupby("g")["log2FC_MLE"].mean()
          .sort_values(ascending=False))
print(f"排序基因数: {len(rnk)}")

# --- 2. GO_BP 基因集 (namespace P) 从 GAF 构建 ---
# GAF 列: 0=db 1=obj_id 4=GO_ID 8=aspect(P/F/C)
gaf = pd.read_csv(ROOT / "gene_association.itag4.0.nr", sep="\t", header=None,
                  usecols=[1, 4, 8], names=["gene", "go", "aspect"],
                  comment="!", dtype=str)
gaf = gaf[gaf["aspect"] == "P"]              # 只要 Biological Process
gaf["g"] = gaf["gene"].map(sv)

# GO 名称 (从 obo)
go_name = {}
with open(ROOT / "data/go-basic.obo") as fh:
    cid = cname = None
    for line in fh:
        line = line.strip()
        if line == "[Term]":
            cid = cname = None
        elif line.startswith("id: GO:"):
            cid = line[4:]
        elif line.startswith("name:"):
            cname = line[6:]
            if cid:
                go_name[cid] = cname

gene_sets = defaultdict(set)
for g, go in zip(gaf["g"], gaf["go"]):
    nm = go_name.get(go, go)
    gene_sets[f"{nm}|{go}"].add(g)
# 过滤过大过小基因集
gene_sets = {k: list(v) for k, v in gene_sets.items() if 10 <= len(v) <= 500}
print(f"GO_BP 基因集数 (10-500): {len(gene_sets)}")

# --- 3. prerank GSEA ---
pre = gp.prerank(rnk=rnk.reset_index().rename(columns={"g": "gene"}),
                 gene_sets=gene_sets, min_size=10, max_size=500,
                 permutation_num=1000, seed=42, threads=8, no_plot=True,
                 outdir=None)
res = pre.res2d.copy()
res.columns = [c.lower().replace(" ", "_") for c in res.columns]
fdr_col = "fdr_q-val" if "fdr_q-val" in res.columns else [c for c in res.columns if "fdr" in c][0]
nes_col = "nes"
res = res.sort_values(fdr_col)
res.to_csv(ROOT / "results/sindica_pathway_gsea.tsv", sep="\t", index=False)

sig = res[res[fdr_col].astype(float) < 0.25]
print(f"\n[OK] results/sindica_pathway_gsea.tsv")
print(f"FDR<0.25 富集通路: {len(sig)}")
print("\n=== 真菌显著推动/抑制的通路 Top15 (按FDR) ===")
for r in res.head(15).itertuples():
    term = getattr(r, "term")[:50]
    print(f"  NES={float(getattr(r,nes_col)):+.2f}  FDR={float(getattr(r,fdr_col)):.3f}  {term}")

# --- 4. 核心可成药靶点落在显著通路的映射 ---
sel = pd.read_csv(ROOT / "results/druggable_target_selection.tsv", sep="\t")
sel["g"] = sel["Solyc_ID"].map(sv)
drug = set(sel[sel["druggable_class"].isin(["DRUGGABLE", "MIXED"])]["g"])
rows = []
for r in sig.itertuples():
    term = getattr(r, "Term") if hasattr(r, "Term") else getattr(r, "term")
    key = term
    members = gene_sets.get(key, [])
    for g in members:
        if g in drug:
            rows.append(dict(pathway=term, nes=float(getattr(r, nes_col)),
                             fdr=float(getattr(r, fdr_col)), Solyc_ID=g))
tg = pd.DataFrame(rows)
if len(tg):
    tg = tg.merge(sel[["g", "abs_ddE", "sind_log2FC", "direction"]].rename(columns={"g": "Solyc_ID"}),
                  on="Solyc_ID", how="left").sort_values("fdr")
    tg.to_csv(ROOT / "results/sindica_pathway_core_targets.tsv", sep="\t", index=False)
    print(f"\n[OK] results/sindica_pathway_core_targets.tsv")
    print(f"可成药靶点落在显著富集通路: {tg['Solyc_ID'].nunique()} 个")
    print("\n=== 通路支撑的可成药靶点 ===")
    for r in tg.drop_duplicates('Solyc_ID').head(10).itertuples():
        print(f"  {r.Solyc_ID}  ddE={r.abs_ddE:.2f} sindLFC={r.sind_log2FC:+.2f}  "
              f"[{r.pathway[:38]} NES={r.nes:+.2f} FDR={r.fdr:.3f}]")
else:
    print("\n无可成药靶点落在 FDR<0.25 通路 (通路层也未显著)")
