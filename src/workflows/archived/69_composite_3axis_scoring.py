"""
69_composite_3axis_scoring.py
建议1 (review) — 三轴硬交集 → 复合评分排序, 取 Top 候选基因集

动机: 硬交集只得 2 基因, 统计功效与下游靶点池太单薄。改为三轴各自打分、标准化求和排序。

评分设计 (每轴 0-1 归一, 缺失=0):
  轴1 响应特异性: score1 = -log10(ΔΔE interaction_FDR), 截断到有限值后 min-max 归一
  轴2 空间(Hub): score2 = min-max( Perc_Identity 归一 * |wilcoxon_score| 归一 的几何均 ), Hub 内才有
  轴3 进化: 加权命中 = 0.4*sweep + 0.4*QTL + 0.3*posSel + 0.15*droughtNet + 0.1*CNV
            (稀有/强信号高权重: sweep/QTL/正选择 > 大集合 droughtNet/CNV), 再 min-max 归一
  composite = score1 + score2 + score3 (等权三轴), 排序取 Top 20

输出 Top 20 作为 Phase 5 候选靶点池 (替代 2 基因硬交集); 2 三轴硬交集基因标记为"最强支持"。

输入: ddE_interaction_GSE39894.tsv, hub_degs_tomato_orthologs.tsv, 文献集+RBH, QTL重叠表
输出: results/composite_3axis_scores.tsv (全部≥1轴), results/composite_top20_candidates.tsv
      drafts/fig2/fig2_composite_score_top20.png
非破坏性。
"""
import os
import re, sys
from datetime import datetime
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(os.environ.get("TOMATO_PROJECT_ROOT", "."))
def sv(s): return re.sub(r"\.\d+$", "", str(s))
def minmax(x):
    x = np.asarray(x, float); lo, hi = np.nanmin(x), np.nanmax(x)
    return (x - lo) / (hi - lo) if hi > lo else np.zeros_like(x)

PAPERS = {
    "sweep": ("data/papers/nph18634_NewPhytol2023_Sweeps/selective_sweeps_genes_Sopen.txt", "sopen", 0.40),
    "posSel": ("data/papers/mec17536_MolEcol2024_DroughtNetwork/positive_selection_core_genes_Schil.txt", "schil", 0.30),
    "droughtNet": ("data/papers/mec17536_MolEcol2024_DroughtNetwork/drought_network_genes_Schil.txt", "schil", 0.15),
    "cnv": ("data/papers/msaf191_MBE2025_CNV/cnv_differentiated_genes_Schil.txt", "schil", 0.10),
}
QTL_WEIGHT = 0.40
LOG = ROOT / f"logs/69_composite_3axis_scoring_{datetime.now():%Y%m%d}.log"


class Tee:
    def __init__(self, *f): self.files = f
    def write(self, o):
        for x in self.files: x.write(o); x.flush()
    def flush(self):
        for x in self.files: x.flush()


def main():
    LOG.parent.mkdir(parents=True, exist_ok=True)
    logf = open(LOG, "w"); sys.stdout = Tee(sys.__stdout__, logf)
    print(f"启动: {datetime.now():%Y-%m-%d %H:%M:%S}")
    print("=" * 72)

    # 轴1 ΔΔE
    dde = pd.read_csv(ROOT/"results/ddE_interaction_GSE39894.tsv", sep="\t")
    dde["g"] = dde["Solyc_ID"].map(sv)
    dde["nlfdr"] = -np.log10(dde["interaction_FDR"].clip(lower=1e-300))
    # 仅对 ΔΔE 显著(FDR<0.1)给轴1分, 其余0 (避免噪声基因刷分)
    dde["score1_raw"] = np.where(dde["interaction_FDR"] < 0.10, dde["nlfdr"], 0)

    # 轴2 Hub
    hub = pd.read_csv(ROOT/"results/hub_degs_tomato_orthologs.tsv", sep="\t")
    hub["g"] = hub["Tomato_Solyc_ID"].map(sv)
    hub["pid_n"] = minmax(hub["Perc_Identity"])
    hub["wil_n"] = minmax(hub["wilcoxon_score"].abs())
    hub["hub_raw"] = np.sqrt(hub["pid_n"] * hub["wil_n"])  # 几何均
    hub_score = hub.groupby("g")["hub_raw"].max()

    # 轴3 进化 (加权命中)
    schil = pd.read_csv(ROOT/"results/Schil_to_Solyc_RBH.tsv", sep="\t")
    sopen = pd.read_csv(ROOT/"results/Sopen_to_Solyc_RBH.tsv", sep="\t")
    smap = dict(zip(schil["Schil_ID"].astype(str), schil["Solyc_ID"].map(sv)))
    omap = dict(zip(sopen["Sopen_ID"].astype(str), sopen["Solyc_ID"].map(sv)))
    evo_w = {}
    for n, (p, k, w) in PAPERS.items():
        ids = [x.strip() for x in open(ROOT/p) if x.strip()]
        mp = omap if k == "sopen" else smap
        for g in {mp[i] for i in ids if i in mp}:
            evo_w[g] = evo_w.get(g, 0) + w
    # QTL 命中
    qtl = pd.read_csv(ROOT/"results/dts_zhang2024_precise_qtl_overlap.tsv", sep="\t")
    qcol = [c for c in qtl.columns if "Status" in c or "Overlap" in c]
    qtl_hit = set()
    if qcol:
        for _, r in qtl.iterrows():
            if "重叠" in str(r[qcol[0]]) and "无" not in str(r[qcol[0]]):
                qtl_hit.add(sv(r["Solyc_ID"]))
    for g in qtl_hit:
        evo_w[g] = evo_w.get(g, 0) + QTL_WEIGHT

    # 汇总所有 ≥1 轴基因
    allg = set(dde.loc[dde.score1_raw > 0, "g"]) | set(hub_score.index) | set(evo_w)
    df = pd.DataFrame({"Solyc_ID": sorted(allg)})
    s1 = dde.groupby("g")["score1_raw"].max()
    df["axis1_ddE"] = df["Solyc_ID"].map(s1).fillna(0)
    df["axis2_hub"] = df["Solyc_ID"].map(hub_score).fillna(0)
    df["axis3_evo"] = df["Solyc_ID"].map(evo_w).fillna(0)
    # 归一化
    df["s1"] = minmax(df["axis1_ddE"])
    df["s2"] = minmax(df["axis2_hub"])
    df["s3"] = minmax(df["axis3_evo"])
    df["composite"] = df["s1"] + df["s2"] + df["s3"]
    df["n_axes_hit"] = (df["axis1_ddE"] > 0).astype(int) + (df["axis2_hub"] > 0).astype(int) + (df["axis3_evo"] > 0).astype(int)
    df = df.sort_values("composite", ascending=False).reset_index(drop=True)

    # 标记硬交集2基因
    hard = {"Solyc01g096470", "Solyc01g011040"}
    df["hard_intersect_core"] = df["Solyc_ID"].isin(hard)
    df.to_csv(ROOT/"results/composite_3axis_scores.tsv", sep="\t", index=False)
    top20 = df.head(20).copy()
    top20.to_csv(ROOT/"results/composite_top20_candidates.tsv", sep="\t", index=False)
    print(f"全部≥1轴基因: {len(df)} | 输出 Top20 候选池")
    print(f"\nTop 20 候选 (composite 排序):")
    print(top20[["Solyc_ID","s1","s2","s3","composite","n_axes_hit","hard_intersect_core"]].round(3).to_string(index=False))
    print(f"\n硬交集2基因在复合排序中的位置:")
    for h in hard:
        r = df.index[df.Solyc_ID == h]
        if len(r): print(f"  {h}: rank {r[0]+1}/{len(df)}, composite={df.loc[r[0],'composite']:.3f}")

    # 图
    fig, ax = plt.subplots(figsize=(10, 8))
    t = top20.iloc[::-1]
    y = np.arange(len(t))
    ax.barh(y, t["s1"], color="#4C72B0", label="axis1 ΔΔE response")
    ax.barh(y, t["s2"], left=t["s1"], color="#55A868", label="axis2 spatial Hub")
    ax.barh(y, t["s3"], left=t["s1"]+t["s2"], color="#C44E52", label="axis3 evolution")
    labels = [f"{'★' if h else ''}{g}" for g,h in zip(t["Solyc_ID"], t["hard_intersect_core"])]
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("composite score (axis1+axis2+axis3, each min-max 0-1)")
    ax.set_title("Top 20 candidate genes by 3-axis composite score\n★ = original hard-intersection core", fontweight="bold")
    ax.legend(loc="lower right")
    plt.tight_layout(); plt.savefig(ROOT/"drafts/fig2/fig2_composite_score_top20.png", dpi=300); plt.close()
    print(f"\n图: drafts/fig2/fig2_composite_score_top20.png")
    print(f"完成: {datetime.now():%Y-%m-%d %H:%M:%S}")
    sys.stdout = sys.__stdout__; logf.close()


if __name__ == "__main__":
    main()
