"""
52_ddE_interaction_GSE39894.py
Section 3N.1 — 基因型×处理交互项 (ΔΔE) 候选提出

重定义逻辑 (2026-06-26):
  抗旱性从未被逐样本测量, 旧 label=物种 与处理 100% 混淆。本脚本不使用任何
  分类器/抗旱标签, 仅用两个真实测量标签 (genotype + treatment), 对每个基因拟合交互项
      ΔΔE = (野生 drought − 野生 control) − (栽培 drought − 栽培 control)
  = 抗旱(野生)基因型在受旱时相对敏感(栽培)基因型特有的转录响应。
  形态学/物种固有差异在 ΔE 做差中被扣除, 构造上规避物种混淆。

⚠ 伪重复修正 (2026-06-26, 自查发现):
  GSE39894 每物种嵌套 3 个 population (lyc: 1421/2618/2789; pimp: 1592/1597/1630),
  每 population 3-4 重复。普通 OLS 把 46 样本当独立重复 → 伪重复 → 低估方差 → p 膨胀
  (初版 p<0.05 占 35.6%, FDR<0.05 达 902, 偏多)。
  本版改用线性混合模型 expr ~ genotype*treat + (1|population), population 作随机效应,
  把真实重复单位还原为 population 层级。混合模型不收敛的基因回退到 population 均值聚合
  OLS (6 pop 均值 × 2 treat), 并标记 method 列。

输入:
  - data/final_clean_matrices/GSE39894_Solyc_expr_log2_clean.csv (4004 基因 × 46 样本, genes×samples)
  - data/GSE39894/sample_metadata.csv (sample_id, species, treatment, ...)
  - results/ml_consensus_dts.tsv (旧 6 基因物种签名, 仅作对照)

输出:
  - results/ddE_interaction_GSE39894.tsv  (全基因: ΔΔE/ΔE_wild/ΔE_cult/interaction_p/FDR/主效应)
  - drafts/fig2/fig2_ddE_volcano.png      (火山图, DPI=300)
  - logs/52_ddE_interaction_GSE39894_<date>.log

非破坏性: 不覆盖任何 3.1-3.7 旧产物。
"""

import os
import os
import random
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

import warnings
from statsmodels.tools.sm_exceptions import ConvergenceWarning

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RANDOM_SEED = 42
random.seed(RANDOM_SEED)
np.random.seed(RANDOM_SEED)
os.environ["PYTHONHASHSEED"] = str(RANDOM_SEED)

ROOT = Path(os.environ.get("TOMATO_PROJECT_ROOT", "."))
EXPR = ROOT / "data/final_clean_matrices/GSE39894_Solyc_expr_log2_clean.csv"
META = ROOT / "data/GSE39894/sample_metadata.csv"
OLD6 = ROOT / "results/ml_consensus_dts.tsv"
OUT_TSV = ROOT / "results/ddE_interaction_GSE39894.tsv"
OUT_FIG = ROOT / "drafts/fig2/fig2_ddE_volcano.png"
LOG = ROOT / f"logs/52_ddE_interaction_GSE39894_{datetime.now():%Y%m%d}.log"


class Tee:
    def __init__(self, *files):
        self.files = files

    def write(self, obj):
        for f in self.files:
            f.write(obj)
            f.flush()

    def flush(self):
        for f in self.files:
            f.flush()


def main():
    OUT_FIG.parent.mkdir(parents=True, exist_ok=True)
    LOG.parent.mkdir(parents=True, exist_ok=True)
    logf = open(LOG, "w")
    sys.stdout = Tee(sys.__stdout__, logf)

    print(f"启动: {datetime.now():%Y-%m-%d %H:%M:%S} | seed={RANDOM_SEED}")
    print("=" * 72)

    # ---- 1. 载入表达矩阵 (genes×samples) → 转置为 samples×genes ----
    expr = pd.read_csv(EXPR, index_col=0)
    expr = expr.T  # samples × genes
    print(f"表达矩阵: {expr.shape[0]} 样本 × {expr.shape[1]} 基因 (转置后)")

    # ---- 2. 载入 metadata, 对齐 ----
    meta = pd.read_csv(META).set_index("sample_id")
    meta = meta.loc[expr.index]  # 对齐顺序
    # 二值编码: genotype (wild=pimp / cult=lyc), treatment (drought / control)
    meta["genotype"] = np.where(
        meta["species"].str.contains("pimpinellifolium"), "wild", "cult"
    )
    meta["treat"] = np.where(meta["treatment"] == "drought", "drought", "control")
    print("\n设计 (genotype × treat):")
    print(pd.crosstab(meta["genotype"], meta["treat"]))

    g = meta["genotype"].values
    t = meta["treat"].values
    g = meta["genotype"].values
    t = meta["treat"].values
    pop = meta["cultivar_or_population"].astype(str).values
    print("\npopulation 嵌套 (真实重复单位):")
    print(pd.crosstab(pop, t))
    # 四组样本掩码
    masks = {
        ("wild", "drought"): (g == "wild") & (t == "drought"),
        ("wild", "control"): (g == "wild") & (t == "control"),
        ("cult", "drought"): (g == "cult") & (t == "drought"),
        ("cult", "control"): (g == "cult") & (t == "control"),
    }

    # ---- 3. 逐基因 混合模型交互项 (population 随机效应) ----
    print(f"\n逐基因拟合 expr ~ genotype*treat + (1|population) ({expr.shape[1]} 基因)...")
    df_meta = pd.DataFrame(
        {"genotype": g, "treat": t, "population": pop}, index=expr.index
    )
    X = expr.values
    genes = expr.columns.values
    # ΔΔE 效应量 (group-means, 不受模型选择影响)
    mean_wd = X[masks[("wild", "drought")]].mean(0)
    mean_wc = X[masks[("wild", "control")]].mean(0)
    mean_cd = X[masks[("cult", "drought")]].mean(0)
    mean_cc = X[masks[("cult", "control")]].mean(0)
    dE_wild = mean_wd - mean_wc
    dE_cult = mean_cd - mean_cc
    ddE = dE_wild - dE_cult  # 交互项效应量 (group-means)

    inter_p = np.full(len(genes), np.nan)
    main_t_p = np.full(len(genes), np.nan)   # treatment 主效应 p
    main_g_p = np.full(len(genes), np.nan)   # genotype 主效应 p
    method = np.array(["mixed"] * len(genes), dtype=object)
    n_fallback = 0

    def get_p(res, want_inter=False, want_treat=False, want_geno=False):
        """从 statsmodels 结果按关键字取交互/主效应 p 值。"""
        idx = res.pvalues.index
        if want_inter:
            keys = [k for k in idx if ":" in k]
        elif want_treat:
            keys = [k for k in idx if "treat" in k.lower() and ":" not in k]
        elif want_geno:
            keys = [k for k in idx if "genotype" in k.lower() and ":" not in k]
        return res.pvalues[keys[0]] if keys else np.nan

    for i, gene in enumerate(genes):
        d = df_meta.copy()
        d["y"] = X[:, i]
        fitted = False
        # 优先: 线性混合模型, population 随机截距
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", ConvergenceWarning)
                warnings.simplefilter("ignore", RuntimeWarning)
                mm = smf.mixedlm(
                    "y ~ C(genotype) * C(treat)", d, groups=d["population"]
                ).fit(reml=False, method="lbfgs")
            if mm.converged:
                inter_p[i] = get_p(mm, want_inter=True)
                main_t_p[i] = get_p(mm, want_treat=True)
                main_g_p[i] = get_p(mm, want_geno=True)
                fitted = True
        except Exception:
            pass
        # 回退: population 均值聚合后 OLS (6 单元 × 2 treat)
        if not fitted or np.isnan(inter_p[i]):
            try:
                agg = d.groupby(["population", "genotype", "treat"])["y"].mean().reset_index()
                res = smf.ols("y ~ C(genotype) * C(treat)", agg).fit()
                inter_p[i] = get_p(res, want_inter=True)
                main_t_p[i] = get_p(res, want_treat=True)
                main_g_p[i] = get_p(res, want_geno=True)
                method[i] = "agg_ols"
                n_fallback += 1
            except Exception:
                method[i] = "failed"
        if (i + 1) % 500 == 0:
            print(f"  ...{i + 1}/{len(genes)}  (fallback so far: {n_fallback})")
    print(f"\n混合模型成功: {(method=='mixed').sum()} | 聚合OLS回退: {n_fallback} | 失败: {(method=='failed').sum()}")

    # ---- 4. BH-FDR ----
    valid = ~np.isnan(inter_p)
    fdr = np.full(len(genes), np.nan)
    fdr[valid] = multipletests(inter_p[valid], method="fdr_bh")[1]

    out = pd.DataFrame({
        "Solyc_ID": genes,
        "ddE": ddE,
        "dE_wild": dE_wild,
        "dE_cult": dE_cult,
        "interaction_p": inter_p,
        "interaction_FDR": fdr,
        "treatment_main_p": main_t_p,
        "genotype_main_p": main_g_p,
        "method": method,
    }).sort_values("interaction_p").reset_index(drop=True)
    out.to_csv(OUT_TSV, sep="\t", index=False)
    print(f"\n输出候选表: {OUT_TSV} ({len(out)} 基因)")

    # ---- 4b. p 值分布健康度 (伪重复修正后应回归正常) ----
    p = out["interaction_p"].dropna()
    hist = np.histogram(p, bins=10, range=(0, 1))[0]
    print("\n--- p 值分布健康度 (修正后) ---")
    print(f"  p 直方图(10 bins 0-1): {hist.tolist()}")
    print(f"  p<0.05 占比: {(p < 0.05).mean():.1%} (随机期望 5%; 初版膨胀为 35.6%)")
    print(f"  p>0.50 占比: {(p > 0.5).mean():.1%} (健康应 ~40-50%; 初版仅 25%)")

    # ---- 5. 显著性汇总 (诚实) ----
    n05 = int((out["interaction_FDR"] < 0.05).sum())
    n10 = int((out["interaction_FDR"] < 0.10).sum())
    nraw = int((out["interaction_p"] < 0.05).sum())
    print("\n--- 交互项显著性 (混合模型, population 随机效应) ---")
    print(f"  interaction_p<0.05 (未校正): {nraw}")
    print(f"  interaction_FDR<0.10        : {n10}")
    print(f"  interaction_FDR<0.05        : {n05}")
    print(f"  (初版伪重复 OLS: p<0.05=1424, FDR<0.05=902 — 应显著下降)")
    print(f"\nTop 10 交互项基因:")
    print(out.head(10).to_string(index=False))

    # ---- 6. 与旧 6 基因对照 ----
    if OLD6.exists():
        old = pd.read_csv(OLD6, sep="\t")
        oldcol = "Solyc_ID" if "Solyc_ID" in old.columns else old.columns[0]
        old_genes = old[oldcol].astype(str).tolist()
        print(f"\n--- 旧 6 基因物种签名在 ΔΔE 排名 (对照, 共 {len(out)} 基因) ---")
        idx = {gname: r for r, gname in enumerate(out["Solyc_ID"], 1)}
        for og in old_genes:
            if og in idx:
                row = out[out["Solyc_ID"] == og].iloc[0]
                print(f"  {og}: ΔΔE_rank={idx[og]}/{len(out)} | ddE={row['ddE']:+.3f} | "
                      f"inter_p={row['interaction_p']:.3g} | FDR={row['interaction_FDR']:.3g}")
            else:
                print(f"  {og}: 不在 GSE39894 4004 背景中")

    # ---- 7. 火山图 ----
    plt.figure(figsize=(7, 6))
    logp = -np.log10(out["interaction_p"].clip(lower=1e-300))
    sig = out["interaction_FDR"] < 0.10
    plt.scatter(out["ddE"][~sig], logp[~sig], s=8, c="lightgray", label="ns (FDR≥0.10)")
    plt.scatter(out["ddE"][sig], logp[sig], s=14, c="crimson", label="FDR<0.10")
    plt.axvline(0, color="k", lw=0.5, ls="--")
    plt.xlabel("ΔΔE = ΔE_wild − ΔE_cult  (interaction effect, log2)")
    plt.ylabel("−log10(interaction p)")
    plt.title("GSE39894 genotype×treatment interaction (ΔΔE)\n"
              f"n={nraw} p<0.05, {n10} FDR<0.10 | honest: n=10-12/group")
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(OUT_FIG, dpi=300)
    plt.close()
    print(f"\n火山图: {OUT_FIG}")

    print("\n" + "=" * 72)
    print(f"完成: {datetime.now():%Y-%m-%d %H:%M:%S}")
    sys.stdout = sys.__stdout__
    logf.close()


if __name__ == "__main__":
    main()
