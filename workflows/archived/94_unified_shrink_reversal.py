#!/usr/bin/env python3
"""
94_unified_shrink_reversal.py
统一 apeGLM 收缩 + 聚焦逆转评分 + 敏感性分析.

目的: 回答"收缩后 Sindica 是否仍是最优正例"——在同口径(全部 apeGLM 收缩)下横向比较.
方法固定为 apeGLM (pydeseq2 lfc_shrink 默认, 自适应经验贝叶斯, 无可调阈值),
对所有可重算候选一视同仁, 避免选择偏倚.

可重算集合 (有count + 有重复 + 对比=干旱下诱抗剂vs干旱, 逻辑可比):
  Sindica     Sind40FC vs Ctrl40FC   ~condition        3v3
  MeJA        WDJA     vs WD          ~condition        3v3
  Ultrasound  UD       vs DR          ~condition        3v3
  ABA_tomato  PA       vs PE          ~condition        4v4
  Ascophyllum SF       vs water       ~timepoint+cond   6v6 (原设计带协变量)

无法纳入(标注沿用旧口径): Melatonin(无重复) Strigolactone(GR24仅常水,对比不可比)
  SA/ABA_ath/BR/Trehalose(无count矩阵)

敏感性分析: 每个可重算候选在 [不收缩 MLE, apeGLM 收缩] 两档下都算方向命中率+加权分,
  证明结论稳健(不靠挑设定).
输出: results/unified_shrink_reversal.tsv
"""
import os
import re, warnings
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
warnings.filterwarnings("ignore")

ROOT = Path(os.environ.get("TOMATO_PROJECT_ROOT", "."))
TARGET = ROOT / "results/correction_target_geneset.tsv"
OUT = ROOT / "results/unified_shrink_reversal.tsv"
sv = lambda s: re.sub(r"\.\d+$", "", str(s))

# (name, count文件, {样本列:组}, trt组名, ctrl组名, design)
JOBS = {
 "Sindica": (
    "data/elicitors/Fungi_Sindica_drought_tomato_GSE279454/GSE279454_Solyc_expr_counts.csv",
    None, "Sind40FC", "Ctrl40FC", "~condition"),  # 列名映射见下
 "MeJA": (
    "data/elicitors/MeJA_drought_tomato_CRA004475/CRA004475_Solyc_expr_counts.csv",
    {"WDJA_rep1":"WDJA","WDJA_rep2":"WDJA","WDJA_rep3":"WDJA",
     "WD_rep1":"WD","WD_rep2":"WD","WD_rep3":"WD"}, "WDJA", "WD", "~condition"),
 "Ultrasound": (
    "data/elicitors/Ultrasound_drought_tomato_GSE303938/GSE303938_Solyc_expr_counts.csv",
    {"UD_1":"UD","UD_2":"UD","UD_3":"UD","DR_1":"DR","DR_2":"DR","DR_3":"DR"},
    "UD", "DR", "~condition"),
 "ABA_tomato": (
    "data/elicitors/ABA_drought_tomato_PRJNA925706/PRJNA925706_Solyc_expr_counts.csv",
    {"PA_1":"PA","PA_2":"PA","PA_3":"PA","PA_4":"PA",
     "PE_1":"PE","PE_2":"PE","PE_3":"PE","PE_4":"PE"}, "PA", "PE", "~condition"),
}
# Sindica 列即 SRR, 单独定义
SIND_GROUPS = {"SRR31025484":"Sind40FC","SRR31025483":"Sind40FC","SRR31025482":"Sind40FC",
               "SRR31025502":"Ctrl40FC","SRR31025501":"Ctrl40FC","SRR31025500":"Ctrl40FC"}
JOBS["Sindica"] = (JOBS["Sindica"][0], SIND_GROUPS, "Sind40FC", "Ctrl40FC", "~condition")
JOBS["Ascophyllum"] = (None, None, "SF", "water", "~timepoint + condition")

# Ascophyllum: featureCounts 格式, 带时间点协变量
ASC_GROUPS = {  # SRR: (condition, timepoint)
  "SRR10550655":("SF","S3"),"SRR10550654":("SF","S3"),"SRR10550643":("SF","S3"),
  "SRR10550632":("SF","S4"),"SRR10550631":("SF","S4"),"SRR10550630":("SF","S4"),
  "SRR10550648":("water","S3"),"SRR10550647":("water","S3"),"SRR10550646":("water","S3"),
  "SRR10550645":("water","S4"),"SRR10550644":("water","S4"),"SRR10550642":("water","S4")}


def load_target():
    t = pd.read_csv(TARGET, sep="\t")
    t["Solyc_ID"] = t["Solyc_ID"].map(sv)
    return dict(zip(t["Solyc_ID"], np.where(t["expected_dir"] == "UP", 1, -1)))


def score(lfc_map, exp_dir):
    """给定 {solyc: log2FC}, 算方向命中 + 幅度加权分."""
    rows = [(g, lfc_map[g], exp_dir[g]) for g in lfc_map if g in exp_dir]
    if len(rows) < 3:
        return None
    df = pd.DataFrame(rows, columns=["solyc", "lfc", "exp"]).dropna()
    df = df[df["lfc"] != 0]
    n = len(df)
    agree = np.sign(df["lfc"]) == df["exp"]
    hit = int(agree.sum())
    p = stats.binomtest(hit, n, 0.5, alternative="greater").pvalue
    a = df["lfc"].abs()
    w = (a * np.where(agree, 1, -1)).sum() / a.sum() if a.sum() else 0.0
    return dict(n=n, hit=hit, frac=hit / n, p=p, weighted=w)


def run_deseq(counts_df, meta, design, contrast, do_shrink):
    """counts_df: 样本×基因(int); meta index=样本; 返回 {solyc: log2FC}.
    强制 ref_level=ctrl, 使 coeff 恒为 <factor>[T.<trt>], 方向 = trt vs ctrl."""
    factor, trt, ctrl = contrast
    sub = counts_df.loc[:, counts_df.sum(0) >= 10]
    dds = DeseqDataSet(counts=sub, metadata=meta, design=design,
                       ref_level=[factor, ctrl], refit_cooks=True, quiet=True)
    dds.deseq2()
    st = DeseqStats(dds, contrast=contrast, quiet=True)
    st.summary()
    if do_shrink:
        # pydeseq2 参照组按字母序自选, coeff 未必是 [T.<trt>].
        # 取该因子实际可用的 coeff; 若它是 [T.<ctrl>] (反方向), shrink 后取负.
        avail = [c for c in st.LFC.columns if c.startswith(f"{factor}[")]
        want = f"{factor}[T.{trt}]"
        flip = False
        if want in avail:
            coeff = want
        else:
            coeff = avail[0]  # [T.<ctrl>], 方向相反
            flip = True
        st.lfc_shrink(coeff=coeff)
        if flip:
            st.results_df["log2FoldChange"] = -st.results_df["log2FoldChange"]
    r = st.results_df.reset_index()
    r.columns.values[0] = "gene"
    r["solyc"] = r["gene"].str.extract(r"(Solyc\d+g\d+)", expand=False)
    r = r.dropna(subset=["solyc", "log2FoldChange"])
    r["abslfc"] = r["log2FoldChange"].abs()
    r = r.sort_values("abslfc", ascending=False).drop_duplicates("solyc")
    return dict(zip(r["solyc"], r["log2FoldChange"]))


def load_simple(cntfile, groups, trt, ctrl):
    """普通 csv (Solyc 行索引, 样本列). 返回 counts_df(样本×基因), meta, contrast.
    列名可为 '组_rep|SRR' 形式, 按 SRR 后缀或精确名匹配 groups 的键."""
    c = pd.read_csv(ROOT / cntfile)
    c = c.set_index(c.columns[0])
    colmap = {}  # 实际列名 -> 组
    for col in c.columns:
        key = col.split("|")[-1] if "|" in col else col  # 取 SRR 部分
        if key in groups:
            colmap[col] = groups[key]
        elif col in groups:
            colmap[col] = groups[col]
    cols = list(colmap)
    sub = c[cols].T.astype(int)
    sub.columns = [sv(g) for g in sub.columns]
    sub = sub.T.groupby(level=0).sum().T  # 合并重复 Solyc (列)
    meta = pd.DataFrame({"condition": [colmap[s] for s in cols]}, index=cols)
    return sub, meta, ["condition", trt, ctrl]


def load_ascophyllum():
    f = ROOT / "data/elicitors/Ascophyllum_seaweed_drought_tomato_PRJNA592330/counts/featureCounts_output.txt"
    df = pd.read_csv(f, sep="\t", comment="#")
    gc = df.columns[0]
    bam = [c for c in df.columns if c.endswith(".sorted.bam")]
    srr = {c: Path(c).name.replace(".sorted.bam", "") for c in bam}
    cnt = df[[gc] + bam].copy()
    cnt.columns = [gc] + [srr[c] for c in bam]
    cnt[gc] = cnt[gc].map(sv)
    cnt = cnt.groupby(gc).sum()
    srrs = [s for s in ASC_GROUPS if s in cnt.columns]
    sub = cnt[srrs].T.astype(int)
    meta = pd.DataFrame({"condition": [ASC_GROUPS[s][0] for s in srrs],
                         "timepoint": [ASC_GROUPS[s][1] for s in srrs]}, index=srrs)
    return sub, meta, ["condition", "SF", "water"]


def main():
    exp_dir = load_target()
    print(f"靶基因集: {len(exp_dir)} 个\n")
    rows = []
    for name, (cntfile, groups, trt, ctrl, design) in JOBS.items():
        try:
            if name == "Ascophyllum":
                sub, meta, contrast = load_ascophyllum()
                design = "~timepoint + condition"
            else:
                sub, meta, contrast = load_simple(cntfile, groups, trt, ctrl)
            for shrink in (False, True):
                lfc = run_deseq(sub, meta, design, contrast, shrink)
                s = score(lfc, exp_dir)
                tag = "apeGLM" if shrink else "MLE"
                if s:
                    rows.append(dict(elicitor=name, method=tag, **s))
                    print(f"[{name:11s}] {tag:7s} 覆盖靶{s['n']:3d} "
                          f"命中{s['hit']:3d}/{s['n']}={s['frac']:.0%} "
                          f"p={s['p']:.2g} 加权{s['weighted']:+.3f}")
        except Exception as e:
            print(f"[{name}] 失败: {e}")
    # Ascophyllum 单独设计, 也走上面(design已在main内改), 但它在JOBS外
    res = pd.DataFrame(rows)
    res.to_csv(OUT, sep="\t", index=False)
    print(f"\n[OK] 写入 {OUT}")
    # 排名总结 (apeGLM 口径)
    ape = res[res["method"] == "apeGLM"].sort_values("weighted", ascending=False)
    print("\n=== apeGLM 收缩后加权分排名 ===")
    for i, r in enumerate(ape.itertuples(), 1):
        print(f"  {i}. {r.elicitor:11s} 加权{r.weighted:+.3f} "
              f"命中{r.frac:.0%} p={r.p:.2g}")


if __name__ == "__main__":
    main()
