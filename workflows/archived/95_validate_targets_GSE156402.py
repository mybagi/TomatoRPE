"""
95_validate_targets_GSE156402.py
GSE156402 (叶脉/脉间干旱时序) 对 143 靶基因做独立方向一致性验证。

背景: 靶基因集 (results/correction_target_geneset.tsv) 的 expected_dir = 野生受旱方向,
      来自 PRJNA800740 (M82) + 野生种对比。本脚本用一套完全独立的番茄干旱数据
      (GSE156402, 栽培种 S.lyc, PNNL) 检验这些期望方向是否可复现。

注意: GSE156402 无诱抗剂处理臂 -> 不进逆转排名, 仅作靶基因方向的独立信号验证。

对比: T1 干旱(D) vs T1 watered(W), 同时间点配对, Vein 与 Intervein 各一套 3v3。
方法: 每组织 CPM 归一 -> log2FC = log2(mean_CPM_D+1) - log2(mean_CPM_W+1)
      -> 符号与 expected_dir 一致率 -> 二项检验 (H0: 一致率=0.5)。
非破坏性。
"""
import os
import sys
from datetime import datetime
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(os.environ.get("TOMATO_PROJECT_ROOT", "."))
XLSX = ROOT / "data/GSE156402_Processed_Data_gene_count_matrix.xlsx"
TARGET = ROOT / "results/correction_target_geneset.tsv"
OUT = ROOT / "results/validate_targets_GSE156402.tsv"
LOG = ROOT / f"logs/95_validate_targets_GSE156402_{datetime.now():%Y%m%d}.log"


class Tee:
    def __init__(self, *f): self.files = f
    def write(self, o):
        for x in self.files: x.write(o); x.flush()
    def flush(self):
        for x in self.files: x.flush()


def cpm(counts):
    """列为样本, 行为基因; 返回 CPM。"""
    lib = counts.sum(axis=0)
    return counts / lib * 1e6


def main():
    LOG.parent.mkdir(parents=True, exist_ok=True)
    logf = open(LOG, "w"); sys.stdout = Tee(sys.__stdout__, logf)
    print(f"启动: {datetime.now():%Y-%m-%d %H:%M:%S}")
    print("=" * 72)

    # ---- 读 GSE156402 (表头在索引2, 数据从索引3) ----
    raw = pd.read_excel(XLSX, header=None)
    hdr = [str(x) for x in raw.iloc[2].tolist()]
    data = raw.iloc[3:].copy()
    data.columns = hdr
    data = data.set_index("gene_id_short")

    sample_cols = [c for c in data.columns if "expression" in c]
    counts = data[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
    counts.columns = [c.replace(".fastq_reads_expression", "") for c in counts.columns]
    print(f"GSE156402: {counts.shape[0]} 基因 x {counts.shape[1]} 样本")

    # ---- 靶基因集 ----
    tgt = pd.read_csv(TARGET, sep="\t").set_index("Solyc_ID")
    print(f"靶基因集: {len(tgt)} 个 (UP={sum(tgt.expected_dir=='UP')}, "
          f"DOWN={sum(tgt.expected_dir=='DOWN')})")

    results = []
    for tissue, pref in [("Vein", "Vein"), ("Intervein", "Intervein")]:
        d_cols = [f"{pref}_T1_D_{i}" for i in (1, 2, 3)]
        w_cols = [f"{pref}_T1_W_{i}" for i in (1, 2, 3)]
        miss = [c for c in d_cols + w_cols if c not in counts.columns]
        if miss:
            print(f"  [{tissue}] 缺列: {miss}"); continue

        sub = counts[d_cols + w_cols]
        cpmv = cpm(sub)
        log2fc = np.log2(cpmv[d_cols].mean(axis=1) + 1) - \
                 np.log2(cpmv[w_cols].mean(axis=1) + 1)

        for g, row in tgt.iterrows():
            if g not in log2fc.index:
                continue
            fc = log2fc[g]
            obs = "UP" if fc > 0 else "DOWN"
            results.append({
                "Solyc_ID": g, "tissue": tissue,
                "expected_dir": row["expected_dir"],
                "GSE156402_log2FC": round(float(fc), 4),
                "GSE156402_dir": obs,
                "concordant": int(obs == row["expected_dir"]),
            })

    res = pd.DataFrame(results)
    res.to_csv(OUT, sep="\t", index=False)
    print(f"\n写出: {OUT}  ({len(res)} 行)")

    # ---- 二项检验 (每组织 + 合并) ----
    print("\n" + "=" * 72)
    print("方向一致性检验 (H0: 一致率=0.5, 二项检验单侧 greater)")
    print("-" * 72)
    for tissue in ["Vein", "Intervein"]:
        sub = res[res.tissue == tissue]
        n, k = len(sub), int(sub.concordant.sum())
        if n == 0: continue
        p = stats.binomtest(k, n, 0.5, alternative="greater").pvalue
        print(f"  {tissue:10s}: {k}/{n} 一致 ({k/n:.1%})  binom p={p:.3e}")
        # 分方向看
        for dd in ["UP", "DOWN"]:
            s2 = sub[sub.expected_dir == dd]
            if len(s2):
                print(f"      expected={dd:4s}: {int(s2.concordant.sum())}/{len(s2)} "
                      f"({s2.concordant.mean():.1%})")

    # 合并 (取两组织平均方向, 每基因一次)
    print("-" * 72)
    piv = res.pivot_table(index="Solyc_ID", columns="tissue",
                          values="GSE156402_log2FC")
    piv["mean_fc"] = piv.mean(axis=1)
    piv = piv.join(tgt["expected_dir"])
    piv["obs"] = np.where(piv["mean_fc"] > 0, "UP", "DOWN")
    piv["conc"] = (piv["obs"] == piv["expected_dir"]).astype(int)
    n, k = len(piv), int(piv["conc"].sum())
    p = stats.binomtest(k, n, 0.5, alternative="greater").pvalue
    print(f"  合并(双组织均值): {k}/{n} 一致 ({k/n:.1%})  binom p={p:.3e}")
    print("=" * 72)
    print(f"完成: {datetime.now():%Y-%m-%d %H:%M:%S}")


if __name__ == "__main__":
    main()
