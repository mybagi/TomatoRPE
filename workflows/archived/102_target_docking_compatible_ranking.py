#!/usr/bin/env python3
"""
102_target_docking_compatible_ranking.py  [5.4.1 定稿: 加对接兼容硬条件]

修补 5.4.1 的下游埋坑: 原遴选漏了"靶点需激活 vs 抑制"判据, 而分子对接
(AutoDock Vina) 天然筛"结合活性口袋的分子"≈抑制剂. 若选 UP 方向靶点
(诱抗剂需激活), 对接筛出的抑制剂会关闭该抗旱酶, 方向反, 下游白做.

硬条件(从源头防坑):
  校正方向 both_DOWN_wild-drought (野生受旱下调, 栽培需下调)
    → 诱抗剂要抑制该靶点 → 对接找抑制剂 → 逻辑自洽, docking_compatible=True
  both_UP_wild-drought (需激活) → 对接筛抑制剂方向反 → docking_compatible=False
    (标注保留, 但不进主对接排名; 如需可走别构激活/激动剂专库另议)

在 docking_compatible=True 的可成药靶点内, 复用 100 的多证据评分+权重敏感性+Borda.

输出: results/target_docking_ranking.tsv  (定稿主对接靶点名单)
"""
import os
import re
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats

ROOT = Path(os.environ.get("TOMATO_PROJECT_ROOT", "."))
sv = lambda s: re.sub(r"\.\d+$", "", str(s))
rank01 = lambda s: s.rank(pct=True)

sel = pd.read_csv(ROOT / "results/druggable_target_selection.tsv", sep="\t")
sel["g"] = sel["Solyc_ID"].map(sv)
comp = pd.read_csv(ROOT / "results/composite_3axis_scores.tsv", sep="\t")
comp["g"] = comp["Solyc_ID"].map(sv)
m = sel.merge(comp[["g", "axis1_ddE", "axis2_hub", "axis3_evo"]], on="g", how="left")
m = m[m["druggable_class"].isin(["DRUGGABLE", "MIXED"])].copy()

# --- 干预类型 + 对接兼容性 ---
m["intervention_type"] = np.where(m["direction"].str.contains("DOWN"),
                                  "inhibit", "activate")
m["docking_compatible"] = m["intervention_type"] == "inhibit"

# 硬条件: 只在 docking_compatible=True 且真菌方向一致的池内排名
pool = m[m["docking_compatible"] & m["sind_agree"].fillna(False)].copy().reset_index(drop=True)
print(f"可成药池 {len(m)} → 对接兼容(DOWN需抑制)且真菌方向一致: {len(pool)}\n")

# --- 多证据评分 (池内秩归一) ---
pool["E_ddE"] = rank01(pool["abs_ddE"].fillna(0))
pool["E_hub"] = rank01(pool["axis2_hub"].fillna(0))
pool["E_evo"] = rank01(pool["axis3_evo"].fillna(0))
sabs = pool["sind_log2FC"].abs().fillna(0)
pool["E_sind"] = rank01(sabs)   # 池内已全部方向一致, 直接用幅度秩
pool["is_hub"] = pool["axis2_hub"].fillna(0) > 0

SCHEMES = {"equal": (1,1,1,1), "hub_low": (1,0.5,1,1),
           "sind_heavy": (1,0.5,1,2), "evo_heavy": (1,0.5,2,1)}
for name, (wa,wh,we,ws) in SCHEMES.items():
    sc = (wa*pool["E_ddE"]+wh*pool["E_hub"]+we*pool["E_evo"]+ws*pool["E_sind"])/(wa+wh+we+ws)
    pool[f"rank_{name}"] = sc.rank(ascending=False).astype(int)

N = len(pool)
pool["borda"] = sum((N - pool[f"rank_{s}"]) for s in SCHEMES)
pool["borda_rank"] = pool["borda"].rank(ascending=False).astype(int)
rank_cols = [f"rank_{s}" for s in SCHEMES]
pool["max_rank"] = pool[rank_cols].max(axis=1)
half = max(3, N // 2)
pool["rank_robust"] = pool["max_rank"] <= half

pool = pool.sort_values("borda_rank")
out_cols = (["Solyc_ID", "direction", "intervention_type", "docking_compatible",
             "druggable_class", "is_hub", "abs_ddE", "sind_log2FC",
             "E_ddE", "E_hub", "E_evo", "E_sind"] + rank_cols +
            ["borda_rank", "max_rank", "rank_robust", "GO_BP_terms"])
pool[out_cols].to_csv(ROOT / "results/target_docking_ranking.tsv", sep="\t", index=False)

print(f"[OK] results/target_docking_ranking.tsv  (rank-robust阈值 Top{half})")
print(f"rank-robust 靶点: {pool['rank_robust'].sum()}\n")
print("=== 对接兼容主靶点 Borda 排名 (全部) ===")
print(f"{'Solyc_ID':16s}{'borda#':>6s} " + " ".join(f"{s[:5]:>6s}" for s in SCHEMES) + "  rob hub  GO")
for r in pool.itertuples():
    ranks = " ".join(f"{getattr(r,f'rank_{s}'):6d}" for s in SCHEMES)
    print(f"{r.Solyc_ID:16s}{r.borda_rank:6d} {ranks}  {'✓' if r.rank_robust else ' ':^3s} "
          f"{'H' if r.is_hub else '-'}  {(r.GO_BP_terms or '')[:32]}")

n=N; hit=n  # 池内全部方向一致
print(f"\n[真菌轴] 对接兼容池 {n} 个全部真菌方向一致(inhibit方向). ")
print("[方向自洽] 全部 DOWN(需抑制)→对接找抑制剂→逻辑闭环, 无 UP 方向反转风险.")
