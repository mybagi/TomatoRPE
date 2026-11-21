#!/usr/bin/env python3
"""
107_L1_fusion_prefilter.py  [5.5 漏斗 L1: 生物学先验 + 快速打分 融合]

融合方向1(RDKit/性质快速打分, 全库无偏)+ 方向3(生物学先验分层)：
全库 45万都打分, 生物学相关分子加权优先, 不排他丢弃; 双通道保底防漏黑马.

L1 融合分 = 快速打分 × 生物学先验权重
  快速打分(0-1): qed(类药) + np_likeness(类天然物,归一) + 性质适配(MW/logP/TPSA/rotB)
                 - 糖苷惩罚(contains_sugar)
  生物学先验权重(≥1): 茄科Solanum来源×1.5 / 其他植物有organisms×1.15 / 无来源×1.0
                      × 抗逆相关NP类别(Terpenoids/Alkaloids/Shikimates&Phenylpropanoids
                        /Amino acids&Peptides[含多胺])×1.2

输出 L2 子集(并集, 去重): 融合分Top N_FUSION ∪ 纯快速打分Top N_PURE
  → data/compounds/L2_candidates.tsv
"""
import os
import zipfile, io
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path(os.environ.get("TOMATO_PROJECT_ROOT", "."))
ZIP = ROOT / "data/compounds/coconut_csv-03-2025.zip"
N_FUSION = 4000   # 融合分通道
N_PURE = 1500     # 纯快速打分通道(保底黑马)

cols = ["identifier","canonical_smiles","standard_inchi_key","name",
        "molecular_weight","alogp","topological_polar_surface_area",
        "rotatable_bond_count","hydrogen_bond_donors","hydrogen_bond_acceptors",
        "qed_drug_likeliness","np_likeness","contains_sugar",
        "np_classifier_pathway","organisms","chemical_super_class"]
with zipfile.ZipFile(ZIP) as z:
    name = [n for n in z.namelist() if n.endswith(".csv")][0]
    with z.open(name) as f:
        df = pd.read_csv(io.TextIOWrapper(f,"utf-8"), usecols=cols, low_memory=False)
print(f"COCONUT 总: {len(df)}")

# --- 类药过滤 + 单片段 ---
df = df.dropna(subset=["canonical_smiles","molecular_weight"])
df = df[~df["canonical_smiles"].str.contains(r"\.", regex=True, na=False)]  # 去多片段/盐
df = df[(df.molecular_weight.between(100,600)) & (df.alogp<=6) &
        (df.hydrogen_bond_donors<=6) & (df.hydrogen_bond_acceptors<=12)]
df = df.drop_duplicates("standard_inchi_key").reset_index(drop=True)
print(f"类药+单片段+去重: {len(df)}")

# --- 快速打分 (0-1) ---
def prop_fit(r):
    s = 0
    s += 1 if 150 <= r.molecular_weight <= 450 else 0.5   # 适中MW
    s += 1 if 0 <= r.alogp <= 4 else 0.5                  # 合理亲脂
    s += 1 if r.topological_polar_surface_area <= 120 else 0.5
    s += 1 if r.rotatable_bond_count <= 8 else 0.5        # 不过柔性
    return s/4
np_norm = (df["np_likeness"] - df["np_likeness"].min()) / (df["np_likeness"].max()-df["np_likeness"].min()+1e-9)
df["fast_score"] = (0.4*df["qed_drug_likeliness"].fillna(0.3)
                    + 0.3*np_norm.fillna(0.5)
                    + 0.3*df.apply(prop_fit, axis=1))
df["fast_score"] -= 0.1*df["contains_sugar"].fillna(False).astype(int)  # 糖苷惩罚

# --- 生物学先验权重 ---
org = df["organisms"].fillna("")
w_org = np.where(org.str.contains("Solanum", case=False), 1.5,
         np.where(org.str.len() > 0, 1.15, 1.0))
STRESS_PW = ["Terpenoids","Alkaloids","Shikimates and Phenylpropanoids","Amino acids and Peptides"]
w_cls = np.where(df["np_classifier_pathway"].fillna("").isin(STRESS_PW), 1.2, 1.0)
df["bio_weight"] = w_org * w_cls
df["fusion_score"] = df["fast_score"] * df["bio_weight"]

# --- 双通道取并集 ---
top_fus = df.nlargest(N_FUSION, "fusion_score")
top_pure = df.nlargest(N_PURE, "fast_score")
sel = pd.concat([top_fus, top_pure]).drop_duplicates("standard_inchi_key")
sel = sel.sort_values("fusion_score", ascending=False)
sel["in_fusion"] = sel["standard_inchi_key"].isin(top_fus["standard_inchi_key"])
sel["in_pure"] = sel["standard_inchi_key"].isin(top_pure["standard_inchi_key"])
sel.to_csv(ROOT/"data/compounds/L2_candidates.tsv", sep="\t", index=False)

print(f"\nL2 子集(并集): {len(sel)}  [融合通道{N_FUSION} ∪ 纯快速{N_PURE}, 重叠{top_fus['standard_inchi_key'].isin(top_pure['standard_inchi_key']).sum()}]")
print(f"  茄科Solanum来源: {sel['organisms'].fillna('').str.contains('Solanum',case=False).sum()}")
print(f"  抗逆相关NP类别: {sel['np_classifier_pathway'].fillna('').isin(STRESS_PW).sum()}")
print("\n=== L2 Top10 ===")
for r in sel.head(10).itertuples():
    print(f"  {r.identifier}  fus={r.fusion_score:.2f} fast={r.fast_score:.2f} w={r.bio_weight:.2f}  "
          f"{str(r.np_classifier_pathway)[:22]}  {str(r.name)[:26]}")
print(f"\n[OK] data/compounds/L2_candidates.tsv")
