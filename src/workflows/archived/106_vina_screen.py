#!/usr/bin/env python3
"""
106_vina_screen.py  [5.5.2 高通量对接]

对 data/compounds/ligands_pdbqt/*.pdbqt 批量对接两个对接兼容靶点
(RPE, FAdesaturase), grid box 取自 results/binding_pockets.tsv.
用 vina python API, 多分子并行(60核上限).
输出: results/vina_screen_{target}.tsv (每分子最佳打分, 按亲和力排序)
"""
import os
import sys
from pathlib import Path
import pandas as pd
from vina import Vina

ROOT = Path(os.environ.get("TOMATO_PROJECT_ROOT", "."))
LIGDIR = ROOT / "data/compounds/L2_pdbqt"
POCKETS = pd.read_csv(ROOT / "results/binding_pockets.tsv", sep="\t")
REC = {"RPE_Solyc03g115820": ROOT / "data/structures/receptors/RPE.pdbqt",
       "FAdesaturase_Solyc03g063110": ROOT / "data/structures/receptors/FAdesaturase.pdbqt"}

TARGET = sys.argv[1] if len(sys.argv) > 1 else "RPE_Solyc03g115820"
EXHAUST = int(sys.argv[2]) if len(sys.argv) > 2 else 16

row = POCKETS[POCKETS["target"] == TARGET].iloc[0]
center = [row.center_x, row.center_y, row.center_z]
box = [row.box_x, row.box_y, row.box_z]
print(f"靶点 {TARGET} | 中心 {center} | 盒 {box}")

v = Vina(sf_name="vina", cpu=0, seed=42, verbosity=0)
v.set_receptor(str(REC[TARGET]))
v.compute_vina_maps(center=center, box_size=box)

ligs = sorted(LIGDIR.glob("*.pdbqt"))
print(f"配体数: {len(ligs)}")
rows = []
for i, lig in enumerate(ligs):
    try:
        v.set_ligand_from_file(str(lig))
        v.dock(exhaustiveness=EXHAUST, n_poses=5)
        e = v.energies(n_poses=1)[0][0]        # 最佳位姿亲和力 kcal/mol
        rows.append(dict(cid=lig.stem, affinity=round(float(e), 2)))
    except Exception as ex:
        rows.append(dict(cid=lig.stem, affinity=None))
    if (i + 1) % 200 == 0:
        print(f"  ...{i+1}/{len(ligs)}")

res = pd.DataFrame(rows).dropna().sort_values("affinity")
out = ROOT / f"results/vina_screen_{TARGET}.tsv"
res.to_csv(out, sep="\t", index=False)
print(f"[OK] {out}  ({len(res)} 个成功)")
print("\n=== Top10 强结合 (最负=最强) ===")
print(res.head(10).to_string(index=False))
