#!/usr/bin/env python3
"""
109_vina_parallel.py  [5.5 L2/L3 对接: 配体级60核并行]

修复 106 的低效(配体串行). 每个 worker 单核 Vina, 60进程并行处理配体.
对 L3_screen_set.tsv 的 2000 精选配体对接指定靶点, exhaust=4(快筛档).
用法: python 109_vina_parallel.py <target> <exhaust> <n_workers>
输出: results/vina_dock_{target}.tsv (按亲和力排序)
"""
import os
import sys, os
from pathlib import Path
from multiprocessing import Pool
import pandas as pd

ROOT = Path(os.environ.get("TOMATO_PROJECT_ROOT", "."))
LIGDIR = ROOT / "data/compounds/L2_pdbqt"
POCKETS = pd.read_csv(ROOT / "results/binding_pockets.tsv", sep="\t")
REC = {"RPE_Solyc03g115820": ROOT/"data/structures/receptors/RPE.pdbqt",
       "FAdesaturase_Solyc03g063110": ROOT/"data/structures/receptors/FAdesaturase.pdbqt"}

TARGET = sys.argv[1] if len(sys.argv) > 1 else "RPE_Solyc03g115820"
EXHAUST = int(sys.argv[2]) if len(sys.argv) > 2 else 4
NWORK = int(sys.argv[3]) if len(sys.argv) > 3 else 60
row = POCKETS[POCKETS.target == TARGET].iloc[0]
CENTER = [row.center_x, row.center_y, row.center_z]
BOX = [row.box_x, row.box_y, row.box_z]

def dock_one(cid):
    from vina import Vina
    lig = LIGDIR / f"{cid}.pdbqt"
    if not lig.exists():
        return (cid, None)
    try:
        v = Vina(sf_name="vina", cpu=1, seed=42, verbosity=0)
        v.set_receptor(str(REC[TARGET]))
        v.compute_vina_maps(center=CENTER, box_size=BOX)
        v.set_ligand_from_file(str(lig))
        v.dock(exhaustiveness=EXHAUST, n_poses=3)
        return (cid, round(float(v.energies(n_poses=1)[0][0]), 2))
    except Exception:
        return (cid, None)

if __name__ == "__main__":
    cids = pd.read_csv(ROOT/"data/compounds/L3_screen_set.tsv", sep="\t")["identifier"].tolist()
    print(f"[{TARGET}] {len(cids)}配体 exhaust={EXHAUST} {NWORK}核并行 中心{CENTER} 盒{BOX}", flush=True)
    with Pool(NWORK) as p:
        res = p.map(dock_one, cids, chunksize=5)
    df = pd.DataFrame([r for r in res if r[1] is not None], columns=["cid","affinity"]).sort_values("affinity")
    out = ROOT / f"results/vina_dock_{TARGET}.tsv"
    df.to_csv(out, sep="\t", index=False)
    print(f"[OK] {out} 成功{len(df)}/{len(cids)}", flush=True)
    print(df.head(15).to_string(index=False), flush=True)
