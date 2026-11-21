#!/usr/bin/env python3
"""
108_gen_L2_ligands.py  [5.5 L2 配体制备: 并行 3D+PDBQT]

对 L2_candidates.tsv 的 4895 分子并行生成 3D(ETKDGv3+MMFF) → meeko PDBQT.
多进程(60核上限). 跳过多片段/失败分子, 如实记录成功率.
输出: data/compounds/L2_pdbqt/{identifier}.pdbqt
"""
import os
import warnings
from pathlib import Path
from multiprocessing import Pool
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger
RDLogger.logger().setLevel(RDLogger.CRITICAL)
warnings.filterwarnings("ignore")

ROOT = Path(os.environ.get("TOMATO_PROJECT_ROOT", "."))
OUTDIR = ROOT / "data/compounds/L2_pdbqt"
OUTDIR.mkdir(parents=True, exist_ok=True)

def one(args):
    cid, smi = args
    out = OUTDIR / f"{cid}.pdbqt"
    if out.exists():
        return (cid, True)
    try:
        from meeko import MoleculePreparation, PDBQTWriterLegacy
        m = Chem.MolFromSmiles(smi)
        if m is None or len(Chem.GetMolFrags(m)) != 1:
            return (cid, False)
        m = Chem.AddHs(m)
        if AllChem.EmbedMolecule(m, AllChem.ETKDGv3()) != 0:
            return (cid, False)
        try: AllChem.MMFFOptimizeMolecule(m)
        except Exception: pass
        prep = MoleculePreparation()
        s = prep.prepare(m)
        pq, ok, _ = PDBQTWriterLegacy.write_string(s[0])
        if ok:
            out.write_text(pq)
            return (cid, True)
    except Exception:
        pass
    return (cid, False)

if __name__ == "__main__":
    df = pd.read_csv(ROOT / "data/compounds/L2_candidates.tsv", sep="\t")
    tasks = list(zip(df["identifier"], df["canonical_smiles"]))
    print(f"L2 配体: {len(tasks)}, 60核并行生成 PDBQT...")
    with Pool(60) as p:
        res = p.map(one, tasks, chunksize=20)
    ok = sum(1 for _, o in res if o)
    print(f"[OK] 成功 {ok}/{len(tasks)} → {OUTDIR}")
    fail = [c for c, o in res if not o]
    pd.Series(fail).to_csv(ROOT/"data/compounds/L2_ligand_failed.txt", index=False, header=False)
    print(f"失败 {len(fail)} (多片段/embed失败, 已记录)")
