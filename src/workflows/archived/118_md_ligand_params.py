#!/usr/bin/env python3
"""
118_md_ligand_params.py  [7.1a MD配体参数化: pose→GROMACS拓扑]

关键(第6节自查发现的隐患): 配体必须用**对接pose坐标**(结合在RPE口袋),
不能从SMILES重新生成构象.
流程: pose PDBQT --(meeko mk_export恢复全氢,保留坐标)--> SDF
     --> RDKit读+加氢补全 --> acpype(GAFF力场) --> GROMACS itp/gro/top
输出: data/md/ligands/{cid}/ 每个候选的 GMX 拓扑
用法: python 118_md_ligand_params.py [cid]  (不给则做全部5个)
"""
import os
import sys, subprocess, os
from pathlib import Path
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger; RDLogger.logger().setLevel(RDLogger.CRITICAL)

ROOT = Path(os.environ.get("TOMATO_PROJECT_ROOT", "."))
POSEDIR = ROOT/"data/structures/poses"
OUTBASE = ROOT/"data/md/ligands"; OUTBASE.mkdir(parents=True, exist_ok=True)

f5 = pd.read_csv(ROOT/"results/final_5_MD_candidates.tsv", sep="\t")
smi_map = dict(zip(f5.cid, f5.smiles))
charge_map = {}
for r in f5.itertuples():
    m = Chem.MolFromSmiles(r.smiles)
    charge_map[r.cid] = Chem.GetFormalCharge(m) if m else 0

targets = [sys.argv[1]] if len(sys.argv) > 1 else list(f5.cid)

for cid in targets:
    outdir = OUTBASE/cid; outdir.mkdir(exist_ok=True)
    pose = POSEDIR/f"{cid}_pose.pdbqt"
    print(f"\n=== {cid} (net charge {charge_map[cid]}) ===")
    # 1) meeko 导出带氢 SDF(保留pose坐标)
    sdf = outdir/f"{cid}.sdf"
    r = subprocess.run(["mk_export.py", str(pose), "-s", str(sdf)],
                       capture_output=True, text=True)
    if not sdf.exists() or sdf.stat().st_size == 0:
        print(f"  meeko导出失败: {r.stderr[:200]}"); continue
    # 2) RDKit 读回, 确认原子/加氢完整
    mol = Chem.MolFromMolFile(str(sdf), removeHs=False)
    if mol is None:
        # SDF可能多构象, 用SDMolSupplier
        mol = next(Chem.SDMolSupplier(str(sdf), removeHs=False))
    if mol is None:
        print("  RDKit读取失败"); continue
    nH = sum(1 for a in mol.GetAtoms() if a.GetSymbol()=="H")
    print(f"  pose→SDF: {mol.GetNumAtoms()}原子(含{nH}H)")
    # 保存为mol给acpype
    molf = outdir/f"{cid}.mol"
    Chem.MolToMolFile(mol, str(molf))
    # 3) acpype GAFF 参数化
    os.chdir(outdir)
    r = subprocess.run(["acpype","-i",str(molf),"-b",cid,"-n",str(charge_map[cid])],
                       capture_output=True, text=True, timeout=300)
    ac = outdir/f"{cid}.acpype"
    gro = ac/f"{cid}_GMX.gro"; itp = ac/f"{cid}_GMX.itp"
    if gro.exists() and itp.exists():
        print(f"  ✓ acpype成功: {gro.name}, {itp.name}")
    else:
        print(f"  ✗ acpype失败: {r.stdout[-300:]} {r.stderr[-200:]}")
    os.chdir(ROOT)

print("\n[完成] 配体拓扑在 data/md/ligands/{cid}/{cid}.acpype/")
