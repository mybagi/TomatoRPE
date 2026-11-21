#!/usr/bin/env python3
"""
105_prepare_ligands.py  [5.5.1 配体库准备 + 预过滤]

COCONUT 70万天然产物全量对接不现实, 先做类药性/性质预过滤缩库:
  - 有效 SMILES + RDKit 可解析
  - 分子量 100-600 (农用小分子合理区间)
  - Lipinski 类药性宽松版 (MW<600, logP<6, HBD<=6, HBA<=12)
  - 去重(InChIKey)
过滤后 → RDKit 生成3D(ETKDG+MMFF优化) → meeko 转 PDBQT.
分批输出到 data/compounds/ligands_pdbqt/ (每个分子一个pdbqt).

输出: data/compounds/coconut_filtered.tsv (过滤后清单)
      data/compounds/ligands_pdbqt/*.pdbqt
"""
import os
import sys, zipfile, io
from pathlib import Path
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, AllChem
from rdkit import RDLogger
RDLogger.DisableAllLogging()

ROOT = Path(os.environ.get("TOMATO_PROJECT_ROOT", "."))
ZIP = ROOT / "data/compounds/coconut_csv-03-2025.zip"
OUTDIR = ROOT / "data/compounds/ligands_pdbqt"
OUTDIR.mkdir(parents=True, exist_ok=True)
MAXN = int(sys.argv[1]) if len(sys.argv) > 1 else 5000   # 预过滤后取前N个做3D(可调)

def find_smiles_col(cols):
    for c in cols:
        if c.lower() in ("canonical_smiles", "smiles", "structure_smiles"):
            return c
    return None

# 读 zip 内 csv
with zipfile.ZipFile(ZIP) as z:
    csvname = [n for n in z.namelist() if n.endswith(".csv")][0]
    print(f"读 {csvname}")
    with z.open(csvname) as f:
        df = pd.read_csv(io.TextIOWrapper(f, "utf-8"), low_memory=False)
print(f"COCONUT 总分子: {len(df)}  列: {list(df.columns)[:8]}...")

smi_col = find_smiles_col(df.columns)
id_col = "identifier" if "identifier" in df.columns else df.columns[0]
print(f"SMILES列={smi_col}  ID列={id_col}")

kept = []
seen = set()
for r in df.itertuples():
    smi = getattr(r, smi_col, None)
    if not isinstance(smi, str) or not smi:
        continue
    m = Chem.MolFromSmiles(smi)
    if m is None:
        continue
    mw = Descriptors.MolWt(m)
    if not (100 <= mw <= 600):
        continue
    logp = Descriptors.MolLogP(m)
    hbd = Lipinski.NumHDonors(m)
    hba = Lipinski.NumHAcceptors(m)
    if logp > 6 or hbd > 6 or hba > 12:
        continue
    ik = Chem.MolToInchiKey(m)
    if ik in seen:
        continue
    seen.add(ik)
    kept.append(dict(cid=getattr(r, id_col), smiles=Chem.MolToSmiles(m),
                     inchikey=ik, mw=round(mw, 1), logp=round(logp, 2),
                     hbd=hbd, hba=hba))

fdf = pd.DataFrame(kept)
fdf.to_csv(ROOT / "data/compounds/coconut_filtered.tsv", sep="\t", index=False)
print(f"预过滤后: {len(fdf)} 个类药天然产物")

# 生成3D + PDBQT (前 MAXN 个)
try:
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    have_meeko = True
except Exception as e:
    print(f"meeko 导入失败: {e}"); have_meeko = False

n_ok = 0
for r in fdf.head(MAXN).itertuples():
    out = OUTDIR / f"{r.cid}.pdbqt"
    if out.exists():
        n_ok += 1; continue
    m = Chem.MolFromSmiles(r.smiles)
    if m is None:
        continue
    m = Chem.AddHs(m)
    if AllChem.EmbedMolecule(m, AllChem.ETKDGv3()) != 0:
        continue
    try:
        AllChem.MMFFOptimizeMolecule(m)
    except Exception:
        pass
    if have_meeko:
        prep = MoleculePreparation()
        setups = prep.prepare(m)
        pdbqt, ok, _ = PDBQTWriterLegacy.write_string(setups[0])
        if ok:
            out.write_text(pdbqt); n_ok += 1
print(f"[OK] 生成 {n_ok} 个配体 PDBQT → {OUTDIR}")
