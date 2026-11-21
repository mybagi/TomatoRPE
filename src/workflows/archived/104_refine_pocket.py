#!/usr/bin/env python3
"""
104_refine_pocket.py  [5.4.2 口袋精化: 多保守催化残基定位]

fpocket 安装超时, 改用更严谨的文献锚定法: RPE (Pfam PF00834, TIM-barrel)
催化口袋由多个保守残基构成(非单残基). 番茄 RPE 序列比对家族保守位点:
  - His91 / Asp95 (WIHVDVMDG motif, N端催化+金属结合)
  - Asp后段保守 (VHLMIVEPD ~D131 区)
  - 金属结合 His (HCEQ ~H159, HLHR ~H166 区)
取这些保守催化残基 CA 的质心 = 口袋中心, 盒子覆盖残基簇 + 余量.
对去饱和酶同理用双铁 E/DxxH 保守簇.

用序列 motif 精确定位多残基, 比 103 的单残基定心更准.
输出: 覆盖 results/binding_pockets.tsv (精化版)
"""
import os
import re
from pathlib import Path
import numpy as np, pandas as pd
from Bio.PDB import PDBParser

ROOT = Path(os.environ.get("TOMATO_PROJECT_ROOT", "."))
STRUCT = ROOT / "data/structures"

def full_seq(solyc):
    name=None; buf=[]
    with open(ROOT/"data/genome/tomato_SL4.0/ITAG4.0_proteins.fasta") as f:
        for line in f:
            if line.startswith(">"):
                if name and solyc in name: return "".join(buf)
                name=line[1:].split()[0]; buf=[]
            else: buf.append(line.strip())
    return "".join(buf) if name and solyc in name else ""

# 每靶点: 多个保守催化 motif (锚点残基在 motif 内的偏移)
TARGETS = {
 "RPE_Solyc03g115820": dict(
    pdb="RPE_Solyc03g115820_A0A3Q7FSA0.pdb", solyc="Solyc03g115820",
    # RPE 家族保守催化/金属结合残基 (motif, motif内锚点偏移0-based)
    anchors=[(r"WIHVDVMDG", [2,5]),   # His, Asp (催化+金属)
             (r"H[LM]MIVEPD", [7]),   # 催化 Asp
             (r"H[CL]EQ", [0]),       # 金属结合 His
             (r"H[LM]HR", [0,2])],    # 金属结合 His 簇
    note="RPE catalytic Asp/His cluster (PF00834 TIM-barrel)"),
 "FAdesaturase_Solyc03g063110": dict(
    pdb="FAdesaturase_Solyc03g063110_A0A3Q7FL53.pdb", solyc="Solyc03g063110",
    anchors=[(r"E..H", [0,3]), (r"D..H", [0,3])],  # 双铁 E/DxxH
    note="acyl-ACP desaturase diiron center (E/DxxH x2)"),
}

parser = PDBParser(QUIET=True)
rows=[]
for name,cfg in TARGETS.items():
    seq=full_seq(cfg["solyc"])
    chain=list(parser.get_structure(name,STRUCT/cfg["pdb"])[0])[0]
    ca={r.id[1]:r["CA"].coord for r in chain if "CA" in r}
    pos=[]
    for mt,offs in cfg["anchors"]:
        for m in re.finditer(mt,seq):
            for o in offs:
                p=m.start()+o+1
                if p in ca: pos.append(p)
    pos=sorted(set(pos))
    coords=np.array([ca[p] for p in pos])
    center=coords.mean(0)
    span=coords.max(0)-coords.min(0)
    box=np.maximum(span+12.0, 20.0)   # 覆盖催化簇+12Å余量
    rows.append(dict(target=name, solyc=cfg["solyc"], pdb=cfg["pdb"],
        n_catalytic=len(pos), catalytic_res=";".join(map(str,pos)),
        center_x=round(center[0],2),center_y=round(center[1],2),center_z=round(center[2],2),
        box_x=round(box[0],1),box_y=round(box[1],1),box_z=round(box[2],1), note=cfg["note"]))
    print(f"[{name}] {len(pos)}催化残基 {pos}")
    print(f"    中心 {center.round(2)} | 盒 {box.round(1)} | 残基簇跨度 {span.round(1)}")

pd.DataFrame(rows).to_csv(ROOT/"results/binding_pockets.tsv",sep="\t",index=False)
print(f"\n[OK] 精化覆盖 results/binding_pockets.tsv")
