#!/usr/bin/env python3
"""
103_define_binding_pocket.py  [5.4.2: 活性口袋定位 → 对接grid box]

为对接兼容主靶点定位活性口袋中心+盒子, 供 5.5 AutoDock Vina grid box.
方法(无fpocket, 用biopython+保守催化残基几何法):
  RPE (ribulose-P 3-epimerase, TIM-barrel): 保守催化 motif D-x-MDG 含催化Asp,
    活性位点由 2个Asp + 金属(Fe/Zn)结合His/Asp 构成. 定位这些残基取质心=口袋中心.
  以催化残基质心为盒中心, 盒子边长取 覆盖残基簇 + 8Å 余量(min 22.5Å).
输出: results/binding_pockets.tsv (每靶点 center_xyz + box_size + 催化残基)
"""
import os
import re
from pathlib import Path
import numpy as np
from Bio.PDB import PDBParser

ROOT = Path(os.environ.get("TOMATO_PROJECT_ROOT", "."))
STRUCT = ROOT / "data/structures"

# 靶点: (pdb文件, 完整蛋白序列取自ITAG, 催化残基定位策略)
def full_seq(solyc):
    name = None; buf = []
    with open(ROOT / "data/genome/tomato_SL4.0/ITAG4.0_proteins.fasta") as f:
        for line in f:
            if line.startswith(">"):
                if name and solyc in name:
                    return "".join(buf)
                name = line[1:].split()[0]; buf = []
            else:
                buf.append(line.strip())
    return "".join(buf) if name and solyc in name else ""

TARGETS = {
    "RPE_Solyc03g115820": {
        "pdb": "RPE_Solyc03g115820_A0A3Q7FSA0.pdb",
        "solyc": "Solyc03g115820",
        # RPE 家族保守催化: DxMDG motif 的 Asp; 另有第二个保守 Asp 及金属结合 His
        "motifs": [r"D.MDG"],           # 催化 Asp (motif起始位)
        "extra_kw": "catalytic Asp/His cluster",
    },
    "FAdesaturase_Solyc03g063110": {
        "pdb": "FAdesaturase_Solyc03g063110_A0A3Q7FL53.pdb",
        "solyc": "Solyc03g063110",
        # 酰基-ACP去饱和酶: 双铁中心, 保守 (D/E)EXXH motif ×2 结合双铁
        "motifs": [r"E..H", r"D..H"],
        "extra_kw": "diiron center (E/DxxH)",
    },
}

parser = PDBParser(QUIET=True)
rows = []
for name, cfg in TARGETS.items():
    seq = full_seq(cfg["solyc"])
    struct = parser.get_structure(name, STRUCT / cfg["pdb"])
    model = struct[0]
    chain = list(model)[0]
    # residue index → CA coord (PDB残基编号从1连续, 对齐序列)
    ca = {}
    for res in chain:
        if "CA" in res:
            ca[res.id[1]] = res["CA"].coord

    # 定位催化残基位点 (序列 motif → 残基号)
    cat_pos = []
    for mt in cfg["motifs"]:
        for m in re.finditer(mt, seq):
            cat_pos.append(m.start() + 1)          # motif 起始残基(催化 D/E)
    cat_pos = sorted(set(p for p in cat_pos if p in ca))
    if not cat_pos:
        print(f"[{name}] 未定位到 motif, 回退用全蛋白质心"); cat_pos = list(ca)

    coords = np.array([ca[p] for p in cat_pos])
    center = coords.mean(axis=0)
    # 盒子: 覆盖催化残基簇范围 + 8Å 余量, 最小 22.5Å
    span = coords.max(0) - coords.min(0)
    box = np.maximum(span + 16.0, 22.5)
    rows.append(dict(target=name, solyc=cfg["solyc"], pdb=cfg["pdb"],
                     n_catalytic=len(cat_pos),
                     catalytic_res=";".join(map(str, cat_pos)),
                     center_x=round(center[0],2), center_y=round(center[1],2),
                     center_z=round(center[2],2),
                     box_x=round(box[0],1), box_y=round(box[1],1), box_z=round(box[2],1),
                     note=cfg["extra_kw"]))
    print(f"[{name}] 催化残基 {cat_pos} | 中心 {center.round(1)} | 盒 {box.round(1)}")

import pandas as pd
pd.DataFrame(rows).to_csv(ROOT / "results/binding_pockets.tsv", sep="\t", index=False)
print(f"\n[OK] results/binding_pockets.tsv ({len(rows)} 靶点)")
