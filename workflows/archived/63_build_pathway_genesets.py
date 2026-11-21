"""
63_build_pathway_genesets.py
Section 4.1 — 构建 8 条真实番茄通路基因集 (本地 GO, 取代 script 50 随机占位)

复用 3NT.3 验证的方法: 番茄 GAF (gene_association.itag4.0.nr) 直接注释 →
  go-basic.obo is_a/part_of 祖先传播 → 每个通路 GO 根节点的基因集 (限 M82 22099 背景)。

8 通路 (对齐 proposal 2.3.2 + 3NT.3 发现):
  防御/胁迫 (预期受旱↑): ABA / JA / Osmotic / ROS / Phenylpropanoid-Suberin
  生长 (预期受旱↓):      Auxin / Cytokinin / CellWall

输出:
  - results/pathway_genesets_tomato.json   (通路→Solyc 列表)
  - results/pathway_genesets_tomato.gmt     (标准 GMT)
  - results/pathway_genesets_summary.tsv    (通路/GO/基因数)
  - logs/63_build_pathway_genesets_<date>.log
非破坏性。
"""
import json
import os
import re
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path

ROOT = Path(os.environ.get("TOMATO_PROJECT_ROOT", "."))
GAF = ROOT / "gene_association.itag4.0.nr"
OBO = ROOT / "data/go-basic.obo"
MAT = ROOT / "data/final_clean_matrices/PRJNA800740_M82_Mechanism.tsv"
OUT_JSON = ROOT / "results/pathway_genesets_tomato.json"
OUT_GMT = ROOT / "results/pathway_genesets_tomato.gmt"
OUT_SUM = ROOT / "results/pathway_genesets_summary.tsv"
LOG = ROOT / f"logs/63_build_pathway_genesets_{datetime.now():%Y%m%d}.log"

# 通路 → GO 根节点 (主选, 已验证覆盖 77-339)
PATHWAYS = {
    "Defense_ABA":            ("GO:0009737", "response to abscisic acid", "defense"),
    "Defense_JA":             ("GO:0009753", "response to jasmonic acid", "defense"),
    "Defense_Osmotic":        ("GO:0006970", "response to osmotic stress", "defense"),
    "Defense_ROS":            ("GO:0006979", "response to oxidative stress", "defense"),
    "Defense_Phenylpropanoid":("GO:0009699", "phenylpropanoid biosynthetic process", "defense"),
    "Growth_Auxin":           ("GO:0009733", "response to auxin", "growth"),
    "Growth_Cytokinin":       ("GO:0009735", "response to cytokinin", "growth"),
    "Growth_CellWall":        ("GO:0042546", "plant-type cell wall biogenesis", "growth"),
}


class Tee:
    def __init__(self, *f): self.files = f
    def write(self, o):
        for x in self.files: x.write(o); x.flush()
    def flush(self):
        for x in self.files: x.flush()


def sv(g): return re.sub(r"\.\d+$", "", g)


def main():
    LOG.parent.mkdir(parents=True, exist_ok=True)
    logf = open(LOG, "w")
    sys.stdout = Tee(sys.__stdout__, logf)
    print(f"启动: {datetime.now():%Y-%m-%d %H:%M:%S}")

    hdr = open(MAT).readline().rstrip("\n").split("\t")
    mat = set(hdr[1:])
    print(f"M82 背景基因: {len(mat)}")

    gene2go = defaultdict(set)
    for line in open(GAF):
        if line.startswith("!"): continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 9: continue
        g = sv(p[1])
        if g in mat:
            gene2go[g].add(p[4])

    parents = defaultdict(set)
    cur = None
    for line in open(OBO):
        line = line.rstrip("\n")
        if line == "[Term]": cur = None
        elif line.startswith("id: GO:"): cur = line[4:]
        elif line.startswith("is_a:") and cur: parents[cur].add(line.split()[1])
        elif line.startswith("relationship: part_of") and cur:
            m = re.search(r"GO:\d+", line)
            if m: parents[cur].add(m.group(0))
    cache = {}
    def anc(t):
        if t in cache: return cache[t]
        seen, st = {t}, list(parents.get(t, ()))
        while st:
            x = st.pop()
            if x not in seen: seen.add(x); st.extend(parents.get(x, ()))
        cache[t] = seen; return seen
    t2g = defaultdict(set)
    for g, gos in gene2go.items():
        e = set()
        for t in gos: e |= anc(t)
        for t in e: t2g[t].add(g)

    genesets = {}
    summ = []
    for name, (go, desc, grp) in PATHWAYS.items():
        genes = sorted(t2g.get(go, set()))
        genesets[name] = genes
        summ.append((name, grp, go, desc, len(genes)))
        print(f"  {name:26} {go} n={len(genes):4} ({grp})")

    OUT_JSON.write_text(json.dumps(genesets, indent=2))
    with open(OUT_GMT, "w") as f:
        for name, (go, desc, grp) in PATHWAYS.items():
            f.write("\t".join([name, f"{go}|{desc}|{grp}"] + genesets[name]) + "\n")
    with open(OUT_SUM, "w") as f:
        f.write("Pathway\tGroup\tGO_ID\tGO_name\tn_genes\n")
        for r in summ:
            f.write("\t".join(map(str, r)) + "\n")
    print(f"\n输出: {OUT_JSON.name} / {OUT_GMT.name} / {OUT_SUM.name}")
    print(f"完成: {datetime.now():%Y-%m-%d %H:%M:%S}")
    sys.stdout = sys.__stdout__
    logf.close()


if __name__ == "__main__":
    main()
