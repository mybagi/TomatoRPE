"""Reproduce selected summaries from the archived analysis tables."""
from __future__ import annotations

import argparse
import json
import statistics
from fractions import Fraction
from pathlib import Path

from prepare_data import read_rows, reference_directory, verify_files


def mean_descending_ranks(scores):
    return [
        Fraction(1 + sum(other > score for other in scores))
        + Fraction(sum(other == score for other in scores) - 1, 2)
        for score in scores
    ]


def borda_ranking(tables):
    source = read_rows(tables / "S4_borda_evidence_ranks_14.tsv")
    expected = read_rows(tables / "S4_borda_tie_corrected_14.tsv")
    fields = ("root_evidence", "stele_evidence", "evolution_evidence", "inoculation_evidence")
    rank_columns = ("rank_equal", "rank_lower_stele", "rank_higher_inoculation", "rank_higher_evolution")
    # Percentile ranks in a pool of 14 have a denominator of at most 28.
    evidence = [
        [Fraction(row[field]).limit_denominator(28) for field in fields]
        for row in source
    ]
    weights = [(2, 2, 2, 2), (2, 1, 2, 2), (2, 1, 2, 4), (2, 1, 4, 2)]
    schemes = [
        mean_descending_ranks([sum(a * b for a, b in zip(row, w)) for row in evidence])
        for w in weights
    ]
    ranks = list(zip(*schemes))
    points = [4 * len(source) - sum(row) for row in ranks]
    positions = mean_descending_ranks(points)
    reference = {row["Solyc_ID"]: row for row in expected}
    for row, rank_vector, points_i, position in zip(source, ranks, points, positions):
        recorded = reference[row["Solyc_ID"]]
        if tuple(Fraction(recorded[k]) for k in rank_columns) != rank_vector:
            raise ValueError(f"Scheme rank mismatch: {row['Solyc_ID']}")
        if (Fraction(recorded["borda_points"]), Fraction(recorded["borda_rank"])) != (points_i, position):
            raise ValueError(f"Borda mismatch: {row['Solyc_ID']}")
    idx = next(i for i, row in enumerate(source) if row["Solyc_ID"] == "Solyc03g115820")
    return {"n": len(source), "rpe_rank": float(positions[idx]), "rpe_points": float(points[idx])}


def transcriptomics(tables):
    targets = read_rows(tables / "Table_S3_corrected_44_targets.tsv")
    leaf = read_rows(tables / "GSE156402_target_44_T1_log2FC.tsv")
    shrinkage = read_rows(tables / "S4_sindica_revised_44_effects.tsv")
    terms = read_rows(tables / "S4_sindica_gsea_848.tsv")
    hits = sum(float(r["log2FC_apeGLM"]) * (1 if r["expected_dir"] == "UP" else -1) > 0 for r in shrinkage)
    photo = next(row for row in terms if row["term"].endswith("|GO:0015979"))
    positive = [float(row["fdr_q-val"]) for row in terms if float(row["nes"]) > 0]
    return {
        "significant_root_interactions": len(read_rows(tables / "Supplementary_Data_1_435_root_interactions.csv", ",")),
        "core": sum(r["branch"] == "Core" for r in targets),
        "expansion": sum(r["branch"] == "Expansion" for r in targets),
        "leaf_loci_with_both_T1_estimates": sum(
            all(r[f"GSE156402_{tissue}_log2FC"] != "NA" for tissue in ("vein", "intervein"))
            for r in leaf
        ),
        "sindica_direction_matches": hits,
        "sindica_direction_total": len(shrinkage),
        "sindica_gene_fdr_below_0_05": sum(float(r["padj"]) < 0.05 for r in shrinkage),
        "photosynthesis_NES": float(photo["nes"]),
        "minimum_positive_gsea_fdr": min(positive),
    }


def docking(data, tables):
    score_pairs = read_rows(tables / "Supplementary_Data_2_Vina_score_mass_pairs.csv", ",")
    pairs = [(float(row["molecular_weight_Da"]), float(row["vina_score_kcal_per_mol"])) for row in score_pairs]
    negatives = [pair for pair in pairs if pair[1] < 0]
    candidate_rows = read_rows(data / "docking/final_elicitor_candidates_clean.tsv")
    qualified = [row for row in candidate_rows if float(row["affinity"]) < -6 and float(row["LE"]) >= 0.32]
    selected = read_rows(data / "docking/final_5_MD_candidates.tsv")
    if not {row["cid"] for row in selected}.issubset(row["cid"] for row in qualified):
        raise ValueError("A chosen ligand does not belong to the qualified pool")
    return {
        "pairs": len(pairs),
        "negative_pairs": len(negatives),
        "all_scores_pearson_r": statistics.correlation(*zip(*pairs)),
        "negative_scores_pearson_r": statistics.correlation(*zip(*negatives)),
        "prefiltered_compounds": len(candidate_rows),
        "qualified_compounds": len(qualified),
        "manually_selected_compounds": len(selected),
    }


def md_summaries(data):
    summaries = {}
    for path in sorted((data / "md_observables").glob("md_metrics_*.tsv")):
        series = read_rows(path)
        system = series[0]["system"]
        if system == "apo":
            summaries[system] = {"frames": len(series)}
            continue
        distances = [
            float(row["fixed_pocket_heavy_minimum_image_distance_nm"])
            for row in series if 50 <= float(row["time_ns"]) <= 100
        ]
        summaries[system] = {
            "frames": len(series),
            "frames_50_to_100_ns": len(distances),
            "frames_below_4_angstrom": sum(d < 0.4 for d in distances),
            "frames_below_5_angstrom": sum(d < 0.5 for d in distances),
        }
    return summaries


def verify_geo_workbook(workbook, tables):
    import numpy as np
    import pandas as pd

    raw = pd.read_excel(workbook, header=None)
    count_table = raw.iloc[3:].copy()
    count_table.columns = raw.iloc[2].astype(str)
    count_table = count_table.set_index("gene_id_short")
    suffix = ".fastq_reads_expression"
    input_columns = [col for col in count_table if col.endswith(suffix)]
    counts = count_table[input_columns].apply(pd.to_numeric, errors="coerce").fillna(0)
    counts.columns = [col.removesuffix(suffix) for col in input_columns]
    expected = pd.read_csv(tables / "GSE156402_target_44_T1_log2FC.tsv", sep="\t").set_index("Solyc_ID")
    for tissue in ("Vein", "Intervein"):
        drought = [f"{tissue}_T1_D_{i}" for i in range(1, 4)]
        watered = [f"{tissue}_T1_W_{i}" for i in range(1, 4)]
        selected = counts[drought + watered]
        cpm = selected / selected.sum() * 1e6
        effect = np.log2(cpm[drought].mean(axis=1) + 1) - np.log2(cpm[watered].mean(axis=1) + 1)
        column = f"GSE156402_{tissue.lower()}_log2FC"
        np.testing.assert_allclose(effect.reindex(expected.index), expected[column], rtol=0, atol=1e-10, equal_nan=True)
    return int(expected.dropna().shape[0])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, help="Directory of extracted reference tables")
    parser.add_argument("--geo-workbook", type=Path, help="Optional GSE156402 supplementary workbook from GEO")
    args = parser.parse_args()
    data = reference_directory(args.data)
    tables = data / "manuscript_tables"
    output = {
        "verified_files": verify_files(data),
        "transcriptomics": transcriptomics(tables),
        "borda": borda_ranking(tables),
        "docking": docking(data, tables),
        "molecular_dynamics": md_summaries(data),
    }
    if args.geo_workbook:
        output["geo_recomputed_paired_loci"] = verify_geo_workbook(args.geo_workbook, tables)
    print(json.dumps(output, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
