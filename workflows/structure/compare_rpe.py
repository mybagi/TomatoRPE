"""Compare full-length tomato RPE predictions in fixed residue windows."""

import csv
import hashlib
import json
from pathlib import Path
import sys
from urllib.request import urlopen

WORK = Path(__file__).resolve().parent
sys.path.insert(0, str(WORK / "analysis_packages"))

import numpy as np
from Bio.Align import PairwiseAligner
from Bio.PDB import MMCIFParser, PDBParser, Superimposer, is_aa
from Bio.SeqUtils import seq1

INPUTS = WORK / "inputs"
MODELS = WORK / "predictions"
ORIGINAL_AF2_SHA256 = "ca9a4760271470c6a88c864e2ef78225c71e5a947de4005182bd01033de2d741"
QUARTET = (89, 91, 122, 233)


def read_model(path):
    parser = MMCIFParser(QUIET=True) if path.suffix.lower() == ".cif" else PDBParser(QUIET=True)
    return next(parser.get_structure(path.stem, str(path)).get_models())


def protein_residues(chain):
    return [res for res in chain if is_aa(res, standard=True)]


def sequence(residues):
    return "".join(seq1(res.get_resname()) for res in residues)


def read_tomato(path, expected=None):
    model = read_model(path)
    chains = [protein_residues(chain) for chain in model]
    chains = [chain for chain in chains if chain]
    if len(chains) != 1 or any(res.id[0] not in (" ", "W") for res in model.get_residues()):
        raise ValueError(f"{path.name}: expected one protein chain without ligand or metal atoms")
    residues = chains[0]
    seq = sequence(residues)
    if len(seq) != 282 or "".join(seq[i - 1] for i in QUARTET) != "HDHD":
        raise ValueError(f"{path.name}: unexpected tomato RPE sequence or quartet")
    if expected is not None and seq != expected:
        raise ValueError(f"{path.name}: the sequence differs from the project's AF2 sequence")
    if any("CA" not in res for res in residues):
        raise ValueError(f"{path.name}: missing C-alpha atoms")
    return residues, seq


def align_reference(tomato_seq, reference):
    aligner = PairwiseAligner(
        mode="global", match_score=2, mismatch_score=-1, open_gap_score=-10, extend_gap_score=-0.5,
        end_gap_score=0,
    )
    alignment = aligner.align(tomato_seq, sequence(reference))[0]
    pairs = [
        (i, j)
        for (t_start, t_end), (r_start, r_end) in zip(*alignment.aligned)
        for i, j in zip(range(t_start, t_end), range(r_start, r_end), strict=True)
    ]
    if pairs != list(zip(range(48, 278), range(230), strict=True)):
        raise ValueError("1RPX alignment does not match the established 230-residue reference mapping")
    if sum(tomato_seq[i] == sequence(reference)[j] for i, j in pairs) != 226:
        raise ValueError("1RPX sequence identity differs from the established reference")
    return pairs


def compare(moving, fixed, pairs):
    superimposer = Superimposer()
    superimposer.set_atoms([fixed[j]["CA"] for i, j in pairs], [moving[i]["CA"] for i, j in pairs])
    rotation, translation = superimposer.rotran
    mapping = dict(pairs)
    squared = []
    for pos in QUARTET:
        lhs = moving[pos - 1]
        rhs = fixed[mapping[pos - 1]]
        if lhs.resname != rhs.resname:
            raise ValueError(f"Catalytic residue {pos} mismatch")
        for atom in lhs:
            if atom.element not in {"H", "D"} and atom.name in rhs:
                displacement = atom.coord @ rotation + translation - rhs[atom.name].coord
                squared.append(displacement @ displacement)
    if len(squared) != 36:
        raise ValueError(f"Expected 36 common quartet heavy atoms, got {len(squared)}")
    return superimposer.rms, float(np.sqrt(np.mean(squared))), len(squared)


def confidence(residues):
    values = np.array([res["CA"].bfactor for res in residues])
    if not np.isfinite(values).all() or values.min() < 0.0 or values.max() > 100.0 or not values.any():
        return None, None
    if values.max() <= 1.0:
        values *= 100.0
    return float(values.mean()), float(values[47:279].mean())


def prepare_inputs():
    INPUTS.mkdir(exist_ok=True)
    MODELS.mkdir(exist_ok=True)
    original = INPUTS / "RPE_AF.pdb"
    if not original.is_file():
        raise FileNotFoundError(f"Provide the archived AlphaFold2 receptor at {original}")
    if hashlib.sha256(original.read_bytes()).hexdigest() != ORIGINAL_AF2_SHA256:
        raise ValueError("AF2 reference differs from the original project file")
    _, seq = read_tomato(original)
    (INPUTS / "RPE_A0A3Q7FSA0.fasta").write_text(
        f">Solyc03g115820|A0A3Q7FSA0|tomato_RPE|282aa\n{seq}\n", encoding="ascii"
    )
    job = [{
        "name": "Tomato_RPE_A0A3Q7FSA0_AF3",
        "modelSeeds": [],
        "sequences": [{"proteinChain": {"sequence": seq, "count": 1}}],
        "dialect": "alphafoldserver",
        "version": 1,
    }]
    (INPUTS / "RPE_AF3_server_job.json").write_text(json.dumps(job, indent=2) + "\n", encoding="ascii")
    crystal = INPUTS / "1RPX.pdb"
    if not crystal.exists():
        with urlopen("https://files.rcsb.org/download/1RPX.pdb", timeout=20) as response:
            data = response.read()
        if data.splitlines()[0][62:66] != b"1RPX":
            raise ValueError("RCSB response does not contain PDB 1RPX")
        crystal.write_bytes(data)
    return original, crystal, seq


def main():
    original, crystal, seq = prepare_inputs()
    af2, _ = read_tomato(original, seq)
    references = []
    for chain in "ABC":
        residues = protein_residues(read_model(crystal)[chain])
        references.append((f"1RPX_{chain}", residues, align_reference(seq, residues)))
    self_pairs = [(i, i) for i in range(47, 279)]
    paths = {"AlphaFold2": original}
    for name in ("AF3", "AF3_no_template", "ESMFold2"):
        options = [p for ext in (".cif", ".pdb") if (p := MODELS / f"RPE_{name}{ext}").exists()]
        if len(options) > 1:
            raise ValueError(f"Keep only the selected coordinate file for {name}")
        if options:
            paths[name] = options[0]
    rows = []
    signatures = {str(crystal.relative_to(WORK)): hashlib.sha256(crystal.read_bytes()).hexdigest()}
    for label, path in paths.items():
        model, _ = read_tomato(path, seq)
        signatures[str(path.relative_to(WORK))] = hashlib.sha256(path.read_bytes()).hexdigest()
        full, core = confidence(model)
        comparisons = references if label == "AlphaFold2" else references + [("AlphaFold2_48-279", af2, self_pairs)]
        for reference, fixed, pairs in comparisons:
            ca_rmsd, q_rmsd, atoms = compare(model, fixed, pairs)
            rows.append([label, reference, len(pairs), atoms, f"{ca_rmsd:.6f}", f"{q_rmsd:.6f}",
                         "" if full is None else f"{full:.6f}", "" if core is None else f"{core:.6f}"])
    baseline = [(round(float(r[4]), 3), round(float(r[5]), 3)) for r in rows[:3]]
    if baseline != [(0.755, 0.313), (0.731, 0.308), (0.726, 0.333)]:
        raise ValueError(f"AF2 reproducibility baseline differs from the prior measured results: {baseline}")
    with (WORK / "RPE_structure_RMSD.tsv").open("w", encoding="utf-8", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(["model", "reference", "paired_Calpha", "quartet_heavy_atoms", "Calpha_RMSD_A",
                         "quartet_RMSD_A", "full_chain_mean_pLDDT", "core_48_279_mean_pLDDT"])
        writer.writerows(rows)
    (WORK / "coordinate_sha256.json").write_text(json.dumps(signatures, indent=2) + "\n", encoding="utf-8")
    print(f"RPE length: {len(seq)}. Completed coordinate models: {', '.join(paths)}.")
    for row in rows:
        print("\t".join(str(value) for value in row))


if __name__ == "__main__":
    main()
