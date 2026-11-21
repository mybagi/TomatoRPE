# TomatoRPE

TomatoRPE combines drought transcriptomics, substrate-guided docking, and molecular dynamics to examine tomato ribulose-5-phosphate 3-epimerase. Use the notebook to reproduce study summaries and inspect simulation measurements.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24425383.svg)](https://doi.org/10.5281/zenodo.24425383) [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/mybagi/TomatoRPE/blob/main/notebooks/TomatoRPE_quickstart.ipynb) [Data and downloads](#data)

![Study workflow for tomato drought transcriptomics and ligand screening](assets/fig1.jpg)

## Repository structure

```text
TomatoRPE/
├── README.md
├── assets/
├── notebooks/TomatoRPE_quickstart.ipynb
└── src/
```

## Run locally

### Install

Use Python 3.11 or later.

```bash
git clone https://github.com/mybagi/TomatoRPE.git
cd TomatoRPE
```

### Prepare data

The preparation command downloads the archive and verifies its checksum before extraction.

```bash
python src/prepare_data.py
```

### Reproduce the summaries

This command reports root interaction counts, target rankings, docking results, and pocket occupancy.

```bash
python src/reproduce.py
```

### Verify a trajectory

The example uses the first zinc-bound replicate of native D-ribulose 5-phosphate. The verifier uses the GROMACS run input and compressed trajectory files for the same system. [System identifiers](assets/systems.tsv)

```bash
python -m pip install -r src/requirements/md.txt
python src/download_trajectory.py --system native_D-Ru5P_with-Zn_rep1 --output trajectories
python src/verify_md_trajectory.py --trajectory-dir trajectories --system native_D-Ru5P_with-Zn_rep1
```

## Data

The [reference archive](https://github.com/mybagi/TomatoRPE/releases/download/reference-v1/tomatorpe-reference-data.zip) contains the docking tables, supplementary analyses, trajectory measurements, structural comparisons, and source dataset accessions. Extract the archive with the preparation command above to browse the tables by analysis.

The [Zenodo record](https://doi.org/10.5281/zenodo.24425383) provides the original trajectories, starting structures, and GROMACS run inputs.

## Method

Drought response contrasts inform tomato gene prioritization. Substrate-guided docking places candidate ligands in the catalytic pocket. The molecular dynamics analysis measures minimum heavy-atom distances between each ligand or substrate and the fixed 31-residue pocket. Pocket occupancy is calculated separately at thresholds of 4 and 5 Å using saved frames from 50 to 100 ns. [Analysis workflows](src/workflows/README.md)

## Contact

[Open a GitHub issue](https://github.com/mybagi/TomatoRPE/issues) for questions about the code or contact [Mingyu Bian](https://github.com/mybagi) and [Yiquan Wang](https://github.com/wyqmath). [Software citation](assets/CITATION.cff)
