# TomatoRPE

TomatoRPE links drought transcriptomics with substrate-guided docking and molecular dynamics of tomato ribulose-5-phosphate 3-epimerase. The Colab notebook reproduces selected study summaries and plots pocket distances or protein displacement.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24425383.svg)](https://doi.org/10.5281/zenodo.24425383) [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/mybagi/TomatoRPE/blob/main/notebooks/TomatoRPE_quickstart.ipynb)

![Workflow from tomato drought transcriptomics to ligand screening](assets/fig1.jpg)

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

Download and check the reference tables.

```bash
python src/prepare_data.py
```

### Reproduce the summaries

```bash
python src/reproduce.py
```

This command reports root interaction counts, candidate rankings, docking summaries, and pocket occupancy.

### Verify one trajectory

The example uses the first zinc-bound replicate of native D-ribulose 5-phosphate. Download that system or use the matching files from Zenodo.

```bash
python src/download_trajectory.py --system native_D-Ru5P_with-Zn_rep1 --output trajectories
python -m pip install -r src/requirements/md.txt
python src/verify_md_trajectory.py --trajectory-dir trajectories --system native_D-Ru5P_with-Zn_rep1
```

## Data

The [reference data package](https://github.com/mybagi/TomatoRPE/releases/download/reference-v1/tomatorpe-reference-data.zip) contains supplementary data 1 to 4, the source accession table and the analysis inputs for the examples. The Zenodo record linked above provides the original GROMACS trajectories and matching files. [System identifiers](assets/systems.tsv)

## Method

The workflow uses drought response contrasts to prioritize tomato genes and screens natural products against a substrate-defined RPE pocket. For each ligand or substrate trajectory, pocket occupancy is the fraction of saved frames below the 4 or 5 Å distance threshold between 50 and 100 ns. [Analysis workflows](src/workflows/README.md)

## Contact

For questions about the code, [open an issue](https://github.com/mybagi/TomatoRPE/issues) or contact [Mingyu Bian](https://github.com/mybagi) and [Yiquan Wang](https://github.com/wyqmath). [Code citation](assets/CITATION.cff)
