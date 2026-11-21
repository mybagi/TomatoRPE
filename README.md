# TomatoRPE

TomatoRPE combines drought transcriptomics, substrate-guided docking, and molecular dynamics to examine tomato ribulose-5-phosphate 3-epimerase. Use Colab to reproduce study summaries and inspect time series for 13 simulation systems.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24425383.svg)](https://doi.org/10.5281/zenodo.24425383) [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/mybagi/TomatoRPE/blob/main/notebooks/TomatoRPE_quickstart.ipynb)

![Study workflow for tomato drought transcriptomics and ligand screening](assets/fig1.jpg)

## Repository structure

```text
TomatoRPE/
├── README.md
├── assets/
├── notebooks/TomatoRPE_quickstart.ipynb
└── src/
    ├── prepare_data.py
    ├── reproduce.py
    ├── download_trajectory.py
    ├── verify_md_trajectory.py
    └── requirements/
```

## Run locally

Use Python 3.11 or later. The preparation command downloads and checks the reference tables. The summary command reports root interaction counts, target rankings, docking statistics, and pocket occupancy.

```bash
git clone https://github.com/mybagi/TomatoRPE.git
cd TomatoRPE
python src/prepare_data.py
python src/reproduce.py
```

### Verify a trajectory

The example checks frame times and recalculates minimum heavy-atom distances to the fixed 31-residue pocket for the first replicate of native D-ribulose 5-phosphate with catalytic zinc. The output counts frames below 4 and 5 Å from 50 to 100 ns. [System identifiers](assets/systems.tsv)

```bash
python -m pip install -r src/requirements/md.txt
python src/download_trajectory.py --system native_D-Ru5P_with-Zn_rep1 --output trajectories
python src/verify_md_trajectory.py --trajectory-dir trajectories --system native_D-Ru5P_with-Zn_rep1
```

## Data

| Material | Access |
| --- | --- |
| Study tables, time series, and source accessions | [GitHub reference archive](https://github.com/mybagi/TomatoRPE/releases/download/reference-v1/tomatorpe-reference-data.zip) |
| Original trajectories, starting structures, and GROMACS run inputs | [Zenodo record](https://doi.org/10.5281/zenodo.24425383) |
| Original transcriptomics, docking, and simulation workflows | [Archived source and input notes](https://github.com/mybagi/TomatoRPE/blob/56c6c2f613b655ad1d47f82c42b9f31a4781a397/src/workflows/README.md) |

## Contact and citation

[Open a GitHub issue](https://github.com/mybagi/TomatoRPE/issues) for code questions. [Software citation](assets/CITATION.cff)
