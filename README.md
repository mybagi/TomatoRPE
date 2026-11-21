# TomatoRPE

This repository contains analysis code and reference tables for tomato drought responses and natural product screening against ribulose-5-phosphate 3-epimerase.

## Reproduce the summaries

Use Python 3.11 or later. Run this command from the repository root.

```bash
python src/reproduce.py
```

## Verify a trajectory

Download the matching GROMACS TPR and XTC files for one system. [Zenodo dataset](https://doi.org/10.5281/zenodo.24425383)

```bash
python -m pip install -r src/requirements/md.txt
python src/verify_md_trajectory.py --trajectory-dir /path/to/downloads --system native_D-Ru5P_with-Zn_rep1
```

The command checks frame times and minimum substrate distances to the fixed protein pocket. [Analysis workflows](src/workflows/README.md) · [Code citation](assets/CITATION.cff)
