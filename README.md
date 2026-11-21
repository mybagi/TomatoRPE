# Tomato drought response and RPE ligand screening

This repository contains analysis code and reference tables for tomato drought responses and natural product screening against ribulose-5-phosphate 3-epimerase.

## Reproduce the analysis summaries

From the repository root, use Python 3.11 or later to reproduce selected transcriptomic, docking, and molecular dynamics summaries.

```bash
python reproduce.py
```

## Verify a trajectory

The Zenodo dataset contains 13 trajectories with starting structures and GROMACS run inputs. Download the matching TPR and XTC files for one system into a single directory. [Zenodo dataset](https://doi.org/10.5281/zenodo.24425383)

```bash
python -m pip install -r requirements-md.txt
python verify_md_trajectory.py --trajectory-dir /path/to/downloads --system native_D-Ru5P_with-Zn_rep1
```

The command checks the frame times and minimum distances between the substrate and the fixed pocket. [Workflow sources and reference tables](workflows/README.md)
