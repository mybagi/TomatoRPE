# Analysis workflows

These scripts document the transcriptomic analyses, docking calculations, and molecular dynamics workflows. Set `TOMATO_PROJECT_ROOT` to the directory containing the original inputs and intermediate files. Source dataset accessions and structural comparison tables are available in the [reference archive](https://github.com/mybagi/TomatoRPE/releases/download/reference-v1/tomatorpe-reference-data.zip).

The archived three-axis ranking omitted quantitative trait locus weights because the parser read the overlap count column. The reference tables retain the results of that run. [Composite scoring script](archived/69_composite_3axis_scoring.py)

The paired molecular dynamics workflow prepares the initial 10 ns runs for ribulose-5-phosphate 3-epimerase with native D-ribulose 5-phosphate, with or without catalytic zinc. Production trajectories extend to 100 ns from the 10 ns checkpoints. The verifier compares frame times and minimum ligand or substrate distances to the fixed protein pocket with the archived measurements. [Zenodo trajectories](https://doi.org/10.5281/zenodo.24425383)
