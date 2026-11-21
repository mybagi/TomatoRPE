# Analysis workflows

The archived scripts cover transcriptomic analysis, docking, and ligand preparation. They require source datasets and intermediate files from the original analyses. Set `TOMATO_PROJECT_ROOT` to the corresponding project directory before running the applicable scripts. The source accessions are listed in the [reference data package](https://github.com/mybagi/TomatoRPE/releases/download/reference-v1/tomatorpe-reference-data.zip).

The archived three-axis ranking omitted quantitative trait locus weights because the parser read the count column. The reference tables reflect the archived run. [Composite scoring script](archived/69_composite_3axis_scoring.py)

The paired molecular dynamics workflow sets up the 10 ns runs and extends them from their checkpoints to 100 ns. The trajectory verifier compares the published trajectories with the archived minimum-distance series. [Published trajectory files](https://doi.org/10.5281/zenodo.24425383)

The structural comparison tables record the coordinate checksums and root mean square deviations for the model comparisons. Both are included in the reference data package.
