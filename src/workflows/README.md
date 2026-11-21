# Analysis workflows

The archived scripts cover transcriptomic analysis, docking, and ligand preparation. They require source datasets and intermediate files from the original analyses. Set `TOMATO_PROJECT_ROOT` to the corresponding project directory before running the applicable scripts. [Source dataset accessions](../../assets/source_accessions.tsv)

The archived three-axis ranking omitted quantitative trait locus weights because the parser read the count column. The bundled analysis tables reflect the archived run. [Composite scoring script](archived/69_composite_3axis_scoring.py)

The paired molecular dynamics workflow sets up the 10 ns runs and extends them from their checkpoints to 100 ns. The trajectory verifier compares the published trajectories with the archived minimum-distance series. [Published trajectory files](https://doi.org/10.5281/zenodo.24425383)

The structural comparison tables record the coordinate checksums and root mean square deviations for the model comparisons. [Coordinate checksums](../../assets/structures/coordinate_sha256.json) · [Structural comparison table](../../assets/structures/RPE_structure_RMSD.tsv)
