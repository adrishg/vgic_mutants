# CaV1.2 structural clustering analysis

`cav12_opm_cluster.py` is an optional exploratory analysis. It uses the OPM
orientation of 8FD7 to define a membrane-frame scaffold, fits only the S1–S3
segments, and retains S4–S6/pore geometry as clustering features. Use one final
or earliest-converged model per trajectory; recycle snapshots are correlated
and are not independent cluster observations.

Inputs are an OPM-framed reference and a TSV manifest with `model_path`.
Generated features and cluster labels should be written to an external output
directory and are not part of the primary manuscript analysis.
