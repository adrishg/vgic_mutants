# Stand-ins for what `screen.py` imports but the checks never use

`pipeline/screen.py` imports three names-only modules at module scope, and one of them pulls in
`gemmi`, which the ColabFold environment on Hive does not have. The geometry itself needs none of
it: the checks take plain coordinate arrays.

So `filter_chain._screen_module()` imports his file unchanged, and only if that import fails for a
missing dependency does it put this directory on the path and try again. His code is still the code
under test; what is replaced is the config and archive plumbing around it, which our loose PDB trees
do not go through anyway.
