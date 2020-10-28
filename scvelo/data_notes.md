# Notes about data

It looks like loompy tells you the number of valid cells just before it creates the .loom file. This number is also reflected in the object created from the .loom file by scvelo (n_obs).
- Maybe best to just specify a large number of cells (in the metadata file for loompy) and hope that it is more than in the dataset.

###KUL01
- .loom file - 3512 x 60662 (3512 cells, 60662 observations)
- processed in loompy with request of 6000 cells

###KUL19 (TUMOR CORE CELLS ONLY)
- .loom file - 4022 x 60662 (4022 cells)
- processed in loompy with request of 4022 cells
