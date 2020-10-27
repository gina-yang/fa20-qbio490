# Running scvelo
### Steps:
1. Follow kallisto installation instructions [here](https://chmi-sops.github.io/mydoc_kallisto.html) to download kallisto. (For Windows, remember to add kallisto to corresponding PATH)
2. [Create loom file with loompy/kallisto.](https://linnarssonlab.org/loompy/kallisto/index.html)
- Getting indexes for non 10X chromium data seems mildly complicated. However, the tutorial implies that non 10X data is not supported with loompy yet anyway.
- Also currently unsure about the *metadata.tab* file's contents. It says that sample IDs go under the "name" column, so running multiple samples means multiple entries in the table. Maybe with a lot of samples (or a lot of metadata columns) it might be easier to use sqlite3
- `loompy fromfq` requires pairs of fastq read files (one filename should have "R1" and the other "R2"). Don't add the index ("I1") file.
3. Use [this guide](https://scvelo.readthedocs.io/VelocityBasics.html) to get velocity analysis going. (see script *crc_scvelo_run.py*)

### Notes:
- Make sure to use Python 3.8 (but not 3.9 because as of Oct 2020 it seems to be badly supported). Anything 3.7 and below can result in using old versions of packages with a bunch of fun bugs.
- If stuff isn't working/weird errors, ensure that packages are updated (ex. with `pip install -U [packagename]` or `conda update [packagename]`) 
- loompy outputs total number of cells in the input files just before creating the loom file. So to include all cells it is probably best to specify a large number of cells in the metadata.tab file.

### Questions/Issues:
- How to remove cells from processing?
- How to change number of genes used in velocity analysis?
