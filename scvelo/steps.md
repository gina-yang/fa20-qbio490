# Running scvelo
### Steps:
1. Follow kallisto installation instructions [here](https://chmi-sops.github.io/mydoc_kallisto.html) to download kallisto. (For Windows, remember to add kallisto to PATH)
2. [Create loom file with loompy/kallisto.](https://linnarssonlab.org/loompy/kallisto/index.html)
- Getting indexes for non 10X chromium data seems mildly complicated. However, the tutorial implies that non 10X data is not supported with loompy yet anyway.
- `loompy fromfq` requires pairs of fastq read files (one filename should have "R1" and the other "R2"). Don't add the index ("I1") file. Additionally, make sure the R1 and R2 file for each sample come one after the other.
3. Use [this guide](https://scvelo.readthedocs.io/VelocityBasics.html) to get velocity analysis going.

### Notes:
- Make sure to use Python 3.8 (but not 3.9 because as of Oct 2020 it seems to be poorly supported). 3.7 and below can result in using old versions of packages with a bunch of fun errors and cross-package compatability issues.
- On that note, if stuff isn't working/weird errors, ensure that packages are updated (ex. with `pip install -U [packagename]` or `conda update [packagename]`) 
- loompy outputs total number of cells in the input files just before creating the loom file. So to include all cells it is probably best to specify a large number of cells in the metadata.tab file.

### Questions/Issues:
- How to remove cells from processing?
- How to change number of genes used in velocity analysis?
