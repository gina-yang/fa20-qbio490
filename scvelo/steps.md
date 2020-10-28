# Running scvelo
### Steps:
1. Follow kallisto installation instructions [here](https://chmi-sops.github.io/mydoc_kallisto.html) to download kallisto. (For Windows, remember to add kallisto to PATH)
2. [Follow the steps here to create a loom file with loompy/kallisto.](https://linnarssonlab.org/loompy/kallisto/index.html)
    * I used the provided pre-built index of the human genome. Getting indexes for non 10X chromium data seems mildly complicated. However, the tutorial implies that non 10X data is not supported with loompy yet anyway.
    * Metadata file: has one row per sample. If the target/expected number of cells is not specified in the metadata file, it will default to 5000.
    * **Important:** the `loompy fromfq` command requires pairs of fastq read files, i.e. one filename should have `R1` (barcodes) and the other `R2` (cDNA sequences). Do not add the  file with`I1` (index/lane info). Additionally, the order the files are presented in the command matters:corresponding R1 and R2 files for a given lane should be listed together. More info about 10x fastq file naming conventions [here.](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/using/fastq-input)
3. Use [this guide](https://scvelo.readthedocs.io/VelocityBasics.html) to get a basic workflow for velocity analysis going.
    * The tutorial data has pre-computed UMAP embedding and colored clusters. To perform clustering and embedding, use `tl.louvain` and `tl.umap`. 
4. More info about functions in the [API](https://scvelo.readthedocs.io/api.html).
5. Summary of workflow:
    * `scv.read` reads a dataset (.loom file obtained from kallisto/loompy)
    * `scv.pp.*` functions for preprocessing
    * `scv.tl.*` analysis tools
    * `scv.pl.*` plotting tools
    * Also `scv.utils*` has useful utilities for data cleaning and merging

### Notes:
- Make sure to use Python 3.8 (but not 3.9 because as of Oct 2020 it seems to be poorly supported). 3.7 and below can result in using old versions of packages with a bunch of fun errors and cross-package compatability issues.
- On that note, if stuff isn't working/weird errors, ensure that packages are updated (ex. with `pip install -U [packagename]` or `conda update [packagename]`) 
- loompy outputs total number of cells in the input files just before creating the loom file. So, to include all cells, it is probably best to specify a large number of cells in the metadata.tab file. (unless you know the number of input cells, in which case put that)

### Questions/Issues:
- How to remove cells from processing?
- How to change number of genes used in velocity analysis?
- What defines a sample in the metadata file? It seems that loompy assumes each R1+R2 file pair represents a sample.
