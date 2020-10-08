# Running scvelo
### Steps:
1. Follow kallisto installation instructions [here](https://chmi-sops.github.io/mydoc_kallisto.html) to download kallisto. (For Windows, remember to add kallisto to corresponding PATH)
2. [Create loom file with loompy/kallisto.](https://linnarssonlab.org/loompy/kallisto/index.html)
3. Use [this guide](https://scvelo.readthedocs.io/VelocityBasics.html) to get velocity analysis going. (see script crc_scvelo_run.py)

### Notes:
- Make sure to use Python 3.8 (but not 3.9 because as of Oct 2020 it seems to be badly supported). Anything 3.7 and below can result in using old versions of packages with a bunch of fun bugs.
- If stuff isn't working/weird errors, ensure that packages are updated (ex. with `pip install -U [packagename]` or `conda update [packagename]`) 
