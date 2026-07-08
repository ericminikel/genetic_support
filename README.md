# Impact of genetic evidence on clinical success

This repository holds the data and source code for the following manuscript:

[Minikel EV, Painter JL, Dong CC, Nelson MR. **Refining the impact of genetic evidence on clinical success.** _Nature_. 2024 May;629(8012):624-629. doi: 10.1038/s41586-024-07316-0. Epub 2024 Apr 17. PMID: 38632401; PMCID: PMC11096124.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11096124/)

## Setup

1. Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or Anaconda.

Now follow one of the following options.

### Step 1: install Conda environment

Here you can use two options below.

#### Option 1: Conda

Create the environment for this project:

   ```bash
   conda env create -f environment.yml
   conda activate gensup
   ```

#### Option 2: `conda-lock`

Only available for `linux-64` and `osx-64`.

1. Install [`pipx`](https://github.com/pypa/pipx).

1. Install `conda-lock` using `pipx`:

   ```bash
   pipx install conda-lock
   ```

1. Create the environment for this project:

   ```bash
   # create conda environment
   conda-lock install --name gensup conda-lock.yml

   # activate environment
   conda activate gensup
   ```

### Step 2: install R packages

Some R packages are not available in conda. Open `R` and install them:

```R
require(remotes)
install_version("lawstat", version = "3.4", upgrade="never", repos = "http://cran.us.r-project.org")
```

## Run the analyses
Here, you can:

+ Run the source code to reproduce the figures from the input datasets. Just say `Rscript `[`src/gensup_analysis.R`](/src/gensup_analysis.R), noting the dependencies at the top of the script. It completes in about 8 minutes on a 2021 MacBook Pro. The script reproduces Figures 1-3, S2-S5, Tables S1-S30, and stats_for_text.txt, all of which you can find in [display_items](/display_items). To run the script in "one target only mode", where drugs with >1 human target are removed, say `Rscript src/gensup_analysis.R --oto` and you'll find the output in [oto](/oto); Figures 1-3 from that version of the analysis become figures S6 - S8 in the manuscript.
+ If you're curious, you can also browse the source code for other scripts that prepared this releasable analytical dataset, in [src](/src). These scripts require inputs that are either too large for GitHub, and/or not approved for public release, thus, you will not be able to successfully run them after cloning the repository; they are provided simply for reference in case you want to see what we did.
+ Browse the input datasets in [data](/data). We have permission from Citeline Pharmaprojects to publicly release the subset of their data that appear here. This includes [data/pp.tsv](/data/pp.tsv), which contains the highest phase reached for all target-indication (T-I) pairs added to Pharmaprojects since 2000.

_Note about dependencies_. This code was written for R 4.2.0 and the following package versions: tidyverse_1.3.1, janitor_2.1.0, binom_1.1-1.1, glue_1.6.2, lawstat_3.4, weights_1.0.4, epitools_0.5-10.1, DescTools_0.99.45, openxlsx_4.2.5, optparse_1.7.1, MASS_7.3-56.

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.

