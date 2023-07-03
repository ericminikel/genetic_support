This repository holds the data and source code for the following manuscript:

[Minikel EV, Painter JL, Dong CL, Nelson MR. **Refining the impact of genetic evidence on clinical success**. _medRxiv_. 2023 Jun 29. 2023.06.23.23291765.](https://doi.org/10.1101/2023.06.23.23291765)

Here, you can:

+ Run the source code to reproduce the figures from the input datasets. Just say `Rscript `[`src/gensup_analysis.R`](/src/gensup_analysis.R), noting the dependencies at the top of the script. It completes in about 8 minutes on a 2021 MacBook Pro. The script reproduces Figures 1-3, S2-S5, Tables S1-S30, and stats_for_text.txt, all of which you can find in [display_items](/display_items). To run the script in "one target only mode", where drugs with >1 human target are removed, say `Rscript src/gensup_analysis.R --oto` and you'll find the output in [oto](/oto); Figures 1-3 from that version of the analysis become figures S6 - S8 in the manuscript.
+ If you're curious, you can also browse the source code for other scripts that prepared this releasable analytical dataset, in [src](/src). These scripts require inputs that are either too large for GitHub, and/or not approved for public release, thus, you will not be able to successfully run them after cloning the repository; they are provided simply for reference in case you want to see what we did.
+ Browse the input datasets in [data](/data). We have permission from Citeline Pharmaprojects to publicly release the subset of their data that appear here. This includes [data/pp.tsv](/data/pp.tsv), which contains the highest phase reached for all target-indication (T-I) pairs added to Pharmaprojects since 2000.

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.

