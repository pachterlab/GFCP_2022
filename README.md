[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6361041.svg)](https://doi.org/10.5281/zenodo.6361041)

# RNA velocity unraveled
This directory contains scripts and code to generate all the figures in the manuscript ["RNA velocity unraveled"](https://www.biorxiv.org/content/10.1101/2022.02.12.480214v1) by Gennady Gorin, Meichen Fang, Tara Chari, and Lior Pachter.

The `notebooks` directory contains the notebooks that generate the figures: 
* Figure 1: `vcy_scvelo_comparison.ipynb`
* Figure 5: `phaseplots_lme.ipynb`
* Figure 6: `embed_neighbors_jaccard_lme.ipynb`
* Figure 7: `embed_stability_lme.ipynb`
* Figure 8: `occup_meas_sim.ipynb`
* Figure 9: `aba_sim.ipynb`
* Figure 10: `embed_stability_sim.ipynb`
* Figure S1: `abcde_sim.ipynb`
* Figure S3: `paired_dataset_comparison.ipynb`
* Figure S4: `embed_neighbors_transf_lme.ipynb`
* Figure S5: `occup_meas_sim_nonorm.ipynb`
* Figures S6-7: `burst_false_positives.ipynb`

The `scripts` directory contains the bash, Python, and R scripts used to generate `kallisto|bustools`, `velocyto`, and `salmon` spliced and unspliced molecule count matrices in the `loom` format.

The `figures` directory contains up-to-date figures for the manuscript.

The raw data for the `paired_dataset_comparison` notebook are available at [the CaltechData repository](https://data.caltech.edu/records/20030). The human forebrain dataset generated was obtained from [this server](http://pklab.med.harvard.edu/velocyto/hgForebrainGlut/hgForebrainGlut.loom), as in the [`velocyto` tutorials](https://github.com/velocyto-team/velocyto-notebooks/blob/master/python/hgForebrainGlutamatergic.ipynb); for convenience we also host it in `notebooks/data`.

The file `vis.py` contains Python code that implements all of the analysis, including processing, simulation, and visualization.
