# Figure Code 

Script | Description
---|---
`project_to_gtdb_space.py` | projects raw abundance tables to GTDB space
`project_to_uhgg_space.py`  | projects raw abundance tables to UHGG space
`calculate_metrics.py` | calculate basic metrics (distances, true/false positives) for a simulation
`unifrac_gtdb_multiproc.py` | calculate UniFrac distances for the data projected onto GTDB space
`unifrac_uhgg_multiproc.py` | calculate UniFrac distances for the data projected onto UHGG space
`aitchison_norms.py` |calculate Aitchison distances for the data
`fp_fn_analysis.py` | calculate correlations between false positive and false negative species
`raw_richness.py` | calculate simulations characteristics
`combine_metrics.py` | combing metrics data obtained for different simulations
`final_figures.py` | plot final figures included in the article
`radar_plots.py` | plot radar plots
`unifrac_aux.py`, `metrics_aux.py` | scripts with auxhiliary functions used by the above scripts

* Simulation for which projections/metrics are calculated is defined by parameter `simu`, taking values `refMet4, refKrak` (described in section Data)
* parameter `space` (or `sp`) in script `calculate_metrics.py` takes values `uhgg` or `gtdb207`, depending on whether the metrics are calculated for data projected onto UHGG or GTDB space

## Pre-calculated metrics 
Folder: [figure_data/analyses](../figure_data/analyses/)

Directory | Description
---|---
`article_results_refMet4_uhgg/`, `article_results_refKrak_uhgg/`, `article_results_refMet4_gtdb207/`, `article_results_refKrak_gtdb207/`| pre-calculated metrics for simulations (as specified in the directory name)
`combined_metrics/`, `confusion_pairs/`, `raw_richness/`| assembledtogether metrics for all simulations for creating final figures

## Figures
Folder: [figure_data/figures](../figure_data/figures/)

Final figures included in the article

# Workflow for final figures

This details how and in what order the different scripts should be run.

## Projections
(need to be run for all possible values of parameter `simu`.)
```
project_to_uhgg_space.py  
project_to_gtdb_space.py
```

## Calculating metrics
(need to be run for all possible values of parameters `simu` and `space`.)
```
calculate_metrics.py
```

(need to be run for all possible values of parameter `simu`.)
```
unifrac_uhgg_multiproc.py  
unifrac_gtdb_multiproc.py
aitchison_norms.py
```

```
combine_metrics.py
fp_fn_analysis.py
raw_richness.py
```


## Creating final figures
```
final_figures.py
radar_plots.py
```
