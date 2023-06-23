# Code and data

## Data 
Folder: [figure_data](../figure_data/)

Directory | Description
---|---
`refMet4_raw/`, `refBioms_raw/`, `refKrak_raw/` | abundance tables in tools native spaces   
`refMet4_to_igc2/`, `refBioms_to_igc2/`, `refKrak_to_igc2/` | abundance data projected onto IGC2/MSP space
`refMet4_to_gtdb207/`, `refBioms_to_gtdb207/`, `refKrak_to_gtdb207/` | abundance data projected onto GTSB space
taxo_indexes | tables for converting between feature spaces, IGC2/MSP and GTDB taxonomies, phylogenetic trees

* `refMet4`, `refBioms` and `refKrak` correspond to RefMet4, RefKrak, and RefBiom simulations
* `supersat_*.sat`, `lost_clades.tsv`, `lost_abundance.tsv` contain combined data for all pipelines (for a given simulation and feature space)

## Code 
This folder:

Script | Description
---|---
`project_to_gtdb_space.py` | projects raw abundance tables to GTDB space
`project_to_igc2_space.py`  | projects raw abundance tables to IGC2/MSP space
`calculate_metrics.py` | calculate basic metrics (distances, true/false positives) for a simulation
`unifrac_gtdb.py` | calculate UniFrac distances for the data projected onto GTDB space
`unifrac_igc2.py` | calculate UniFrac distances for the data projected onto IGC2/MSP space
`combine_metrics.py` | combing metrics data obtained for different simulations
`final_figures.py` | plot final figures included in the article
`unifrac_aux.py`, `metrics_aux.py` | scripts with auxhiliary functions used by the above scripts

* Simulation for which projections/metrics are calculated is defined by parameter `simu`, taking values `refMet4, refBioms, refKrak` (described in section Data)
* parameter `space` (or `sp`) in script `calculate_metrics.py` takes values `igc2` or `grdb207`, depending on whethe rteh matrics are calculate for data projected onto IGC2/MSP or GTDB space

## Pre-calculated metrics 
Folder: [figure_data/analyses](../figure_data/analyses/)

Directory | Description
---|---
`article_results_refMet4_igc2/`, `article_results_refBioms_igc2/`, `article_results_refKrak_igc2/`, `article_results_refMet4_gtdb207/`, `article_results_refBioms_gtdb207/`, `article_results_refKrak_gtdb207/`| pre-calculated metrics for simulations (as specified in the directory name)
`combined_metrics/ `| assembledtogether metrics for all simulations for creating final figures

## Figures
Folder: [figure_data/figures](../figure_data/figures/)

Final figures included in the article

