# Figure Data 

## Data folder

About the [data](./data) folder

Directory | Description
---|---
`refMet4_raw/`, `refKrak_raw/` | abundance tables in tools native spaces   
`refMet4_to_uhgg/`, `refKrak_to_uhgg/` | abundance data projected onto UHGG space
`refMet4_to_gtdb207/`, `refKrak_to_gtdb207/` | abundance data projected onto GTDB space
taxo_indexes | tables for converting between feature spaces, UHGG and GTDB taxonomies, phylogenetic trees

* `refMet4` and `refKrak` correspond to RefMet4 and RefKrak simulations
* `supersat_*.sat`, `lost_clades.tsv`, `lost_abundance.tsv` contain combined data for all pipelines (for a given simulation and feature space)
