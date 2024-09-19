# Codes for publication on pipeline benchmarking

This github repository is the central point for all codes and data for GMT Science publication on microbiome pipeline benchmarking.

Item|description
--|--
[scitq](https://github.com/gmtsciencedev/scitq)|A task distribution system used to generate and analyse the samples
[scitq-examples](https://github.com/gmtsciencedev/scitq-examples)|The code used to generate and analyse the samples
[containers](https://hub.docker.com/u/gmtscience)|The containers used in the code mentionned above (except for [CAMISIM container](https://hub.docker.com/r/cami/camisim)). They can all be used with or without scitq.
[container source](https://github.com/gmtsciencedev/bioit-dockers)|The source code of above mentionned containers
[PRJNA987980](https://www.ncbi.nlm.nih.gov/bioproject/987980)|The simulated sample datasets
[PRJEB6070](https://www.ebi.ac.uk/ena/browser/view/PRJEB6070) [PRJEB7774](https://www.ebi.ac.uk/ena/browser/view/PRJEB7774)|The original samples from which the simulation was derived|
[UHGG](http://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0.1/)|The MGNIFY genomes upon which the simulation was based on|

## Locally in this repository

- the figure code available in the [figure_code](./figure_code) folder
- the figure data available in the [figure_data](./figure_data) folder
- the reference, the true species composition of sample, the one that was used for simulation is available in the [reference](./reference/) folder
- the different feature space (a.k.a. reference catalogs) projections in [feature_space_projections](./feature_space_projections/)



## Processing steps

This work consists in a sample simulation based upon measurements of real samples, which were then estimated in their species abundance. The estimations were then compared to the simulation source. All the required tools are described above but we wanted to clarify the different steps and when to use which tool to reproduce this work or parts of this work.

Most of these steps require a significant amount of computational power. It may be possible to use a very powerful server, but the method proposed here is to split the work in small tasks and distribute them on a cluster, either of physical servers or on cloud instances, using an orchestration solution. The orchestration solution proposed is an open source in house solution, [scitq](https://github.com/gmtsciencedev/scitq), that allow life cycle management of cloud instances (automatic creation and destruction) in different cloud environments (OVH (or any Openstack compatible cluster) and Azure), but that can also be used on a physical cluster.

### Initial measurements (reference)

The initial measurements are obtained through 2 different tools: [MetaPhlAn](https://github.com/biobakery/MetaPhlAn) and [Kraken](https://github.com/DerrickWood/kraken2), and some advices are proposed here to ease their use.

The initial data was fetched from EBI, using [PRJEB6070](https://www.ebi.ac.uk/ena/browser/view/PRJEB6070) and [PRJEB7774](https://www.ebi.ac.uk/ena/browser/view/PRJEB7774).

A sample code to download and filter the data from a bioproject accession is explained in detail at the bottom of [scitq README](https://github.com/gmtsciencedev/scitq/blob/main/README.md).

Some python scripts to distribute the tasks on scitq are proposed for [MetaPhlAn](https://github.com/gmtsciencedev/scitq-examples/tree/main/metaphlan4) and [Kraken](https://github.com/gmtsciencedev/scitq-examples/tree/main/kraken2). 

For Kraken, non default options were used: the [UHGG MGNIFY database](http://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0.1/kraken2_db_uhgg_v2.0.1/) was used, the files were deposited in a folder, `krakenmgnify`. A tar archive was then created with `tar czf krakenmgnify2.0.1.tgz krakenmgnify/`, this archived was then uploaded to an S3 bucket, `s3://rnd/resource/krakenmgnify2.0.1.tgz`, and the script was used with the following options: 

```bash
python ~/scitq-examples/kraken2/scitq_kraken2.py --bracken --database krakenmgnify --fastq --download --batch refKrak s3://rnd/data/raw/refKrak s3://rnd/resource/krakenmgnify2.0.1.tgz s3://rnd/refKrak/rawresult
```

For MetaPhlAn, the script was used with default options (and resource specified with the script), but results needed projections on the MGNIFY feature space. The projection used is available here: [feature_space_projections/metaphlan4_to_mgnify2.0.1.tsv](./feature_space_projections/metaphlan4_to_mgnify2.0.1.tsv).

The method used to obtain this projection is:

- extraction of genomes from MetaPhlAn database using the script [feature_space_projections/code/metaphlan_genome_extractor.py](./feature_space_projections/code/metaphlan_genome_extractor.py), 
- projection of genomes on MGNIFY using Kraken (same as above but without `--bracken` option).


All results once projected to MGNIFY are proposed in [reference](./reference/).

### Simulation

The simulation makes extensive use of [CAMISIM](https://github.com/CAMI-challenge/CAMISIM). This required the species composition of the samples to generate ([reference](./reference/)), the genomes specified in these files downloaded from [MGNIFY database](http://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0.1/). Next, the script specified in [scitq-example/camisim](https://github.com/gmtsciencedev/scitq-examples/tree/main/camisim) must be used, as explained in the link. 

Note that the simulation requires extensive resources and time, which is why samples are also directly proposed in [PRJNA987980](https://www.ncbi.nlm.nih.gov/bioproject/987980).

### Estimation

Estimation makes use of [MetaPhlAn](https://github.com/biobakery/MetaPhlAn) version 3 and 4, [Kraken](https://github.com/DerrickWood/kraken2) version 2, [mOTUs](https://github.com/motu-tool/mOTUs) version 3. 

These tools are used in scitq using the following scripts:

- [scitq MetaPhlAn](https://github.com/gmtsciencedev/scitq-examples/tree/main/metaphlan4) (a specific options exist to specify the MetaPhlAn version to use, `--metaphlan-version` which should be 3.1.0 or 4.0.6),
- [scitq Kraken](https://github.com/gmtsciencedev/scitq-examples/tree/main/kraken2),
- [scitq mOTUs](https://github.com/gmtsciencedev/scitq-examples/tree/main/motus3)

Details are given in each link. The obtained (raw) results are given in [figure_data/data](./figure_data/data/) `...raw` folders.

### Projection

The raw results obtained above were then projected either to UHGG or to GTDB. The projection tables are available in [feature_space_projections](./feature_space_projections/). And the transposed results are proposed in [figure_data/data](./figure_data/data/).

The projection tables were obtained :

- in the case of GTDB, using Kraken on reference genomes extracted from tool catalogs in the same manner as explained for MetaPhlAn in [initial mesurements](#initial-measurements-reference), except that GTDB database was used instead of UHGG database, 
- in the case of UHGG, reads were extracted from reference genomes, and submitted to Biomscope, and then when it was unsufficient, reference genomes were projected to UHGG (as explained for MetaPhlAn in [initial mesurements](#initial-measurements-reference)), and reads from the corresponding UHGG genomes were submitted to Biomscope.

### Comparison

The comparison and figure code, and all intermediate figure data are explained in detail in [figure_code](./figure_code/) [README](./figure_code/README.md).

## How to cite

To be completed.