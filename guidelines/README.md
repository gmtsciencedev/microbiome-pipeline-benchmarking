# Guidelines

This is a practical guide applying the principles described into the publication abstract.

## use or create realistic simulations tailored to your biological context (BC)

### Using the provided simulated data

The most simple option, which applies only for human gut microbiome context, is to use the provided simulated data. We propose a skeleton code here to launch a custom taxonomic profiler on these data: [scitq_custom_profiler_launcher.py](./scitq_custom_profiler_launcher.py)

### Creating your own realisatic simulated data

- find a study with real samples matching your use case,
- create a taxonomic profile with an adapted profiler (maybe using your own, but ideally a reference profiler, like MetaPhlAn 4)
- find a genome collection adapted to your BC (for instance looking at MGnify collections)
- project your profiles to this collection (use Kraken2, a database should be provided with the MGnify collection, you can use https://github.com/gmtsciencedev/scitq-examples/tree/main/kraken2 if you want to use scitq to do your Kraken projection), you will want to use the annotate_kraken2.py script provided in [feature_space_projections/code](../feature_space_projections/code/) folder.
- use CAMISIM to do the simulation (again using scitq, this can be done with this resource: https://github.com/gmtsciencedev/scitq-examples/tree/main/camisim)

## identify a common feature space suited to your BC and independent of the catalogues used by the profilers, 

We recommand using a curated feature space, like GTDB (avoid NCBI with is too open allowing for various substrains and candidate species), or even better a narrow feature space adapted to your BC (like a MGnify collection).

- project profiles to this feature space (use Kraken2, a database is provided either for GTDB or MGnify, you can use https://github.com/gmtsciencedev/scitq-examples/tree/main/kraken2 if you want to use scitq to do your Kraken projection), you will want to use the annotate_kraken2.py script provided in [feature_space_projections/code](../feature_space_projections/code/) folder.

## apply a comprehensive set of metrics covering accuracy (sensitivity/precision), overall representativity (richness/Shannon), and quantification (Unifrac and/or Aitchison distance).

Use the code provided in the [figure_code](../figure_code/) folder to apply the metrics proposed in our article.