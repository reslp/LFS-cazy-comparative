# Run a comparative genomic analysis of 83 ascomycete genomes including lichen-fungal symbionts

For the pipeline to run you will have to change a few things:

1. Adjust the `data/cluster-config.yaml` file to fit your own cluster enviroment (if you work on a HPC cluster).
3. Place input files in the correct folder. The pipeline needs several input files corresponding to the individual genomes such as GFF3 files, predicted proteins and transcripts, and annotation files produced by funannotate. Currently these files are not included in the repository. When the manuscript and the associated data gets published, download of this data will be automated and included here. 

## Run the pipeline:

Different steps of the pipeline can be run individually to generate results:
The commands here assume the pipeline is run an SGE cluster (in this case Sauron, the HPC cluster of the University of Graz). The commands below have to be executed in the base directory of the repository.

1. To calculate gene and site concordance factors and plot phylogenetic trees. The phylogeny itself was calculated using phylociraptor (https://github.com/reslp/phylociraptor):
	
```
./submit.sh -t sge -c data/cluster_config-sauron.yaml -s "-r phylogeny"
```

2. To calculate ancestral CAZyme family size and create plots:

```
./submit.sh -t sge -c data/cluster_config-sauron.yaml -s "-r reconstruct_ancestral_states"
```

3. To calculate PCA and phylogenetically corrected PCA for CAZyme sets, as well as % differences between different taxonomic groups:

```
./submit.sh -t sge -c data/cluster_config-sauron.yaml -s "-r geneset_similarity"
```

4. To perform subcellular location prediction, data mining of cazy.org and phylogenetic characterization of CAZyme families run: 

```
./submit.sh -t sge -c data/cluster_config-sauron.yaml -s "-r characterize_cazymes"
```

5. To calculate genome overview statistics run:

```
./submit.sh -t sge -c data/cluster_config-sauron.yaml -s "-r statistics"
```

6. To run sugar transporter orthology analysis:

```
./submit.sh -t sge -c data/cluster_config-sauron.yaml -s "-r characterize_transporters"
```

7. To run peroxidase orthology analysis:
```
./submit.sh -t sge -c data/cluster_config-sauron.yaml -s "-r characterize_peroxidases"
```

8. To run gene family expansion analysis with CAFE:

```
./submit.sh -t sge -c data/cluster_config-sauron.yaml -s "-r infer_gene_family_evolution"
```




 
