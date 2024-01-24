| Filename      | Description |
| ----------- | ----------- |
| palmdb_rdrp_seqs.fa       | Unique RdRP amino acid sequences from https://github.com/rcedgar/palmdb/blob/main/2021-03-14/uniques.fa.gz. Sequences that were not unique after reverse translation to comma-free code (due to ambiguous amino acid annotation) were 'merged' and given a representative ID. |
| palmdb_clustered_t2g.txt       | Transcript to gene file. Virus "transcript" IDs with identical taxonomy annotation (from phylum to species) were assigned to a single representative virus "gene" ID.  |
| ID_to_taxonomy_mapping.csv   | Maps virus IDs to their species operational taxonomic units (sOTUs) (taxonomies). This mapping was originally released here: https://github.com/rcedgar/palmdb/blob/main/2021-03-14/u_tax.tsv. This file also contains the mapping of all virus IDs to their representative virus IDs.         |

The code used to modify the original PalmDB files can be found [here](https://github.com/pachterlab/LSCHWCP_2023/tree/main/Notebooks/create_optimized_palmdb).
