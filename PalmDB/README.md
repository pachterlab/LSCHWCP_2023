| Filename      | Description |
| ----------- | ----------- |
| palmdb_rdrp_seqs.fa       | Unique RdRP amino acid sequences from https://github.com/rcedgar/palmdb/blob/main/2021-03-14/uniques.fa.gz. IDs with sequences that were not unique after reverse translation to comma-free code (due to ambiguous amino acid annotation) were merged. |
| palmdb_clustered_t2g.txt       | Transcript to gene file. Virus "transcript" IDs with identical taxonomy annotation (from phylum to species) were assigned the same virus "gene" ID.  |
| ID_to_taxonomy_mapping.csv   | Maps virus "transcript" IDs to their representative virus "gene" ID and complete taxonomy. Taxonomies were originally released here: https://github.com/rcedgar/palmdb/blob/main/2021-03-14/u_tax.tsv.         |

The code used to modify the original PalmDB files can be found [here](https://github.com/pachterlab/LSCHWCP_2023/tree/main/Notebooks/create_optimized_palmdb).
