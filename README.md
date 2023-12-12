# Efficient and accurate detection of viral sequences at single-cell resolution reveals novel viruses perturbing host gene expression

Data, code, and figures generated in the manuscript
```
Laura Luebbert, Delaney K Sullivan, Maria Carilli, Kristjan Eldjarn Hjorleifsson, Alexander Viloria Winnett, Tara Chari, Lior Pachter (2023). Efficient and accurate detection of viral sequences at single-cell resolution reveals novel viruses perturbing host gene expression. bioRxiv 2023.12.11.571168; doi: https://doi.org/10.1101/2023.12.11.571168
```
Read the article here: [https://www.biorxiv.org/content/10.1101/2023.12.11.571168](https://www.biorxiv.org/content/10.1101/2023.12.11.571168)

The [Notebooks](https://github.com/pachterlab/LSCHWCP_2023/tree/main/Notebooks) folder contains all analyses that were performed within the scope of this manuscript, from the raw data to the final figure, in immediately executable Google Colab notebooks. 

Large datasets are stored on [Caltech Data](https://data.caltech.edu/records/sh33z-hrx98?token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6IjlhNDNkZWVkLTRiODYtNDIwMS1hNTcwLTYyNDZhOGYwZjU3YyIsImRhdGEiOnt9LCJyYW5kb20iOiI3YTU1MDY5MjEzY2Y0ZmMyNjVlODMyYTZlOWQ4MTUxMCJ9.RkUlR18JUioegjOX_7m89ngFcatseZGRLZaadwc8X0GgzCxztvnkNc6rUMT8ozAta2LEcpwhdOq33QOH9Slj7g).

[Click here](https://htmlpreview.github.io/?https://github.com/pachterlab/LSCHWCP_2023/blob/main/krona_plot.html) to view the interactive Krona plot showing all viruses expressed above the QC threshold in macaque cells that passed quality control broken down by animal, timepoint, taxonomy, and fraction of positive cells occupied by each virus.

The [precomputed_refs](https://github.com/pachterlab/LSCHWCP_2023/tree/main/precomputed_refs) folder contains precomputed reference indices for the detection of viral RNA in sequencing data (through alignment to the [optimized PalmDB](https://github.com/pachterlab/LSCHWCP_2023/tree/main/PalmDB)) with masked human (or mouse) genome and transcriptome.
