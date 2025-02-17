# Efficient and accurate detection of viral sequences at single-cell resolution reveals novel viruses perturbing host gene expression

This repository contains data, code, and figures generated for the manuscript:
```
Laura Luebbert, Delaney K Sullivan, Maria Carilli, Kristján Eldjárn Hjörleifsson, Alexander Viloria Winnett, Tara Chari, Lior Pachter (2023). [Efficient and accurate detection of viral sequences at single-cell resolution reveals novel viruses perturbing host gene expression](https://www.biorxiv.org/content/10.1101/2023.12.11.571168). bioRxiv 2023.12.11.571168; doi: https://doi.org/10.1101/2023.12.11.571168
```
The preprint is posted on the bioRxiv: [https://www.biorxiv.org/content/10.1101/2023.12.11.571168](https://www.biorxiv.org/content/10.1101/2023.12.11.571168)

> 💡 **General tutorials with example data can be found on the [kallisto bustools](https://kallisto.readthedocs.io/en/latest/) website:**  
> - [Detecting viral sequences in bulk RNA sequencing data](https://kallisto.readthedocs.io/en/latest/translated/notebooks/virus_detection_bulk.html)  
> - [Detecting viral sequences in single-cell RNA sequencing data](https://kallisto.readthedocs.io/en/latest/translated/notebooks/virus_detection_sc.html)

When interpreting the presence of RdRP-like sequences / virus IDs, keep in mind that there will likely be many RdRP-like sequences introduced by contamination of laboratory reagents. A (non-comprehensive) list of **virus IDs observed in blank sequencing data** is available [here](https://github.com/pachterlab/LSCHWCP_2023/blob/main/viruses_in_blank_reagents/total_raw_count_per_virus_id_in_laboratory_reagents.csv).  

The [Notebooks](https://github.com/pachterlab/LSCHWCP_2023/tree/main/Notebooks) folder contains notebooks to reproduce all of our analyses, starting with pre-processing of the raw data all the way to final figure generation. The notebooks are organized by figure (based on the bioRxiv preprint) and immediately executable via Google Colab.  

Since the figure order was updated between the bioRxiv preprint and the subsequent publication of the manuscript in _Nature Biotechnology_, the [Notebooks_Nature_Biotech](https://github.com/pachterlab/LSCHWCP_2023/tree/main/Notebooks_Nature_Biotech) folder links to the appropriate notebooks based on the figure numbering in the _Nature Biotechnology_ version.  

Large intermediary files that are generated/used in these notebooks are stored on Caltech Data and can be accessed under the DOIs [10.22002/krqmp-5hy81](https://data.caltech.edu/records/krqmp-5hy81) and [10.22002/k7xqw-88d74](https://data.caltech.edu/records/k7xqw-88d74).

[Click here](https://htmlpreview.github.io/?https://github.com/pachterlab/LSCHWCP_2023/blob/main/krona_plot.html) to view the interactive Krona plot showing all viruses expressed above the QC threshold in macaque cells that passed quality control, broken down by animal, timepoint, taxonomy, and fraction of positive cells occupied by each virus. [Code to reproduce the Krona plot](https://github.com/pachterlab/LSCHWCP_2023/tree/main/Notebooks/krona_plot)

The [precomputed_refs](https://github.com/pachterlab/LSCHWCP_2023/tree/main/precomputed_refs) folder contains precomputed reference indices for the detection of viral RNA in sequencing data (through alignment to the [optimized PalmDB](https://github.com/pachterlab/LSCHWCP_2023/tree/main/PalmDB)) and with masked human (or mouse) genome **and** transcriptome.

A description of kallisto, bustools, and kb-python including tutorials for their use can be found here: [https://www.nature.com/articles/s41596-024-01057-0](https://www.nature.com/articles/s41596-024-01057-0)

<br>
</br>

```bash
# 1. Install kb-python (optional: install gget to fetch the host genome and transcriptome)
pip install kb-python gget

# 2. Download optimized PalmDB reference files
!wget https://raw.githubusercontent.com/pachterlab/LSCHWCP_2023/main/PalmDB/palmdb_rdrp_seqs.fa
!wget https://raw.githubusercontent.com/pachterlab/LSCHWCP_2023/main/PalmDB/palmdb_clustered_t2g.txt

# 3. Create reference index (+ optional masking of the host, here human, genome using the D-list)
# Single-thread runtime: 1.5 h; Max RAM: 4.4 GB; Size of generated index: 593 MB
# Without D-list: Single-thread runtime: 3.5 min; Max RAM: 3.9 GB; Size of generated index: 592 MB
kb ref \
    --aa \
    -k 55 \
    --d-list $(gget ref --ftp -w dna,cdna homo_sapiens) \
    -i index.idx --workflow custom \
    palmdb_rdrp_seqs.fa
    
# 4. Align sequencing reads
# Single-thread runtime: 1.5 min / 1 million sequences; Max RAM: 2.1 GB
kb count \
    --aa \
    -k 55 \
    -i index.idx -g palmdb_clustered_t2g.txt \
    --parity single \
    -x default \
    $USER_DATA.fastq.gz
```
  
![Overview_v3_noCode](https://github.com/pachterlab/LSCHWCP_2023/assets/56094636/e5cc1c24-3ce3-47cc-893b-93efc5a7329f)



