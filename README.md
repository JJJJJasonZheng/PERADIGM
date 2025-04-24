### PERADIGM: Phenotype Embedding Similarity-based Rare Disease Gene Mapping

This repository contains the R code used for the analysis described in the paper:

PERADIGM: Phenotype Embedding Similarity-based Rare Disease Gene Mapping
https://www.biorxiv.org/content/10.1101/2025.04.01.646670v1

The code consists of two main files:

1. main.R: The primary script that includes data loading, pre-processing, analysis, and output generation.
2. function.R: A collection of functions supporting the analysis, including regression testing, embedding computations, similarity scoring, carrier processing, and permutation testing.

# =============================================================================
# SETUP: Load Data Files Using Relative Paths
# =============================================================================
# Adjust your repository structure as follows:
# ├── data
# │   ├── R_doc
# │   │   ├── hesin_diag_all_new.RData
# │   │   ├── eid_all.RData
# │   │   ├── cov_adjust.RData
# │   │   └── IC_hesin_500k.csv
# │   ├── icd_related
# │   │   └── ICD10_mapping.csv
# │   ├── generate_all_gene_pos
# │   │   └── gene_info.RData
# │   ├── embedding
# │   │   └── hesin_icd10_descrip_embed.txt
# │   ├── wes470K_fam.txt
# │   └── hesin_diag.txt          # redundant diagnosis file
