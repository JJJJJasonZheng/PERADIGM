# PERADIGM  
**Phenotype Embedding Similarity-based Rare Disease Gene Mapping**

This repository contains the R code supporting the analysis described in the paper:  
**[PERADIGM: Phenotype Embedding Similarity-based Rare Disease Gene Mapping](https://www.biorxiv.org/content/10.1101/2025.04.01.646670v1)**

## Overview

PERADIGM is a framework that integrates phenotype embedding and patient similarity to identify rare disease-associated genes using large-scale biobank data. This repository includes code to replicate the key analyses and figures from the study.

### Contents

- `main.R`: Main script for the analysis, including:
  - Data loading and preprocessing  
  - Running phenotype-gene association tests  
  - Generating similarity matrices and embeddings  
  - Outputting statistical results  

- `function.R`: Contains all helper functions for:
  - Embedding computation  
  - Similarity scoring  
  - Regression-based testing  
  - Carrier/control selection  

---

## 📁 Repository Structure

Place your data files using the following directory structure:


```text
data/
├── R_doc/
│   ├── hesin_diag_all_new.RData
│   ├── eid_all.RData
│   ├── cov_adjust.RData
│   └── IC_hesin_500k.csv
├── icd_related/
│   └── ICD10_mapping.csv
├── generate_all_gene_pos/
│   └── gene_info.RData
├── embedding/
│   └── hesin_icd10_descrip_embed.txt
└── hesin_diag.txt   # Optional/redundant diagnosis file
```

## 🔧 Getting Started

To reproduce the analysis:

1. Ensure R and required packages are installed.
2. Place the data files in the correct subfolders as shown above.
3. Run `main.R` to initiate the pipeline.

---
