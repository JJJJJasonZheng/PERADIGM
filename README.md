# PERADIGM  
**Phenotype Embedding Similarity-based Rare Disease Gene Mapping**

This repository contains the R and Python code supporting the analysis described in the paper:  
**[PERADIGM: Phenotype Embedding Similarity-based Rare Disease Gene Mapping](https://www.biorxiv.org/content/10.1101/2025.04.01.646670v2)**

---

## Overview

PERADIGM is a framework that integrates phenotype embedding and patient similarity to identify rare disease-associated genes using large-scale biobank data.  
This repository includes code to replicate the key analyses and figures from the study, as well as scripts to generate custom ICD-10 code embeddings from medical text.

---

## Contents

- `main.R`: Main R script for the analysis pipeline, including:
  - Data loading and preprocessing  
  - Phenotype-gene association tests  
  - Generating similarity matrices and embeddings  
  - Outputting statistical results  

- `function.R`: R helper functions for:
  - Embedding computation  
  - Similarity scoring  
  - Regression-based testing  
  - Carrier/control selection  

- `train_icd10_embeddings.py`: **Python script to train ICD-10 code embeddings** using Word2Vec on medical descriptions.
  - Use this if you wish to regenerate or customize the ICD-10 code embeddings used in downstream analysis.

---

## 📁 Repository Structure

Your files should be organized as follows:

```text
data/
├── R_doc/
│   ├── hesin_diag_all_new.RData
│   ├── eid_all.RData
│   ├── cov_adjust.RData
│   └── IC_hesin_200k.csv
├── icd_related/
│   └── ICD10_mapping.csv
├── generate_all_gene_pos/
│   └── gene_info.RData
├── embedding/
│   └── hesin_icd10_descrip_embed.txt
├── icd_200k_codes.csv
├── hesin_200k_descrip_icd10.txt
└── hesin_diag.txt   # Optional/redundant diagnosis file

train_icd10_embeddings.py
main.R
function.R
requirements.txt
README.md
