# Retrieve command-line arguments
args <- commandArgs(trailingOnly = TRUE)
phenotype <- args[1]
compare_icd <- args[2]
# Example usage:
# phenotype = "Q612;Q613"
# compare_icd = "Q612;Q613;Q446;I710;I711;I712;I720;I721;I722;I723;I724;I725;I254"
# Split the ICD codes by semicolon into a vector
compare_icd <- strsplit(compare_icd, ";")[[1]]

# Load required libraries
library(data.table)
library(dplyr)
library(lsa)
library(scales)

# Option to disable scientific notation
options(scipen = 999)

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
#
# Load main RData files
load("data/R_doc/hesin_diag_all_new.RData")  # loads hesin_diag_all
load("data/R_doc/eid_all.RData")             # loads eid_all
load("data/R_doc/cov_adjust.RData")          # loads cov_adjust

# Filter covariate adjustment data for European individuals
cov_adjust <- cov_adjust[cov_adjust$eur == 1, ]

# Load the ICD10 mapping file
ICD10_mapping <- fread("data/icd_related/ICD10_mapping.csv")

# Load gene information
load("data/generate_all_gene_pos/gene_info.RData")  # loads gene_info

# Load the family data and get unique IDs
eid_500k <- fread("data/wes470K_fam.txt")
eid_500k <- unique(eid_500k$eid)

# Subset diagnosis data from the main dataset using eid_all list
hesin_diag_200k <- hesin_diag_all[hesin_diag_all$eid %in% eid_all, ]

# =============================================================================
# SOURCE HELPER FUNCTIONS
# =============================================================================
# Make sure the script is saved in a "script" folder within the repository.
source("script/simi_embed_function.R")

# =============================================================================
# Load and prepare additional datasets
# =============================================================================
# Load the IC risk file
IC_500k <- read.csv("data/R_doc/IC_hesin_500k.csv", row.names = 1)
colnames(IC_500k) <- c("icd10", "Freq", "Prob", "IC")
IC_500k_rare <- IC_500k[IC_500k$Prob <= 0.01, ]

# Load the redundant diagnosis file and subset its relevant columns.
hesin_diag_redundant <- fread("data/hesin_diag.txt")
hesin_diag_200k_redundant <- hesin_diag_redundant[, c(1, 7)]
hesin_diag_200k_redundant <- hesin_diag_200k_redundant[hesin_diag_200k_redundant$eid %in% eid_all, ]
hesin_diag_200k_redundant <- hesin_diag_200k_redundant[hesin_diag_200k_redundant$diag_icd10 != "", ]

# Load disease description embedding data
embed_icd_descrip <- fread("data/embedding/hesin_icd10_descrip_embed.txt")
embed_icd_descrip <- as.data.frame(embed_icd_descrip)
rownames(embed_icd_descrip) <- embed_icd_descrip$coding
embed_icd_descrip <- embed_icd_descrip[, -1]
embed_use <- embed_icd_descrip

# =============================================================================
# DATA PRE-PROCESSING
# =============================================================================
# Keep only records with rare IC codes in the diagnosis file
hesin_use <- hesin_diag_200k_redundant[hesin_diag_200k_redundant$diag_icd10 %in% IC_500k_rare$icd10, ]
eid_use <- intersect(eid_all, cov_adjust$eid)
flag_use <- "200k"
hesin_use <- hesin_use[hesin_use$eid %in% cov_adjust$eid, ]
IC_500k_rare <- IC_500k_rare[IC_500k_rare$icd10 %in% rownames(embed_use), ]
eid_use <- intersect(eid_use, unique(hesin_use$eid))

# =============================================================================
# DISEASE PATIENT EMBEDDING AND WEIGHTING
# =============================================================================
# Process the phenotype(s)
phenotype_list <- strsplit(phenotype, ";")[[1]]
disease_id <- unique(hesin_use[hesin_use$diag_icd10 %in% phenotype_list, ]$eid)
disease_non_id <- setdiff(eid_use, disease_id)

# Calculate significance between disease and non-disease groups
# (regres_sig_phe must be defined in your sourced script or another file)
compare_disease <- regres_sig_phe(disease_id, disease_non_id, hesin_use, cov_adjust, ICD10_mapping)
# Set adjusted p-values for the phenotype (i.e. target disease) to zero
compare_disease[compare_disease$icd10 %in% phenotype_list, ]$pvalue_adjusted <- 0
compare_disease$log10Pvalue <- -log10(compare_disease$pvalue_adjusted)

# Apply sigmoid weighting (assign_weights_sigmoid is assumed to be defined)
compare_disease$weight <- sapply(compare_disease$log10Pvalue, assign_weights_sigmoid)
compare_disease_merge <- merge(compare_disease, IC_500k_rare, by = "icd10", all = TRUE)
compare_disease_merge <- compare_disease_merge[, c("icd10", "weight", "IC")]
compare_disease_merge[is.na(compare_disease_merge)] <- 0.2
compare_disease_merge$weight_combine <- compare_disease_merge$weight * compare_disease_merge$IC
disease_weight_sigPhe <- compare_disease_merge[, c("icd10", "weight_combine")]
colnames(disease_weight_sigPhe) <- c("icd10", "weight")
rownames(disease_weight_sigPhe) <- disease_weight_sigPhe$icd10

# =============================================================================
# CALCULATE AVERAGE EMBEDDING FOR DISEASE PATIENTS (DESCRIPTION EMBEDDING)
# =============================================================================
# Create a matrix to store embeddings (extra column for ID will be removed later)
embedding_disease <- data.frame(matrix(nrow = length(disease_id), ncol = ncol(embed_use) + 1))
embedding_disease[, 1] <- disease_id
rownames(embedding_disease) <- embedding_disease[, 1]
embedding_disease <- embedding_disease[, -1]

# Subset the diagnosis records for disease patients
disease_hesin <- hesin_use[hesin_use$eid %in% disease_id, ]
# Keep only entries with ICD codes that are in the compare_icd list
disease_hesin <- disease_hesin[disease_hesin$diag_icd10 %in% compare_icd, ]

# Loop over disease patient IDs and compute average embedding with weights
for (i in 1:nrow(embedding_disease)) {
  eid_temp <- rownames(embedding_disease)[i]
  embed_temp <- get_avg_embed_weight(eid_temp, embed_use, disease_hesin, weight_table = disease_weight_sigPhe)
  embedding_disease[i, ] <- embed_temp
  if (i %% 50 == 0) { cat(i, " ") }
}
embedding_disease <- na.omit(embedding_disease)
print("embed disease patients done!")

# =============================================================================
# COMPUTE SIMILARITY SCORES FOR THE WHOLE COHORT
# =============================================================================
disease_result <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(disease_result) <- c("eid", "simi_score")
for (i in 1:length(eid_use)) {
  tmp_id <- eid_use[i]
  tmp <- compare_group_disease(id = tmp_id,
                               hesin_diag = hesin_use,
                               eid_whole = eid_use,
                               embed_icd = embed_use,
                               embedding_disease = embedding_disease,
                               weight_sigPhe = disease_weight_sigPhe)
  disease_result <- rbind(disease_result, tmp)
  if (i %% 50 == 0) { cat(i, " ") }
}

# Save the computed risk score results in a relative results directory.
output_file <- paste0("data/riskscore_summary/riskscore/", phenotype,
                      "_riskscore_eur_weighted_IC_rare_all_diseaseCompare.RData")
save(disease_result, file = output_file)
print("riskscore_eur finished!")

# =============================================================================
# GENE RISK SCORE CALCULATION
# =============================================================================
# Reload the previously saved disease_result for further processing
filename <- paste0("data/riskscore_summary/riskscore/", phenotype,
                   "_riskscore_eur_weighted_IC_rare_all_diseaseCompare.RData")
load(filename)

# Identify individuals present in both the redundant diagnosis and the selected cohort
eid_hesin <- unique(hesin_diag_200k_redundant$eid)
eid_inter <- intersect(eid_hesin, eid_use)  # use IDs available in both datasets
disease_result <- disease_result[disease_result$eid %in% eid_inter, ]
all_simi <- disease_result$simi_score

# Initialize a data frame to store gene risk scores and related statistics
gene_risk_pvalue <- data.frame(matrix(NA, nrow = 0, ncol = 4))
colnames(gene_risk_pvalue) <- c("gene", "riskscore", "p_value", "carriers_num")
count <- 0

# Loop through each gene in gene_info, calculate risk scores for carriers
for (gene in gene_info$gene) {
  chr <- gene_info[gene_info$gene == gene, ]$chr_list
  if (length(chr) > 1) {
    chr <- chr[1]
  }
  carriers <- get_lof_carriers(gene, chr, hesin_diag_200k_redundant, eid_use)
  carriers <- intersect(carriers, eid_inter)
  if (length(carriers) >= 1) {
    carriers_simi <- disease_result[disease_result$eid %in% carriers, ]$simi_score
    riskscore <- mean(carriers_simi)
    wilcox_result <- wilcox.test(carriers_simi, all_simi, alternative = "greater")
    p_value <- wilcox_result$p.value
    tmp <- c(gene, riskscore, p_value, length(carriers))
    gene_risk_pvalue <- rbind(gene_risk_pvalue, tmp)
    count <- count + 1
    if (count %% 100 == 0) { cat(count, " ") }
  }
}
colnames(gene_risk_pvalue) <- c("gene", "riskscore", "p_value", "carriers_num")

# Save the gene-level risk scores using a relative path
outputfile <- paste0("data/riskscore_summary/riskscore/", phenotype,
                     "_gene_riskscore_weighted_IC_rare_all_diseaseCompare.RData")
save(gene_risk_pvalue, file = outputfile)
print("riskscore_gene finished!")

# =============================================================================
# BOOTSTRAP ANALYSIS TO COMPUTE GENE-RISK STATISTICS
# =============================================================================
# Load necessary files for bootstrap analysis
load(paste0("data/riskscore_summary/riskscore/", phenotype, "_riskscore_eur_weighted_IC_rare_all_diseaseCompare.RData"))
load(paste0("data/riskscore_summary/riskscore/", phenotype, "_gene_riskscore_weighted_IC_rare_all_diseaseCompare.RData"))
load("data/gene_carriers_num.RData")   # Ensure this file is in the data/ folder

# If merging with gene_carriers_num is needed, do so here.
gene_risk_pvalue_num <- gene_risk_pvalue
disease_result$simi_score <- as.numeric(disease_result$simi_score)
n_bootstrap <- 10000
gene_risk_pvalue_num$p_value_bootstrap <- NA

# Loop through each gene to perform bootstrap sampling and calculate p-values
for (j in 1:nrow(gene_risk_pvalue_num)) {
  gene <- gene_risk_pvalue_num$gene[j]
  carrier_num <- gene_risk_pvalue_num$carriers_num[j]
  bootstrap_means <- numeric(n_bootstrap)
  
  for (i in 1:n_bootstrap) {
    bootstrap_sample <- sample(disease_result$simi_score, size = carrier_num, replace = TRUE)
    bootstrap_means[i] <- mean(bootstrap_sample, na.rm = TRUE)
  }
  
  observed_mean <- as.numeric(gene_risk_pvalue_num$riskscore[j])
  bootstrap_mean <- mean(bootstrap_means)
  bootstrap_sd <- sd(bootstrap_means)
  z_score <- (observed_mean - bootstrap_mean) / bootstrap_sd
  p_value <- 1 - pnorm(z_score)
  gene_risk_pvalue_num$p_value_bootstrap[j] <- p_value
  
  if (j %% 100 == 0) { cat(j, " ") }
}

# Save the final bootstrap results using a relative file path
outputfile <- paste0("data/riskscore_summary/riskscore/", phenotype,
                     "_gene_riskscore_weighted_IC_rare_all_diseaseCompare_bootstrap.RData")
save(gene_risk_pvalue_num, file = outputfile)