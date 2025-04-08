library(lsa)
library(data.table)
library(dplyr)

# Define a base directory for carrier files (modify this to match your repository structure)
carrier_base_dir <- "data/phe_sig_pipeline/file/carriers/LoFDmis"

regres_sig_phe <- function(carry_eid, non_carry_eid, hesin_diag, cov_adjust, ICD10_mapping){
  icd10_carry <- unique(hesin_diag[hesin_diag$eid %in% carry_eid, "diag_icd10"])
  print(length(icd10_carry))
  dat <- cov_adjust[, c('eid','age_recruit','sex')]
  dat$icd <- 0
  dat$var <- 0
  dat <- dat[dat$eid %in% c(carry_eid, non_carry_eid), ]
  dat[dat$eid %in% carry_eid, "var"] <- 1
  model_result <- data.frame(matrix(ncol = 12))
  colnames(model_result) <- c('icd10','carry_number','carry_proportion', 'non-carry_proportion',
                              'Intercept', 'sexMale', 'age', 'variant',
                              'p_Intercept', 'p_sexMale', 'p_age', 'p_variant')
  flag <- 0
  for(i in icd10_carry){
    flag <- flag + 1
    cat(flag, '')
    if(flag %% 100 == 0){ print('') }
    eid_temp <- hesin_diag[hesin_diag$diag_icd10 == i, "eid"]
    temp <- dat
    temp[temp$eid %in% eid_temp, "icd"] <- 1
    model <- glm(formula = icd ~ sex + age_recruit + var, data = temp, family = binomial())
    summary_temp <- summary(model)
    prop_carry <- nrow(temp[temp$icd == 1 & temp$var != 0, ]) / nrow(temp[temp$var != 0, ])
    num_carry <- nrow(temp[temp$icd == 1 & temp$var != 0, ])
    prop_non_carry <- nrow(temp[temp$icd == 1 & temp$var == 0, ]) / nrow(temp[temp$var == 0, ])
    vec_temp <- c(i, num_carry, prop_carry, prop_non_carry,
                  as.numeric(summary_temp$coefficients[, 1]),
                  as.numeric(summary_temp$coefficients[, 4]))
    model_result <- rbind(model_result, vec_temp)
  }
  model_result <- model_result[-1, ]
  model_result$p_variant <- as.numeric(model_result$p_variant)
  compare <- model_result %>% dplyr::select(icd10,
                                            pheno_num_in_carrier = `carry_number`,
                                            pheno_prop_in_carrier = `carry_proportion`,
                                            pheno_prop_in_noncarrier = `non-carry_proportion`,
                                            pvalue_variants = `p_variant`,
                                            effsize_allele = variant)
  compare$OR <- exp(as.numeric(compare$effsize_allele))
  compare$ICD_diag <- ICD10_mapping$meaning[match(compare$icd10, ICD10_mapping$coding)]
  compare <- compare %>% mutate(across(2:7, as.numeric))
  pvalue_adjusted <- p.adjust(compare$pvalue_variants, method = "bonferroni")
  compare$pvalue_adjusted <- pvalue_adjusted
  
  return(compare)
}

regres_sig_phe_SE <- function(carry_eid, non_carry_eid, hesin_diag, cov_adjust, ICD10_mapping){
  icd10_carry <- unique(hesin_diag[hesin_diag$eid %in% carry_eid, "diag_icd10"])
  print(length(icd10_carry))
  dat <- cov_adjust[, c('eid','age_recruit','sex')]
  dat$icd <- 0
  dat$var <- 0
  dat <- dat[dat$eid %in% c(carry_eid, non_carry_eid), ]
  dat[dat$eid %in% carry_eid, "var"] <- 1
  model_result <- data.frame(matrix(ncol = 16))
  colnames(model_result) <- c('icd10','carry_number','carry_proportion', 'non-carry_proportion',
                              'Intercept', 'sexMale', 'age', 'variant',
                              'p_Intercept', 'p_sexMale', 'p_age', 'p_variant',
                              'se_Intercept', 'se_sexMale','se_age', 'se_variant')
  flag <- 0
  for(i in icd10_carry){
    flag <- flag + 1
    cat(flag, '')
    if(flag %% 100 == 0){ print('') }
    eid_temp <- hesin_diag[hesin_diag$diag_icd10 == i, "eid"]
    temp <- dat
    temp[temp$eid %in% eid_temp, "icd"] <- 1
    model <- glm(formula = icd ~ sex + age_recruit + var, data = temp, family = binomial())
    summary_temp <- summary(model)
    prop_carry <- nrow(temp[temp$icd == 1 & temp$var != 0, ]) / nrow(temp[temp$var != 0, ])
    num_carry <- nrow(temp[temp$icd == 1 & temp$var != 0, ])
    prop_non_carry <- nrow(temp[temp$icd == 1 & temp$var == 0, ]) / nrow(temp[temp$var == 0, ])
    vec_temp <- c(i, num_carry, prop_carry, prop_non_carry,
                  as.numeric(summary_temp$coefficients[, 1]),
                  as.numeric(summary_temp$coefficients[, 4]),
                  as.numeric(summary_temp$coefficients[, 2]))
    model_result <- rbind(model_result, vec_temp)
  }
  model_result <- model_result[-1, ]
  model_result$p_variant <- as.numeric(model_result$p_variant)
  compare <- model_result %>% dplyr::select(icd10,
                                            pheno_num_in_carrier = `carry_number`,
                                            pheno_prop_in_carrier = `carry_proportion`,
                                            pheno_prop_in_noncarrier = `non-carry_proportion`,
                                            pvalue_variants = `p_variant`,
                                            effsize_allele = variant,
                                            se_variants = `se_variant`)
  compare$OR <- exp(as.numeric(compare$effsize_allele))
  compare$ICD_diag <- ICD10_mapping$meaning[match(compare$icd10, ICD10_mapping$coding)]
  compare <- compare %>% mutate(across(2:8, as.numeric))
  pvalue_adjusted <- p.adjust(compare$pvalue_variants, method = "bonferroni")
  compare$pvalue_adjusted <- pvalue_adjusted
  
  return(compare)
}

assign_weights <- function(weight) {
  # Transform weight (e.g., -log10 p-value) into a significance category
  if(weight > 10) {
    return(2)       # Very significant
  } else if(weight > 2) {
    return(1)       # Moderately significant
  } else if(weight > 1) {
    return(0.7)     # Less moderately significant
  } else {
    return(0.2)     # Slightly significant
  }
}

assign_weights_logPvalue <- function(weight){
  if(weight == Inf){
    return(300)
  } else {
    return(weight + 0.01)
  }
}

assign_weights_sigmoid <- function(weight){
  1 / (1 + exp(-weight))
}

get_avg_embed <- function(eid1, embed = embed_icd_redundant, hesin_table = hesin_diag_200k_redundant, disease = c('Q612','Q613')){
  eid1 <- hesin_table[hesin_table$eid == eid1, "diag_icd10"]
  eid_embed <- rep(0, ncol(embed))
  for(i in seq_along(eid1)){
    embed_temp <- as.numeric(embed[eid1[i], ])
    eid_embed <- eid_embed + embed_temp
  }
  eid_embed <- eid_embed / length(eid1)
  return(eid_embed)
}

get_avg_embed_weight <- function(eid1, embed = embed_icd_redundant,
                                 hesin_table = hesin_500k_redundant,
                                 disease = c('Q612','Q613'),
                                 weight_table){
  eid1 <- hesin_table[hesin_table$eid == eid1, "diag_icd10"]
  if(length(eid1) == 0){ return(9) }
  eid_embed <- rep(0, ncol(embed))
  sum_weight <- 0
  for(i in seq_along(eid1)){
    if(eid1[i] %in% rownames(embed)){
      embed_temp <- as.numeric(embed[eid1[i], ])
      weight_temp <- as.numeric(weight_table[eid1[i], 2])
      if(is.na(weight_temp)){ weight_temp <- min(weight_table[, 2]) }
      eid_embed <- eid_embed + embed_temp * weight_temp
      sum_weight <- sum_weight + weight_temp
    }
  }
  eid_embed <- eid_embed / sum_weight
  return(eid_embed)
}

get_avg_embed_weight_IC <- function(eid1, embed = embed_icd_redundant,
                                    hesin_table = hesin_500k_redundant,
                                    weight_table,
                                    IC_table){
  eid1 <- hesin_table[hesin_table$eid == eid1, "diag_icd10"]
  if(length(eid1) == 0){ return(9) }
  eid_embed <- rep(0, ncol(embed))
  sum_weight <- 0
  for(i in seq_along(eid1)){
    if(eid1[i] %in% intersect(rownames(embed), rownames(IC_table))){
      embed_temp <- as.numeric(embed[eid1[i], ])
      weight_temp <- as.numeric(weight_table[eid1[i], 2])
      IC_tmp <- as.numeric(IC_table[eid1[i], 'IC'])
      if(is.na(weight_temp)){ weight_temp <- min(weight_table[, 2]) }
      eid_embed <- eid_embed + embed_temp * weight_temp * IC_tmp
      sum_weight <- sum_weight + weight_temp * IC_tmp
    }
  }
  eid_embed <- eid_embed / sum_weight
  return(eid_embed)
}

generate_dataset <- function(disease_id, non_disease_id, hesin_table = hesin_500k_redundant, cov_adjust){
  set.seed(1234)
  cov_adjust <- as.data.frame(cov_adjust)
  rownames(cov_adjust) <- cov_adjust$eid
  compare <- sample(disease_id, length(disease_id) / 4)
  treat <- setdiff(disease_id, compare)
  control <- sample(non_disease_id, length(treat))
  cov_train <- cov_adjust[as.character(union(treat, control)), ]
  cov_train <- cov_train[, c(1, 2, 3)]
  train_id <- union(disease_id, control)
  hesin_train <- hesin_table[hesin_table$eid %in% train_id, ]
  cov_train$embed_simi <- 0
  cov_train$disease <- 0
  cov_train[cov_train$eid %in% disease_id, "disease"] <- 1
  out <- list(compare = compare, cov_train = cov_train, hesin_train = hesin_train)
  return(out)
}

comapre_simi <- function(embedding_carriers, embedding_patients){
  simi_score <- data.frame(matrix(nrow = nrow(embedding_carriers), ncol = 2))
  colnames(simi_score) <- c('eid', 'simi_score')
  simi_score$eid <- rownames(embedding_carriers)
  for(i in seq_len(nrow(simi_score))){
    eid1 <- simi_score$eid[i]
    eid1_embed <- as.numeric(embedding_carriers[as.character(eid1), ])
    score <- 0
    for(j in seq_len(nrow(embedding_patients))){
      eid2 <- rownames(embedding_patients)[j]
      eid2_embed <- as.numeric(embedding_patients[as.character(eid2), ])
      if(is.na(as.numeric(cosine(eid1_embed, eid2_embed)))){
        cat(i, j, 'NA', '\n')
      }
      score <- score + as.numeric(cosine(eid1_embed, eid2_embed))
    }
    score <- score / nrow(embedding_patients)
    simi_score$simi_score[i] <- score
  }
  return(simi_score)
}

compare_carrier_disease <- function(gene, chr,
                                    hesin_diag,
                                    eid_whole,
                                    embed_icd,
                                    embedding_disease,
                                    weight_sigPhe){
  # Build the relative file path using the carrier_base_dir
  pathname <- paste0(carrier_base_dir, "/chr", chr, "/Carrier_chr", chr, "_", gene, ".txt")
  if(file.exists(pathname)){
    lofDmis_carriers <- fread(pathname)
    lof_carriers <- unique(lofDmis_carriers[lofDmis_carriers$anno == 'LoF', "sample"])
    lof_carriers <- intersect(lof_carriers, hesin_diag$eid)
    non_lof_carriers <- intersect(setdiff(eid_whole, unique(lofDmis_carriers$sample)), hesin_diag$eid)
    cat(gene, 'LoF carriers number is:', length(lof_carriers), '\n')
  
    # Embed for carriers
    embedding_gene <- data.frame(matrix(nrow = length(lof_carriers), ncol = ncol(embed_icd) + 1))
    embedding_gene[, 1] <- lof_carriers
    rownames(embedding_gene) <- embedding_gene[, 1]
    embedding_gene <- embedding_gene[, -1]
    
    hesin <- hesin_diag[hesin_diag$eid %in% lof_carriers, ]
    for(i in seq_len(nrow(embedding_gene))){
      eid_temp <- rownames(embedding_gene)[i]
      embed_temp <- get_avg_embed_weight(eid_temp, embed_icd, hesin, weight_table = weight_sigPhe)
      embedding_gene[i, ] <- embed_temp
      if(i %% 50 == 0){ cat(i, ' ') }
    }
    simi_score_gene <- comapre_simi(embedding_carriers = embedding_gene,
                                    embedding_patients = embedding_disease)
    cat('compare', gene, 'carriers done! \n')
    return(simi_score_gene)
  } else {
    simi_score_gene <- data.frame(simi_score = 0)
    return(simi_score_gene)
  }
}

get_info <- function(gene){
  result <- getBM(attributes = c('hgnc_symbol', 'chromosome_name','start_position', 'end_position'),
                  filters = 'hgnc_symbol', values = gene, mart = mart)
  chr <- result$chromosome_name[nrow(result)]
  start <- result$start_position[nrow(result)]
  end <- result$end_position[nrow(result)]
  info <- c(chr, start, end)
  return(info)
}

compare_group_disease <- function(id, hesin_diag = hesin_diag_200k,
                                  eid_whole = eid_all,
                                  embed_icd = embed_icd_redundant,
                                  embedding_disease,
                                  weight_sigPhe = PKD_weight_sigPhe){
  # Embed for the given IDs (carriers)
  embedding_gene <- data.frame(matrix(nrow = length(id), ncol = ncol(embed_icd) + 1))
  embedding_gene[, 1] <- id
  rownames(embedding_gene) <- embedding_gene[, 1]
  embedding_gene <- embedding_gene[, -1]
  
  hesin <- hesin_diag[hesin_diag$eid %in% id, ]
  for(i in seq_len(nrow(embedding_gene))){
    eid_temp <- rownames(embedding_gene)[i]
    embed_temp <- get_avg_embed_weight(eid_temp, embed_icd, hesin, weight_table = weight_sigPhe)
    embedding_gene[i, ] <- embed_temp
    if(i %% 50 == 0){ cat(i, ' ') }
  }
  simi_score_gene <- comapre_simi(embedding_carriers = embedding_gene,
                                  embedding_patients = embedding_disease)
  return(simi_score_gene)
}

compare_group_disease_avg <- function(id, hesin_diag = hesin_diag_200k,
                                  eid_whole = eid_all,
                                  embed_icd = embed_icd_redundant,
                                  embedding_disease){
  # Embed for the given IDs using the unweighted average embed function
  embedding_gene <- data.frame(matrix(nrow = length(id), ncol = ncol(embed_icd) + 1))
  embedding_gene[, 1] <- id
  rownames(embedding_gene) <- embedding_gene[, 1]
  embedding_gene <- embedding_gene[, -1]
  
  hesin <- hesin_diag[hesin_diag$eid %in% id, ]
  for(i in seq_len(nrow(embedding_gene))){
    eid_temp <- rownames(embedding_gene)[i]
    embed_temp <- get_avg_embed(eid_temp, embed_icd, hesin)
    embedding_gene[i, ] <- embed_temp
    if(i %% 50 == 0){ cat(i, ' ') }
  }
  simi_score_gene <- comapre_simi(embedding_carriers = embedding_gene,
                                  embedding_patients = embedding_disease)
  return(simi_score_gene)
}

compare_noncarrier_disease <- function(gene, chr,
                                       hesin_diag = hesin_diag_200k,
                                       eid_whole = eid_all,
                                       embed_icd = embed_icd_redundant,
                                       embedding_disease,
                                       weight_sigPhe = PKD_weight_sigPhe,
                                       disease_id = PKD_id){
  pathname <- paste0(carrier_base_dir, "/chr", chr, "/Carrier_chr", chr, "_", gene, ".txt")
  if(file.exists(pathname)){
    lofDmis_carriers <- fread(pathname)
    lof_carriers <- unique(lofDmis_carriers[lofDmis_carriers$anno == 'LoF', "sample"])
    lof_carriers <- intersect(lof_carriers, hesin_diag$eid)
    non_lof_carriers <- intersect(setdiff(eid_whole, unique(lofDmis_carriers$sample)), hesin_diag$eid)
    if(length(lof_carriers) > 40) {
      size1 <- 3 * length(lof_carriers)
    } else {
      size1 <- 80
    }
    non_lof_carriers_random <- sample(non_lof_carriers, size = size1, replace = FALSE)
    
    cat(gene, 'random non-carriers number is:', length(non_lof_carriers_random), '\n')
    
    embedding_gene <- data.frame(matrix(nrow = length(non_lof_carriers_random), ncol = ncol(embed_icd) + 1))
    embedding_gene[, 1] <- non_lof_carriers_random
    rownames(embedding_gene) <- embedding_gene[, 1]
    embedding_gene <- embedding_gene[, -1]
    
    hesin <- hesin_diag[hesin_diag$eid %in% non_lof_carriers_random, ]
    for(i in seq_len(nrow(embedding_gene))){
      eid_temp <- rownames(embedding_gene)[i]
      embed_temp <- get_avg_embed_weight(eid_temp, embed_icd, hesin, weight_table = weight_sigPhe)
      embedding_gene[i, ] <- embed_temp
      if(i %% 50 == 0){ cat(i, ' ') }
    }
    simi_score_gene <- comapre_simi(embedding_carriers = embedding_gene,
                                    embedding_patients = embedding_disease)
    cat('compare', gene, 'non-carriers done! \n')
    return(simi_score_gene)
  } else {
    simi_score_gene <- data.frame(simi_score = 0)
    return(simi_score_gene)
  }
}


get_lof_carriers <- function(gene, chr, hesin_diag, eid_whole){
  pathname <- paste0(carrier_base_dir, "/chr", chr, "/Carrier_chr", chr, "_", gene, ".txt")
  if(file.exists(pathname)){
    lofDmis_carriers <- fread(pathname)
    lof_carriers <- unique(lofDmis_carriers[lofDmis_carriers$anno == 'LoF', "sample"])
    lof_carriers <- intersect(lof_carriers, hesin_diag$eid)
    lof_carriers <- intersect(lof_carriers, eid_whole)
    non_lof_carriers <- intersect(setdiff(eid_whole, unique(lofDmis_carriers$sample)), hesin_diag$eid)
    return(lof_carriers)
  }
  return(c())
}

get_non_carriers <- function(gene, chr, hesin_diag, eid_whole){
  pathname <- paste0(carrier_base_dir, "/chr", chr, "/Carrier_chr", chr, "_", gene, ".txt")
  if(file.exists(pathname)){
    lofDmis_carriers <- fread(pathname)
    lof_carriers <- unique(lofDmis_carriers[lofDmis_carriers$anno == 'LoF', "sample"])
    lof_carriers <- intersect(lof_carriers, hesin_diag$eid)
    non_lof_carriers <- intersect(setdiff(eid_whole, unique(lofDmis_carriers$sample)), hesin_diag$eid)
    return(non_lof_carriers)
  }
  return(c())
}