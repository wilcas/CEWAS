#!/usr/bin/Rscript
library(glmnet)
library(R.matlab)
library(data.table)
library(DBI)
library(RSQLite)
library(dplyr)
library(doParallel)

###############
# GLOBAL VARS
###############
rnaseq_lib_cluster <- list(
  `11/11/13` = "base",
  `11/26/13` = "base",
  `17` = "A",
  `18` = "A",
  `10` = "A",
  `8/28/13` = "A",
  `14` = "B",
  `7` = "B",
  `11` = "B",
  `6` = "B",
  `16` = "B",
  `25` = "B",
  `5` = "C",
  `1` = "C",
  `3` = "C",
  `2` = "D",
  `28` = "D",
  `8` = "D",
  `12` = "D",
  `4` = "D",
  `20` = "E",
  `24` = "E",
  `26` = "E",
  `27` = "E",
  `9` = "E",
  `13` = "F",
  `15` = "F",
  `22` = "F",
  `21` = "F",
  `23` = "F",
  `10/15/13` = "G",
  `19` = "G",
  `10/9/13` = "H"
)
###############
# Utilities
###############
match_order <- function(ids_list) {
  shared <- Reduce(intersect, ids_list)
  res <- lapply(ids_list, function(x) na.omit(match(shared, x)))
  return(res)
}
##############
# Data Loading
##############

load_methylation <- function(fname) {
  pred_methy <- fread(fname)
  pred_methy$FID <- gsub("^([A-Z].*-[0-9]*)-.*", "\\1", pred_methy$FID)
  pred_methy$IID <- gsub("^([A-Z].*-[0-9]*)-.*", "\\1", pred_methy$FID)
  pred_methy$FID <- gsub("-", "_", pred_methy$FID)
  pred_methy$IID <- gsub("-", "_", pred_methy$IID)
  # replace with expression IDs
  mapping <- read.csv("/zfs3/scratch/saram_lab/commonMind/data/clinical/clinical.csv")
  pred_methy$IID <- mapping$DLPFC_RNA_Sequencing_Sample_ID[match(pred_methy$IID, mapping$Genotyping_Sample_ID)]
  out_methy <- as.matrix(pred_methy[, -c(1, 2)])
  rownames(out_methy) <- pred_methy$IID
  return(out_methy)
}

load_expression_remove_confounds <- function(fname) {
  clinical <- read.csv("/zfs3/scratch/saram_lab/commonMind/data/clinical/clinical.csv")
  clinical$Age_of_Death <- ifelse(clinical$Age_of_Death == "90+", 90, as.numeric(clinical$Age_of_Death))
  experimental <- read.csv("/zfs3/scratch/saram_lab/commonMind/data/expression/metaData.csv")
  experimental$LIB_BATCH <- c(unlist(rnaseq_lib_cluster[experimental$DLPFC_RNA_Sequencing_Library_Batch[1:239]]), NA, unlist(rnaseq_lib_cluster[experimental$DLPFC_RNA_Sequencing_Library_Batch[240:613]]))
  ancestry <- read.delim("/zfs3/scratch/saram_lab/commonMind/data/genotypeImputed/ancestry.tsv")
  all_covariates <- merge(clinical, ancestry, by = "Genotyping_Sample_ID")
  all_covariates <- merge(all_covariates, experimental, by = "DLPFC_RNA_Sequencing_Sample_ID")
  expression <- fread("/zfs3/scratch/saram_lab/commonMind/data/expression/CMC_MSSM-Penn-Pitt_DLPFC_mRNA_IlluminaHiSeq2500_gene-adjustedNoSVA-dataNormalization-noAncestry-adjustedLogCPM.tsv")
  expr <- as.matrix(expression[, -c(1)])
  rownames(expr) <- expression$GeneFeature
  colnames(expr) <- colnames(expression)[-c(1)]
  expr <- t(expr)
  all_covariates <- all_covariates %>% filter(LIB_BATCH != "H", LIB_BATCH != "base") # < 5% of total after filtering
  all_covariates <- all_covariates[as.character(all_covariates$DLPFC_RNA_Sequencing_Sample_ID) %in% rownames(expr), ]
  expr <- expr[as.character(all_covariates$DLPFC_RNA_Sequencing_Sample_ID), ]
  fit <- lm(expr ~ Dx + Institution + Gender + Age_of_Death + PMI_hrs + pH + DLPFC_RNA_isolation_RIN + DLPFC_RNA_isolation_RIN^2 + EV.1 + EV.2 + EV.3 + EV.4 + EV.5 + LIB_BATCH, data = all_covariates)
  expr <- resid(fit)
  print(summary(fit)[[1]])
  return(expr)
}
get_probes_in_1MB <- function(gene, probe_pos, gene_pos) {
  to_consider <- probe_pos[probe_pos$chr == gene_pos$chr[gene_pos$gene == gene], ]
  cur_gene_pos <- gene_pos$pos[gene_pos$gene == gene]
  return(to_consider$cpg[abs(to_consider$pos - cur_gene_pos) <= 500000])
}
write_to_predictdb <- function(db, res, gene, genename) {
  coef_mat <- coef(res, s = "lambda.min")
  cur_weights <- data.frame(
    rsid = rownames(coef_mat)[-c(1)],
    gene = gene,
    genename = genename,
    weight = as.numeric(coef_mat)[-c(1)],
    ref_allele = "A",
    eff_allele = "C"
  )
  cur_extra <- data.frame(
    gene = gene,
    genename = genename,
    pred.perf.R2 = res$glmnet.fit$dev.ratio[which(res$glmnet.fit$lambda == res$lambda.min)],
    n.snps.in.model = sum(as.numeric(coef_mat[-c(1)]) != 0),
    pred.perf.pval = 0.0,
    pred.perf.qval = 0.0
  )
  dbWriteTable(db, "weights", cur_weights, append = T)
  dbWriteTable(db, "extra", cur_extra, append = T)
}

main <- function() {
  # Load methy_data
  methy_data <- load_methylation("/zfs3/scratch/william.casazza/CMC_alt_fmt/robust_across_all_values_methy_rosmap_nonan.txt")
  expr_data <- load_expression_remove_confounds("/zfs3/scratch/saram_lab/commonMind/data/expression/CMC_MSSM-Penn-Pitt_DLPFC_mRNA_IlluminaHiSeq2500_gene-adjustedNoSVA-dataNormalization-noAncestry-adjustedLogCPM.tsv")
  # match IDs
  idx <- match_order(list(e_ids = rownames(expr_data), m_ids = rownames(methy_data)))
  expr_data <- expr_data[idx$e_ids, ]
  methy_data <- methy_data[idx$m_ids, ]
  # Load and process mapping
  tss_map <- fread("/zfs3/users/william.casazza/william.casazza/hg19_info/tss_map_hg19_filtered.txt")
  colnames(tss_map) <- c("gene", "chr", "pos", "genename")
  methy_map <- fread("~/jaffe_mQTL_coord.csv")
  registerDoParallel(4)
  # output
  db <- dbConnect(SQLite(), "cmc_imputed_methy_expr_model_par.db")
  for (cur_gene in colnames(expr_data)) {
    in_table <- tbl(db, "extra") %>%
      filter(gene == cur_gene) %>%
      count() %>%
      collect()
    if (in_table < 1) {
      cpgs <- get_probes_in_1MB(cur_gene, methy_map, tss_map)
      cpgs <- cpgs[cpgs %in% colnames(methy_data)]
      if (sum(cpgs %in% colnames(methy_data)) > 1) {
        X <- methy_data[, cpgs]
        y <- expr_data[, cur_gene]
        res <- cv.glmnet(X, y, parallel = TRUE)
        genename <- tss_map$genename[tss_map$gene == cur_gene]
        print(cur_gene)
        print(genename)
        write_to_predictdb(db, res, cur_gene, genename)
      } else {
        genename <- tss_map$genename[tss_map$gene == cur_gene]
        print(cur_gene)
        print(genename)
        print(sprintf("gene %s, %s has no nearby imputed methylation", cur_gene, genename))
      }
    }
  }
  dbDisconnect(db)
}

main()
