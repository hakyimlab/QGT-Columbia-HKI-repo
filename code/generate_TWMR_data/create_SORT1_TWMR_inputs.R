# --------------------------------------------------------------------------------------
# ENVIRONMENT
library(tidyverse)
library(Matrix)

# --------------------------------------------------------------------------------------
# DEFINE FILE PATHS

cov_fp <- "/Users/owenmelia/data/QGT-Columbia-HKI/models/SORT1_cov.txt.gz"
eqtl_fp <- "/Users/owenmelia/data/QGT-Columbia-HKI/models/SORT1_eqtl.csv"
cad_gwas_fp <- "/Users/owenmelia/data/QGT-Columbia-HKI/data/spredixcan/imputed_CARDIoGRAM_C4D_CAD_ADDITIVE.txt"
weights_fp <- "/Users/owenmelia/data/QGT-Columbia-HKI/models/SORT1_liver_locus_predictive_snps.csv"

out_ld <- "/Users/owenmelia/data/QGT-Columbia-HKI/repos/TWMR-master/ENSG00000134243.ld"
out_matrix <- "/Users/owenmelia/data/QGT-Columbia-HKI/repos/TWMR-master/ENSG00000134243.matrix"
out_dir <- "/Users/owenmelia/data/QGT-Columbia-HKI/repos/TWMR-master/"
K_SORT1 <- "ENSG00000134243"

# --------------------------------------------------------------------------------------
# LOAD INPUT DATA

cov_df <- read.table(cov_fp, header=TRUE)
eqtl_df <- read.csv(eqtl_fp)
cad_gwas <- read.table(cad_gwas_fp, header=TRUE, sep="\t", nrows = 2000000)
weights <- read.csv(weights_fp)
'%!in%' <- function(x,y)!('%in%'(x,y))
# --------------------------------------------------------------------------------------
# DEFINE FUNCTIONS TO GENERATE TABLES

keep_snps <- weights$varID
keep_genes <- c("ENSG00000134222.16", "ENSG00000143126.7", "ENSG00000134243.11")
eqtl_to_matrix <- function(eqtl, trait_effects, keep_snps, keep_genes){
  eqtl <- eqtl %>% filter(variant_id %in% keep_snps & gene_id %in% keep_genes)
  eqtl$gene_id <- gsub("\\..*", "", eqtl$gene_id)
  out_df <- (eqtl %>% 
               select(c(variant_id, gene_id, slope)) %>% 
               pivot_wider(names_from = gene_id, values_from= slope))
  out_df$GENES <- out_df$variant_id
  out_df$variant_id <- NULL
  # print(head(out_df))
  trait_effects <- (trait_effects %>% select(c(panel_variant_id, effect_size)))

  out_df <- left_join(out_df, trait_effects, by=c('GENES'='panel_variant_id')) %>% rename('BETA_GWAS'='effect_size')
  print(out_df)
  out_df <- drop_na(out_df)
  return(out_df)
}
SORT1_matrix <- eqtl_to_matrix(eqtl_df, cad_gwas, keep_snps, keep_genes)

makeSymm <- function(m) {
  # Thanks to StackOverflow user Ben Bolker
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  return(m)
}

create_ld <- function(cov_df, snp_lst){
  cov_df <- cov_df %>% filter(RSID1 %in% snp_lst & RSID2 %in% snp_lst)
  cov_df <- cov_df %>% pivot_wider(names_from='RSID1', values_from='VALUE') %>% select(-c(GENE))
  # cov_df <- cov_df[match(snp_lst, cov_df$RSID2),]
  print(head(cov_df))
  cov_rows <- cov_df$RSID2
  # cov_df <- cov_df[snp_lst]
  cov_df$RSID2 <- NULL
  cov_matrix <- cov_df %>% data.matrix
  rownames(cov_matrix) <- cov_rows
  cov_matrix <- makeSymm(cov_matrix)

  if (all(rownames(cov_matrix) == colnames(cov_matrix))){
    return(cov_matrix)
  }
  else{
    print("ERROR")
    return(NULL)
  }
}

SORT1_ld <- create_ld(cov_df, SORT1_matrix$GENES)

# Reorder SNPs
reorder_matrix <- function(df, snp_lst){
  df <- df[match(snp_lst, df$GENES),]
  
  df <- df[, c(4, 1:3, 5)]
  return(df)
}
SORT1_matrix <- reorder_matrix(SORT1_matrix, rownames(SORT1_ld))


fix_singularity <- function(df){
  return(df %>% select(-c(ENSG00000134222)))
}

# SORT1_matrix <- fix_singularity(SORT1_matrix)
# --------------------------------------------------------------------------------------
# WRITE DATA

writer <- function(df, fp, c=TRUE){
  print(paste0("Writing to ", fp))
  write.table(df, fp, sep=" ", quote=FALSE, row.names = FALSE, col.names = c)
}

writer(SORT1_matrix, out_matrix)

writer(SORT1_ld, out_ld, c=FALSE)


write_permutations <- function(dir, m, ld ){
  gene_names <- names(m)[str_detect(names(m), 'ENS.*')]
  for ( gene in gene_names ){
    cols <- c("GENES", gene, setdiff(gene_names, gene), "BETA_GWAS")
    m_i <- m[, cols]
    fp_matrix <- file.path(dir, paste(gene, 'matrix', sep="."))
    fp_ld <- file.path(dir, paste(gene, 'ld', sep="."))
    writer(m_i, fp_matrix)
    writer(ld, fp_ld, c=FALSE)
  }
}

write_permutations(out_dir, SORT1_matrix, SORT1_ld)

# --------------------------------------------------------------------------------------
# TEST

source("~/projects/QGT-Columbia-HKI/code/TWMR_script.R")
DIR <- "~/data/QGT-Columbia-HKI/repos/TWMR-master"
analyseMR(K_SORT1, DIR)
