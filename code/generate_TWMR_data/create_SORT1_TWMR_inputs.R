# --------------------------------------------------------------------------------------
# ENVIRONMENT
library(tidyverse)
library(Matrix)

# --------------------------------------------------------------------------------------
# DEFINE FILE PATHS

cov_fp <- "/Users/owenmelia/data/QGT-Columbia-HKI/models/SORT1_cov.txt.gz"
eqtl_fp <- "/Users/owenmelia/data/QGT-Columbia-HKI/models/SORT1_eqtl.csv"
cad_gwas_fp <- "/Users/owenmelia/data/QGT-Columbia-HKI/data/spredixcan/imputed_CARDIoGRAM_C4D_CAD_ADDITIVE.txt"

out_ld <- "/Users/owenmelia/data/QGT-Columbia-HKI/repos/TWMR-master/ENSG00000134243.ld"
out_matrix <- "/Users/owenmelia/data/QGT-Columbia-HKI/repos/TWMR-master/ENSG00000134243.matrix"

# --------------------------------------------------------------------------------------
# LOAD INPUT DATA

cov_df <- read.table(cov_fp, header=TRUE)
weights_df <- read.csv(weights_fp)
cad_gwas <- read.table(cad_gwas_fp, header=TRUE, sep="\t", nrows = 2000000)

# --------------------------------------------------------------------------------------
# DEFINE FUNCTIONS TO GENERATE TABLES

weights_to_matrix <- function(weights, expression_effects, trait_effects, gene_id="ENSG00000134243"){
  out_df <- (weights %>% distinct(varID, .keep_all=TRUE) %>% select(c('GENES'='varID', 'weight'='weight')))
  out_df[[gene_id]] <- out_df$weight
  out_df$weight <- NULL
  trait_effects <- (trait_effects %>% select(c(panel_variant_id, effect_size)))

  out_df <- left_join(out_df, trait_effects, by=c('GENES'='panel_variant_id')) %>% rename('BETA_GWAS'='effect_size')
  return(out_df)
}
SORT1_matrix <- weights_to_matrix(weights_df, NULL, cad_gwas)

SORT1_matrix %>% drop_na %>% select(GENES)

create_ld <- function(cov_df, snp_lst){
  cov_df %>% pivot_wider(names_from='RSID1', values_from='VALUE')
}
aa <- create_ld(cov_df, SORT1_matrix$GENES)



# --------------------------------------------------------------------------------------
# WRITE DATA

writer <- function(df, fp){
  print(paste0("Writing to ", fp))
  write.table(df, fp, sep=" ", quote=FALSE)
}

writer(SORT1_matrix, out_matrix)

writer(SORT1_ld, out_ld)
