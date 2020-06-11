import pandas as pd
import sqlite3
# -----------------------------------------------------------------------------
# SET FILEPATHS

weights_fp = "/Users/owenmelia/data/QGT-Columbia-HKI/models/SORT1_locus_predictive_snps.csv"
out_prefix = "/Users/owenmelia/data/QGT-Columbia-HKI/models/SORT1_"
gwas_fp = "/Users/owenmelia/data/QGT-Columbia-HKI/data/spredixcan/imputed_CARDIoGRAM_C4D_CAD_ADDITIVE.txt"
query_out_fp = "/Users/owenmelia/projects/QGT-Columbia-HKI/code/generate_TWMR_data/BQ_eqtl_query.txt"

db_out = out_prefix + "model.db"
txt_out = out_prefix + "model.txt"

# -----------------------------------------------------------------------------
# DEFINE CONSTANTS

K_GENE = "ENSG00000000"
K_QUERY = 'SELECT * FROM `gtex-awg-im.GTEx_v8_eQTL.Whole_Blood_allpairs` WHERE substr(gene_id,1,15) in ("ENSG00000134243", "ENSG00000143106", "ENSG00000134222", "ENSG00000143126", "ENSG00000031698") AND variant_id in ({})'

# -----------------------------------------------------------------------------
# GENERATE MODEL

conn = sqlite3.connect(db_out)
weights_df = pd.read_csv(weights_fp)
weights_df = weights_df.drop('tissue', axis=1).drop_duplicates(subset=['varID'])

gwas_df = pd.read_csv(gwas_fp, sep="\t", nrows=2000000)

weights_df = weights_df[weights_df['varID'].isin(gwas_df['panel_variant_id'])]

weights_df['gene'] = [K_GENE] * weights_df.shape[0]

extra_df = pd.DataFrame({'gene':[K_GENE],
                         'gene_id': [K_GENE],
                         'gene_name': [K_GENE]})

# -----------------------------------------------------------------------------
# WRITE

weights_df.to_sql('weights', conn)

extra_df.to_sql('extra', conn)

weights_df.to_csv(txt_out, sep="\t", index=False)

# -----------------------------------------------------------------------------
# GENERATE BQ QUERY

snp_lst = ['"{}"'.format(x) for x in weights_df.varID.values]
snp_str = ", ".join(snp_lst)

with open(query_out_fp, 'w') as f:
    f.write(K_QUERY.format(snp_str))