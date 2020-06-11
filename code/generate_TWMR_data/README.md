# Generating TWMR data for the SORT1 gene

1. Haky gave me a list of genes in the SORT1 locus. Here is the list of genes: `PSMA5, PSRC1, CELSR2, SORT1, SARS`. This translates to ensemblIDs `"ENSG00000134243", "ENSG00000143106", "ENSG00000134222", "ENSG00000143126", "ENSG00000031698"`
1. Query BigQuery for prediction snps in these genes. The query script is at `BQ_prediction_SNP_query.txt`. 
1. Filter out SNPs for which we don't have GWAS results. Then put the resulting SNPs in a predictDB database. This code is in `weights_to_model.py` 
1. Run Alvaro's `covariance_for_model.py` script. I did this on Washington. The code is at `generate_covariance.sh`. This uses a reference panel from 1000G as the source of the LD. 
1. 