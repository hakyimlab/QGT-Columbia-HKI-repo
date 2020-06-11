python /vol/bmd/meliao/software/summary-gwas-imputation/src/covariance_for_model.py \
-parquet_genotype_folder /vol/bmd/meliao/data/reference_panel_1000G \
-parquet_genotype_pattern "chr(.*).variants.parquet" \
-model_db /vol/bmd/meliao/data/QGT/SORT1_model.db \
--parquet_metadata /vol/bmd/meliao/data/reference_panel_1000G/variant_metadata.parquet \
-parsimony 2 \
-output /vol/bmd/meliao/data/QGT/SORT1_cov.txt