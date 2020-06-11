import logging
import os
import sys
import argparse
import pandas as pd
from timeit import default_timer as timer

def load_gwas(fp):
    logging.info("Loading GWAS")
    df = pd.read_csv(fp, sep="\t", usecols=['zscore',
                                              'panel_variant_id',
                                              'position',
                                              'chromosome'])
    df = df.rename(mapper={"panel_variant_id": "id"}, axis=1)
    df['chromosome'] = df['chromosome'].str.lstrip('chr').astype(int)
    return df

def _load_ld_regions(fp, name, start, end, chromosome):
    rename_dd = {name: 'region_id',
                 start: 'start',
                 end: 'end',
                 chromosome: 'chromosome'}
    region_df = pd.read_csv(fp)
    region_df = region_df.rename(rename_dd, axis=1)
    region_df['chr'] = region_df.chromosome.str.lstrip('chr').astype(int)
    return region_df


def annotate_gwas_from_regions(df, region_fp, name, start, end):
    logging.info("Annotating GWAS with LD Regions")
    region_df = _load_ld_regions(region_fp, name, start, end, 'chromosome')

    df['region_id'] = ["NA"] * df.shape[0]
    for chr, region_df_chr in region_df.groupby('chr'):
        logging.info("Annotating chromosome {}".format(chr))
        df_chr = df.loc[df['chromosome'] == chr]
        for region_ in region_df_chr.itertuples():
            df_chr.loc[(df_chr['position'] < region_.end) &
                       (df['position'] >= region_.start), 'region_id'] = region_.region_id
    return df


def write_torus_GWAS(df, fp):
    df = df[['id', 'region_id', 'zscore']]
    logging.info("Writing output")
    df.to_csv(fp, sep="\t", header=False, index=False)


def run(args):
    start = timer()
    if os.path.isfile(args.output_fp):
        raise ValueError("Output filepath already exists")
    gwas_df = load_gwas(args.input_gwas)
    annot_df = annotate_gwas_from_regions(gwas_df,
                                          args.input_ld_regions,
                                          args.name_key,
                                          args.start_key,
                                          args.end_key)
    write_torus_GWAS(annot_df, args.output_fp)
    logging.info("Finished in {:.2f} seconds".format(timer() - start))



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-input_gwas')
    parser.add_argument('-input_ld_regions')
    parser.add_argument('-output_fp')
    parser.add_argument("--name_key", default="region")
    parser.add_argument("--start_key", default="start")
    parser.add_argument("--end_key", default="end")

    args = parser.parse_args()

    logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        stream=sys.stdout,
                        level=logging.INFO)

    run(args)