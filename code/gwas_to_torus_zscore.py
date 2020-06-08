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


def annotate_gwas_from_regions(df, region_fp):
    logging.info("Annotating GWAS with LD Regions")

    region_df = pd.read_csv(region_fp, engine='python', sep="\s+")
    region_df['chr'] = region_df['chr'].str.lstrip('chr').astype(int)

    df['region_id'] = ["-1"] * df.shape[0]
    for indx, region_ in enumerate(region_df.itertuples()):
        loc_str = "Loc" + str(indx + 1)
        df.loc[(df['chromosome'] == region_.chr) &
               (df['position'] < region_.stop) &
               (df['position'] >= region_.start), 'region_id'] = loc_str
    return df


def write_torus_GWAS(df, fp):
    df = df[['id', 'region_id', 'zscore']]
    logging.info("Writing output")
    df.to_csv(fp, sep="\t", compression='gzip', header=False, index=False)


def run(args):
    start = timer()
    if os.path.isfile(args.output_fp):
        raise ValueError("Output filepath already exists")
    gwas_df = load_gwas(args.input_gwas)
    annot_df = annotate_gwas_from_regions(gwas_df, args.input_ld_regions)
    write_torus_GWAS(annot_df, args.output_fp)
    logging.info("Finished in {:.2f} seconds".format(timer() - start))



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-input_gwas')
    parser.add_argument('-input_ld_regions')
    parser.add_argument('-output_fp')

    args = parser.parse_args()

    logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        stream=sys.stdout,
                        level=logging.INFO)

    run(args)