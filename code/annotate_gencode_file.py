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


def annotate_gencode_from_regions(df, fp, name, start, end):
    logging.info("Annotating Gencode with LD Regions")
    region_df = _load_ld_regions(fp, name, start, end, chromosome='chromosome')
    df['region_id'] = ["NA"] * df.shape[0]
    for indx, region_ in enumerate(region_df.itertuples()):
        loc_str = region_.region_id
        df.loc[(df['chr'] == region_.chr) &
               (df['start_location'] < region_.end) &
               (df['start_location'] >= region_.start), 'region_id'] = loc_str
    return df


def write_file(df, fp):
    logging.info("Writing output")
    df.to_csv(fp, sep="\t", compression='gzip', index=False)


def load_gencode(fp):
    logging.info("Loading gencode table")
    gencode_df = pd.read_csv(fp, sep = "\t")
    autosomal_chr = { 'chr' + str(x) for x in range(1,23) }
    gencode_df = gencode_df[gencode_df['chromosome'].isin(autosomal_chr)]
    gencode_df['chr'] = gencode_df['chromosome'].str.lstrip('chr').astype(int)
    return gencode_df


def run(args):
    start = timer()
    if os.path.isfile(args.output_fp):
        raise ValueError("Output filepath already exists")
    gencode_df = load_gencode(args.input_gencode)

    annot_df = annotate_gencode_from_regions(gencode_df,
                                             args.input_ld_regions,
                                             args.name_key,
                                             args.start_key,
                                             args.end_key)
    write_file(annot_df, args.output_fp)
    logging.info("Finished in {:.2f} seconds".format(timer() - start))



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-input_gencode')
    parser.add_argument('-input_ld_regions', help="Filepath for LD Block file")
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