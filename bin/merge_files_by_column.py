#!/usr/bin/env python3

import argparse
import os
import pandas

def parse_args():
    """ Parse commandline arguments """
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-i", "--indir", default=None, type=str,
                        required=True, help="The directory containing all files to combine")
    parser.add_argument("-s", "--suffix", default=None, type=str,
                        required=True, help="suffix of the files to merge")
    parser.add_argument("-o", "--outfile", default=None, type=str,
                        required=True, help="The output file path")

    args = parser.parse_args()

    return args

def find_files_in_dir(in_dir, suffix):
    return [in_dir + '/' + f for f in os.listdir(in_dir) if os.path.isfile(in_dir + '/' + f) and f.endswith(suffix)]

def merge_files_by_column(in_files, out_file):
    df_merged = None
    for in_file in in_files:
        in_file_df = pandas.read_csv(in_file, delimiter='\t')

        if df_merged is None:
            df_merged = in_file_df
        else:
            df_merged = df_merged.astype({'gene': 'string'})
            in_file_df = in_file_df.astype({'gene' : 'string'})

            df_merged = df_merged.merge(in_file_df, how='outer', on='gene')

    df_merged = df_merged.fillna(0)
    df_merged.to_csv(out_file, sep='\t', index=False)

def main():
    """ Main subroutine """
    args = parse_args()

    in_files = find_files_in_dir(args.indir, args.suffix)

    merge_files_by_column(in_files, args.outfile)

if __name__ == '__main__':
    main()
