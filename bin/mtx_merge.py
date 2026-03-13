#!/usr/bin/env python3

import argparse
import os
import pandas as pd

LSUFFIX='_left'
RSUFFIX='_right'

def parse_args():
    """ Parse commandline arguments """
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i",
        "--in_dir",
        default=None,
        type=str,
        required=True,
        help="The input directory containing all matrices to merge"
    )

    parser.add_argument(
        "-x",
        "--in_extension",
        default=".mtx",
        type=str,
        required=True,
        help="The file extension for matrices to merge"
    )

    parser.add_argument(
        "-o",
        "--out_file",
        default="out.mtx",
        type=str,
        required=False,
        help="The name of the resulting matrix"
    )

    return parser.parse_args()

def get_mtx_list(in_dir, mtx_ext):
    mtx_list = []

    for mtx in os.listdir(os.fsencode(in_dir)):
        mtx_name = os.fsdecode(mtx)

        if mtx_name.endswith(mtx_ext):
            mtx_list.append('/'.join([in_dir,mtx_name]))

    return mtx_list

def merge_matrices(in_mtx):
    final_mtx = None
    for mtx in in_mtx:
        print(mtx)
        mtx_df = pd.read_csv(mtx, delimiter="\t", header=0, index_col=0).transpose()

        # This means the matrix is empty so we can skip it
        if len(mtx_df.columns) <= 1 and 'count' in mtx_df.columns:
            continue

        if final_mtx is None:
            final_mtx = mtx_df

        else:
            final_mtx = final_mtx.join(mtx_df, how='outer', lsuffix=LSUFFIX, rsuffix=RSUFFIX)

            # We do expect duplicate feature names, we just need to combine them
            if final_mtx.columns.str.contains(LSUFFIX).any():
                dupe_bc_cols = final_mtx.columns[final_mtx.columns.str.contains(LSUFFIX)].str.replace(LSUFFIX, '')

                # Iterate through all duplicated columns and sum them
                for dupe_bc_col in dupe_bc_cols:
                    bc_cols = final_mtx.columns[final_mtx.columns.str.contains(dupe_bc_col)]
                    final_mtx[dupe_bc_col] = final_mtx[bc_cols].sum(axis=1)
                    final_mtx = final_mtx.drop(columns = bc_cols)

    return final_mtx.transpose().fillna(value=0.0)

def main():
    args = parse_args()

    mtx_list = get_mtx_list(args.in_dir, args.in_extension)

    final_mtx = merge_matrices(mtx_list)
    final_mtx.to_csv(args.out_file, sep='\t')

if __name__ == '__main__':
    main()
