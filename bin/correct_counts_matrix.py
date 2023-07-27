#!/usr/bin/env python3

"""
This script acts to split up and correct feature counts file that used multi
overlap settings
"""

import argparse
import pandas
import cProfile
import timeit

def parse_args():
    """ Argument Parser """
    arg_parser = argparse.ArgumentParser()

    arg_parser.add_argument('-i', '--input_file')
    arg_parser.add_argument('-o', '--output_file')

    return arg_parser.parse_args()

def split_overlap_features(args):
    """ Will split the overlapping features """
    # Read in file
    print("Reading in file...")
    in_features = pandas.read_csv(args.input_file, delimiter='\t')

    # Split file into two tables, multi_feature_rows and solo_feature_rows
    print("Splitting file into two separate matrices...")
    multi_feature_rows, solo_feature_rows = get_matrices(in_features)

    # Iterate over multi_feature_rows, and add their counts to the individual
    # solo_feature_rows
    print("Iterating over the multiple features...")
    for multi_feature_row in multi_feature_rows.itertuples():
        solo_feature_rows = process_features(multi_feature_row, solo_feature_rows)

    # Output file
    solo_feature_rows.to_csv(args.output_file, sep='\t', index=False)

def get_matrices(in_features):
    multi_feature_rows = in_features[in_features.gene.str.contains(',')]
    multi_feature_rows = multi_feature_rows[~multi_feature_rows.loc[:, multi_feature_rows.columns != 'gene'].eq(0, axis=0).all(1)]
    empty_multi_features = multi_feature_rows[multi_feature_rows.loc[:, multi_feature_rows.columns != 'gene'].eq(0, axis=0).all(1)]

    solo_feature_rows = in_features[~in_features.gene.str.contains(',')].reset_index(drop = True)
    solo_feature_rows = merge_features(empty_multi_features.gene.tolist(), solo_feature_rows)

    return multi_feature_rows, solo_feature_rows

def merge_features(gene_features, solo_feature_rows):
    for gene_feature in gene_features:
        for gene in gene_feature.split(','):
            if not solo_feature_rows['gene'].eq(gene).any():
                new_row = [gene].append([0] * (len(solo_feature_rows.columns)-1))
                solo_feature_rows.loc[len(df)] = new_row
    return solo_feature_rows


def process_features(multi_feature_row, solo_feature_rows):
    # Split the features
    features = multi_feature_row.gene.split(',')

    for feature in features:
        # Iterate over the fields in multi_feature_row
        gene_row = solo_feature_rows[solo_feature_rows.gene == feature]

        # Gene row should only have 0 or 1 rows
        if not len(gene_row):
            # This means we don't have any hits, so we need to add it
            gene_row_idx = len(solo_feature_rows)

            solo_feature_rows.loc[gene_row_idx,:] = 0
            solo_feature_rows.loc[gene_row_idx, 'gene'] = feature

        elif len(gene_row) == 1:
            # We have a single hit, lets set the gene_index to the
            # correct record
            gene_row_idx = int(gene_row.index.values)

        elif len(gene_row) > 1:
            # We had multiple hits, this should not happen as each
            # gene should only have one row dedicated to it, so
            # something must've gone wrong
            print("Too many matches found!")

        solo_feature_rows = update_counts(multi_feature_row, solo_feature_rows, features, gene_row_idx)

    return solo_feature_rows

def update_counts(multi_feature_row, solo_feature_rows, features, gene_row_idx):
    for field in multi_feature_row._fields:
        # We want to skip gene because we have that, and Index because
        # that's namedtuple specific
        if field == 'gene' or field == 'Index':
            continue
        else:
            # Lets get the count for the feature. Per FeatureCounts
            # documentation each feature on the row adds one to the
            # count, e.g. if a row has 3 features, the count will be 3
            solo_feature_rows.loc[gene_row_idx,field] += int(getattr(multi_feature_row, field)) / len(features)
    return solo_feature_rows

def main():
    """ Main subroutine """
    args = parse_args()
    split_overlap_features(args)

if __name__ == '__main__':
    profiler = cProfile.Profile()
    profiler.enable()

    main()

    profiler.disable()
    profiler.dump_stats("run.stats")
