#!/usr/bin/env python3

import sys
import pysam
import argparse

def parse_args():
    """ Read in command line arguments """
    parser = argparse.ArgumentParser()

    parser.add_argument("--infile", "-i", default=None, type=str,
                        help='nanopore infile bam file')
    parser.add_argument("--featurefile", "-f", default=None, type=str,
                        help='name for feature CORE file')

    parser.add_argument("--outfile", "-o", default=None, type=str,
                        help='name for output file')
    
    args = parser.parse_args()
    return args

def tag_bam_file(input_file, feature_file, output_file):
    """ Produce a bam file with the tags produced by featureCounts """
    samfile = pysam.AlignmentFile(input_file, "rb")
    featurein = open(feature_file, 'rt')
    outfile = pysam.AlignmentFile(output_file, "wb", template=samfile)


    for read in samfile:
        read_name, read_status, feature_count, gene_assignment = featurein.readline().strip('\n').split('\t')

        # The orders should match, on the off chance they don't, lets error out
        if read_name != read.query_name:
            print('{r1} does not match {r2}'.format(r1=read_name, r2=read.query_name))
            sys.exit(1)
        
        read.tags += [('XS', read_status)]

        # This will mimic the tagging that feature counts would do
        if not 'Unassigned' in read_status:
            read.tags += [('XN', feature_count)]
            read.tags += [('XT', gene_assignment)]
    
        outfile.write(read)
    
    samfile.close()
    featurein.close()
    outfile.close()

def main():
    ''' Main Subroutine '''
    args = parse_args()

    tag_bam_file(args.infile, args.featurefile, args.outfile)

if __name__ == '__main__':
    main()
