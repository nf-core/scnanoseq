#!/usr/bin/env python3

""" This file will reformat the whitelist produced by umi_tools. """

import argparse

def parse_args():
    """ Parse the commandline arguments """

    parser = argparse.ArgumentParser()

    parser.add_argument("--infile", "-i", default=None, type=str,
                        required=True, help="The input whitelist")
    parser.add_argument("--outfile", "-o", default=None, type=str,
                        required=True, help="The output file")

    args = parser.parse_args()
    return args

def reformat_whitelist(infile, outfile):
    """ This will take in the whitelist that was created by umi_tools and
        reformat it into a singular column

    Args:
        infile (str): The input file
        outfile (str): The output file

    Returns: None
    """
    outlines = get_outlines(infile)
    print(outlines)
    write_outfile(outfile, outlines)

def get_outlines(infile):
    """ This will read the infile and reformat it so that the alternative
        barcodes are their own rows

    Args:
        infile (str): The input file

    Returns:
        outlines (list): The list of reformatted lines
    """
    outlines = {'primary': [], 'secondary': []}

    with open(infile, 'r', encoding = "utf-8") as f_in:
        for line in f_in.readlines():
            line = line.strip('\n')

            prim_bc, alt_bcs, prim_bc_count, alt_bc_counts = line.split('\t')

            if prim_bc:
                outlines['primary'].append((prim_bc, '', prim_bc_count, ''))

            if alt_bcs:
                for alt_bc, alt_count in zip(alt_bcs.split(','),
                                     alt_bc_counts.split(',')):
                    outlines['secondary'].append((alt_bc, '', alt_count, ''))

    return outlines

def write_outfile(outfile, outlines):
    """ Write out the outlines to the output file """

    with open(outfile, 'w', encoding = "utf-8") as f_out_all_bcs, \
        open(outfile + '.bc_only', 'w', encoding = "utf-8") as f_out_main_bcs:
        for key, value in outlines.items():
            for barcode_info in value:
                f_out_main_bcs.write(barcode_info[0] + '\n')
                f_out_all_bcs.write('\t'.join(barcode_info) + '\n')

def main():
    """ Main Subroutine """
    args = parse_args()

    reformat_whitelist(args.infile, args.outfile)

if __name__ == '__main__':
    main()
