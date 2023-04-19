#!/usr/bin/env python3

""" Given a fastq file and a blaze output file, this will extract the barcode
    and umi and place them in the header of the fastq, as well as stripping
    them from teh read.
"""

import argparse
import regex as re
from Bio import SeqIO
from Bio.Seq import Seq

def parse_args():
    """ Parse the commandline arguments """

    arg_parser = argparse.ArgumentParser()

    arg_parser.add_argument('-i', '--input_file', required=True, type=str,
                             help="The input fastq file")
    arg_parser.add_argument('-b', '--barcode_file', required=True, type=str,
                            help="The file containing the readname and barcode")
    arg_parser.add_argument('-o', '--output_file', required=True, type=str,
                            help="The output fastq")

    args = arg_parser.parse_args()
    return args

def read_barcode_list(barcode_file):
    """ TODO
    """
    barcode_list = {}

    with open(barcode_file, 'r', encoding = "utf-8") as bc_in:
        for bc_line in bc_in.readlines():
            read_name, bc, bc_score = bc_line.split(',')

            # Not all reads had barcodes, ignore those that don't
            if bc:
                barcode_list[read_name] = bc

    return barcode_list


def extract_barcode(input_file, barcode_file, output):
    """ TODO
    """

    barcode_list = read_barcode_list(barcode_file)
    
    with open(f"{output}.R1.fastq", 'w', encoding = "utf-8") as r1_out, \
            open(f"{output}.R2.fastq", 'w', encoding = "utf-8") as r2_out:

        for record in SeqIO.parse(input_file, "fastq"):
            sequence = str(record.seq)
            qualities = ''.join([chr(score + 33) for score in
                                 record.letter_annotations['phred_quality']])

            if record.id in barcode_list:
                bc_index, sequence, qualities = find_seq_indices(barcode_list[record.id], sequence, qualities)
                
                # Strip the primer, bc, umi, and poly-T
                
                # BC_LENGTH = 16
                # UMI_LENGTH = 12
                # POLYT_LENGTH = 10
                print(sequence)
                print(qualities)
                print()

                bc = sequence[bc_index:bc_index+16]
                print(bc)
                bc_quals = qualities[bc_index:bc_index+16]
                print(bc_quals)

                umi = sequence[bc_index+16:bc_index+16+12]
                print(umi)
                umi_quals = qualities[bc_index+16:bc_index+16+12]
                print(umi_quals)
                
                no_bc = sequence[bc_index+16+12+10:]
                no_bc_quals = qualities[bc_index+16+12+10:]
                print(no_bc)
                print(no_bc_quals)

                r1_out.write(f"@{record.id}\n{bc}{umi}\n+\n{bc_quals}{umi_quals}\n")
                r2_out.write(f"@{record.id}\n{no_bc}\n+\n{no_bc_quals}\n")

def find_seq_indices(barcode, sequence, qualities):
    # See if the barcode is in the same direction as the read
    index = sequence.find(barcode)

    # If barcode not found, check the reverse complement
    if index <= 0:
        sequence = str(Seq(sequence).reverse_complement())
        qualities = qualities[::-1]

        index = sequence.find(barcode)

    return index, sequence, qualities

def main():
    """ Main subroutine """

    args = parse_args()
    extract_barcode(args.input_file, args.barcode_file, args.output_file)

if __name__ == '__main__':
    main()
