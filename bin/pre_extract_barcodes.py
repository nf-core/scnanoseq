#!/usr/bin/env python3

""" Given a fastq file, this script will find all sequences that match the
    regex passed in from the command line and output the sequences in cell
    ranger format, i.e. R1 contains the barcode and UMI, while R2 contains
    the sequence without the barcode and UMI.
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
    arg_parser.add_argument('-r', '--regex', required=True, type=str,
                            help="The regex for the barcode and umi")
    arg_parser.add_argument('-o', '--output_file', required=True, type=str,
                            help="The output prefix for the two fastqs")

    args = arg_parser.parse_args()
    return args

def pre_extract_fastq(input_file, regex, output_prefix):
    """ Using the regex, find all the sequences from the fastq that have the
        regex.

    Args:
        input_file (str): The input fastq file
        regex (str): The barcode and umi regex
        output_file (str): The output file prefix to output sequences to

    Returns: None
    """

    print(regex)

    with open(f"{output_prefix}.R1.fastq", 'w', encoding = "utf-8") as r1_out, \
        open(f"{output_prefix}.R2.fastq", 'w', encoding = "utf-8") as r2_out:

        for record in SeqIO.parse(input_file, "fastq"):

            sequence = str(record.seq)
            qualities = ''.join([chr(score + 33) for score in
                                 record.letter_annotations['phred_quality']])

            # Check if there's a match
            regex_match, sequence, qualities = find_regex_sequence(regex,
                                                                   sequence,
                                                                   qualities)

            if regex_match is not None:
                # Find the actual start location of the bc and umi

                matched_seq = sequence[regex_match.start():regex_match.end()]

                # Extract the barcode and umi from the sequence
                bc_umi, bc_umi_quals = get_bc_umi_info(matched_seq,
                                                       regex_match,
                                                       qualities)

                no_bc = sequence[regex_match.end() + 1:]
                no_bc_quals = qualities[regex_match.end() + 1:]

                # Output the barcode and umi into a file
                r1_out.write(f"@{record.id}\n{bc_umi}\n+\n{bc_umi_quals}\n")

                r2_out.write(f"@{record.id}\n{no_bc}\n+\n{no_bc_quals}\n")

def get_bc_umi_info(seq_substring, regex_match, seq_qualities):
    """ Extract the barcode and umi (and their qualities) from the base
        sequence """

    bc_umi = ''
    bc_umi_quals = ''

    for group_name in regex_match.groupdict():
        if 'discard' not in group_name:
            seq = regex_match[group_name]

            start_idx = seq_substring.find(seq)
            end_idx = start_idx + len(seq)

            bc_umi += seq
            bc_umi_quals += seq_qualities[start_idx:end_idx]

        else:
            pass

    return bc_umi, bc_umi_quals


def find_regex_sequence(regex, sequence, qualities):
    """ Find the regex in the sequence """

    reg_ex = re.compile(regex)
    r_match = reg_ex.search(sequence)

    # If there is no match, check the reverse complement of the sequence
    if r_match is None:
        sequence = str(Seq(sequence).reverse_complement())
        qualities = qualities[::-1] # Also reverse the quality string

        r_match = reg_ex.search(sequence)

    return r_match, sequence, qualities


def main():
    """ Main subroutine """

    args = parse_args()
    pre_extract_fastq(args.input_file, args.regex, args.output_file)

if __name__ == '__main__':
    main()
