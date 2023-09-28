#!/usr/bin/env python3

""" Given a fastq file and a blaze output file, this will extract the barcode
    and umi and place them in the header of the fastq, as well as stripping
    them from teh read.
"""

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
import gzip

def parse_args():
    """ Parse the commandline arguments """

    arg_parser = argparse.ArgumentParser()

    arg_parser.add_argument('-i', '--input_file', required=True, type=str,
                             help="The input fastq file")
    arg_parser.add_argument('-b', '--barcode_file', required=True, type=str,
                            help="The file containing the readname and barcode")
    arg_parser.add_argument('-o', '--output_file', required=True, type=str,
                            help="The output fastq")
    arg_parser.add_argument('-f', '--barcode-format', required=False, type=str,
                            help="The barcode/umi format (Options: cellranger)")

    args = arg_parser.parse_args()
    return args

def read_barcode_list(barcode_file):
    """ Read in the blaze putative barcode list into a dict in memory. The
        file contains the read name, predicted barcode, and average bc score

    Args:
        barcode_file (str): The path to the blaze putative barcode list

    Returns:
        barcode_list (dict): The dictionary containing the barcodes. Key is
            read identifier and value is the barcode
    """
    barcode_list = {}

    with open(barcode_file, 'r', encoding = "utf-8") as bc_in:
        for bc_line in bc_in.readlines():
            read_name, barcode, _ = bc_line.split(',')

            # Not all reads had barcodes, ignore those that don't
            if barcode:
                barcode_list[read_name] = barcode

    return barcode_list


def extract_barcode(input_file, barcode_file, output, bc_format):
    """ This will use the information stored in the barcode file to extract
        the predicted barcode and umi out of the main read (along with
        qualities) and place them in a separate fastq. This somewhat mimics
        the way that 10x provides single cell files (minus the index files)

    Args:
        input_file (str): The path to the input fastq
        barcode_file (str): The path to the blaze putative barcode list
        output (str): The output prefix
    """

    barcode_list = read_barcode_list(barcode_file)

    # We will open two fastqs to output to.
    # R1 will contain the barcode + umi
    # R2 contains the actual read

    with gzip.open(f"{output}.R1.fastq.gz", 'wt') as r1_out, \
            gzip.open(f"{output}.R2.fastq.gz", 'wt') as r2_out, \
            gzip.open(input_file, 'rt') as fastq_in:

        for record in SeqIO.parse(fastq_in, "fastq"):
            orig_seq = str(record.seq)
            orig_quals = ''.join([chr(score + 33) for score in
                                 record.letter_annotations['phred_quality']])

            # NOTE: Reads that do not have a predicted barcode will be filtered
            #   out at this point
            if record.id in barcode_list:
                bc_index, seq, quals = find_seq_indices(barcode_list[record.id],
                                                        orig_seq,
                                                        orig_quals)

                # If bc_index is < 0, the barcode was not found
                if bc_index >= 0:
                    read_info = {}

                    # Strip the primer, bc, umi, and poly-T
                    if bc_format == "cellranger":
                        read_info = strip_read_cellranger(bc_index, seq, quals)

                    if read_info:
                        r1_out.write('\n'.join(["@" + record.id,
                                                read_info["r1_read"],
                                                "+",
                                                read_info["r1_qual"],
                                                '']))

                        r2_out.write('\n'.join(["@" + record.id,
                                                read_info["r2_read"],
                                                "+",
                                                read_info["r2_qual"],
                                                '']))

def find_seq_indices(barcode, sequence, qualities):
    """ Find the location in the read where the predictoed barcode exists. If
        it cannot be found, reverse-complement the read to find i.

    Args:
        barcode (str): The predicted barcode for the read
        sequence (str): The original read from the fastq
        qualities (str): The original qualities from the fastq

    Return:
        index (int): The start index in the read where the barcode was found
        sequence (str): The sequence. We return this in case we had to
            reverse-complement it to find the barcode
        qualities (str): The qualities of the sequence. We return this in case
            we had to reverse-complement it to find the barcode.

    """
    # See if the barcode is in the same direction as the read
    index = sequence.find(barcode)

    # If barcode not found, check the reverse complement
    if index < 0:
        sequence = str(Seq(sequence).reverse_complement())
        qualities = qualities[::-1]

        index = sequence.find(barcode)

    return index, sequence, qualities

def strip_read_cellranger(bc_index, seq, quals):
    """ Strip the bc and umi from a read, and convert it from a single read
        format to paired read. This function is used for when the barcode is in
        the 10X format, so we expect that the read would look this:
        {bc}{umi}{polyT}{read}
        where bc is 16 base pairs, umi is 12 base pairs, and the poly T is 10
        base pairs

    Args:
        bc_index (int): The location of the barcode in the read
        seq (str): The full sequence
        quals (str): The full sequence qualities

    Returns:
        read_info (dict): The dictionary containing the read in paired format.
            The keys are listed below:
            * r1_read (str): The barcode and umi sequences
            * r1_qual (str): The qualities of the barcode and umi sequences
            * r2_read (str): The read without the barcode, umi, or polyT
            * r2_qual (str): The quality of the read

    """
    read_info = {}

    # TODO: Need to store these somewhere like log_file rather than be hardcoded
    bc_length = 16
    umi_length = 12
    polyt_length = 10

    read_info["r1_read"] = seq[bc_index:bc_index+bc_length+umi_length]
    read_info["r1_qual"] = quals[bc_index:bc_index+bc_length+umi_length]

    read_info["r2_read"] = seq[bc_index+bc_length+umi_length+polyt_length:]
    read_info["r2_qual"] = quals[bc_index+bc_length+umi_length+polyt_length:]

    return read_info

def main():
    """ Main subroutine """

    args = parse_args()
    extract_barcode(args.input_file,
                    args.barcode_file,
                    args.output_file,
                    args.barcode_format)

if __name__ == '__main__':
    main()
