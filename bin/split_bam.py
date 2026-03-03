#!/usr/bin/env python3

import argparse
import pysam

def parse_args():
    """ Parse commandline arguments"""

    parser = argparse.ArgumentParser(description="Split a BAM file into multiple BAM files based on read groups.")

    parser.add_argument(
        "-i",
        "--input",
        default=None,
        type=str,
        required=True,
        help="Input BAM file"
    )

    parser.add_argument(
        "-o",
        "--output",
        default=None,
        type=str,
        required=True,
        help="Output BAM file"
    )

    parser.add_argument(
        "-c",
        "--contig_file",
        default=None,
        type=str,
        required=True,
        help="File containing contig names to include in the output BAM file."
    )

    return parser.parse_args()

def split_bam_by_contig(input_bam, output_bam, contigs_file):
    """ Split a BAM file into multiple BAM files based on contig names.

    Args:
        input_bam (str): Path to the input BAM file.
        output_bam (str): Path to the output BAM file.
        contigs_file (str): Path to the file containing contig names to include in the output BAM file.

    Returns:
        None
    """

    # Read the contig names from the contigs file
    contigs_to_include = set()
    with open(contigs_file, 'r') as f:
        contigs_to_include = set(line.strip() for line in f)

    # Open the input BAM file
    bam_in = pysam.AlignmentFile(input_bam, "rb")

    # Create a new BAM file for output
    bam_out = pysam.AlignmentFile(output_bam, "wb", template=bam_in)

    # Iterate through each read in the input BAM file
    for read in bam_in.fetch():
        # Check if the read's reference name is in the set of contigs to include
        if read.reference_name in contigs_to_include:
            bam_out.write(read)

    # Close the BAM files
    bam_in.close()
    bam_out.close()

def main():
    args = parse_args()

    split_bam_by_contig(args.input, args.output, args.contig_file)

if __name__ == "__main__":
    main()

