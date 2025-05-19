#!/usr/bin/env python3

""" This script will iterate over a bam and tag the reads with the barcode,
    barcode quality, umi, umi quality that are obtained from the read1 fastq
"""

import argparse
import pysam

BC_TAG = "CR"
BC_QUAL_TAG = "CY"
CORRECTED_TAG = "CB"

UMI_TAG = "UR"
UMI_QUAL_TAG = "UY"


class UmiBcRead:
    """This class holds and parses out the barcode and umi from a fastq
    read"""

    def __init__(self, read_name, bc, bc_qual, umi, umi_qual, corrected_bc):
        self.read_name = read_name
        self.bc = bc
        self.bc_qual = bc_qual
        self.umi = umi
        self.umi_qual = umi_qual
        self.corrected_bc = corrected_bc

def parse_args():
    """Parse the commandline arguments"""

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-b",
        "--in_bam",
        default=None,
        type=str,
        required=True,
        help="The input bam file"
    )

    parser.add_argument(
        "-i",
        "--in_bc_info",
        default=None,
        type=str,
        required=True,
        help="The input tsv containing bc info"
    )

    parser.add_argument(
        "-o",
        "--out_bam",
        default=None,
        type=str,
        required=True,
        help="The output bam file"
    )

    args = parser.parse_args()

    return args


def tag_bams(in_bam, bc_info, out_bam):
    """This will add tags to the read based on the information from the R1
        fastq

    Args:
        in_bam (str): The input bam
        in_fastq (str): The input fastq containing only the barcoding
        out_bam (str): The output bam containing the various tags

    Returns: None
    """

    with pysam.AlignmentFile(in_bam, "rb") as fh_in_bam, pysam.AlignmentFile(
        out_bam, "wb", template=fh_in_bam
    ) as fh_out_bam:

        umi_bc_infos = read_bc_info(bc_info)

        for read in fh_in_bam:
            parsed_read_name = read.query_name.split("_")[0]

            if parsed_read_name in umi_bc_infos:
                umi_bc_info = umi_bc_infos[parsed_read_name]

                # Add the barcode
                if not read.has_tag(BC_TAG):
                    read.tags += [(BC_TAG, umi_bc_info.bc)]

                # Add the barcode quality
                if not read.has_tag(BC_QUAL_TAG):
                    read.tags += [(BC_QUAL_TAG, umi_bc_info.bc_qual)]

                # Add the umi
                if not read.has_tag(UMI_TAG):
                    read.tags += [(UMI_TAG, umi_bc_info.umi)]

                # Add the umi quality
                if not read.has_tag(UMI_QUAL_TAG):
                    read.tags += [(UMI_QUAL_TAG, umi_bc_info.umi_qual)]

                # Add the corrected_bc tag
                if not read.has_tag(CORRECTED_TAG):
                    read.tags += [(CORRECTED_TAG, umi_bc_info.corrected_bc)]

                read.query_name = "_".join([parsed_read_name, umi_bc_info.bc, umi_bc_info.umi])

                fh_out_bam.write(read)


def read_bc_info(bc_info):
    """This will read in the input barcode information

    Args:
        bc_info (str): The input barcode info containing sequence and qualities for a reads
            barcode and umi

    Returns:
        bc_info_dict (dict): The dictionary that contains the quick barcode and umi information
            for each read. The key is the sequence name while the value is the class created from
            that fastq
    """

    bc_info_dict = {}

    with open(bc_info) as fh_bc_info:
        for entry in fh_bc_info:
            # Sequence names can contain other information and this other
            #   information can use the underscore as a delimiter
            entry = entry.strip('\n')
            read_id, bc, bc_info, umi, umi_qual, corrected_bc = entry.split('\t')

            if corrected_bc:
                bc_info_dict[read_id] = UmiBcRead(
                    read_id, bc, bc_info, umi, umi_qual, corrected_bc
                )

    return bc_info_dict


def main():
    """Main subroutine"""
    args = parse_args()
    tag_bams(args.in_bam, args.in_bc_info, args.out_bam)


if __name__ == "__main__":
    main()
