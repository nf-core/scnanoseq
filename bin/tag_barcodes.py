#!/usr/bin/env python3

""" This script will iterate over a bam and tag the reads with the barcode,
    barcode quality, umi, umi quality that are obtained from the read1 fastq
"""

import argparse
import pysam

BC_TAG = 'CR'
BC_QUAL_TAG = 'CY'

UMI_TAG = 'UR'
UMI_QUAL_TAG = 'UY'

class UmiBcRead:
    """ This class holds and parses out the barcode and umi from a fastq
        read """
    def __init__(self, read_name, sequence, qualities, bc_length, umi_length):
        self.read_name = read_name
        self.sequence = sequence
        self.qualities = qualities
        self.bc_length = bc_length
        self.umi_length = umi_length

    def get_bc(self):
        """ Gets the barcode out of the read sequence """
        return self.sequence[0:self.bc_length]

    def get_umi(self):
        """ Gets the umi out of the read name sequence"""
        return self.sequence[self.bc_length:]

    def get_bc_qual(self):
        """ Gets the barcode qualities """
        return self.get_qual(self.get_bc())

    def get_umi_qual(self):
        """ Gets the umi qualities """
        return self.get_qual(self.get_umi())

    def get_qual(self, src_seq):
        """ Generic function that will get the corresponding quality for
            a given sequence """

        # This finds where the input sequence starts
        strt_idx = self.sequence.index(src_seq)

        # We use the start index gathered from above to get the
        #   correct quality string we need
        quality = self.qualities[strt_idx:(strt_idx+len(src_seq))]

        return quality

def parse_args():
    """ Parse the commandline arguments """

    parser = argparse.ArgumentParser()

    parser.add_argument("-b", "--in_bam", default=None, type=str,
                        required=True, help="The input bam file")

    parser.add_argument("-f", "--in_fastq", default=None, type=str,
                        required=True, help="The input R1 fastq file")

    parser.add_argument("-o", "--out_bam", default=None, type=str,
                        required=True, help="The output bam file")

    parser.add_argument("-bl", "--barcode_length", default=None, type=int,
                        required=True, help="The length of the barcode")

    parser.add_argument("-ul", "--umi_length", default=None, type=int,
                        required=True, help="The length of the umi")

    args = parser.parse_args()

    return args

def tag_bams(in_bam, in_fastq, out_bam, bc_length, umi_length):
    """ This will add tags to the read based on the information from the R1
        fastq

    Args:
        in_bam (str): The input bam
        in_fastq (str): The input fastq containing only the barcoding
        out_bam (str): The output bam containing the various tags
        bc_length (int): The length of the barcode
        umi_length (int): The length of the umi

    Returns: None
    """

    with pysam.AlignmentFile(in_bam, "rb") as fh_in_bam, \
        pysam.AlignmentFile(out_bam, "wb", template=fh_in_bam) as fh_out_bam:

        # Read the R1 fastq into a dictionary
        umi_bc_infos = read_r1_fastq_info(in_fastq, bc_length, umi_length)

        for read in fh_in_bam:
            parsed_read_name = read.query_name.split('_')[0]

            umi_bc_info = umi_bc_infos[parsed_read_name]

            # Add the barcode
            if not read.has_tag(BC_TAG):
                read.tags += [(BC_TAG, umi_bc_info.get_bc())]

            # Add the barcode quality
            if not read.has_tag(BC_QUAL_TAG):
                read.tags += [(BC_QUAL_TAG, umi_bc_info.get_bc_qual())]

            # Add the umi
            if not read.has_tag(UMI_TAG):
                read.tags += [(UMI_TAG, umi_bc_info.get_umi())]

            # Add the umi quality
            if not read.has_tag(UMI_QUAL_TAG):
                read.tags += [(UMI_QUAL_TAG, umi_bc_info.get_umi_qual())]

            read.query_name = '_'.join([parsed_read_name,
                                        umi_bc_info.get_bc(),
                                        umi_bc_info.get_umi()])

            fh_out_bam.write(read)

def read_r1_fastq_info(in_fastq, bc_length, umi_length):
    """ This will read in the input fastq and return the information as a
        list

    Args:
        in_fastq (str): The input fastq containing only the barcode sequences
            and qualities
        bc_length (int): The length of the barcode
        umi_length (int): The lenght of the umi

    Returns:
        r1_fastq_info (dict): The dictionary that contains the quick barcode
            and umi information for each read. The key is the sequence name
            while the value is the class created from that fastq
    """

    r1_fastq_info = {}

    with pysam.FastxFile(in_fastq) as fh_in_fastq:
        for entry in fh_in_fastq:
            # Sequence names can contain other information and this other
            #   information can use the underscore as a delimiter
            parsed_read_name = entry.name.split('_')[0]

            r1_fastq_info[parsed_read_name] = UmiBcRead(entry.name,
                                                        entry.sequence,
                                                        entry.quality,
                                                        bc_length,
                                                        umi_length)

    return r1_fastq_info

def main():
    """ Main subroutine """
    args = parse_args()
    tag_bams(args.in_bam, args.in_fastq, args.out_bam, args.barcode_length,
             args.umi_length)

if __name__ == '__main__':
    main()
