#!/usr/bin/env python3

""" Given a fastq file and a blaze output file, this will extract the barcode
    and umi and place them in the header of the fastq, as well as stripping
    them from teh read.
"""

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
import gzip

import threading
import queue

OUTPUT_LOCK = threading.Lock()

def parse_args():
    """Parse the commandline arguments"""

    arg_parser = argparse.ArgumentParser()

    arg_parser.add_argument("-i", "--input_file", required=True, type=str, help="The input fastq file")
    arg_parser.add_argument(
        "-b", "--barcode_file", required=True, type=str, help="The file containing the readname and barcode"
    )
    arg_parser.add_argument("-o", "--output_file", required=True, type=str, help="The output fastq")
    arg_parser.add_argument(
        "-f", "--barcode_format", required=True, type=str, help="The barcode/umi format (Options: cellranger)"
    )
    arg_parser.add_argument(
        "-t", "--threads", type=int, help="The number of threads to use")

    args = arg_parser.parse_args()
    return args

def extract_bc_umi(bc_format, bc_queue, bc_out, r2_out):
    """ Will extract the barcode and umi from a read and write them out to a file """
    while True:
        read_id, barcode, orig_seq, orig_quals = bc_queue.get()
        #print(read_id)
        bc_index, seq, quals = find_seq_indices(barcode, orig_seq, orig_quals)

        # If bc_index is < 0, the barcode was not found
        if bc_index >= 0:
            read_info = {}

            # Strip the primer, bc, umi, and poly-T
            if bc_format in ("10X_3v3", "10X_3v4"):
                read_info = strip_read_10X_3v3_3v4(bc_index, seq, quals)

            elif bc_format in ("10X_5v2"):
                read_info = strip_read_10X_5v2(bc_index, seq, quals)

            elif bc_format in ("10X_5v3"):
                read_info = strip_read_10X_5v3(bc_index, seq, quals)

            if read_info:
                bc_output = "\t".join(
                    [read_id,
                    read_info["bc"],
                    read_info["bc_qual"],
                    read_info["umi"],
                    read_info["umi_qual"]]
                    ) + "\n"

                fq_output = "\n".join(
                    ["@" + read_id,
                    read_info["r2_read"],
                    "+",
                    read_info["r2_qual"],
                    ""])

                with OUTPUT_LOCK:
                    bc_out.write(bc_output)
                    r2_out.write(fq_output)

        bc_queue.task_done()

def extract_barcode(input_file, barcode_file, output, bc_format, threads):
    """ Reads in a fastq and BLAZE putative bc file and strip the bc and umi from read """
    bc_queue = queue.Queue(5000)

    bc_out = open(f"{output}.putative_bc_umi.tsv", "w", encoding="utf-8")
    bc_out.write("read_id\tbc\tbc_qual\tumi\tumi_qual\n")

    r2_out = open(f"{output}.fastq", "w", encoding="utf-8")

    # start worker threads
    for i in range(threads):
        t = threading.Thread(target=extract_bc_umi, args=(bc_format, bc_queue, bc_out, r2_out))
        t.daemon = True
        t.start()

    # read file
    fastq_in = gzip.open(input_file, "rt") if '.gz' in input_file else open(input_file, "r", encoding="utf-8")

    with open(barcode_file, "r", encoding="utf-8") as bc_in:

        for bc_umi_line,record in zip(bc_in,SeqIO.parse(fastq_in, "fastq")):
            bc_umi_line = bc_umi_line.strip()

            _, barcode, _, _, _, _, _ = bc_umi_line.split(",")
            orig_seq = str(record.seq)
            orig_quals = "".join([chr(score + 33) for score in record.letter_annotations["phred_quality"]])
            if barcode:
                bc_queue.put((record.id, barcode, orig_seq, orig_quals))

    fastq_in.close()

    bc_queue.join()

    bc_out.close()
    r2_out.close()

def find_seq_indices(barcode, sequence, qualities):
    """Find the location in the read where the predictoed barcode exists. If
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

def strip_read_10X_3v3_3v4(bc_index, seq, quals):
    bc_length = 16
    umi_length = 12
    polyt_length = 10

    return strip_read_10X(bc_index, seq, quals, bc_length, umi_length, polyt_length)

def strip_read_10X_5v3(bc_index, seq, quals):
    bc_length = 16
    umi_length = 12
    polyt_length = 10

    return strip_read_10X(bc_index, seq, quals, bc_length, umi_length, polyt_length)

def strip_read_10X_5v2(bc_index, seq, quals):
    bc_length = 16
    umi_length = 10
    polyt_length = 10

    return strip_read_10X(bc_index, seq, quals, bc_length, umi_length, polyt_length)

def strip_read_10X(bc_index, seq, quals, bc_length, umi_length, polyt_length):
    """Strip the bc and umi from a read, and convert it from a single read
        format to paired read. This function is used for when the barcode is in
        the 10X format, so we expect that the read would look this:
        {bc}{umi}{polyT}{read}

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

    read_info["bc"] = seq[bc_index : bc_index + bc_length ]
    read_info["bc_qual"] = quals[bc_index: bc_index + bc_length ]
    read_info["umi"] = seq[bc_index + bc_length : bc_index + bc_length + umi_length ]
    read_info["umi_qual"] = quals[bc_index + bc_length : bc_index + bc_length + umi_length ]

    read_info["r2_read"] = seq[bc_index + bc_length + umi_length + polyt_length :]
    read_info["r2_qual"] = quals[bc_index + bc_length + umi_length + polyt_length :]

    return read_info

def main():
    """Main subroutine"""

    args = parse_args()
    extract_barcode(args.input_file, args.barcode_file, args.output_file, args.barcode_format, args.threads)


if __name__ == "__main__":
    main()
