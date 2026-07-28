#!/usr/bin/env python3

""" Given a flexiplex-demultiplexed fastq, split the barcode and umi back out
    of the read name into a separate synthetic fastq. lr-kallisto reads the
    barcode and umi positionally (see the -x string), so they need to live in a
    read of their own rather than in the read name.

    Flexiplex writes read names as:
        @{barcode}_{umi}#{original_read_id}_{strand}{n}of{m}
    with the adapter, barcode, umi and poly-T already stripped from the
    sequence.
"""

import argparse
import gzip
import sys

# Flexiplex does not keep the per-base barcode/umi qualities, so a constant
# high quality is used. bustools only uses the barcode sequence when
# correcting against the on-list.
DUMMY_QUAL = "I"


def parse_args():
    """Parse the commandline arguments"""

    arg_parser = argparse.ArgumentParser()

    arg_parser.add_argument("-i", "--input_file", required=True, type=str, help="The flexiplex fastq file")
    arg_parser.add_argument("-o", "--output_prefix", required=True, type=str, help="Prefix for the output files")
    arg_parser.add_argument("-b", "--bc_length", required=True, type=int, help="The expected barcode length")
    arg_parser.add_argument("-u", "--umi_length", required=True, type=int, help="The expected umi length")

    return arg_parser.parse_args()


def open_fastq(path):
    """Open a fastq for reading, transparently handling gzip"""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r", encoding="utf-8")


def parse_read_name(header):
    """Pull the barcode, umi and original read id out of a flexiplex read name

    Args:
        header (str): The fastq header line, including the leading '@'

    Returns:
        tuple: (barcode, umi, read_id), or None if the name is not in the
            flexiplex format
    """
    # Drop the leading '@' and anything after the first whitespace. Flexiplex
    # appends the CB:Z:/UB:Z: tags as a comment, which we do not need.
    name = header[1:].split(None, 1)[0]

    bc_umi, sep, read_id = name.partition("#")
    if not sep:
        return None

    barcode, sep, umi = bc_umi.partition("_")
    if not sep:
        return None

    return barcode, umi, read_id


def split_bc_umi(input_file, output_prefix, bc_length, umi_length):
    """Split a flexiplex fastq into a barcode+umi fastq and a cdna fastq

    The two fastqs are written in lockstep from a single pass, which is what
    kallisto requires: it reads the files positionally and never matches read
    names between them.
    """
    written = 0
    skipped_unparsed = 0
    skipped_length = 0
    skipped_empty = 0

    # compresslevel=1 keeps this io-bound rather than cpu-bound; the files are
    # intermediates that only kallisto reads.
    with open_fastq(input_file) as fastq_in, gzip.open(
        f"{output_prefix}.bc_umi.fastq.gz", "wt", compresslevel=1
    ) as bc_out, gzip.open(f"{output_prefix}.cdna.fastq.gz", "wt", compresslevel=1) as cdna_out:

        while True:
            header = fastq_in.readline()
            if not header:
                break
            seq = fastq_in.readline().rstrip("\n")
            fastq_in.readline()
            qual = fastq_in.readline().rstrip("\n")

            parsed = parse_read_name(header.rstrip("\n"))
            if parsed is None:
                skipped_unparsed += 1
                continue

            barcode, umi, read_id = parsed

            if len(barcode) != bc_length or len(umi) != umi_length:
                skipped_length += 1
                continue

            if not seq:
                skipped_empty += 1
                continue

            bc_umi = barcode + umi
            bc_out.write(f"@{read_id}\n{bc_umi}\n+\n{DUMMY_QUAL * len(bc_umi)}\n")
            cdna_out.write(f"@{read_id}\n{seq}\n+\n{qual}\n")
            written += 1

    # seurat_qc.R only needs a total read count out of a flagstat file, and the
    # lr-kallisto path never produces a bam. Emit the count in flagstat format
    # so the QC_SCRNA subworkflow can be reused as-is.
    with open(f"{output_prefix}.flagstat", "w", encoding="utf-8") as flagstat_out:
        flagstat_out.write(f"{written} + 0 in total (QC-passed reads + QC-failed reads)\n")

    print(
        f"wrote {written} reads; skipped {skipped_unparsed} (unparseable read name), "
        f"{skipped_length} (unexpected barcode/umi length), {skipped_empty} (empty sequence)",
        file=sys.stderr,
    )

    if written == 0:
        sys.exit("No reads could be split. Are these reads demultiplexed by flexiplex?")


def main():
    """Main subroutine"""

    args = parse_args()
    split_bc_umi(args.input_file, args.output_prefix, args.bc_length, args.umi_length)


if __name__ == "__main__":
    main()
