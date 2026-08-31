#!/usr/bin/env python3

""" Build a kallisto/bustools transcript-to-gene mapping from a GTF.

    The mapping is a headerless tsv of transcript_id, gene_id, gene_name, in
    the order the transcripts appear in the GTF. Column 1 must match the
    sequence names in the kallisto index, which is why the same GTF that
    gffread extracted the transcript sequences from has to be used here.
"""

import argparse
import gzip
import re
import sys

ATTRIBUTE_RE = re.compile(r'(\S+)\s+"([^"]*)"')


def parse_args():
    """Parse the commandline arguments"""

    arg_parser = argparse.ArgumentParser()

    arg_parser.add_argument("-i", "--input_file", required=True, type=str, help="The input gtf")
    arg_parser.add_argument("-o", "--output_file", required=True, type=str, help="The output t2g tsv")

    return arg_parser.parse_args()


def open_gtf(path):
    """Open a gtf for reading, transparently handling gzip"""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r", encoding="utf-8")


def parse_attributes(field):
    """Parse the key "value" pairs of a GTF attribute column into a dict"""
    return dict(ATTRIBUTE_RE.findall(field))


def gtf_to_t2g(input_file, output_file):
    """Write a transcript-to-gene mapping for every transcript in the gtf"""
    seen = set()
    written = 0

    with open_gtf(input_file) as gtf_in, open(output_file, "w", encoding="utf-8") as t2g_out:
        for line in gtf_in:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "transcript":
                continue

            attributes = parse_attributes(fields[8])
            transcript_id = attributes.get("transcript_id")
            gene_id = attributes.get("gene_id")

            if not transcript_id or not gene_id:
                continue

            if transcript_id in seen:
                continue
            seen.add(transcript_id)

            gene_name = attributes.get("gene_name", gene_id)
            t2g_out.write(f"{transcript_id}\t{gene_id}\t{gene_name}\n")
            written += 1

    print(f"wrote {written} transcript-to-gene entries", file=sys.stderr)

    if written == 0:
        sys.exit("No transcripts found in the gtf. Does it contain 'transcript' feature rows?")


def main():
    """Main subroutine"""

    args = parse_args()
    gtf_to_t2g(args.input_file, args.output_file)


if __name__ == "__main__":
    main()
