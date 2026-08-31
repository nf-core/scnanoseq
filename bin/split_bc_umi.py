#!/usr/bin/env python3

""" Given a flexiplex-demultiplexed fastq, split the barcode and umi back out
    of the read name into a separate synthetic fastq. lr-kallisto reads the
    barcode and umi positionally (see the -x string), so they need to live in a
    read of their own rather than in the read name.

    Flexiplex writes read names as:
        @{barcode}_{umi}#{original_read_id}_{strand}{n}of{m}
    with the adapter, barcode, umi and poly-T already stripped from the
    sequence, followed by a tab separated comment carrying CB/CR/UB/UR and the
    XB tag the assign module derives.

    By default the barcode comes from the read name, which is the barcode
    flexiplex assigned against the known list. With --barcode_tag it comes from
    that tag in the comment instead, which is what produces an all-droplet
    matrix: XB holds the assigned barcode for a called cell and the uncorrected
    one for a droplet below the knee, so both land in a single barcode space.
    The umi still comes from the read name either way -- flexiplex gives every
    read that keeps a barcode region a full width umi.

    XB is the RAW, uncorrected barcode for any droplet below the knee, so
    sequencing error alone mints a fresh droplet per read and nothing downstream
    bounds the column count -- on a full sample that is millions of near
    singletons against a few thousand real cells, and the lr-kallisto EM cost
    scales with it. --allowlist restricts the run to a given set of barcodes,
    which is how the all-droplet pass is kept to a sane width.
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
    arg_parser.add_argument(
        "-t",
        "--barcode_tag",
        required=False,
        default=None,
        type=str,
        help="Take the barcode from this fastq comment tag (e.g. 'XB') rather than from the read name",
    )
    arg_parser.add_argument(
        "-l",
        "--allowlist",
        required=False,
        default=None,
        type=str,
        help="Only keep reads whose barcode appears in this file. Either a bare list or a "
        "`barcode<TAB>count` file such as the flexiplex barcode counts",
    )
    arg_parser.add_argument(
        "-m",
        "--allowlist_min_reads",
        required=False,
        default=0,
        type=int,
        help="Only load allowlist barcodes seen at least this many times. Requires a "
        "`barcode<TAB>count` allowlist. 0 keeps every barcode in the file",
    )

    return arg_parser.parse_args()


def open_fastq(path):
    """Open a fastq for reading, transparently handling gzip"""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r", encoding="utf-8")


def read_allowlist(path, min_reads=0):
    """Read an allowlist into a set, optionally thresholding on a read count

    The file is either a bare list of barcodes or the `barcode<TAB>count` form
    the flexiplex discovery/merge steps produce. Blank lines are ignored and
    anything after the barcode is dropped, so both shapes work.

    min_reads > 0 keeps only the barcodes seen at least that many times, which
    is what bounds the all-droplet barcode space: the counts file lists every
    barcode flexiplex ever saw, the overwhelming majority of them single-read
    artefacts of sequencing error.

    The filtering is done here rather than with an awk pre-step so the module
    depends on nothing but python -- the container carries no guaranteed awk.
    """
    allowed = set()
    with open(path, "r", encoding="utf-8") as allowlist_in:
        for line_number, line in enumerate(allowlist_in, start=1):
            fields = line.split()
            if not fields:
                continue

            if min_reads > 0:
                if len(fields) < 2:
                    sys.exit(
                        f"--allowlist_min_reads needs a `barcode<TAB>count` allowlist, but "
                        f"{path} line {line_number} has no count column."
                    )
                try:
                    count = int(fields[1])
                except ValueError:
                    sys.exit(f"Could not read a count from {path} line {line_number}: {line.strip()!r}")
                if count < min_reads:
                    continue

            allowed.add(fields[0])

    if not allowed:
        sys.exit(
            f"The allowlist {path} yielded no barcodes"
            + (f" at --allowlist_min_reads {min_reads}." if min_reads > 0 else ".")
        )

    print(f"allowlist: {len(allowed)} barcodes from {path} (min reads {min_reads})", file=sys.stderr)
    return allowed


def parse_header(header, barcode_tag=None):
    """Pull the barcode, umi and original read id out of a flexiplex read name

    Args:
        header (str): The fastq header line, including the leading '@'
        barcode_tag (str): Take the barcode from this comment tag instead of
            from the read name. None reads it from the name.

    Returns:
        tuple: (barcode, umi, read_id), or None if the name is not in the
            flexiplex format. barcode is None when barcode_tag was asked for
            but the header does not carry it, which the caller counts apart
            from a malformed name.
    """
    # Drop the leading '@' and anything after the first whitespace. Flexiplex
    # appends the CB:Z:/CR:Z:/UB:Z:/UR:Z: tags, and the assign module appends
    # XB:Z:, as a tab separated comment.
    name = header[1:].split(None, 1)[0]

    bc_umi, sep, read_id = name.partition("#")
    if not sep:
        return None

    barcode, sep, umi = bc_umi.partition("_")
    if not sep:
        return None

    if barcode_tag is not None:
        prefix = f"{barcode_tag}:Z:"
        barcode = None
        for field in header.split("\t"):
            if field.startswith(prefix):
                barcode = field[len(prefix) :].strip()
                break

    return barcode, umi, read_id


def split_bc_umi(input_file, output_prefix, bc_length, umi_length, barcode_tag=None, allowlist=None):
    """Split a flexiplex fastq into a barcode+umi fastq and a cdna fastq

    The two fastqs are written in lockstep from a single pass, which is what
    kallisto requires: it reads the files positionally and never matches read
    names between them.
    """
    written = 0
    skipped_unparsed = 0
    skipped_untagged = 0
    skipped_unassigned = 0
    skipped_length = 0
    skipped_offlist = 0
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

            parsed = parse_header(header.rstrip("\n"), barcode_tag)
            if parsed is None:
                skipped_unparsed += 1
                continue

            barcode, umi, read_id = parsed

            # A read demultiplexed by something other than the flexiplex assign
            # module carries no XB. Count these apart from malformed headers so
            # the summary says which of the two went wrong.
            if barcode is None:
                skipped_untagged += 1
                continue

            # flexiplex runs with -a, so it reports the reads it could not assign
            # too, writing "-" in place of the barcode. kallisto has nothing to
            # place these against; count them apart from malformed headers so the
            # summary stays diagnostic. XB never carries "-": the assign module
            # drops the reads with no barcode region at all.
            if barcode == "-":
                skipped_unassigned += 1
                continue

            if len(barcode) != bc_length or len(umi) != umi_length:
                skipped_length += 1
                continue

            # Bound the barcode space. Checked after the length test so a
            # malformed barcode is still reported as malformed rather than
            # disappearing into the off-list count.
            if allowlist is not None and barcode not in allowlist:
                skipped_offlist += 1
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

    source = f"the {barcode_tag} tag" if barcode_tag else "the read name"
    untagged = f"{skipped_untagged} (no {barcode_tag} tag), " if barcode_tag else ""
    offlist = f"{skipped_offlist} (barcode not on the allowlist), " if allowlist is not None else ""
    print(
        f"wrote {written} reads, barcode from {source}; "
        f"skipped {skipped_unparsed} (unparseable read name), "
        f"{untagged}"
        f"{skipped_unassigned} (no barcode assigned), "
        f"{skipped_length} (unexpected barcode/umi length), "
        f"{offlist}"
        f"{skipped_empty} (empty sequence)",
        file=sys.stderr,
    )

    if written == 0:
        sys.exit("No reads could be split. Are these reads demultiplexed by flexiplex?")


def main():
    """Main subroutine"""

    args = parse_args()
    # An empty --barcode_tag means the same thing as not passing one at all.
    barcode_tag = args.barcode_tag or None
    allowlist = read_allowlist(args.allowlist, args.allowlist_min_reads) if args.allowlist else None
    split_bc_umi(
        args.input_file, args.output_prefix, args.bc_length, args.umi_length, barcode_tag, allowlist
    )


if __name__ == "__main__":
    main()
