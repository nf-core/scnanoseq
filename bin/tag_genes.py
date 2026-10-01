#!/usr/bin/env python3

""" Annotate alignments with a gene assignment taken from GTF gene-body overlap.

    Three tags are written per alignment:

      GX:Z  gene_id of the maximum-overlap gene, or "__no_gene" if nothing overlaps
      GN:Z  gene_name of that gene (same sentinel when nothing overlaps)
      GS:Z  assignment status: "unique" | "ambiguous" | "none"

    GX always carries the winning gene even when the assignment is ambiguous, so
    a strict arm (skip ambiguous and none) and a lenient max-overlap-wins arm
    (skip none only) can both be driven off the same tags through umi_tools
    --skip-tags-regex on GS. Changing that policy therefore never requires the
    bam to be tagged again.

    Assignment is strand-agnostic: the pipeline aligns with minimap2 -un, so read
    strand carries no reliable orientation information.

    "Ambiguous" means a second gene also covers at least AMBIG_FRAC of the read's
    aligned length, mirroring the way featureCounts refuses to assign a read that
    straddles two features.

    umi_tools --per-gene DROPS reads whose gene tag is absent or empty, or whose
    --assigned-status-tag matches --skip-tags-regex, so the GS distribution
    written to the summary is the hard ceiling on what per-gene grouping can act
    on. The remainder has to be deduplicated positionally and merged back.
"""

import argparse
import collections
import gzip
import re
import sys

import pysam

ATTRIBUTE_RE = re.compile(r'(\S+)\s+"([^"]*)"')

GENE_TAG = "GX"
GENE_NAME_TAG = "GN"
STATUS_TAG = "GS"

NO_GENE = "__no_gene"

# A second gene covering at least this fraction of the aligned length makes the
# assignment ambiguous.
AMBIG_FRAC = 0.10

# Genes are registered in every bin they span so a read's bins yield all
# candidates without scanning every gene on the contig.
BIN_SIZE = 100_000


def parse_args():
    """Parse the commandline arguments"""

    arg_parser = argparse.ArgumentParser()

    arg_parser.add_argument("-b", "--in_bam", required=True, type=str, help="The input bam file")
    arg_parser.add_argument("-g", "--gtf", required=True, type=str, help="The gtf to take genes from")
    arg_parser.add_argument("-o", "--out_bam", required=True, type=str, help="The output tagged bam file")
    arg_parser.add_argument("-s", "--summary", required=True, type=str, help="The output tsv summarising GS")
    arg_parser.add_argument("-t", "--threads", default=4, type=int, help="Threads for bam (de)compression")

    return arg_parser.parse_args()


def open_gtf(path):
    """Open a gtf for reading, transparently handling gzip"""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r", encoding="utf-8")


def parse_attributes(field):
    """Parse the key "value" pairs of a GTF attribute column into a dict"""
    return dict(ATTRIBUTE_RE.findall(field))


def load_genes(gtf_path):
    """Return {contig: (starts, ends, ids, names)} with the genes sorted by start"""

    per_contig = collections.defaultdict(list)

    with open_gtf(gtf_path) as gtf_in:
        for line in gtf_in:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue

            attributes = parse_attributes(fields[8])
            gene_id = attributes.get("gene_id")
            if not gene_id:
                continue

            gene_name = attributes.get("gene_name", gene_id)

            # The gtf is 1-based inclusive, convert to 0-based half-open.
            per_contig[fields[0]].append((int(fields[3]) - 1, int(fields[4]), gene_id, gene_name))

    gene_tables = {}
    for contig, genes in per_contig.items():
        genes.sort()
        gene_tables[contig] = (
            [gene[0] for gene in genes],
            [gene[1] for gene in genes],
            [gene[2] for gene in genes],
            [gene[3] for gene in genes],
        )

    return gene_tables


def build_index(starts, ends):
    """Bin gene indices so an overlap lookup does not scan every gene on a contig"""

    bins = collections.defaultdict(list)

    for gene_idx, (start, end) in enumerate(zip(starts, ends)):
        for gene_bin in range(start // BIN_SIZE, end // BIN_SIZE + 1):
            bins[gene_bin].append(gene_idx)

    return bins


def assign(read, genes, bins):
    """Return the (GX, GN, GS) values for one alignment"""

    starts, ends, ids, names = genes

    # Aligned reference intervals, so introns do not count towards the overlap.
    blocks = read.get_blocks()
    if not blocks:
        return NO_GENE, NO_GENE, "none"

    aligned_len = sum(end - start for start, end in blocks)
    if aligned_len == 0:
        return NO_GENE, NO_GENE, "none"

    candidates = set()
    for block_start, block_end in blocks:
        for block_bin in range(block_start // BIN_SIZE, (block_end - 1) // BIN_SIZE + 1):
            candidates.update(bins.get(block_bin, ()))

    overlaps = {}
    for gene_idx in candidates:
        gene_start, gene_end = starts[gene_idx], ends[gene_idx]
        total = 0
        for block_start, block_end in blocks:
            low = block_start if block_start > gene_start else gene_start
            high = block_end if block_end < gene_end else gene_end
            if high > low:
                total += high - low
        if total:
            overlaps[gene_idx] = total

    if not overlaps:
        return NO_GENE, NO_GENE, "none"

    # Ties break on gene_id so the assignment does not depend on set ordering.
    ranked = sorted(overlaps.items(), key=lambda item: (-item[1], ids[item[0]]))
    best_idx = ranked[0][0]

    status = "unique"
    if len(ranked) > 1 and ranked[1][1] >= AMBIG_FRAC * aligned_len:
        status = "ambiguous"

    return ids[best_idx], names[best_idx], status


def tag_genes(in_bam, gtf, out_bam, summary, threads):
    """Write a copy of the bam with every alignment carrying GX, GN and GS"""

    gene_tables = load_genes(gtf)
    if not gene_tables:
        sys.exit("No genes found in the gtf. Does it contain 'gene' feature rows?")

    indexes = {contig: build_index(table[0], table[1]) for contig, table in gene_tables.items()}
    print(f"loaded genes on {len(gene_tables)} contigs", file=sys.stderr)

    counts = collections.Counter()

    with pysam.AlignmentFile(in_bam, "rb", threads=threads) as bam_in, pysam.AlignmentFile(
        out_bam, "wb", template=bam_in, threads=threads
    ) as bam_out:
        for read in bam_in:
            contig = read.reference_name

            if read.is_unmapped or contig not in gene_tables:
                gene_id, gene_name, status = NO_GENE, NO_GENE, "none"
            else:
                gene_id, gene_name, status = assign(read, gene_tables[contig], indexes[contig])

            read.set_tag(GENE_TAG, gene_id, value_type="Z")
            read.set_tag(GENE_NAME_TAG, gene_name, value_type="Z")
            read.set_tag(STATUS_TAG, status, value_type="Z")
            bam_out.write(read)

            counts[(contig, status, read.is_supplementary or read.is_secondary)] += 1

    with open(summary, "w", encoding="utf-8") as summary_out:
        summary_out.write("contig\tgs\tnon_primary\tcount\n")
        for (contig, status, non_primary), count in sorted(
            counts.items(), key=lambda item: (str(item[0][0]), item[0][1], item[0][2])
        ):
            summary_out.write(f"{contig}\t{status}\t{int(non_primary)}\t{count}\n")

    total = sum(counts.values())
    print(f"tagged {total} alignments", file=sys.stderr)

    by_status = collections.Counter()
    for (_, status, _), count in counts.items():
        by_status[status] += count

    for status, count in by_status.most_common():
        pct = 100 * count / total if total else 0
        print(f"  {status:<10} {count:>12}  ({pct:5.2f}%)", file=sys.stderr)


def main():
    """Main subroutine"""

    args = parse_args()
    tag_genes(args.in_bam, args.gtf, args.out_bam, args.summary, args.threads)


if __name__ == "__main__":
    main()
