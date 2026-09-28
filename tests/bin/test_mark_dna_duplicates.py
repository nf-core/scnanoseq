#!/usr/bin/env python3

""" Unit test for bin/mark_dna_duplicates.py on a toy bam with one case per rule.

    Not run by nf-test (which only collects *.nf.test). Run it in the module's container,
    which has pysam and numpy but no pytest:

        apptainer exec <mulled-v2-bb96c7354781ab52d8e69ccff89587598dc87fea image> \
            python -m unittest tests/bin/test_mark_dna_duplicates.py
"""

import os
import random
import subprocess
import sys
import tempfile
import unittest

import pysam

# the tool's output is read before it is indexed; htslib would warn about that each time
pysam.set_verbosity(0)

SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "bin", "mark_dna_duplicates.py")

ME_RC = "CTGTCTCTTATACACATCT"
CLIP_B = ME_RC + "CCGAGCCCACGAGAC" + "ACGTACGTAC"   # read-through into the Read2 adapter (A-B)
CLIP_A = ME_RC + "GACGCTGCCGACGA" + "ACGTACGTAC"    # read-through into the Read1 adapter (A-A)
RC = str.maketrans("ACGT", "TGCA")

# name, contig, start, aligned length, reverse, 3' clip (read orientation), XB, base quality
PRIMARIES = [
    # R1: same 5' and strand, far ends 200 bp apart
    ("r1a", "chr1", 1000, 300, False, "", "AAA", 30),
    ("r1b", "chr1", 1000, 500, False, "", "AAA", 10),
    # R1 with a read-through member: it is kept although its base qualities are lower
    ("r1c", "chr1", 3000, 300, False, "", "AAA", 40),
    ("r1d", "chr1", 3000, 300, False, CLIP_B, "AAA", 20),
    # R2: 5' +5 bp, far end +10 bp, both read through
    ("r2a", "chr1", 5000, 400, False, CLIP_B, "AAA", 20),
    ("r2b", "chr1", 5005, 405, False, CLIP_B, "AAA", 30),
    # not R2: 5' +15 bp
    ("r2c", "chr1", 8000, 400, False, CLIP_B, "AAA", 20),
    ("r2d", "chr1", 8015, 395, False, CLIP_B, "AAA", 30),
    # R3: one raw read written out once per orientation
    ("x_+1of1", "chr1", 11000, 300, False, "", "AAA", 20),
    ("x_-1of1", "chr1", 11000, 300, True, "", "AAA", 30),
    # not R3: concatemer pieces at different loci
    ("y_+1of2", "chr1", 14000, 300, False, "", "AAA", 20),
    ("y_+2of2", "chr1", 20000, 300, False, "", "AAA", 30),
    # R4: an A-A fragment and the same fragment read from the other end
    ("aa1", "chr1", 25000, 400, False, CLIP_A, "AAA", 30),
    ("aa2", "chr1", 25003, 402, True, CLIP_B, "AAA", 20),
    # different barcode: never merged
    ("bcx", "chr1", 30000, 300, False, "", "AAA", 30),
    ("bcy", "chr1", 30000, 300, False, "", "CCC", 20),
    # chrM: R1 and R3 apply, R2 and R4 do not
    ("m1a", "chrM", 100, 300, False, "", "AAA", 30),
    ("m1b", "chrM", 100, 300, False, "", "AAA", 10),
    ("m2a", "chrM", 2000, 400, False, CLIP_B, "AAA", 20),
    ("m2b", "chrM", 2005, 405, False, CLIP_B, "AAA", 30),
    ("maa1", "chrM", 5000, 400, False, CLIP_A, "AAA", 30),
    ("maa2", "chrM", 5003, 402, True, CLIP_B, "AAA", 20),
    ("mx_+1of1", "chrM", 8000, 300, False, "", "AAA", 20),
    ("mx_-1of1", "chrM", 8000, 300, True, "", "AAA", 30),
]
# name, contig, start, flag: the secondary/supplementary records of one duplicate (r1b)
# and one representative (r1a)
OTHERS = [
    ("r1b", "chr1", 60000, 2048),
    ("r1b", "chr1", 70000, 256),
    ("r1a", "chr1", 71000, 256),
]

DUPS_DEFAULT = {"r1b", "r1c", "r2a", "x_+1of1", "aa2", "m1b", "mx_+1of1"}
DUPS_NO_EXACT_ONLY = DUPS_DEFAULT | {"m2a", "maa2"}


def make_record(header, name, contig, start, length, reverse, clip, xb, q, flag=0):
    rng = random.Random(name + str(flag))
    body = "".join(rng.choice("ACGT") for _ in range(length))
    r = pysam.AlignedSegment(header)
    r.query_name = name
    r.reference_name = contig
    r.reference_start = start
    r.mapping_quality = 60
    r.flag = flag | (16 if reverse else 0)
    if reverse:
        # the 3' end of a reverse read is its leftmost stored base
        stored_clip = clip.translate(RC)[::-1]
        r.query_sequence = stored_clip + body
        r.cigartuples = ([(4, len(clip))] if clip else []) + [(0, length)]
    else:
        r.query_sequence = body + clip
        r.cigartuples = [(0, length)] + ([(4, len(clip))] if clip else [])
    r.query_qualities = pysam.qualitystring_to_array(chr(q + 33) * len(r.query_sequence))
    r.set_tag("XB", xb)
    return r


def build_bam(path):
    header = pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "unsorted"},
        "SQ": [{"SN": "chr1", "LN": 100000}, {"SN": "chrM", "LN": 16569}],
    })
    unsorted = path + ".unsorted.bam"
    with pysam.AlignmentFile(unsorted, "wb", header=header) as out:
        for name, contig, start, length, rev, clip, xb, q in PRIMARIES:
            out.write(make_record(header, name, contig, start, length, rev, clip, xb, q))
        for name, contig, start, flag in OTHERS:
            out.write(make_record(header, name, contig, start, 200, False, "", "AAA", 30, flag))
        u = pysam.AlignedSegment(header)
        u.query_name = "unmapped1"
        u.flag = 4
        u.query_sequence = "ACGT" * 50
        u.query_qualities = pysam.qualitystring_to_array("I" * 200)
        out.write(u)
    pysam.sort("-o", path, unsorted)
    pysam.index(path)
    os.remove(unsorted)


def read_summary(path):
    with open(path) as fh:
        next(fh)
        return dict(line.rstrip("\n").split("\t") for line in fh)


class MarkDnaDuplicatesTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.inp = os.path.join(cls.tmp.name, "in.bam")
        build_bam(cls.inp)
        cls.default = cls.run_tool("default")
        cls.no_exact_only = cls.run_tool("noexact", "--exact-only-contig", "chrX")

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    @classmethod
    def run_tool(cls, tag, *extra):
        d = cls.tmp.name
        out = {k: os.path.join(d, f"{tag}.{k}") for k in ("bam", "metrics", "summary")}
        subprocess.run([sys.executable, SCRIPT, "-i", cls.inp, "-o", out["bam"], "-m", out["metrics"],
                        "-s", out["summary"], "--barcode-tag", "XB", "--d5", "10", "--d3", "20",
                        "-t", "2", "--tmpdir", d, *extra], check=True, capture_output=True)
        return out

    def records(self, bam):
        with pysam.AlignmentFile(bam, check_sq=False) as fh:
            return list(fh.fetch(until_eof=True))

    def flagged_primaries(self, bam):
        return {r.query_name for r in self.records(bam)
                if r.is_duplicate and not (r.is_secondary or r.is_supplementary)}

    def test_rules_by_default(self):
        self.assertEqual(self.flagged_primaries(self.default["bam"]), DUPS_DEFAULT)

    def test_exact_only_contig_flag(self):
        # with chrM not in the exact-only set, R2 and R4 merge its pairs too
        self.assertEqual(self.flagged_primaries(self.no_exact_only["bam"]), DUPS_NO_EXACT_ONLY)

    def test_propagation_to_secondary_and_supplementary(self):
        others = {(r.reference_start, r.flag & 0x900): r.is_duplicate
                  for r in self.records(self.default["bam"]) if r.is_secondary or r.is_supplementary}
        self.assertEqual(others, {(60000, 0x800): True, (70000, 0x100): True, (71000, 0x100): False})

    def test_duplicate_set_tags(self):
        recs = {r.query_name: r for r in self.records(self.default["bam"])
                if not (r.is_secondary or r.is_supplementary or r.is_unmapped)}
        self.assertEqual((recs["r1a"].get_tag("DS"), recs["r1b"].get_tag("DS")), (2, 2))
        self.assertEqual(recs["r1a"].get_tag("DI"), recs["r1b"].get_tag("DI"))
        self.assertNotEqual(recs["r1a"].get_tag("DI"), recs["r2a"].get_tag("DI"))
        self.assertFalse(recs["bcx"].has_tag("DS"))
        self.assertFalse(recs["y_+1of2"].has_tag("DS"))

    def test_summary(self):
        s = read_summary(self.default["summary"])
        self.assertEqual(s["primary_reads"], str(len(PRIMARIES)))
        self.assertEqual(s["duplicates"], str(len(DUPS_DEFAULT)))
        self.assertEqual(s["dup_by_key_5prime_strand"], "3")
        self.assertEqual(s["dup_added_by_same_raw_read"], "2")
        self.assertEqual(s["dup_added_by_5prime_jitter"], "1")
        self.assertEqual(s["dup_added_by_aa_opposite_end"], "1")
        self.assertEqual(s["exact_only_contigs_primary_reads"], "8")
        self.assertEqual(s["exact_only_contigs_duplicates"], "2")
        self.assertEqual(s["nuclear_duplicate_rate"], f"{5 / 16:.6f}")
        self.assertEqual(s["flagged_secondary"], "1")
        self.assertEqual(s["flagged_supplementary"], "1")
        self.assertEqual(s["unmapped"], "1")
        self.assertEqual(s["param_exact_only_contigs"], "chrM,MT")
        s2 = read_summary(self.no_exact_only["summary"])
        self.assertEqual(s2["exact_only_contigs_primary_reads"], "0")
        self.assertEqual(s2["nuclear_duplicate_rate"], f"{len(DUPS_NO_EXACT_ONLY) / len(PRIMARIES):.6f}")

    def test_metrics(self):
        with open(self.default["metrics"]) as fh:
            lines = fh.read().splitlines()
        i = next(k for k, line in enumerate(lines) if line.startswith("LIBRARY"))
        row = dict(zip(lines[i].split("\t"), lines[i + 1].split("\t")))
        self.assertEqual(row["UNPAIRED_READS_EXAMINED"], str(len(PRIMARIES)))
        self.assertEqual(row["UNPAIRED_READ_DUPLICATES"], str(len(DUPS_DEFAULT)))

    def test_output_sorted_and_not_indexed(self):
        bam = self.default["bam"]
        self.assertFalse(os.path.exists(bam + ".bai"))
        recs = self.records(bam)
        self.assertEqual(len(recs), len(PRIMARIES) + len(OTHERS) + 1)
        keys = [(r.reference_id if r.reference_id >= 0 else 1 << 30, r.reference_start) for r in recs]
        self.assertEqual(keys, sorted(keys))
        pysam.index(bam)
        with pysam.AlignmentFile(bam) as fh:
            self.assertEqual(fh.header.to_dict()["HD"]["SO"], "coordinate")
            self.assertEqual(sum(1 for _ in fh.fetch("chrM")), 8)


if __name__ == "__main__":
    unittest.main()
