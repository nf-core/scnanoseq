#!/usr/bin/env python3
"""Deduplication QC for the scnanoseq DNA branch.

Measures what Picard MarkDuplicates actually keys on in an ONT single-cell DNA
BAM, and what alternative grouping keys would do to the duplicate rate. Written
to answer whether the soft-clip problem that broke `umi_tools` on the cDNA
branch also affects the DNA branch, and whether the per-gene/per-contig fix
applied there should be ported.

MarkDuplicates keys single-end ("fragment") reads on

    (library, barcode, refIndex, 5' UNCLIPPED coordinate, strand)

i.e. getUnclippedStart() for forward reads and getUnclippedEnd() for reverse.
The other end of the fragment and its length are not part of the key. This
script reconstructs that key, checks it reproduces the flag-1024 marks in the
BAM, and only then reports the comparisons that depend on it.

Windows are processed one at a time and reads are attributed to the window
containing their reference_start, so no read is counted twice. Duplicate groups
therefore stay intact within a window, except for the rare reverse-strand pair
that shares a 5' anchor while its reference_starts straddle a window boundary.

Usage:
    dedup_qc.py --bam in.bam --auto-windows --window-size 400000 --out report.txt
    dedup_qc.py --bam in.bam --region chr1:20000000-20300000
"""

import argparse
import collections
import json
import random
import sys

import pysam

# Aligned-length strata, for checking whether short alignments behave differently
LENGTH_STRATA = [(0, 200, "<200bp"), (200, 1000, "200bp-1kb"),
                 (1000, 5000, "1-5kb"), (5000, 10 ** 9, ">5kb")]

# Bin sizes for the "what if we grouped by a feature instead of a position" ladder
BIN_LADDER = [10, 100, 1000, 10000, 100000]

# Two reads are treated as plausibly the same molecule when their 3' ends agree
# this closely. Used for the jitter measurement only.
SAME_MOLECULE_3P_TOL = 20

# Both-end clustering tolerances reported as (5' tol, 3' tol)
CLUSTER_TOLERANCES = [(0, 20), (2, 20), (5, 20), (10, 20),
                      (5, 50), (10, 50), (10, 200), (20, 200)]


def stratum(length):
    for lo, hi, label in LENGTH_STRATA:
        if lo <= length < hi:
            return label
    return LENGTH_STRATA[-1][2]


class Quantiles:
    """Exact quantiles from an integer-keyed histogram."""

    def __init__(self):
        self.hist = collections.Counter()

    def add(self, value):
        self.hist[value] += 1

    @property
    def n(self):
        return sum(self.hist.values())

    def summary(self):
        total = self.n
        if not total:
            return {}
        targets = [("p10", .10), ("p25", .25), ("median", .50),
                   ("p75", .75), ("p90", .90), ("p99", .99)]
        out, seen, pending = {}, 0, list(targets)
        keys = sorted(self.hist)
        for k in keys:
            seen += self.hist[k]
            while pending and seen >= pending[0][1] * total:
                out[pending.pop(0)[0]] = k
        for name, _ in pending:
            out[name] = keys[-1]
        out.update(n=total, min=keys[0], max=keys[-1],
                   mean=round(sum(k * c for k, c in self.hist.items()) / total, 2))
        return out

    def share_above(self, threshold):
        total = self.n
        if not total:
            return 0.0
        return 100.0 * sum(c for k, c in self.hist.items() if k > threshold) / total

    def share_equal(self, value):
        total = self.n
        return 100.0 * self.hist.get(value, 0) / total if total else 0.0


def read_records(bam, chrom, start, end, stats):
    """Return per-primary-read tuples for one window, accumulating flag counts.

    Each record is (cb, refid, pos5_unclipped, pos5_clipped, pos3_clipped,
    is_reverse, is_duplicate, ref_len).
    """
    recs = []
    for read in bam.fetch(chrom, start, end):
        if read.reference_start is None or not (start <= read.reference_start < end):
            continue
        stats["records"] += 1
        if read.is_secondary:
            stats["secondary"] += 1
            stats["secondary_flagged_dup"] += read.is_duplicate
            continue
        if read.is_supplementary:
            stats["supplementary"] += 1
            stats["supplementary_flagged_dup"] += read.is_duplicate
            continue
        cigar = read.cigartuples
        if not cigar:
            continue
        stats["primary"] += 1
        stats["primary_dup"] += read.is_duplicate
        stats["reverse"] += read.is_reverse
        if read.has_tag("DT"):
            stats["dt_" + str(read.get_tag("DT"))] += 1
        for tag in ("UB", "UR", "UY"):
            if read.has_tag(tag):
                stats["umi_tag_" + tag] += 1

        lead = cigar[0][1] if cigar[0][0] == 4 else 0
        trail = cigar[-1][1] if cigar[-1][0] == 4 else 0
        # The read's own 5' end is the left end in reference coordinates for a
        # forward alignment and the right end for a reverse one. Picard anchors
        # on that end, so that is the clip length that can move the key.
        if read.is_reverse:
            clip5, clip3 = trail, lead
            pos5_unclipped = read.reference_end + trail
            pos5_clipped = read.reference_end
            pos3_clipped = read.reference_start
        else:
            clip5, clip3 = lead, trail
            pos5_unclipped = read.reference_start - lead
            pos5_clipped = read.reference_start
            pos3_clipped = read.reference_end

        ref_len = read.reference_length or 0
        strat = stratum(ref_len)
        stats["clip5"].add(clip5)
        stats["clip3"].add(clip3)
        stats["clip5_by_len"][strat].add(clip5)
        stats["aligned_len"].add(ref_len)
        cb = read.get_tag("CB") if read.has_tag("CB") else "NA"
        stats["per_cb"][cb][0] += 1
        stats["per_cb"][cb][1] += 1 if read.is_duplicate else 0
        recs.append((sys.intern(cb), read.reference_id, pos5_unclipped,
                     pos5_clipped, pos3_clipped, read.is_reverse,
                     read.is_duplicate, ref_len))
    return recs


def count_dups(recs, keyfn):
    """Duplicates = members beyond the first in every group of an equality key."""
    groups = collections.Counter(keyfn(r) for r in recs)
    return sum(c - 1 for c in groups.values() if c > 1)


def cluster_dups(recs, tol5, tol3):
    """Greedy clustering on the 5' anchor within tol5, requiring the 3' end to
    agree within tol3. Models a jitter-tolerant, both-end-aware deduplicator."""
    by_cell = collections.defaultdict(list)
    for i, r in enumerate(recs):
        by_cell[(r[0], r[1], r[5])].append(i)
    dups = 0
    for members in by_cell.values():
        if len(members) < 2:
            continue
        members.sort(key=lambda i: recs[i][2])
        used = [False] * len(members)
        for a in range(len(members)):
            if used[a]:
                continue
            used[a] = True
            size = 1
            anchor = recs[members[a]]
            for b in range(a + 1, len(members)):
                if used[b]:
                    continue
                cand = recs[members[b]]
                if cand[2] - anchor[2] > tol5:
                    break
                if abs(cand[4] - anchor[4]) <= tol3:
                    used[b] = True
                    size += 1
            dups += size - 1
    return dups


def hamming(a, b):
    return sum(1 for x, y in zip(a, b) if x != y) if len(a) == len(b) else -1


def analyse_window(recs, acc):
    """Accumulate every key comparison for one window."""
    acc["n_primary"] += len(recs)
    acc["picard_marked"] += sum(1 for r in recs if r[6])

    # --- The Picard-equivalent key, and validation against the actual flags ---
    groups = collections.defaultdict(list)
    for i, r in enumerate(recs):
        groups[(r[0], r[1], r[2], r[5])].append(i)
    in_multi = set()
    for members in groups.values():
        if len(members) > 1:
            in_multi.update(members)
            acc["picard_key_dups"] += len(members) - 1
            acc["dup_set_sizes"].add(len(members))
            ends = [recs[i][4] for i in members]
            lens = [max(1, recs[i][7]) for i in members]
            acc["set_3p_spread"].add(max(ends) - min(ends))
            acc["set_len_ratio"].add(int(round(10 * max(lens) / min(lens))))
    acc["picard_marked_in_key_group"] += sum(1 for i, r in enumerate(recs)
                                             if r[6] and i in in_multi)

    # --- Alternative equality keys ---
    for name, keyfn in (
        ("picard_equivalent", lambda r: (r[0], r[1], r[2], r[5])),
        ("clipped_5p_anchor", lambda r: (r[0], r[1], r[3], r[5])),
        ("both_ends_exact", lambda r: (r[0], r[1], r[2], r[4], r[5])),
        ("no_barcode", lambda r: (r[1], r[2], r[5])),
    ):
        acc["keys"][name] += count_dups(recs, keyfn)
    for size in BIN_LADDER:
        acc["keys"]["bin_%d" % size] += count_dups(
            recs, lambda r, s=size: (r[0], r[1], r[2] // s))
    # Whole window as one bin: a lower bound on what --per-contig would collapse
    acc["keys"]["whole_window"] += count_dups(recs, lambda r: (r[0], r[1]))

    # --- Jitter-tolerant clustering ---
    for tol5, tol3 in CLUSTER_TOLERANCES:
        acc["cluster"]["5p%d_3p%d" % (tol5, tol3)] += cluster_dups(recs, tol5, tol3)

    # --- 5'-anchor offset histogram for pairs that look like the same molecule ---
    by_cell = collections.defaultdict(list)
    for i, r in enumerate(recs):
        by_cell[(r[0], r[1], r[5])].append(i)
    for members in by_cell.values():
        if len(members) < 2:
            continue
        members.sort(key=lambda i: recs[i][2])
        for a in range(len(members) - 1):
            for b in range(a + 1, min(a + 6, len(members))):
                d5 = recs[members[b]][2] - recs[members[a]][2]
                if d5 > 200:
                    break
                if abs(recs[members[b]][4] - recs[members[a]][4]) <= SAME_MOLECULE_3P_TOL:
                    acc["jitter"][d5] += 1
                    acc["jitter_by_len"][stratum(recs[members[a]][7])][d5] += 1

    # --- Barcode Hamming check on coordinate-identical read pairs ---
    coord = collections.defaultdict(list)
    for i, r in enumerate(recs):
        coord[(r[1], r[2], r[4], r[5])].append(i)
    for members in coord.values():
        if len(members) < 2:
            continue
        cbs = [recs[i][0] for i in members]
        for a in range(len(cbs)):
            for b in range(a + 1, len(cbs)):
                acc["coord_pairs"] += 1
                if cbs[a] == cbs[b]:
                    acc["coord_pairs_same_cb"] += 1
                else:
                    acc["cb_hamming"][hamming(cbs[a], cbs[b])] += 1


def build_windows(bam, args):
    if args.region:
        out = []
        for spec in args.region:
            chrom, _, span = spec.partition(":")
            if span:
                lo, _, hi = span.partition("-")
                out.append((chrom, int(lo.replace(",", "")), int(hi.replace(",", ""))))
            else:
                out.append((chrom, 0, bam.get_reference_length(chrom)))
        return out
    names = set(bam.references)
    wanted = ["chr%s" % c for c in list(range(1, 23)) + ["X"]]
    wanted = [c for c in wanted if c in names]
    if not wanted:
        wanted = [c for c in [str(x) for x in list(range(1, 23)) + ["X"]] if c in names]
    if not wanted:
        sys.exit("no standard chromosomes found in the BAM header; pass --region")
    out = []
    for chrom in wanted:
        length = bam.get_reference_length(chrom)
        mid = length // 2
        start = max(0, mid - args.window_size // 2)
        out.append((chrom, start, min(length, start + args.window_size)))
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bam", required=True)
    p.add_argument("--region", action="append",
                   help="chrom:start-end (repeatable). Default: one window per chromosome.")
    p.add_argument("--auto-windows", action="store_true",
                   help="one window per standard chromosome, centred on its midpoint")
    p.add_argument("--window-size", type=int, default=400000)
    p.add_argument("--out", help="write the report here as well as to stdout")
    p.add_argument("--json", help="write all measurements as JSON")
    p.add_argument("--min-key-agreement", type=float, default=99.9,
                   help="fail if the reconstructed Picard key explains fewer than "
                        "this percent of flagged duplicates (default: 99.9)")
    p.add_argument("--seed", type=int, default=1)
    args = p.parse_args()

    bam = pysam.AlignmentFile(args.bam, "rb")
    windows = build_windows(bam, args)

    stats = collections.defaultdict(int)
    stats["clip5"] = Quantiles()
    stats["clip3"] = Quantiles()
    stats["aligned_len"] = Quantiles()
    stats["clip5_by_len"] = collections.defaultdict(Quantiles)
    stats["per_cb"] = collections.defaultdict(lambda: [0, 0])

    acc = collections.defaultdict(int)
    acc["keys"] = collections.Counter()
    acc["cluster"] = collections.Counter()
    acc["jitter"] = collections.Counter()
    acc["jitter_by_len"] = collections.defaultdict(collections.Counter)
    acc["cb_hamming"] = collections.Counter()
    acc["dup_set_sizes"] = Quantiles()
    acc["set_3p_spread"] = Quantiles()
    acc["set_len_ratio"] = Quantiles()

    for chrom, start, end in windows:
        recs = read_records(bam, chrom, start, end, stats)
        print("  %s:%s-%s  %s primary reads" % (chrom, format(start, ","),
                                                format(end, ","), format(len(recs), ",")),
              file=sys.stderr)
        analyse_window(recs, acc)
        del recs

    report = []
    out = report.append
    n = acc["n_primary"]
    if not n:
        sys.exit("no primary reads found in the requested windows")

    def pct(v, tot=None):
        return "%.2f%%" % (100.0 * v / (tot or n))

    out("# Deduplication QC -- %s" % args.bam)
    out("")
    out("windows: %d" % len(windows))
    for chrom, start, end in windows:
        out("  %s:%s-%s" % (chrom, format(start, ","), format(end, ",")))
    out("")

    out("## 1. Record composition")
    out("")
    total_records = stats["records"]
    out("records in windows      : %10s" % format(total_records, ","))
    for key in ("primary", "secondary", "supplementary"):
        out("  %-20s: %10s  (%s)" % (key, format(stats[key], ","),
                                     pct(stats[key], total_records)))
    out("primary flagged dup     : %10s  (%s)" % (format(stats["primary_dup"], ","),
                                                  pct(stats["primary_dup"])))
    out("reverse strand          : %10s  (%s)" % (format(stats["reverse"], ","),
                                                  pct(stats["reverse"])))
    out("secondary flagged dup   : %10s   <- Picard does not propagate the flag to these"
        % format(stats["secondary_flagged_dup"], ","))
    out("supplementary flagged   : %10s   <- nor to these; filter with -F 0xD04, not -F 0x400"
        % format(stats["supplementary_flagged_dup"], ","))
    dt = dict((k[3:], v) for k, v in stats.items()
              if isinstance(k, str) and k.startswith("dt_"))
    out("DT tag values           : %s" % (dt or "none"))
    umi = dict((k[8:], v) for k, v in stats.items()
               if isinstance(k, str) and k.startswith("umi_tag_"))
    out("UMI tags present        : %s"
        % (umi or "NONE -- position is the only molecular identifier"))
    out("")

    out("## 2. Soft-clip geometry")
    out("")
    out("Picard keys single-end reads on the 5' UNCLIPPED coordinate, so only the")
    out("clip at the read's own 5' end can move the key.")
    out("")
    out("%-34s%8s%8s%8s%10s%9s%9s" % ("", "median", "p90", "p99", "mean", ">20bp", "==0"))
    for key, label in (("clip5", "5' end  (Picard keys on this)"),
                       ("clip3", "3' end  (Picard ignores this)")):
        q, s = stats[key].summary(), stats[key]
        out("%-34s%8d%8d%8d%10s%8.2f%%%8.2f%%"
            % (label, q["median"], q["p90"], q["p99"], q["mean"],
               s.share_above(20), s.share_equal(0)))
    out("")
    out("5' clip by aligned length:")
    for _, _, label in LENGTH_STRATA:
        s = stats["clip5_by_len"].get(label)
        if s and s.n:
            out("  %-12s n=%9s  median=%5d  ==0 %6.2f%%  >20bp %6.2f%%"
                % (label, format(s.n, ","), s.summary()["median"],
                   s.share_equal(0), s.share_above(20)))
    out("")
    q = stats["aligned_len"].summary()
    out("aligned length: median=%s p75=%s p90=%s p99=%s max=%s"
        % tuple(format(q[k], ",") for k in ("median", "p75", "p90", "p99", "max")))
    out("")

    out("## 3. Validation of the reconstructed Picard key")
    out("")
    marked = acc["picard_marked"]
    agreement = 100.0 * acc["picard_marked_in_key_group"] / marked if marked else 0.0
    out("primary reads                          : %10s" % format(n, ","))
    out("flagged duplicates in the BAM          : %10s  (%s)" % (format(marked, ","), pct(marked)))
    out("duplicates implied by the rebuilt key  : %10s" % format(acc["picard_key_dups"], ","))
    out("flagged duplicates inside a key group  : %10s  (%.3f%%)"
        % (format(acc["picard_marked_in_key_group"], ","), agreement))
    out("")

    out("## 4. Duplicate rate under alternative keys")
    out("")
    out("%-52s%11s%9s%12s" % ("key", "dups", "rate", "kept"))
    labels = [
        ("picard_equivalent", "(CB, ref, 5' unclipped, strand)  = current"),
        ("clipped_5p_anchor", "(CB, ref, 5' CLIPPED, strand)"),
        ("both_ends_exact", "(CB, ref, 5' unclipped, 3', strand)"),
        ("no_barcode", "(ref, 5' unclipped, strand)  -- no CB"),
    ]
    labels += [("bin_%d" % s, "(CB, ref, 5' // %s bp bin)  -- feature-style" % format(s, ","))
               for s in BIN_LADDER]
    labels += [("whole_window", "(CB, ref)  -- whole window as one feature")]
    for key, label in labels:
        d = acc["keys"][key]
        out("%-52s%11s%9s%12s" % (label, format(d, ","), pct(d), format(n - d, ",")))
    out("")
    out("A real --per-contig grouping spans a whole chromosome, so it would collapse")
    out("at least as much as the whole-window row above.")
    out("")

    out("## 5. Jitter-tolerant, both-end-aware clustering")
    out("")
    base = acc["keys"]["picard_equivalent"]
    out("%-28s%11s%9s%13s" % ("tolerance", "dups", "rate", "vs current"))
    for tol5, tol3 in CLUSTER_TOLERANCES:
        d = acc["cluster"]["5p%d_3p%d" % (tol5, tol3)]
        out("%-28s%11s%9s%+13s" % ("5p +-%d, 3p +-%d" % (tol5, tol3),
                                   format(d, ","), pct(d), format(d - base, ",")))
    out("")
    jit = acc["jitter"]
    tot = sum(jit.values())
    out("5'-anchor offset between read pairs whose 3' ends agree within %d bp (n=%s):"
        % (SAME_MOLECULE_3P_TOL, format(tot, ",")))
    if tot:
        for lo, hi, label in ((0, 0, "0  (caught by Picard)"), (1, 10, "1-10  (MISSED)"),
                              (11, 50, "11-50"), (51, 200, "51-200")):
            v = sum(c for d, c in jit.items() if lo <= d <= hi)
            out("  offset %-26s%10s  (%6.2f%%)" % (label, format(v, ","), 100.0 * v / tot))
        out("  head of the distribution: "
            + ", ".join("%dbp=%s" % (d, format(jit[d], ",")) for d in sorted(jit)[:8]))
        out("")
        out("  by aligned length (share of same-molecule pairs missed at 1-10 bp):")
        for _, _, label in LENGTH_STRATA:
            h = acc["jitter_by_len"].get(label)
            if h and sum(h.values()):
                sub = sum(h.values())
                miss = sum(c for d, c in h.items() if 1 <= d <= 10)
                out("    %-12s n=%9s  missed %6.2f%%"
                    % (label, format(sub, ","), 100.0 * miss / sub))
    out("")

    out("## 6. Over-dedup: how different are the molecules inside a duplicate set")
    out("")
    for key, label in (("dup_set_sizes", "set size"),
                       ("set_3p_spread", "3'-end spread within set (bp)")):
        q = acc[key].summary()
        if q:
            out("%-32s median=%d p90=%d p99=%d max=%s"
                % (label, q["median"], q["p90"], q["p99"], format(q["max"], ",")))
    spread = acc["set_3p_spread"]
    nsets = spread.n
    if nsets:
        for thr in (20, 100, 500, 1000):
            out("  sets whose 3' ends differ by >%5d bp: %6.2f%% of %s sets"
                % (thr, spread.share_above(thr), format(nsets, ",")))
        ratio = acc["set_len_ratio"]  # stored as round(10 * ratio)
        for thr in (1.5, 2, 5, 10):
            out("  sets with aligned-length ratio >%4sx  : %6.2f%%"
                % (thr, ratio.share_above(int(thr * 10))))
    out("")

    out("## 7. Barcode integrity")
    out("")
    cp, same = acc["coord_pairs"], acc["coord_pairs_same_cb"]
    out("read pairs at identical (ref, 5', 3', strand): %10s" % format(cp, ","))
    if cp:
        out("  same CB   (called duplicates)             : %10s  (%s)"
            % (format(same, ","), pct(same, cp)))
        out("  different CB (kept as distinct cells)     : %10s  (%s)"
            % (format(cp - same, ","), pct(cp - same, cp)))
    hd = acc["cb_hamming"]
    random.seed(args.seed)
    cbs = [cb for cb in stats["per_cb"] if cb != "NA"]
    bg = collections.Counter()
    if len(cbs) > 2:
        for _ in range(200000):
            a, b = random.sample(cbs, 2)
            bg[hamming(a, b)] += 1
    hdt, bgt = sum(hd.values()), sum(bg.values())
    out("Hamming distance of the differing CBs (n=%s) vs a random-pair background:"
        % format(hdt, ","))
    for d in range(1, 5):
        o = 100.0 * hd.get(d, 0) / hdt if hdt else 0.0
        e = 100.0 * bg.get(d, 0) / bgt if bgt else 0.0
        flag = "  <- barcode errors would show up here" if d == 1 else ""
        out("  d%d: observed %7.4f%% (%s)   background %7.4f%%%s"
            % (d, o, format(hd.get(d, 0), ","), e, flag))
    out("")
    percb = stats["per_cb"]
    out("distinct CB: %s   reads without CB: %s"
        % (format(len(percb), ","), format(percb.get("NA", [0, 0])[0], ",")))
    out("duplicate rate by per-cell read depth in the sampled windows:")
    for lo, hi, label in ((1, 10, "1-10"), (11, 50, "11-50"), (51, 200, "51-200"),
                          (201, 1000, "201-1000"), (1001, 10 ** 9, ">1000")):
        sel = [(t, d) for cb, (t, d) in percb.items() if lo <= t <= hi]
        if sel:
            t = sum(x for x, _ in sel)
            d = sum(x for _, x in sel)
            out("  %9s: %6s cells, %10s reads, dup rate %5.2f%%"
                % (label, format(len(sel), ","), format(t, ","), 100.0 * d / t))

    text = "\n".join(report)
    print(text)
    if args.out:
        with open(args.out, "w") as fh:
            fh.write(text + "\n")
    if args.json:
        def plain(v):
            if isinstance(v, Quantiles):
                return v.summary()
            if isinstance(v, (collections.Counter, dict, collections.defaultdict)):
                return dict((str(k), plain(x)) for k, x in v.items())
            return v
        with open(args.json, "w") as fh:
            json.dump({"bam": args.bam,
                       "windows": ["%s:%d-%d" % w for w in windows],
                       "composition": dict((k, plain(v)) for k, v in stats.items()),
                       "analysis": dict((k, plain(v)) for k, v in acc.items())},
                      fh, indent=2, default=str)

    # With nothing flagged there is nothing to validate against -- that is the
    # case for a --skip_dedup bam, not a broken key.
    if marked and agreement < args.min_key_agreement:
        sys.exit("\nFAIL: the reconstructed Picard key explains only %.3f%% of flagged "
                 "duplicates (need %s%%). Every comparison above depends on that key "
                 "being right, so the numbers are not usable."
                 % (agreement, args.min_key_agreement))


if __name__ == "__main__":
    main()
