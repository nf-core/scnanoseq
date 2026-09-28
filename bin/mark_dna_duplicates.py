#!/usr/bin/env python3

""" Mark PCR/optical duplicates in single-cell ONT DNA (Tn5, no UMI) alignments.

    Replaces Picard MarkDuplicates for the DNA branch. Picard keys single-end reads on
    (barcode, 5' unclipped coordinate, strand). Measured against a heterozygous-SNP
    ground truth on three libraries, that key misses most same-molecule reads for two
    reasons it cannot see:

      - one raw ONT read is often two records. Flexiplex reports a read once per
        barcode it finds, so a read carrying the barcode pattern at both ends (mostly
        A-A fragments, barcode adapter on both sides) is written twice, once per
        orientation (`..._+1of1` / `..._-1of1`). The two map to the same locus on
        opposite strands, so no positional key merges them;
      - basecalling/alignment jitter moves the 5' anchor by a few bases.

    Within one barcode, two primary reads are one molecule when any of these hold:

      R1  same 5' unclipped coordinate and strand (Picard's key). Applied whatever the
          far end does: 5'-identical reads with far ends >1 kb apart are still 95-98%
          one molecule.
      R2  both read through into the far adapter, 5' anchors within --d5 and
          adapter-anchored far ends within --d3 (greedy, anchored on the family's first
          member so a family's extent stays bounded). 1-10 bp of 5' jitter is <=4%
          distinct molecules; 11-20 bp is 12-52%, hence --d5 10.
      R3  same raw read (read name without flexiplex's `_+NofM` suffix).
      R4  an A-A fragment (Read1/barcode adapter at the 3' end too) read from the other
          end: opposite strand, fragment interval within --d5 on both ends.

    On the reference sequences given by --exact-only-contig (default chrM and MT) only
    R1 and R3 are applied. Mitochondrial coverage is so dense that chance 5'/far-end
    coincidences between distinct molecules dominate the tolerance rules (on one
    library chrM made 14.4M of 18.3M R2 merges and 101k of 109k R4 merges), and with
    no heterozygous sites there is nothing to validate them against.

    One representative per molecule stays unflagged: a read-through read if there is
    one, then the highest sum of base qualities >= 15 (Picard's default score).
    Duplicate status is propagated to the secondary and supplementary records of every
    duplicate primary, so `-F 0x400` alone gives a clean bam. Members of a duplicate set
    carry DS:i (set size) and DI:i (set index, unique within a reference sequence).

    The far end of a read is only used when the read reads through into the adapter,
    i.e. when its aligned 3' end is a Tn5 insertion site. That is decided by finding
    the Nextera mosaic end at the start of the 3' soft clip.

    Also writes a Picard-format DuplicationMetrics file, so MultiQC's picard module
    keeps reporting the rate, and a per-rule summary.
"""

import argparse
import array
import collections
import hashlib
import multiprocessing as mp
import os
import re
import shutil
import sys
import tempfile

import numpy as np
import pysam

ME_RC = "CTGTCTCTTATACACATCT"
ME_SEED = ME_RC[:12]
B_TAIL = "CCGAGCCCACGAGAC"       # Nextera Read2 rc: an A-B fragment
A_TAIL = "GACGCTGCCGACGA"        # Nextera Read1 rc: an A-A fragment
CLIP_WINDOW = 80
MAX_ED_ME = 4
MAX_ED_TAIL = 3
MAX_ED_FULL = 7

CLS_NONE, CLS_OTHER, CLS_ME, CLS_B, CLS_A = 0, 1, 2, 3, 4
CLASS_NAMES = ["none", "other", "ME", "B", "A"]

NAME_SUFFIX_RE = re.compile(r"_[+-]\d+of\d+$")
RC = str.maketrans("ACGTNacgtn", "TGCANtgcan")

class Infix:
    """Myers/Hyyro bit-vector search: best edit distance of a short pattern anywhere in a text."""

    def __init__(self, pattern):
        self.m = len(pattern)
        self.peq = collections.defaultdict(int)
        for i, c in enumerate(pattern):
            self.peq[c] |= 1 << i
        self.full = (1 << self.m) - 1
        self.high = 1 << (self.m - 1)

    def search(self, text):
        """(edit distance, start offset in text) of the best match."""
        m, full, high, peq = self.m, self.full, self.high, self.peq
        pv, mv, score = full, 0, m
        best, best_end = m + 1, -1
        for j, ch in enumerate(text):
            eq = peq.get(ch, 0)
            xv = eq | mv
            xh = (((eq & pv) + pv) ^ pv) | eq
            ph = mv | (~(xh | pv) & full)
            mh = pv & xh
            if ph & high:
                score += 1
            elif mh & high:
                score -= 1
            ph = (ph << 1) & full
            mh = (mh << 1) & full
            pv = mh | (~(xv | ph) & full)
            mv = ph & xv
            if score < best:
                best, best_end = score, j
        return best, max(0, best_end - m + 1)


ME_SEARCH = Infix(ME_RC)
B_SEARCH = Infix(B_TAIL[:12])
A_SEARCH = Infix(A_TAIL[:12])
ADAPT_B_SEARCH = Infix(ME_RC + B_TAIL)
ADAPT_A_SEARCH = Infix(ME_RC + A_TAIL)


def classify_clip3(clip):
    """(class, offset of the mosaic end) for a 3' soft clip given in read orientation.

    Fast path: an exact 12-mer of the mosaic end near the start of the clip, then an
    exact 8-mer of either tail right after it. Only reads that miss an exact seed pay
    for the bit-vector edit-distance search.
    """
    window = clip[:CLIP_WINDOW]
    if len(window) < 12:
        return CLS_NONE, -1
    off = window.find(ME_SEED, 0, 40)
    if off < 0:
        ed, off = ME_SEARCH.search(window)
        if ed > MAX_ED_ME:
            # a noisy mosaic end can still be recognised with its tail attached
            eb, ob = ADAPT_B_SEARCH.search(window)
            ea, oa = ADAPT_A_SEARCH.search(window)
            if min(eb, ea) > MAX_ED_FULL:
                return CLS_OTHER, -1
            if eb == ea:
                return CLS_ME, min(ob, oa)
            return (CLS_B, ob) if eb < ea else (CLS_A, oa)
    tail = window[off + len(ME_RC) - 2: off + len(ME_RC) + 18]
    if len(tail) < 10:
        return CLS_ME, off
    if tail.find(B_TAIL[:8], 0, 12) >= 0:
        return CLS_B, off
    if tail.find(A_TAIL[:8], 0, 12) >= 0:
        return CLS_A, off
    eb, _ = B_SEARCH.search(tail)
    ea, _ = A_SEARCH.search(tail)
    if eb <= MAX_ED_TAIL and eb < ea:
        return CLS_B, off
    if ea <= MAX_ED_TAIL and ea < eb:
        return CLS_A, off
    return CLS_ME, off


def name_hash(name):
    return int.from_bytes(hashlib.blake2b(name.encode(), digest_size=8).digest(), "little")


def raw_read_id(name):
    return NAME_SUFFIX_RE.sub("", name)


class DSU:
    """Union-find on a compact int64 array. A root is always the lowest index of its
    set, so every parent index is below its child's and roots resolve in one pass."""

    def __init__(self, n):
        self.p = array.array("q", range(n))

    def find(self, x):
        p = self.p
        while p[x] != x:
            p[x] = p[p[x]]
            x = p[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return False
        if ra < rb:
            self.p[rb] = ra
        else:
            self.p[ra] = rb
        return True

    def labels(self):
        p = self.p
        for i in range(len(p)):
            p[i] = p[p[i]]
        return np.frombuffer(p, dtype=np.int64)


def as_array(idx):
    """numpy index vector -> compact array.array, iterated without a Python list."""
    out = array.array("q")
    out.frombytes(np.ascontiguousarray(idx, dtype=np.int64).tobytes())
    return out


# ---------------------------------------------------------------------------------------
# Pass 1: decide molecules on one reference sequence
# ---------------------------------------------------------------------------------------

def decide_contig(job):
    bam_path, contig, args_d, outdir = job
    d5, d3, bc_tag = args_d["d5"], args_d["d3"], args_d["barcode_tag"]
    exact_only = contig in args_d["exact_only"]
    stats = collections.Counter()

    qhash = array.array("Q")
    bc = array.array("q")
    rev = array.array("b")
    p5u = array.array("q")
    far = array.array("q")
    lo = array.array("q")
    hi = array.array("q")
    rt = array.array("b")
    cls = array.array("b")
    rawh = array.array("Q")
    score = array.array("q")
    bc_ids = {}

    with pysam.AlignmentFile(bam_path) as bam:
        for r in bam.fetch(contig):
            if r.is_unmapped or r.is_secondary or r.is_supplementary:
                continue
            cig = r.cigartuples
            lead = cig[0][1] if cig[0][0] == 4 else 0
            trail = cig[-1][1] if cig[-1][0] == 4 else 0
            seq = r.query_sequence or ""
            if r.is_reverse:
                pos5u = r.reference_end + trail
                a5 = r.reference_end
                p3 = r.reference_start
                clip3 = seq[:lead].translate(RC)[::-1] if lead else ""
            else:
                pos5u = r.reference_start - lead
                a5 = r.reference_start
                p3 = r.reference_end
                clip3 = seq[len(seq) - trail:] if trail else ""
            c, off = classify_clip3(clip3)
            through = c >= CLS_ME and not r.has_tag("SA")
            if through:
                f = p3 - off if r.is_reverse else p3 + off
            else:
                f = p3
            b = r.get_tag(bc_tag) if r.has_tag(bc_tag) else ""
            if b not in bc_ids:
                bc_ids[b] = len(bc_ids)
            q = r.query_qualities
            if q is not None:
                qa = np.frombuffer(q, dtype=np.uint8)
                s = int(qa[qa >= 15].sum())
            else:
                s = 0
            qhash.append(name_hash(r.query_name))
            bc.append(bc_ids[b])
            rev.append(r.is_reverse)
            p5u.append(pos5u)
            far.append(f)
            lo.append(min(a5, f))
            hi.append(max(a5, f))
            rt.append(through)
            cls.append(c)
            rawh.append(name_hash(raw_read_id(r.query_name)))
            score.append(s)
            stats["class_" + CLASS_NAMES[c]] += 1

    n = len(qhash)
    stats["primary"] = n
    out = os.path.join(outdir, f"decide.{contig}.npz")
    if n == 0:
        np.savez(out, dup=np.zeros(0, np.uint64), mem=np.zeros(0, np.uint64),
                 di=np.zeros(0, np.int64), ds=np.zeros(0, np.int64))
        return contig, out, stats

    bc_a = np.frombuffer(bc, dtype=np.int64)
    rev_a = np.frombuffer(rev, dtype=np.int8)
    p5u_a = np.frombuffer(p5u, dtype=np.int64)
    far_a = np.frombuffer(far, dtype=np.int64)
    lo_a = np.frombuffer(lo, dtype=np.int64)
    hi_a = np.frombuffer(hi, dtype=np.int64)
    rt_a = np.frombuffer(rt, dtype=np.int8).astype(bool)
    cls_a = np.frombuffer(cls, dtype=np.int8)
    raw_a = np.frombuffer(rawh, dtype=np.uint64)
    score_a = np.frombuffer(score, dtype=np.int64)

    dsu = DSU(n)

    def link_equal(order, keys):
        """Union consecutive members of `order` whose key columns are all equal."""
        same = np.ones(len(order) - 1, dtype=bool)
        for k in keys:
            ks = k[order]
            same &= ks[1:] == ks[:-1]
        merged = 0
        for x in np.flatnonzero(same):
            merged += dsu.union(int(order[x]), int(order[x + 1]))
        return merged

    # R1: Picard's key
    o1 = np.lexsort((p5u_a, rev_a, bc_a))
    stats["dup_key"] = link_equal(o1, (bc_a, rev_a, p5u_a)) if n > 1 else 0

    # R3: records of one raw read under the same barcode. The intervals must overlap:
    # flexiplex also splits concatemers into `_+1of2`/`_+2of2` segments, which are
    # different molecules at different loci.
    merged = 0
    if n > 1:
        o3 = np.lexsort((lo_a, bc_a, raw_a))
        same = (raw_a[o3][1:] == raw_a[o3][:-1]) & (bc_a[o3][1:] == bc_a[o3][:-1])
        for x in np.flatnonzero(same).tolist():
            i, j = int(o3[x]), int(o3[x + 1])
            if min(hi_a[i], hi_a[j]) > max(lo_a[i], lo_a[j]):
                merged += dsu.union(i, j)
    stats["dup_raw_read"] = merged

    # R2: 5' jitter with the far end agreeing, read-through reads only. R2 and R4 both
    # iterate the read-through reads, so an empty set switches both off.
    merged = 0
    idx_rt = np.zeros(0, dtype=np.int64) if exact_only else np.flatnonzero(rt_a)
    o2 = idx_rt[np.lexsort((p5u_a[idx_rt], rev_a[idx_rt], bc_a[idx_rt]))]
    heads = []
    prev_group = None
    for i in as_array(o2):
        g = (bc[i], rev[i])
        if g != prev_group:
            heads = []
            prev_group = g
        placed = False
        pi, fi = p5u[i], far[i]
        for h in reversed(heads):
            if pi - h[0] > d5:
                break
            if abs(fi - h[1]) <= d3:
                merged += dsu.union(i, h[2])
                placed = True
                break
        if not placed:
            heads.append((pi, fi, i))
    stats["dup_jitter"] = merged

    # R4: A-A fragments read from the other end
    merged = 0
    o4 = as_array(idx_rt[np.lexsort((lo_a[idx_rt], bc_a[idx_rt]))])
    n4 = len(o4)
    for k in np.flatnonzero(cls_a[np.frombuffer(o4, dtype=np.int64)] == CLS_A).tolist():
        i = o4[k]
        for step in (1, -1):
            t = k + step
            while 0 <= t < n4:
                j = o4[t]
                if bc[j] != bc[i] or abs(lo[j] - lo[i]) > d5:
                    break
                if rev[j] != rev[i] and abs(hi[j] - hi[i]) <= d5:
                    merged += dsu.union(i, j)
                t += step
    stats["dup_aa"] = merged

    lab = dsu.labels()
    del dsu
    # representative: read-through first, then base-quality score, then lowest index
    order = np.lexsort((np.arange(n), -score_a, -rt_a.astype(np.int64), lab))
    first = np.ones(n, dtype=bool)
    first[1:] = lab[order][1:] != lab[order][:-1]
    is_rep = np.zeros(n, dtype=bool)
    is_rep[order[first]] = True
    uniq, inv, counts = np.unique(lab, return_inverse=True, return_counts=True)
    size = counts[inv]
    hashes = np.frombuffer(qhash, dtype=np.uint64)
    dup_h = np.sort(hashes[~is_rep])
    mem = size > 1
    mo = np.argsort(hashes[mem])
    np.savez(out, dup=dup_h, mem=hashes[mem][mo], di=inv[mem][mo].astype(np.int64),
             ds=size[mem][mo].astype(np.int64))
    stats["duplicates"] = int((~is_rep).sum())
    stats["molecules"] = int(len(uniq))
    if exact_only:
        stats["primary_exact_only"] = n
        stats["duplicates_exact_only"] = stats["duplicates"]
    return contig, out, stats


# ---------------------------------------------------------------------------------------
# Pass 2: write one reference sequence with flags and tags
# ---------------------------------------------------------------------------------------

def write_contig(job):
    bam_path, contig, dup_path, npz_path, out_path, header = job
    dup = np.load(dup_path, mmap_mode="r")
    own = np.load(npz_path) if npz_path else None
    mem, di, ds = (own["mem"], own["di"], own["ds"]) if own is not None else (None, None, None)
    stats = collections.Counter()
    with pysam.AlignmentFile(bam_path) as bam, \
            pysam.AlignmentFile(out_path, "wb", header=header) as out:
        it = bam.fetch(contig) if contig != "*" else bam.fetch("*")
        batch = []

        def flush():
            if not batch:
                return
            h = np.fromiter((name_hash(r.query_name) for r in batch), dtype=np.uint64, count=len(batch))
            if len(dup):
                k = np.searchsorted(dup, h)
                k[k >= len(dup)] = 0
                isdup = dup[k] == h
            else:
                isdup = np.zeros(len(batch), dtype=bool)
            if mem is not None and len(mem):
                km = np.searchsorted(mem, h)
                km[km >= len(mem)] = 0
                ismem = mem[km] == h
            else:
                ismem = np.zeros(len(batch), dtype=bool)
            for x, r in enumerate(batch):
                primary = not (r.is_secondary or r.is_supplementary)
                d = bool(isdup[x]) and not r.is_unmapped
                r.is_duplicate = d
                if d:
                    stats["flag_primary" if primary else
                          ("flag_secondary" if r.is_secondary else "flag_supplementary")] += 1
                if primary and ismem[x] and not r.is_unmapped:
                    r.set_tag("DI", int(di[km[x]]), "i")
                    r.set_tag("DS", int(ds[km[x]]), "i")
                if r.is_unmapped:
                    stats["unmapped"] += 1
                elif not primary:
                    stats["secondary_or_supplementary"] += 1
                out.write(r)
            batch.clear()

        for r in it:
            if contig == "*" and r.reference_id >= 0:
                # placed reads were already written with their reference sequence
                continue
            batch.append(r)
            if len(batch) >= 20000:
                flush()
        flush()
    return contig, out_path, stats


def write_metrics(path, input_name, argv, st):
    examined = st["primary"]
    dups = st["duplicates"]
    pct = dups / examined if examined else 0.0
    with open(path, "w") as fh:
        fh.write("## htsjdk.samtools.metrics.StringHeader\n")
        fh.write(f"# MarkDuplicates --INPUT {input_name} (emulated by mark_dna_duplicates.py: {' '.join(argv)})\n")
        fh.write("## htsjdk.samtools.metrics.StringHeader\n")
        fh.write("# Started on: n/a\n\n")
        fh.write("## METRICS CLASS\tpicard.sam.DuplicationMetrics\n")
        fh.write("LIBRARY\tUNPAIRED_READS_EXAMINED\tREAD_PAIRS_EXAMINED\tSECONDARY_OR_SUPPLEMENTARY_RDS\t"
                 "UNMAPPED_READS\tUNPAIRED_READ_DUPLICATES\tREAD_PAIR_DUPLICATES\t"
                 "READ_PAIR_OPTICAL_DUPLICATES\tPERCENT_DUPLICATION\tESTIMATED_LIBRARY_SIZE\n")
        fh.write(f"Unknown Library\t{examined}\t0\t{st['secondary_or_supplementary']}\t{st['unmapped']}\t"
                 f"{dups}\t0\t0\t{pct:.6f}\t\n\n")


def write_summary(path, st, args):
    nuc_primary = st["primary"] - st["primary_exact_only"]
    nuc_dups = st["duplicates"] - st["duplicates_exact_only"]
    rows = [
        ("primary_reads", st["primary"]),
        ("molecules", st["molecules"]),
        ("duplicates", st["duplicates"]),
        ("duplicate_rate", f"{st['duplicates'] / st['primary']:.6f}" if st["primary"] else "0"),
        # the rate without the exact-only (mitochondrial) contigs, which otherwise
        # dominate the headline in libraries with a high chrM fraction
        ("nuclear_duplicate_rate", f"{nuc_dups / nuc_primary:.6f}" if nuc_primary else "0"),
        ("exact_only_contigs_primary_reads", st["primary_exact_only"]),
        ("exact_only_contigs_duplicates", st["duplicates_exact_only"]),
        ("dup_by_key_5prime_strand", st["dup_key"]),
        ("dup_added_by_same_raw_read", st["dup_raw_read"]),
        ("dup_added_by_5prime_jitter", st["dup_jitter"]),
        ("dup_added_by_aa_opposite_end", st["dup_aa"]),
        ("flagged_secondary", st["flag_secondary"]),
        ("flagged_supplementary", st["flag_supplementary"]),
        ("secondary_or_supplementary", st["secondary_or_supplementary"]),
        ("unmapped", st["unmapped"]),
    ] + [(f"clip3_class_{c}", st["class_" + c]) for c in CLASS_NAMES] + [
        ("param_barcode_tag", args.barcode_tag), ("param_d5", args.d5), ("param_d3", args.d3),
        ("param_exact_only_contigs", ",".join(args.exact_only_contig))]
    with open(path, "w") as fh:
        fh.write("metric\tvalue\n")
        for k, v in rows:
            fh.write(f"{k}\t{v}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--input", required=True, help="coordinate-sorted, indexed bam")
    ap.add_argument("-o", "--output", required=True,
                    help="output bam (coordinate-sorted, not indexed: index it with samtools index -@)")
    ap.add_argument("-m", "--metrics", required=True, help="Picard-format DuplicationMetrics file")
    ap.add_argument("-s", "--summary", required=True, help="per-rule summary tsv")
    ap.add_argument("--barcode-tag", default="XB")
    ap.add_argument("--d5", type=int, default=10, help="5' anchor tolerance for read-through reads (bp)")
    ap.add_argument("--d3", type=int, default=20, help="far-end tolerance for read-through reads (bp)")
    ap.add_argument("--exact-only-contig", action="append", metavar="NAME",
                    help="reference sequence on which only R1 and R3 apply (repeatable; "
                         "default: chrM and MT)")
    ap.add_argument("-t", "--threads", type=int, default=1)
    ap.add_argument("--tmpdir", default=".")
    args = ap.parse_args()
    if args.exact_only_contig is None:
        args.exact_only_contig = ["chrM", "MT"]

    tmp = tempfile.mkdtemp(prefix="markdup_dna.", dir=args.tmpdir)
    with pysam.AlignmentFile(args.input) as bam:
        if not bam.has_index():
            sys.exit(f"{args.input} has no index; the tool fetches one reference sequence at a time")
        hdr = bam.header.to_dict()
        if hdr.get("HD", {}).get("SO") != "coordinate":
            sys.exit(f"{args.input} is not coordinate-sorted")
        refs = list(bam.references)
        mapped = {s.contig: s.mapped for s in bam.get_index_statistics()}
    work = [c for c in refs if mapped.get(c, 0) > 0]
    work.sort(key=lambda c: -mapped[c])     # largest first, for load balance
    args_d = {"d5": args.d5, "d3": args.d3, "barcode_tag": args.barcode_tag,
              "exact_only": frozenset(args.exact_only_contig)}

    total = collections.Counter()
    npz = {}
    with mp.Pool(max(1, args.threads)) as pool:
        for contig, path, st in pool.imap_unordered(decide_contig, [(args.input, c, args_d, tmp) for c in work]):
            npz[contig] = path
            total.update(st)
            print(f"[decide] {contig}: {st['primary']:,} primary, {st['duplicates']:,} duplicates", file=sys.stderr, flush=True)

    dup_all = np.sort(np.concatenate([np.load(p)["dup"] for p in npz.values()] or [np.zeros(0, np.uint64)]))
    dup_path = os.path.join(tmp, "dup_hashes.npy")
    np.save(dup_path, dup_all)
    del dup_all

    pg = {"ID": "mark_dna_duplicates", "PN": "mark_dna_duplicates.py", "CL": " ".join(sys.argv)}
    prev = [p["ID"] for p in hdr.get("PG", [])]
    if prev:
        pg["PP"] = prev[-1]
    ids = set(prev)
    while pg["ID"] in ids:
        pg["ID"] += ".1"
    hdr.setdefault("PG", []).append(pg)
    jobs = [(args.input, c, dup_path, npz.get(c), os.path.join(tmp, f"out.{i:05d}.bam"), hdr)
            for i, c in enumerate(refs) if mapped.get(c, 0) > 0]
    jobs.append((args.input, "*", dup_path, None, os.path.join(tmp, "out.99999.unmapped.bam"), hdr))
    parts = {}
    with mp.Pool(max(1, args.threads)) as pool:
        for contig, path, st in pool.imap_unordered(write_contig, sorted(jobs, key=lambda j: -mapped.get(j[1], 0))):
            parts[path] = True
            total.update(st)
    ordered = [j[4] for j in jobs]
    # Parts are concatenated in header order with the unmapped reads last, so the result
    # is coordinate-sorted. It is not indexed here: a single-threaded index took over a
    # third of the run on a full library, and samtools index -@ does it downstream.
    pysam.cat("-o", args.output, *ordered)

    write_metrics(args.metrics, os.path.basename(args.input), sys.argv, total)
    write_summary(args.summary, total, args)
    shutil.rmtree(tmp, ignore_errors=True)
    print(f"primary {total['primary']:,}; duplicates {total['duplicates']:,} "
          f"({total['duplicates'] / max(total['primary'], 1):.2%}); "
          f"propagated to {total['flag_secondary']:,} secondary and {total['flag_supplementary']:,} supplementary",
          file=sys.stderr)


if __name__ == "__main__":
    main()
