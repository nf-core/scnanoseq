# DNA deduplication: does the cDNA soft-clip problem apply here?

A rendered version of this write-up, with a diagram of the read geometry behind §2, is at
<https://claude.ai/code/artifact/02851b03-fcd0-43e9-a93a-5f6d4351b82b>.

## The question

The cDNA branch had a deduplication bug: `umi_tools` bundles reads on a **soft-clip-corrected**
position, and on ONT cDNA the clip length tracks untrimmed adapter and polyA rather than a real
fragment boundary, so members of one PCR family landed in different position groups even when they
aligned at an identical coordinate. Duplicate removal was recovering roughly a seventh of the
duplicates actually present. The fix replaced the position key with a **feature** key — `--per-gene`
on a new `GX` tag for genome alignments, `--per-contig` for the transcriptome.

The DNA branch deduplicates with `picard MarkDuplicates` instead
(`subworkflows/local/align_deduplicate_dna/main.nf`). Two questions follow:

1. Does MarkDuplicates suffer the same soft-clip problem?
2. Should the feature-grouping fix be ported to it?

**No, and no.** MarkDuplicates does apply the same soft-clip correction, but on ONT DNA the variable
adapter clip sits at the end it ignores. And porting the feature-grouping fix would be destructive,
because the DNA library has no UMI — so the position is not a flawed proxy for molecular identity,
it is the _only_ identifier there is.

## Reproducing

```bash
python3 analysis/dna_dedup/dedup_qc.py \
    --bam <sample>.bam --auto-windows --window-size 800000 \
    --out report.txt --json report.json
```

`--auto-windows` places one window at the midpoint of each of chr1–22 and chrX. The script rebuilds
the MarkDuplicates key from the CIGAR and flags, and **exits non-zero unless that key explains
≥99.9% of the flag-1024 marks already in the BAM** — every comparison below depends on the key being
right, so it is checked rather than assumed.

Numbers below are from `B2233535_dedup_qc.txt` in this directory:
`SPLONGGET_pipeline_runs/old/B2233535/.../dna/bam/dedup/B2233535.bam`, 23 × 800 kb windows,
**3,991,957 primary reads**, key agreement **99.998%**. Duplicate rate 22.72%, sets of median size 2
(max 9), all `DT:LB` — no optical duplicates, as expected for ONT.

Caveat: B2233535 was the only DNA sample available when this was written, so the figures describe one
library. The two structural conclusions (§2 and §3) follow from library chemistry and the absence of
a UMI and will hold for any run of this assay; the percentages in §4–§6 should be re-measured on a
new sample before being quoted for it. Now that the duplication metrics are published per run
(`<sample>/dna/qc/dedup/*.metrics.txt`), the duplicate rate at least is available without re-running
this script.

## 1. MarkDuplicates does key on soft-clip-corrected coordinates

For single-end ("fragment") reads the key is

```
(library, barcode, refIndex, 5' UNCLIPPED coordinate, strand)
```

`getUnclippedStart()` for forward reads, `getUnclippedEnd()` for reverse — **the identical
correction `umi_tools` applies**. The fragment's other end and its length are not part of the key.
So the mechanism that broke the cDNA branch is present here.

## 2. But the variable adapter clip is on the end MarkDuplicates ignores

This is a Tn5-tagmented 10x Multiome ATAC library. Flexiplex strips the 5' barcode adapter, and the
untrimmed Nextera read-through stays at the 3' end. Resolving each clip to the read's own
orientation (the left end in reference coordinates for a forward alignment, the right end for a
reverse one):

| soft clip                                | median | p90 | p99 | mean | share > 20 bp | share == 0 |
| ---------------------------------------- | ------ | --- | --- | ---- | ------------- | ---------- |
| **5' end — MarkDuplicates keys on this** | 0      | 0   | 266 | 17.9 | 2.5%          | **93.8%**  |
| 3' end — MarkDuplicates ignores this     | 69     | 92  | 357 | 88.0 | 98.4%         | 0.7%       |

The two ends could hardly be more different: essentially every read carries an adapter clip, and
essentially none of it lands on the anchor. Re-deriving the key with a **clipped** 5' anchor instead
of the unclipped one moves the duplicate rate 22.72% → 22.77%, a difference of **0.05%**. The share
of zero-length 5' clips only drifts from 94.6% to 90.4% across aligned-length strata from under
200 bp to over 5 kb, so this is not an artefact of pooling short and long reads.

**The cDNA failure mode does not transfer.** The soft-clip correction is real but lands on a clean
coordinate — the Tn5 insertion site.

## 3. Porting the feature-grouping fix would destroy the data

`umi_tools --per-gene` works because the UMI still separates molecules _within_ a feature; the
feature only replaces the unreliable position. The DNA branch has no such fallback:

- flexiplex is configured `-b '????????????????'` with **no `-u` UMI pattern**
  (`conf/modules.config`, the `DEMULTIPLEX_FLEXIPLEX_DNA` blocks) — 10x Multiome ATAC has no UMI by
  design;
- the BAM carries `CB` but **no `UB`/`UR`/`UY` tag** (confirmed across all 4.0 M reads);
- read names are `BARCODE_#readid`, with the UMI field empty.

So position _is_ the molecular identifier. Replacing it with a feature collapses every read in a
(cell, feature) group to one:

| grouping key                                   | duplicate rate | reads kept of 3,991,957 |
| ---------------------------------------------- | -------------- | ----------------------- |
| **(CB, ref, 5' unclipped, strand) — current**  | **22.72%**     | **3,084,864**           |
| (CB, ref, 5' // 10 bp bin)                     | 24.11%         | 3,029,355               |
| (CB, ref, 5' // 100 bp bin)                    | 27.63%         | 2,888,797               |
| (CB, ref, 5' // 1 kb bin)                      | 36.98%         | 2,515,890               |
| (CB, ref, 5' // 10 kb bin)                     | 58.31%         | 1,664,288               |
| (CB, ref, 5' // 100 kb bin — gene-sized)       | 82.78%         | 687,337                 |
| (CB, ref) — whole 800 kb window as one feature | 95.96%         | 161,263                 |

A gene-sized feature would discard ~83% of reads, and a real `--per-contig` grouping spans an entire
chromosome, so it would collapse at least as much as the 95.96% row. **Do not port the fix.**

## 4. The DNA-appropriate variant is jitter tolerance, and it is worth about 1%

ONT indel and alignment wobble does move the tagmentation-site coordinate. Taking read pairs that
share `(CB, contig, strand)` and whose 3' ends agree within 20 bp — i.e. that look like the same
molecule — and histogramming their 5'-anchor offset (n = 1,152,914):

| 5'-anchor offset                | share     |
| ------------------------------- | --------- |
| 0 bp — caught by MarkDuplicates | 91.55%    |
| **1–10 bp — missed**            | **5.97%** |
| 11–50 bp                        | 0.96%     |
| 51–200 bp                       | 1.52%     |

The distribution decays steeply over the first few bases (1 bp: 31,579; 2 bp: 14,349; 3 bp: 7,487),
the classic ONT signature, and is flat across aligned-length strata (5.6–6.3% missed), so it is
basecall/alignment noise rather than a property of a particular read class.

Clustering with tolerance at both ends instead of exact matching at one:

| tolerance                       | duplicate rate | vs current |
| ------------------------------- | -------------- | ---------- |
| 5' ±5, 3' ±50                   | 22.82%         | +4,065     |
| 5' ±10, 3' ±50                  | 22.95%         | +9,169     |
| 5' ±10, 3' ±200                 | 23.35%         | +25,102    |
| 5' ±0, 3' ±20 (strict both-end) | 21.48%         | −49,574    |

So jitter costs about **+0.5% to +1.0% relative** in missed duplicates (0.1–0.2% of reads).

## 5. There is an opposite bias of similar size

Because only the 5' coordinate is keyed, fragment length is ignored entirely — and aligned length
here is broad (median 347 bp, p75 1,033, p90 2,533, max 13,489). Within the 722,595 duplicate sets:

- **6.76%** contain members whose 3' ends differ by more than 20 bp; 5.07% by more than 100 bp;
- **3.08%** have a max/min aligned-length ratio above 2× (1.22% above 5×, 0.55% above 10×).

Requiring both ends to agree exactly gives 20.29% versus 22.72%, so roughly **2% of marked
duplicates are probably distinct molecules** that happen to share a tagmentation site within one
cell. Collapsing a 13 kb read into a 200 bp one because they start at the same base is a real loss
of information on a long-read platform.

**This bias points the opposite way to §4 and is of comparable magnitude, so the two largely
cancel.** Recommendation: leave the algorithm alone. A jitter-tolerant deduplicator would need to
become both-end-aware at the same time to avoid trading one bias for a larger one, and the net
correction is around 1% of reads — not worth a custom implementation and the risk it carries.

## 6. What is working correctly

**`--BARCODE_TAG CB` is essential and effective.** Dropping the barcode from the key gives a 50.18%
duplicate rate instead of 22.72%: 19.74% of coordinate-identical read pairs come from _different_
cells sharing a Tn5 insertion site, and the barcode correctly keeps them apart. Do not remove it.

**Barcode correction is clean.** Across all 243,746 coordinate-identical read pairs with differing
`CB`, there are **zero** at Hamming distance 1. An uncorrected barcode error would split a real
duplicate family and show up exactly here; it does not happen. (d2–d4 are enriched ~20–40× over a
random-pair background but amount to 1,135 pairs, i.e. 0.5% of the differing-CB pairs — negligible.)

**Per-cell duplicate rate rises with depth** — 19.1% at 11–50 reads per cell in the sampled windows,
22.6% at 201–1,000, 23.1% above 1,000 — the direction PCR saturation predicts.

## 7. One caveat for downstream consumers

MarkDuplicates does **not** propagate the duplicate flag to secondary or supplementary records of a
duplicate read. Verified across the full sample: of 387,019 secondary and 199,278 supplementary
records, **zero** carry flag 1024. minimap2 runs with `--secondary=yes`, so those records are
12.8% of the BAM.

The DNA branch marks duplicates without removing them, by design. Anyone filtering the published BAM
should therefore use

```bash
samtools view -b -F 0xD04 <sample>.bam    # unmapped + secondary + supplementary + duplicate
```

Filtering on `-F 0x400` alone leaves the secondary and supplementary records of removed primaries
behind as orphans.

## Summary

|                                                           | verdict                                                                                             |
| --------------------------------------------------------- | --------------------------------------------------------------------------------------------------- |
| Does MarkDuplicates key on soft-clip-corrected positions? | Yes, same as `umi_tools`                                                                            |
| Does that hurt the DNA branch?                            | No — 0.05% effect; the adapter clip is on the ignored 3' end                                        |
| Should the per-gene/per-contig fix be ported?             | **No** — no UMI exists, so it would discard 83–96% of reads                                         |
| Is anything else wrong with the dedup?                    | Two ~1% biases in opposite directions; not worth fixing                                             |
| Barcode handling                                          | Correct and essential                                                                               |
| Action taken                                              | Publish the duplication metrics so this is visible per run without re-deriving it from a 585 GB BAM |
