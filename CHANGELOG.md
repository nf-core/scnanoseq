# nf-core/scnanoseq: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### Enhancements

- Added support for processing single-cell/nuclei DNA samples alongside cDNA samples in the same pipeline run, via a new optional `type` column (`dna`/`cdna`) in the input samplesheet
- Added `flexiplex/discovery`, `flexiplex/filter` and `flexiplex/assign` modules for barcode extraction, filtering and assignment
- Added `demultiplex_flexiplex` and `demultiplex_blaze` subworkflows, replacing the previous inline demultiplexing steps in the main workflow
- Added `align_deduplicate_dna` subworkflow (`minimap2` alignment, `picard MarkDuplicates`, `BAM_SORT_STATS_SAMTOOLS`) for DNA sample processing
- Added `demux_tool_cdna`/`demux_tool_dna` parameters to select `flexiplex` or `blaze` for cDNA demultiplexing (DNA currently supports `flexiplex` only)
- Added `lrkallisto` as a third `--quantifier` option: an alignment-free path that pseudoaligns the demultiplexed FASTQ with lr-kallisto (`kallisto bus --long`), collapses UMIs with `bustools`, and quantifies genes and transcripts with `kallisto quant-tcc --long`. Adds a `test_lrkallisto` profile and the `gffread` module for transcript sequence extraction
- Bumped IsoQuant to v3.13.1
- Deduplication now groups UMIs by feature rather than by alignment position. Genome alignments are tagged with a gene by gene-body overlap (new `tag_genes` and `split_gene_status` modules, `GX`/`GN`/`GS` tags) and grouped with `umi_tools --per-gene`; reads with an ambiguous gene call or none are deduplicated positionally and merged back, so no reads are dropped. Transcriptome alignments use `umi_tools --per-contig`, where a contig is a transcript. **This substantially changes absolute UMI counts** — on the data this was developed against, total duplicate removal rose from 6.5% to ~34%, i.e. previous counts were inflated by roughly 1.6x. Relative expression is essentially unchanged (Spearman 0.998, no genes lost), so clustering and ratio-based analyses are unaffected, but UMIs per cell, library complexity and saturation estimates will differ from earlier releases. Revert with `--dedup_per_gene=false` / `--dedup_per_contig=false`
- `umi_tools dedup` now always takes the cell barcode from the corrected `CB` tag. Previously this applied only when `--demux_tool_cdna flexiplex` was set; the `blaze` path fell back to reading the raw, uncorrected barcode out of the read name, which fragmented cells and made the two demultiplexing paths incomparable
- Bumped flexiplex to v1.02.7 and run it with `-a true`, so it reports every read rather than only the reads it could assign. Reads from a droplet below the knee are no longer dropped for want of a _known_ barcode, and the uncorrected barcode comes back in `CR`, which is what makes such a droplet distinguishable from a read with no barcode region at all. v1.02.7 also fixes an out-of-bounds read in flexiplex's UMI handling that only triggers on UMI-less protocols, i.e. the DNA path
- Removed the `flexiformatter` module. It existed to move the barcode and UMI from the read name into BAM tags, which is no longer needed: flexiplex writes `CB`/`CR`/`UB`/`UR` into the FASTQ comment itself, and `minimap2 -y` now carries them onto every alignment. `-y` is applied only when demultiplexing with flexiplex, since a basecaller comment on any other FASTQ would be copied in verbatim and produce unparseable records
- Deduplication now groups on `XB`/`UB` rather than `CB`/`UR`. `XB` is a derived tag holding the corrected barcode for a called cell and the uncorrected one for a droplet flexiplex could not match to the known list, so uncalled droplets stay separate instead of pooling into a single `-` cell. `UB` is flexiplex's corrected UMI; `UR` holds the raw one. The `blaze` path writes both tags too, so the two demultiplexers stay comparable. DNA `picard MarkDuplicates` keys on `XB` for the same reason
- Reads in which flexiplex found no barcode region at all are now discarded during assignment rather than carried through the pipeline. They belong to no droplet, called or empty, and flexiplex gives them a one-character UMI placeholder instead of a real UMI. Dropping them at the FASTQ stage keeps them out of alignment entirely and means every alignment downstream has a real `CR`, a real `XB` and a full-width UMI
- Added a second IsoQuant run over the same alignments, published to `isoquant_all_droplets/`, grouping on `XB` so that droplets below the knee are quantified alongside the called cells. The existing `isoquant/` run is unchanged and still groups on `CB`

### Parameter changes

- `--whitelist` replaced by `--whitelist_dna` and `--whitelist_cdna`
- Bundled whitelists in `assets/` changed from `.zip` to `.gz`
- Added `--split_amount` parameter to control FASTQ/barcode splitting for parallelization
- Added `--kallisto_kmer_size` (default: `63`), `--kallisto_threshold` (default: `0.8`), `--kallisto_platform` (default: `ONT`) and `--kallisto_dlist` (default: `true`) for the `lrkallisto` quantifier
- Added `--dedup_per_gene` (default: `true`, genome alignments) and `--dedup_per_contig` (default: `true`, transcriptome alignments) to control UMI grouping. `--dedup_per_gene` requires `--gtf`

### Fixes

- Fixed the three `flexiplex` modules declaring `conda "${moduleDir}/environment.yml"` when the file on disk is `environment.yaml`, which left them unresolvable under `-profile conda`
- `split_gene_status` no longer emits an empty BAM. `umi_tools` aborts on one, and a contig with no annotation produces exactly that — 50 of 85 contigs on a real GRCh38 run
- Fixed `split_amount` parameter type coercion so it is read as an integer under the Nextflow v2 strict syntax parser
- Fixed `SAMTOOLS_INDEX_DEDUP` in the `dedup_umis` subworkflow indexing the `umi_tools` output channel unconditionally, which left the `--dedup_tool picard` branch indexing a channel that was never populated
- `--quantifier` now defaults to `null` rather than an empty string. `--quantifier` has been optional since DNA support was added, since DNA samples are aligned and deduplicated but never quantified, but the empty default was still checked against the parameter's own regex and failed it. Every DNA-only run therefore aborted at parameter validation with `"" does not match regular expression`, `-profile test_dna` included. The samplesheet-aware check that errors only when a `cdna` row is present and no quantifier was given is now what actually reports the problem. Added `tests/dna.nf.test` so the DNA-only path is covered by nf-test
- Removed four module includes from the `process_longread_scrna` subworkflow that were never invoked (`PICARD_MARKDUPLICATES`, `SAMTOOLS_FLAGSTAT_DEDUP`, `SAMTOOLS_INDEX_DEDUP`, `SAMTOOLS_FILTER_DEDUP`), along with the `conf/modules.config` selector for `SAMTOOLS_FILTER_DEDUP`, which matched no process. The underlying modules are all still used elsewhere; only the dead aliases are gone
- DNA samtools QC files are no longer named with a doubled type suffix. The `BAM_STATS_SAMTOOLS` prefix appended a literal `.dna` on top of the `meta.type` interpolation, so a DNA sample produced `<id>.dna.dna.flagstat`. They are now `<id>.dna.flagstat`/`.idxstats`/`.stats`. File contents are unchanged; MultiQC sample labels lose the duplicate suffix too
- Fixed two `tests/.nftignore` patterns that had silently stopped matching. `**/genome/isoquant/*.mtx` and its `isoquant_all_droplets` counterpart were added to keep the count matrices out of content checking, but the published path later gained two nested `<sample>.<type>/` levels, so the matrices were being md5-compared again. Also added the DNA dedup BAM and its index, which are not byte-reproducible between runs -- `samtools sort` ties are resolved differently across runs -- and whose cDNA equivalents were already ignored. Both files are still checked by name

## v1.3.0 [2026-06-26]

### Credits

Special thanks to a new contributor to scnanoseq:

- [Nick Youngblut](https://github.com/nick-youngblut)

### Enhancements

- [#94](https://github.com/nf-core/scnanoseq/issues/94) Strict syntax conversion: converted entire workflow to strict syntax and reorganized into `<name>/main.nf` directory structure
- [#93](https://github.com/nf-core/scnanoseq/issues/93) and [#55](https://github.com/nf-core/scnanoseq/issues/55) Upgraded IsoQuant from v3.6.1 to v3.13.0 and removed chromosome-splitting logic in the IsoQuant subworkflow due to improvements in IsoQuant; IsoQuant now processes all chromosomes in a single invocation
- [#65](https://github.com/nf-core/scnanoseq/pull/65) and [#61](https://github.com/nf-core/scnanoseq/issues/61) Replaced NanoFilt with Chopper for read filtering, with gzip compression of intermediate split FASTQ files to reduce disk usage
- [#87](https://github.com/nf-core/scnanoseq/issues/87) Added `--skip_blaze_demux` parameter to allow skipping BLAZE demultiplexing
- Added `SPLIT_SEQ` module using seqkit for splitting FASTQ files, replacing the previous split approach
- Updated `CAT_FASTQ` to the nf-core module which now supports compressed inputs
- Moved gzip compression to pre-extraction step to minimize uncompressed FASTQ footprint in work directories
- Upgraded nf-core template from 3.2.1 to 3.5.1

### Parameter changes

- `--skip_fastqc` default changed from `false` to `true` (disabled by default due to runtime issues with long-read data)
- `--skip_fastq_nanocomp` default changed from `false` to `true`
- `--skip_bam_nanocomp` default changed from `false` to `true`
- [#99](https://github.com/nf-core/scnanoseq/issues/99) `--skip_toulligqc` default changed from `false` to `true`
- `--skip_blaze_demux` added (default: `true`)

### Fixes

- [#83](https://github.com/nf-core/scnanoseq/issues/83) Isoquant count matrix headers are now compatible with Seurat (removed leading `#` from header line). Note that lastest version of IsoQuant also resolves this issue.

### Software dependencies

| Dependency | Old version | New version |
| ---------- | ----------- | ----------- |
| `IsoQuant` | 3.6.1       | 3.13.0      |
| `Chopper`  | -           | 0.10.0      |
| `NanoFilt` | 2.8.0       | removed     |

## v1.2.2 [2026-03-09]

### Fixes

- Fixed error that would occur when a large amount of transcripts were used to split a bam by replacing samtools_view with a custom script.

## v1.2.1 [2025-09-02]

### Fixes

- Fixed issue where `-resume` would not always cache IsoQuant steps, resulting in silently skipping chromosomes

## v1.2.0 [2025-06-09]

### Enhancements

- Fixed issue with linting.yml preventing automatic template PRs
- Upgraded nf-core template to 3.2.1
- Upgraded nanocomp nf-core module (no version change)
- (https://github.com/nf-core/scnanoseq/issues/44) All output files produced by isoquant are now produced in the results file
- (https://github.com/nf-core/scnanoseq/issues/45) Reference files are now accepted in .zip format
- (https://github.com/nf-core/scnanoseq/issues/47) BLAZE scripts has been removed from the repo so the actual published code can be used
- (https://github.com/nf-core/scnanoseq/issues/47) Added new whitelists for 10X 3v4 and 10X 5v3
- (https://github.com/nf-core/scnanoseq/issues/56) Fixed an error where using `--skip_dedup` would end the pipeline early
- (https://github.com/nf-core/scnanoseq/issues/58) Fixed UMI length for 5 prime chemistries
- Fixed an error caused by `--skip_qc` and `--skip_seurat`
- Seurat process now places the seurat object to pipeline outputs
- No longer output uncorrected correct barcodes
- Updated metro diagram

### Software dependencies

| Dependency | Old version | New version |
| ---------- | ----------- | ----------- |
| `BLAZE`    | 2.2.0       | 2.5.1       |

## v1.1.0 [2025-03-18]

### Enhancements

- Inputs for IsoQuant are split on chromosome to allow for faster processing
- The read counts QC metric is now able to leverage NanoPlot counts if FastQC is skipped
- Added `oarfish` as an option for quantification
- Added `picard markdupes` as an option for deduplication

### Fixes

- The 'Post Trim Read QC' and 'Post Extract Read QC' nodes on the metro diagram have been placed in correct locations (closes issue #36)
- The BLAZE process in the example config has been corrected to use cpus instead of `--threads`

### Software dependencies

| Dependency | Old version | New version |
| ---------- | ----------- | ----------- |
| `IsoQuant` | 3.5.0       | 3.6.1       |
| `MultiQC`  | 1.25        | 1.25.1      |

## v1.0.0 [2024-10-07]

Initial release of nf-core/scnanoseq, created with the [nf-core](https://nf-co.re/) template.
