# nf-core/scnanoseq: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
