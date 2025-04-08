# nf-core/scnanoseq: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.2.0 [TBD]

### Enhancements

- (https://github.com/nf-core/scnanoseq/issues/44) All output files produced by isoquant are now produced in the results file
- (https://github.com/nf-core/scnanoseq/issues/45) Reference files are now accepted in .zip format
- (https://github.com/nf-core/scnanoseq/issues/47) BLAZE scripts has been removed from the repo so the actual published code can be used
- (https://github.com/nf-core/scnanoseq/issues/47) Added new whitelists for 10X 3v4 and 10X 5v3
- Fixed an error caused by `--skip_qc` and `--skip_seurat`

### Software dependencies

| Dependency | Old version | New version |
| ---------- | ----------- | ----------- |
| `BLAZE`    | 2.2.0       | 2.5.1       |

## v1.1.0 [2024-03-18]

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
