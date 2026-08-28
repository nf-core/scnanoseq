# nf-core/scnanoseq: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

TODO: Should here be added which output is cDNA/DNA specific?
TODO: Go over entire output section and remove/add flexiplex etc where needed.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Barcode Calling](#barcode-calling)
  - [Flexiplex](#flexiplex) - Barcode caller
  - [BLAZE](#blaze) - Barcode caller
- [Alignment](#alignment)
  - [Minimap2](#minimap2) - Long read alignment
- [Alignment Post-processing](#alignment-post-processing)
  - [Samtools](#samtools) - Sort and index alignments and make alignment qc
  - [Barcode Tagging Blaze](#barcode-tagging-blaze) - Barcode tagging with quality metrics and barcode information
  - [Barcode Tags Flexiplex](#barcode-tags-flexiplex) - The barcode and UMI tags flexiplex puts on every alignment
  - [UMI-tools Dedup](#umi-tools-dedup) - UMI-based Read deduplication
  - [Picard MarkDuplicates](#picard-markduplicates) - Read deduplication
- [Feature-Barcode Quantification](#feature-barcode-quantification)\*
  - [IsoQuant](#isoquant) - Feature-barcode quantification (gene and transcript level)
  - [oarfish](#oarfish) - Feature-barcode quantification (transcript-level only)
  - [Seurat](#seurat) - Feature-barcode matrix QC
- [Other steps](#other-steps)\*
  - [UCSC](#ucsc) - Annotation BED file
- [Quality Control](#quality-control)
  - [FastQC](#fastqc) - FASTQ QC
  - [NanoComp](#nanocomp) - Long Read FASTQ QC
  - [NanoPlot](#nanoplot) - Long Read FASTQ QC
  - [ToulligQC](#toulligqc) - Long Read FASTQ QC
  - [RSeQC](#rseqc) - Various RNA-seq QC metrics\*
  - [Read Counts](#read-counts) - Read Counts QC
  - [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

\* Indicates RNA only output

\*\* Indicates DNA only output

## Barcode Calling

### Flexiplex

<details markdown="1">
<summary>Output files</summary>

- `<sample_identifier>/`
  - `flexiplex/`
    - `*.barcodes_counts.txt` : This is a file containing each barcode and the counts of how many reads support it.
    - `*.known_barcodes` : This file is a list of all "true" barcodes and the counts associated to it in the sample. Can be used as whitelist for downstream tools.

</details>

[Flexiplex](https://github.com/DavidsonGroup/flexiplex/) is a fast, multithreaded, and user-configurable demultiplexer. Given a set of reads as either FASTQ or FASTA, it will demultiplex and/or identify a sequence of interest, reporting matching reads and read-barcode assignment. Flexiplex works in two modes: (i) when one or more sequences of interest are known, such as barcodes, and (ii) discovery mode—when only the sequence which flanks the region of interest is known.

### BLAZE

<details markdown="1">
<summary>Output files</summary>

- `<sample_identifier>/`
  - `blaze/`
    - `blaze/*.bc_count.txt` : This is a file containing each barcode and the counts of how many reads support it.
    - `blaze/*.knee_plot.png` : The knee plot detailing the ranking of each barcode.
    - `blaze/*.putative_bc.csv` : This file contains the naively detected barcode for each read.
    - `blaze/*.whitelist.csv` : This is a list of the "true" barcodes detected for a sample. The length of the file should roughly match the expected amount of cells that is expected for the sample.

</details>

![BLAZE - knee plot](images/blaze.png)

[BLAZE](https://github.com/shimlab/BLAZE) enables the accurate identification of barcodes and UMIs from Nanopore reads. The files produced by BLAZE can be used to assess the quality of the barcode calling and the data.

The knee plot (an example is listed above) that is provided by BLAZE shows all barcodes detected in a sample, ranked from highest to lowest read count. The "cliff-and-knee" shape (similar to the image above) is indicative of good quality. Deviations from this shape can be indicative of concerns with the data, such as low barcode counts. The `*.bc_count.txt` file can be used to accompany this figure to show every barcode and its abundance in the dataset.

## Alignment

### Minimap2

<details markdown="1">
<summary>Output files</summary>

- `<sample_identifier>/`
  - `genome`
    - `bam/`
      - `original/`
        - `*.sorted.bam` : The genome mapped and sorted bam.
        - `*.sorted.bam.bai` : The bam index for the genome mapped and sorted bam.
  - `transcriptome`
    - `bam/`
      - `original/`
      - `*.sorted.bam` : The transcriptome mapped and sorted bam.
      - `*.sorted.bam.bai` : The bam index for the transcriptome mapped and sorted bam.

</details>

[Minimap2](https://github.com/lh3/minimap2) is a versatile sequence alignment program that aligns DNA or mRNA sequences against a large reference database. Minimap2 is optimized for large, noisy reads making it a staple for alignment of nanopore reads.

## Alignment Post-processing

### Samtools

<details markdown="1">
<summary>Output files</summary>

- `<sample_identifier>/`
  - `genome/`
    - `bam/`
      - `mapped_only/`
        - `*.sorted.bam` : The genome aligned bam containing only reads that were able to be mapped.
        - `*.sorted.bam.bai` : The genome aligned bam index for the bam containing only reads that were able to be mapped.
    - `qc/`
      - `samtools/`
        - `minimap/`
          - `*.minimap.flagstat` : The flagstat file for the genome aligned bam obtained from minimap.
          - `*.minimap.idxstats` : The idxstats file for the genome aligned bam obtained from minimap.
          - `*.minimap.stats` : The stats file for the genome aligned bam obtained from minimap.
        - `mapped_only/`
          - `*.mapped_only.flagstat` : The flagstat file for the genome aligned bam containing only mapped reads.
          - `*.mapped_only.idxstats` : The idxstats file for the genome aligned bam containing only mapped reads.
          - `*.mapped_only.stats` : The stats file for the genome aligned bam containing only mapped reads.
        - `dedup/`
          - `*.dedup.flagstat` : The flagstat file for the genome aligned bam containing deduplicated umis.
  - `transcriptome/`
    - `bam/`
      - `mapped_only/`
        - `*.sorted.bam` : The transcriptome aligned bam containing only reads that were able to be mapped.
        - `*.sorted.bam.bai` : The transcriptome aligned bam index for the bam containing only reads that were able to be mapped.
    - `qc/`
      - `samtools/`
        - `minimap/`
          - `*.minimap.flagstat` : The flagstat file for the transcriptome aligned bam obtained from minimap.
          - `*.minimap.idxstats` : The idxstats file for the transcriptome aligned bam obtained from minimap.
          - `*.minimap.stats` : The stats file for the transcriptome aligned bam obtained from minimap.
        - `mapped_only/`
          - `*.mapped_only.flagstat` : The flagstat file for the transcriptome aligned bam containing only mapped reads.
          - `*.mapped_only.idxstats` : The idxstats file for the transcriptome aligned bam containing only mapped reads.
          - `*.mapped_only.stats` : The stats file for the transcriptome aligned bam containing only mapped reads.
        - `dedup/`
          - `*.dedup.flagstat` : The flagstat file for the transcriptome aligned bam containing deduplicated umis.

</details>

![MultiQC - samtools idxstats](images/samtools_idxstats.png)
![MultiQC - samtools stats](images/samtools_stats.png)

[Samtools](https://www.htslib.org/) is a suite of programs for reading, writing, editing, indexing, and viewing files that are in SAM, BAM, or CRAM format

### Barcode Tagging Blaze

<details markdown="1">
<summary>Output files</summary>

- `<sample_identifier>/`
  - `genome/`
    - `bam/`
      - `barcode_tagged/`
        - `*.tagged.bam` : The genome aligned bam containing tagged barcode and UMI metadata.
  - `transcriptome/`
    - `bam/`
      - `barcode_tagged/`
        - `*.tagged.bam` : The transcriptome aligned bam containing tagged barcode and UMI metadata.

</details>

Barcode tagging is a custom script which adds metadata to the BAM files with commonly used single-cell tags which can be useful for custom down stream analysis (e.g.: subsetting BAMs based on cell barcodes). Specifically the following tags are added:

```
barcode tag = "CR"
barcode quality tag = "CY"
corrected barcode tag = "CB"
UMI tag = "UR"
UMI quality tag = "UY"
```

Note that barcodes are corrected with the custom script, `correct_barcodes.py`.

### Barcode Tags Flexiplex

No separate tagging step is needed on the flexiplex path. Flexiplex writes the barcode
and UMI into the FASTQ header comment, and the pipeline runs `minimap2 -y`, which copies
that comment onto every alignment. The tags are:

```
CB   corrected cell barcode, or "-" when no known barcode matched
CR   cell barcode as observed in the read
UB   corrected UMI
UR   UMI as observed in the read
XB   CB when a barcode was called, otherwise CR
```

`XB` is derived by the pipeline rather than written by flexiplex. It is what
deduplication and the all-droplet quantification group on, so that a droplet which fell
below the knee — real barcode, just not on the known list — stays distinct instead of
pooling with every other unassigned read under `CB:Z:-`.

Since flexiplex is run with `-a true` it reports every read, not only the ones it could
match to the known barcode list, which is what makes `CR` — and therefore the uncalled
droplets — available at all. Reads where flexiplex found no barcode region whatsoever
are dropped during assignment: they belong to no droplet, called or empty, and carry a
one-character UMI placeholder rather than a real UMI. Every alignment downstream
therefore has a real `CR`, a real `XB` and a full-width UMI. Barcodes are corrected
during the flexiplex run itself and are not post-corrected.

### Gene assignment

<details markdown="1">
<summary>Output files</summary>

- `<sample>/<type>/genome/qc/gene_assignment/`
  - `*.gene_assignment.tsv` : The distribution of gene assignment statuses, per contig.

</details>

When `--dedup_per_gene` is enabled (the default for genome alignments), every alignment is tagged with the gene whose body it overlaps most before deduplication:

```
gene id tag     = "GX"
gene name tag   = "GN"
gene status tag = "GS"    # unique | ambiguous | none
```

Assignment is by gene-body overlap and is strand-agnostic. `GX` always holds the winning gene, even when a second gene also covers a meaningful share of the read, in which case `GS` is `ambiguous`. Only `unique` alignments are grouped by gene; the rest are deduplicated by position and merged back, so the `GS` distribution in the summary tsv is the ceiling on what gene grouping can act on.

### UMI-tools Dedup

<details markdown="1">
<summary>Output files</summary>

- `<sample_identifier>/`
  - `genome/`
    - `bam/`
      - `dedup_umitools/`
        - `*.dedup.bam` : The genome aligned bam containing corrected barcodes and deduplicated umis.
        - `*.dedup.bam.bai` : The genome aligned bam index for the bam containing corrected barcodes and deduplicated umis.
  - `transcriptome/`
    - `bam/`
      - `dedup_umitools/`
        - `*.dedup.bam` : The transcriptome aligned bam containing corrected barcodes and deduplicated umis.
        - `*.dedup.bam.bai` : The transcriptome aligned bam index for the bam containing corrected barcodes and deduplicated umis.

</details>

[UMI-Tools](https://umi-tools.readthedocs.io/en/latest/reference/dedup.html) deduplicate reads based on the mapping co-ordinate and the UMI attached to the read. The identification of duplicate reads is performed in an error-aware manner by building networks of related UMIs.

Users should note that `oarfish` requires input reads to be deduplicated. As a result, the `skip_dedup` option is only applicable to `IsoQuant`. By default, `scnanoseq` will perform deduplication for IsoQuant unless the `skip_dedup` option is explicitly enabled, while deduplication will always be executed for `oarfish` quantification.

### Picard MarkDuplicates

<details markdown="1">
<summary>Output files</summary>

- `<sample_identifier>/`
  - `genome/`
    - `bam/`
      - `dedup_picard/`
        - `*.dedup.bam` : The genome aligned bam containing corrected barcodes and deduplicated umis.
        - `*.dedup.bam.bai` : The genome aligned bam index for the bam containing corrected barcodes and deduplicated umis.
  - `transcriptome/`
    - `bam/`
      - `dedup_picard/`
        - `*.dedup.bam` : The transcriptome aligned bam containing corrected barcodes and deduplicated umis.
        - `*.dedup.bam.bai` : The transcriptome aligned bam index for the bam containing corrected barcodes and deduplicated umis.
  - `qc/`
    - `dedup/`
      - `*.metrics.txt` : The MarkDuplicates duplication metrics, also summarised in the MultiQC report. \*\*

</details>

[Picard MarkDuplicates](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard) locates and tags duplicate reads in a BAM or SAM file.

Users should note that `oarfish` requires input reads to be deduplicated. As a result, the `skip_dedup` option is only applicable to `IsoQuant`. By default, `scnanoseq` will perform deduplication for IsoQuant unless the `skip_dedup` option is explicitly enabled, while deduplication will always be executed for `oarfish` quantification.

## Feature-Barcode Quantification

### IsoQuant

<details markdown="1">
<summary>Output files</summary>

- `<sample_identifier>/`
  - `genome/`
    - `isoquant/`
      - `*.gene_counts.tsv` : The feature-barcode matrix from gene quantification.
      - `*.transcript_counts.tsv` : The feature-barcode matrix from transcript quantification.
    - `isoquant_all_droplets/` (flexiplex only)
      - `*.all_droplets.gene_counts.tsv` : As above, over every droplet rather than the called cells.
      - `*.all_droplets.transcript_counts.tsv` : As above, over every droplet rather than the called cells.

</details>

[IsoQuant](https://github.com/ablab/IsoQuant) is a tool for the genome-based analysis of long RNA reads, such as PacBio or Oxford Nanopores. IsoQuant allows to reconstruct and quantify transcript models with high precision and decent recall. If the reference annotation is given, IsoQuant also assigns reads to the annotated isoforms based on their intron and exon structure. IsoQuant further performs annotated gene, isoform, exon and intron quantification. The outputs of IsoQuant can be important for downstream analysis with tools specialized in single-cell/nuclei analysis (e.g.: `Seurat`).

In order to assist with the performance of IsoQuant, the inputs are split by chromosome to add a further degree of parallelization.

It should also be noted that IsoQuant can only accurately perform quantification on a **genome** aligned bam, and will produce both gene and transcript level matrices

When demultiplexing with flexiplex, IsoQuant is run twice off the same alignments:

- `isoquant/` groups on `CB`, so it covers the cells flexiplex matched to the known
  barcode list. This is the matrix to use for ordinary analysis, and the one Seurat QC
  is run against. Reads that matched no known barcode collect in a single `-` column,
  which should be dropped.
- `isoquant_all_droplets/` groups on `XB`, which falls back to the uncorrected barcode.
  Droplets below the knee called by `flexiplex-filter` therefore appear here under their
  own barcode alongside the called cells, which is what makes ambient/empty-droplet
  estimation possible. It has no `-` column: reads with no barcode region at all were
  dropped during assignment.

`oarfish` groups on `CB` too, so its transcript matrix gains the same droppable `-`
column.

### oarfish

<details markdown="1">
<summary>Output files</summary>

- `<sample_identifier>/`
  - `transcriptome/`
    - `oarfish/`
      - `barcodes.tsv.gz`
      - `features.tsv.gz`
      - `matrix.mtx.gz`

</details>

[oarfish](https://github.com/COMBINE-lab/oarfish) is a program, written in Rust (https://www.rust-lang.org/), for quantifying transcript-level expression from long-read (i.e. Oxford nanopore cDNA and direct RNA and PacBio) sequencing technologies. oarfish requires a sample of sequencing reads aligned to the transcriptome (currently not to the genome). It handles multi-mapping reads through the use of probabilistic allocation via an expectation-maximization (EM) algorithm.

It should also be noted that oarfish can only accurately perform quantification on a **transcript** aligned bam, and will only produce transcript level matrices. It's also recommended to ensure that the `--save_transcript_secondary_alignment` is enabled to produce the most accurate oarfish results (true by default for `oarfish` quantification). Notably, this can lead to much higher number of reads reported as aligned, however, this is expected behavior when secondary alignments are included in the analysis.

### lr-kallisto

<details markdown="1">
<summary>Output files</summary>

- `reference/`
  - `lrkallisto/`
    - `kallisto.idx` : The k=63 kallisto index.
    - `*.t2g.tsv` : The transcript-to-gene mapping.
    - `transcripts.fasta` : The transcript sequences extracted from the genome and GTF by `gffread`.
- `<sample_identifier>/`
  - `cdna/`
    - `lrkallisto/`
      - `gene/` : Gene-level feature-barcode matrix (`barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz`).
      - `transcript/` : Transcript-level feature-barcode matrix (`barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz`).
      - `quant/` : The raw `kallisto quant-tcc` output, including TPM-normalised matrices.
      - `counts_unfiltered/` : The transcript compatibility count (TCC) matrix and its equivalence classes.
      - `bus/` : The BUS file, equivalence class map, transcript names and `run_info.json` pseudoalignment statistics.
    - `lrkallisto_all_droplets/` : The same five outputs, over every droplet rather than the called cells.

</details>

[lr-kallisto](https://kallisto.readthedocs.io/en/latest/lr/pseudoalignment.html) is the long-read mode of [kallisto](https://github.com/pachterlab/kallisto). Rather than aligning reads, it pseudoaligns them against an index built with a longer k-mer than short-read kallisto uses (63 rather than 31), then quantifies transcript abundances with an expectation-maximization algorithm adapted to long-read error profiles. Because it does not align, this branch of the pipeline produces no BAM file and therefore no alignment-derived QC.

lr-kallisto is run twice off the same reads and the same index, mirroring the two
IsoQuant passes described above. It reads the barcode out of the flexiplex read, so this
quantifier requires `--demux_tool_cdna flexiplex` and both passes always run:

- `lrkallisto/` takes the barcode from the flexiplex read name, which is the one matched
  to the known barcode list, and `bustools correct` then drops the reads that matched
  nothing. This is the matrix to use for ordinary analysis, and the one Seurat QC is run
  against.
- `lrkallisto_all_droplets/` takes the barcode from the `XB` tag, which falls back to the
  uncorrected barcode. Droplets below the knee called by `flexiplex-filter` therefore
  appear here under their own barcode alongside the called cells, which is what makes
  ambient/empty-droplet estimation possible. Barcodes are not corrected in this pass:
  `XB` is already the final droplet identity, and correcting would discard the very
  droplets the matrix exists to keep. Sequencing errors in the below-knee barcodes
  therefore inflate the column count, which on a full sample can make this matrix much
  wider than `lrkallisto/`.

Barcodes and UMIs are taken from the flexiplex read names and written into a synthetic barcode read, so that kallisto can locate them positionally. UMI deduplication is performed by `bustools` during counting rather than by UMI-tools or Picard, so `--skip_dedup` and `--dedup_tool` do not apply to this path.

Both gene and transcript matrices come from a single pseudoalignment pass. `kallisto bus --long` writes the BUS file, `bustools` corrects barcodes against the flexiplex known-barcode list and collapses UMIs into a transcript compatibility count matrix, and `kallisto quant-tcc --long` resolves that into transcript abundances which are also aggregated to genes. Because gene counts are derived from the equivalence classes rather than from reads assigned to a single gene, reads that are compatible with transcripts of more than one gene are distributed by the EM rather than discarded.

`bus/run_info.json` is worth checking after a run: a low `p_pseudoaligned` usually indicates a barcode geometry or strandedness mismatch rather than poor data.

### Seurat

<details markdown="1">
<summary>Output files</summary>

- `<sample_identifier>/`
  - `genome/`
    - `qc/`
      - `gene/`
        - `*.csv`: A file containing statistics about the isoquant generated cell-read distribution for genes.
        - `*.png`: A series of qc images to determine the quality of the isoquant generated gene quantification.
      - `transcript/`
        - `*.csv`: A file containing statistics about the isoquant generated cell-read distribution for transcripts.
        - `*.png`: A series of qc images to determine the quality of the isoquant generated transcript quantification.
  - `transcriptome/`
    - `qc/`
      - `transcript/`
        - `*.csv`: A file containing statistics about the oarfish generated cell-read distribution for transcript.
        - `*.png`: A series of qc images to determine the quality of the oarfish generated transcript quantification.

</details>

[Seurat](https://satijalab.org/seurat/) is an R package designed for QC, analysis, and exploration of single-cell RNA-seq data.

## Other steps

### UCSC

<details markdown="1">
<summary>Output files</summary>

- `ucsc/`
  - `*.annotation.bed`: BED file format from input GTF
  - `*.annotation.genepred`: genepred file format from input GTF

</details>

[`ucsc-gtftogenepred` and `ucsc-genepredtobed`](https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/) are stand-alone applications developed by UCSC which, together, converts a GTF file the BED file format.

## Quality Control

### FastQC

<details markdown="1">
<summary>Output files</summary>

- `<sample_identifier>/`
  - `qc/`
    - `fastqc/`
      - `pre_trim/` and `post_trim/` and `post_extract/`
        - `*_fastqc.html`: FastQC report containing quality metrics.
        - `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

![MultiQC - FastQC sequence counts plot](images/mqc_fastqc_counts.png)

![MultiQC - FastQC mean quality scores plot](images/mqc_fastqc_quality.png)

![MultiQC - FastQC adapter content plot](images/mqc_fastqc_adapter.png)

:::note
The FastQC plots displayed in the MultiQC report shows _untrimmed_ reads. They may contain adapter sequence and potentially regions with low quality.
:::

### NanoComp

<details markdown="1">
<summary>Output files</summary>

- `batch_qcs/`
  - `nanocomp/`
    - `fastq/`
      - `NanoComp_*.log`: This is the log file detailing the nanocomp run.
      - `NanoComp-report.html` - This is browser-viewable report that contains all the figures in a single location.
      - `*.html`: Nanocomp outputs all the figures in the report as individual files that can be inspected separately.
      - `NanoStats.txt`: This file contains quality control statistics about the dataset.
  - `genome`
    - `nanocomp/`
      - `bam/`
        - `NanoComp_*.log`: This is the log file detailing the nanocomp run.
        - `NanoComp-report.html` - This is browser-viewable report that contains all the figures in a single location.
        - `*.html`: Nanocomp outputs all the figures in the report as individual files that can be inspected separately.
        - `NanoStats.txt`: This file contains quality control statistics about the dataset.
  - `transcriptome`
    - `nanocomp/`
      - `bam/`
        - `NanoComp_*.log`: This is the log file detailing the nanocomp run.
        - `NanoComp-report.html` - This is browser-viewable report that contains all the figures in a single location.
        - `*.html`: Nanocomp outputs all the figures in the report as individual files that can be inspected separately.
        - `NanoStats.txt`: This file contains quality control statistics about the dataset.

</details>

![Nanocomp](images/nanocomp_1.png)
![Nanocomp](images/nanocomp_2.png)

[NanoComp](https://github.com/wdecoster/nanocomp) compares multiple runs of long read sequencing data and alignments. It creates violin plots or box plots of length, quality and percent identity and creates dynamic, overlaying read length histograms and a cumulative yield plot

**Note**: Please note that `NanoComp` is enabled by default for FASTQ files but **disabled** for BAM files due to its high-memory usage. Users who are interesting in executing `NanoComp` for BAMs should set `skip_bam_nanocomp` to `false`.

### NanoPlot

<details markdown="1">
<summary>Output files</summary>

- `<sample_identifier>/`
  - `qc/`
    - `nanoplot/`
      - `pre_trim/` and `post_trim/` and `post_extract/`
        - `NanoPlot_*.log`: This is the log file detailing the nanoplot run
        - `NanoPlot-report.html` - This is browser-viewable report that contains all the figures in a single location.
        - `*.html`: Nanoplot outputs all the figures in the report as individual files that can be inspected separately.
        - `NanoStats.txt`: This file contains quality control statistics about the dataset.
        - `NanoStats_post_filtering.txt`: If any filtering metrics are used for nanoplot this will contain the differences. This is produced by default and should contain no differences from `NanoStats.txt` if the process was unmodified

</details>

![Nanoplot](images/nanoplot_1.png)
![Nanoplot](images/nanoplot_2.png)

[NanoPlot](https://github.com/wdecoster/NanoPlot) is a plotting tool for long read sequencing data and alignments.

### ToulligQC

<details markdown="1">
<summary>Output files</summary>

- `<sample_identifier>/`
  - `qc/`
    - `toulligqc/`
      - `pre_trim/` and `post_trim/` and `post_extract/`
        - `<sample_identifier>ToulligQC-report-<date>/`
          - `report.html`: This is browser-viewable report that contains all the figures in a single location.
          - `report.data`: A log file containing information about ToulligQC execution, environment variables and full statistics
          - `images/*`: This is folder containing all the individual images produced by ToulligQC

</details>

![ToulligQC](images/toulligqc_1.png)
![ToulligQC](images/toulligqc_2.png)

[ToulligQC](https://github.com/GenomiqueENS/toulligQC) is a post sequencing QC tool for Oxford Nanopore sequencers.

### RSeQC

<details markdown="1">
<summary>Output files</summary>

- `<sample_identifier>/`
  - `qc/`
    - `rseqc/`
      - `*.read_distribution.txt`: This file contains statistics noting the type of reads located within the dataset

</details>

![RSeQC](images/rseqc.png)

[RSeQC](https://rseqc.sourceforge.net/) package provides a number of useful modules that can comprehensively evaluate high throughput sequence data especially RNA-seq data.

### Read Counts

<details markdown="1">
<summary>Output files</summary>

- `batch_qcs/`
  - `read_counts/`
    - `read_counts.csv`: This file contains the read counts for each sample at various points in the pipeline. Each row is a different sample, and the columns are the amount of reads the sample contained at that point in the pipeline.

</details>

![Read Counts](images/read_counts.png)

Since flexiplex is run with `-a true` it passes through every read in which it found a
barcode region, whether or not that barcode matched the known list. `extracted_read_counts`
therefore no longer drops to the reads that were given a known barcode; it is
`trimmed_read_counts` minus the reads with no barcode region at all, which are discarded
during assignment. Read the barcode-calling rate from `corrected_read_counts`, which still
counts only reads assigned a known barcode.

This is a custom script written using BASH scripting. Its purpose is to report the amount of reads that are filtered out at steps in the pipeline that will result in filtered reads, such as barcode detection, barcode correction, alignment, etc. Elevated levels of filtering can be indicative of quality concerns.

For performance, this step parses the read counts from the output of either FastQC or NanoPlot rather than computing it. If the options `--skip_fastqc` and `--skip_nanoplot` or `--skip_qc` is used, this file will not be produced.

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
