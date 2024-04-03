# nf-core/scnanoseq: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Preprocessing](#preprocessing)
  - [Nanofilt](#nanofilt) - Read Quality Filtering and Trimming
- [Barcode Calling](#barcode-calling)
  - [BLAZE](#blaze) - Barcode caller
- [Alignment](#alignment)
  - [Minimap2](#minimap2) - Long read alignment
- [Alignment Post-processing](#alignment-post-processing)
  - [Samtools](#samtools) - Sort and index alignments and make alignment qc
  - [Barcode Tagging](#barcode-tagging) - Barcode tagging with quality metrics and barcode information
  - [Barcode Correction](#barcode-correction) - Barcode whitelist correction
  - [UMI Deduplication](#umi-deduplication) - UMI-based deduplication
- [Feature-Barcode Quantification](#feature-barcode-quantification)
  - [IsoQuant](#isoquant) - Feature-barcode quantification (gene and transcript level)
  - [Seurat](#seurat) - Feature-barcode matrix QC
- [Other steps](#other-steps)
  - [UCSC](#ucsc) - Annotation BED file
- [Quality Control](#quality-control)
  - [FastQC](#fastqc) - Fastq QC
  - [Nanocomp](#nanocomp) - Long Read Fastq QC
  - [Nanoplot](#nanoplot) - Long Read Fastq QC
  - [RSeQC](#rseqc) - Various RNA-seq QC metrics
  - [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

## Preprocessing
### Nanofilt
<details markdown="1">
<summary>Output files</summary>

- `<sample_identifier>/`
  - `fastq/`
    - `trimmed_nanofilt/`
      - `*_filtered.fastq.gz`: The post-trimmed fastq. By default this will be mostly quality trimmed.

</details>

[Nanofilt](https://github.com/wdecoster/nanocomp) is a tool used for filtering and trimming of long read sequencing data.

## Barcode Calling
### BLAZE
<details markdown="1">
<summary>Output files</summary>

- `<sample_identifier>/`
  - `blaze/`
    - `blaze/*.bc_count.txt` : This is a file containing each barcode and the counts of how many reads support it.
    - `blaze/*.knee_plot.png` : The knee plot detailing the ranking of each barcode.
    - `blaze/*.putative_bc.csv` : This file contains the naively detected barcode for each read.
    - `blaze/*.whitelist.csv` : This is the detected "true" barcodes for the dataset.

</details>

![BLAZE - knee plot](images/blaze.png)

[BLAZE](https://github.com/shimlab/BLAZE) enables the accurate identification of barcodes and UMIs from Nanopore reads.

## Alignment
### Minimap2

<details markdown="1">
<summary>Output files</summary>

- `<sample_identifier>/`
  - `bam/`
    - `original/`
      - `*.sorted.bam` : The mapped and sorted bam.
      - `*.sorted.bam.bai` : The bam index for the mapped and sorted bam.
</details>

[Minimap2](https://github.com/lh3/minimap2) is a versatile sequence alignment program that aligns DNA or mRNA sequences against a large reference database. Minimap2 is optimized for large, noisy reads making it a staple for alignment of nanopore reads

## Alignment Post-processing
### Samtools
<details markdown="1">
<summary>Output files</summary>

- `<sample_identifier>/`
  - `bam/`
    - `mapped_only/`
      - `*.sorted.bam` : The bam contaning only reads that were able to be mapped.
      - `*.sorted.bam.bai` : The bam index for the bam containing only reads that were able to be mapped.
  - `qc/`
    - `samtools/`
      - `minimap/`
        - `*.minimap.flagstat` : The flagstat file for the bam obtained from minimap.
        - `*.minimap.idxstats` : The idxstats file for the bam obtained from minimap.
        - `*.minimap.stats` : The stats file for the bam obtained from minimap.
      - `mapped_only/`
        - `*.mapped_only.flagstat` : The flagstat file for the bam containing only mapped reads.
        - `*.mapped_only.idxstats` : The idxstats file for the bam containing only mapped reads.
        - `*.mapped_only.stats` : The stats file for the bam containing only mapped reads.
      - `corrected/`
        - `*.corrected.flagstat` : The flagstat file for the bam containing corrected barcodes.
        - `*.corrected.idxstats` : The idxstat file for the bam containing corrected barcodes.
        - `*.corrected.stats` : The stat file for the bam containing corrected barcodes.
      - `dedup/`
        - `*.dedup.flagstat` : The flagstat file for the bam containing deduplicated umis.
        - `*.dedup.idxstats` : The idxstats file for the bam containing deduplicated umis.
        - `*.dedup.stats` : The stats file for the bam containing deduplicated umis.
</details>

![MultiQC - samtools idxstats](images/samtools_idxstats.png)
![MultiQC - samtools stats](images/samtools_stats.png)

[Samtools](https://www.htslib.org/) is a suite of programs for reading, writing, editing, indexing, and viewing files that are in SAM, BAM, or CRAM format

### Barcode Tagging
<details markdown="1">
<summary>Output files</summary>

- `<sample_identifier>/`
  - `bam/`
    - `barcode_tagged/`
      - `*.tagged.bam` : The bam containing tagged barcode and UMI metadata.

</details>

Barcode tagging is a custom script which adds metadata to the BAM files with commonly used single-cell tags which can be useful for custom down stream analysis (e.g.: subsetting BAMs based on cell barcodes). Specifically the following tags are added: 

```
barcode tag = "CR"
barcode quality tag = "CY"
UMI tag = "UR"
UMI quality tag = "UY"
```

Please see "Barcode Correction" below for metadata added post-correction.

### Barcode Correction
<details markdown="1">
<summary>Output files</summary>

- `<sample_identifier>/`
  - `bam/`
    - `corrected/`
      - `*.corrected.bam` : The bam containing corrected barcodes.
      - `*.corected.bam.bai` : The bam index for the bam containing corrected barcodes.

</details>

Barcode correction is a custom script that uses the whitelist generated by BLAZE in order to correct barcodes that are not on the whitelist into a whitelisted barcode. During this step, an additional BAM tag is added, `CB`, to indicate a barcode sequence that is error-corected.

### UMI Deduplication
<details markdown="1">
<summary>Output files</summary>

- `<sample_identifier>/`
  - `bam/`
    - `dedup/`
      - `*.dedup.bam` : The bam containing corrected barcodes and deduplicated umis.
      - `*.dedup.bam.bai` : The bam index for the bam containing corrected barcodes and deduplicated umis.

</details>

![MultiQC - umitools dedup](images/umitools_dedup.png)

[UMI-Tools](https://umi-tools.readthedocs.io/en/latest/reference/dedup.html) deduplicate reads based on the mapping co-ordinate and the UMI attached to the read. The identification of duplicate reads is performed in an error-aware manner by building networks of related UMIs

## Feature-Barcode Quantification
### IsoQuant
<details markdown="1">
<summary>Output files</summary>

- `<sample_identifier>/`
  - `isoquant/`
    - `*.gene_counts.tsv` : The feature-barcode matrix from gene quantification.
    - `*.transcript_counts.tsv` : The feature-barcode matrix from transcript quantification.

</details>

[IsoQuant](https://github.com/ablab/IsoQuant) is a tool for the genome-based analysis of long RNA reads, such as PacBio or Oxford Nanopores. IsoQuant allows to reconstruct and quantify transcript models with high precision and decent recall. If the reference annotation is given, IsoQuant also assigns reads to the annotated isoforms based on their intron and exon structure. IsoQuant further performs annotated gene, isoform, exon and intron quantification

### Seurat
<details markdown="1">
<summary>Output files</summary>

- `<sample_identifier>/`
  - `qc/`
    - `gene/`
      - `*.csv`: A file containing statistics about the cell-read distribution for genes.
      - `*.png`: A series of qc images to determine the quality of the gene quantification.
    - `transcript/`
      - `*.csv`: A file containing statistics about the cell-read distribution for transcript.
      - `*.png`: A series of qc images to determine the quality of the transcript quantification.
</details>

![MultiQC - seurat](images/seurat.png)

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
      - `pre_trim/`
        - `*_fastqc.html`: FastQC report containing quality metrics.
        - `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.
      - `post_trim/`
        - `*_fastqc.html`: FastQC report containing quality metrics.
        - `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.
      - `post_extract/`
        - `*_fastqc.html`: FastQC report containing quality metrics.
        - `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

![MultiQC - FastQC sequence counts plot](images/mqc_fastqc_counts.png)

![MultiQC - FastQC mean quality scores plot](images/mqc_fastqc_quality.png)

![MultiQC - FastQC adapter content plot](images/mqc_fastqc_adapter.png)

> **NB:** The FastQC plots displayed in the MultiQC report shows _untrimmed_ reads. They may contain adapter sequence and potentially regions with low quality.

### Nanocomp

<details markdown="1">
<summary>Output files</summary>

- `batch_qcs/`
  - `nanocomp/`
    - `fastq/` and `bam/`
      - `NanoComp_*.log`: This is the log file detailing the nanocomp run.
      - `NanoComp-report.html` - This is browser-viewable report that contains all the figures in a single location.
      - `*.html`: Nanocomp outputs all the figures in the report as individual files that can be inspected separately.
      - `NanoStats.txt`: This file contains quality control statistics about the dataset.

</details>

![Nanocomp](images/nanocomp_1.png)
![Nanocomp](images/nanocomp_2.png)

[Nanocomp](https://github.com/wdecoster/nanocomp) compares multiple runs of long read sequencing data and alignments. It creates violin plots or box plots of length, quality and percent identity and creates dynamic, overlaying read length histograms and a cumulative yield plot

### Nanoplot
<details markdown="1">
<summary>Output files</summary>

- `<sample_identifier>/`
  - `qc/`
    - `nanoplot/`
      - `pre_trim/` and `post_trim/` and `post_extract`
        - `NanoPlot_*.log`: This is the log file detailing the nanoplot run
        - `NanoPlot-report.html` - This is browser-viewable report that contains all the figures in a single location.
        - `*.html`: Nanoplot outputs all the figures in the report as individual files that can be inspected separately.
        - `NanoStats.txt`: This file contains quality control statistics about the dataset.
        - `NanoStats_post_filtering.txt`: If any filtering metrics are used for nanoplot this will contain the differences. This is produced by default and should contain no differences from `NanoStats.txt` if the process was unmodified

</details>

![Nanoplot](images/nanoplot_1.png)
![Nanoplot](images/nanoplot_2.png)

[Nanoplot](https://github.com/wdecoster/NanoPlot) is a plotting tool for long read sequencing data and alignments.

### RSeQC

<details markdown="1">
<summary>Output files</summary>

- `<sample_identifier>/`
  - `qc/`
    - `rseqc/`
      - `*.read_distribution.txt`: This file contains statisitics noting the type of reads located within the dataset

</details>

![RSeQC](images/rseqc.png)

[RSeQC](https://rseqc.sourceforge.net/) package provides a number of useful modules that can comprehensively evaluate high throughput sequence data especially RNA-seq data

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

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
