# nf-core/scnanoseq: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/scnanoseq/usage](https://nf-co.re/scnanoseq/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 4 columns, and a header row as shown in the examples below.

```bash
--input '[path to samplesheet file]'
```

### Full samplesheet

The example `samplesheet.csv` below contains a single FASTQ file per biological replicate with sample specific cell counts.

```csv title="samplesheet.csv"
sample,fastq,cell_count,type
CONTROL_REP1,AEG588A1_S1.fastq.gz,5000,cdna
CONTROL_REP2,AEG588A2_S1.fastq.gz,6000,cdna
TREATMENT_REP1,AEG588A4_S1.fastq.gz,5500,cdna
TREATMENT_REP2,AEG588A5_S1.fastq.gz,6000,cdna
CONTROL_REP1,AEG588A1_S1.fastq.gz,5000,dna
CONTROL_REP2,AEG588A2_S1.fastq.gz,6000,dna
TREATMENT_REP1,AEG588A4_S1.fastq.gz,5500,dna
TREATMENT_REP2,AEG588A5_S1.fastq.gz,6000,dna
```

| Column       | Description                                                                                                                                                                            |
| ------------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`     | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (`_`). |
| `fastq`      | Full path to FastQ file for Oxford Nanopore. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                                    |
| `cell_count` | Expected number of cells/nuclei. This value is used by the barcode calling tool (BLAZE and/or Flexiplex) as a baseline when determining an acceptable number of detected barcodes.     |
| `type`       | An optional column specifiying whether the sample is DNA or cDNA. If omitted, the default `cdna` is used.                                                                              |

Note: DNA samples are only compatible with `flexiplex` demultiplexing.

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

### Multiple runs of the same sample

The `sample` identifiers have to be the same when you have re-sequenced the same sample more than once e.g. to increase sequencing depth. The pipeline will concatenate the raw reads before performing any downstream analysis. Below is an example for the same sample sequenced across replicates 1 and 4 (`REP1` and `REP4` respectively):

```csv title="samplesheet.csv"
sample,fastq,cell_count,type
CONTROL_REP1,AEG588A1_S1.fastq.gz,5000,cdna
CONTROL_REP1,AEG588A1_S2.fastq.gz,5000,cdna
CONTROL_REP2,AEG588A2_S1.fastq.gz,2000,cdna
CONTROL_REP3,AEG588A3_S1.fastq.gz,7500,cdna
CONTROL_REP4,AEG588A4_S1.fastq.gz,9000,cdna
CONTROL_REP4,AEG588A4_S2.fastq.gz,9000,cdna
CONTROL_REP4,AEG588A4_S3.fastq.gz,9000,cdna
```

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/scnanoseq \
  --input ./samplesheet.csv \
  --outdir ./results \
  --genome_fasta /path/to/genome.fa \
  --transcript_fasta /path/to/transcriptome.fa \
  --gtf /path/to/file.gtf \
  --quantifier "isoquant,oarfish" \
  --demux_tool_cdna flexiplex \
  --demux_tool_dna flexiplex \
  --barcode_format 10X_3v3 \
  -profile <docker/singularity/institute>
```

Please note that while the above command specifies both transcriptome and genome fasta files, only one is needed for the pipeline and is dependent on which quantifier you wish to use. Isoquant requires a genome fasta, while oarfish requires a transcript fasta. Furthermore, if you have any DNA samples, the `genome_fasta` is required.
Additionally, for the `quantifier` parameter in the above command, we've listed the quantifiers as a comma-delimited string. It is possible to only use one quantifier, and can be accomplished by just providing the name of the quantifying tool you wish to run as a single value, i.e. providing `oarfish` if you only wish to run `oarfish`.

### Quantifiers

Three quantifiers are available through `--quantifier`, and any combination of them can be requested as a comma-delimited string:

| Quantifier   | Reference required             | Levels             | Notes                                                                             |
| ------------ | ------------------------------ | ------------------ | --------------------------------------------------------------------------------- |
| `isoquant`   | `genome_fasta` + `gtf`         | gene + transcript  | Quantifies from a genome alignment                                                |
| `oarfish`    | `transcript_fasta`             | transcript         | Quantifies from a transcriptome alignment; always deduplicates                    |
| `lrkallisto` | `genome_fasta` + `gtf`         | gene + transcript  | Pseudoaligns straight from the FASTQ; no alignment step                            |

#### lr-kallisto

[lr-kallisto](https://kallisto.readthedocs.io/en/latest/lr/pseudoalignment.html) is the long-read mode of kallisto (`kallisto --long`). Unlike the other two quantifiers, it does not align reads: it pseudoaligns them against a k=63 index and quantifies with a long-read expectation-maximisation step. That makes it considerably faster, at the cost of producing no BAM file and therefore no alignment-based QC for this branch of the pipeline.

Because lr-kallisto locates the cell barcode and UMI positionally rather than by name, the pipeline moves them out of the flexiplex read name into a synthetic barcode read before handing the FASTQs to `kallisto bus`. This reuses the barcodes that flexiplex already corrected against the whitelist, so cell barcodes are directly comparable with the other quantifiers. UMI deduplication is performed by `bustools` while counting, so `--dedup_tool` and `--skip_dedup` have no effect on this path.

This means lr-kallisto has some extra requirements:

- `--genome_fasta` and `--gtf` must both be provided. The transcript sequences are extracted from them with `gffread`, the transcript-to-gene mapping is derived from the same GTF, and the genome is used as the kallisto d-list so that reads originating outside the transcriptome are discarded rather than misassigned. Deriving both from one GTF is what guarantees the index sequence names and the mapping agree.
- `--demux_tool_cdna` must be `flexiplex`; the `blaze` path keeps barcodes in a separate per-read table rather than in the read name.
- `--custom_flexiplex_barcode_cdna` cannot be used, because the barcode and UMI widths are needed in order to locate them. Use one of the `--barcode_format` presets instead.

The k-mer length, mapping threshold, platform and d-list can be tuned with `--kallisto_kmer_size` (default 63), `--kallisto_threshold` (default 0.8), `--kallisto_platform` (`ONT` or `PacBio`, default `ONT`) and `--kallisto_dlist` (default `true`).

> [!WARNING]
> To turn the d-list off, write `--kallisto_dlist=false` with an equals sign. `params.kallisto_dlist` defaults to a boolean, so Nextflow treats the space-separated `--kallisto_dlist false` as a bare flag and drops the `false` as a positional argument: the d-list stays **on** and the only sign is a `positional argument` warning in the log. The same applies to any other boolean parameter you want to set to `false` on the command line.

> [!IMPORTANT]
> Check `p_pseudoaligned` in `bus/run_info.json` after your first run. The two defaults below are both what lr-kallisto recommends, and on high-error reads each one costs sensitivity independently:
>
> | | `--kallisto_dlist=false` | `--kallisto_dlist=true` (default) |
> | ------------------------------------ | ------------------------ | --------------------------------- |
> | `--kallisto_kmer_size 31` | 20.3% | 4.1% |
> | `--kallisto_kmer_size 63` (default) | 4.5% | **3.4%** |
>
> Measured on the chr21 `test_lrkallisto` data, where minimap2 aligns 21.7% of the same reads to the same transcriptome.
>
> - **k-mer length.** Pseudoalignment needs *exact* k-mer matches, so k interacts strongly with per-base error rate. A 63-mer survives only if all 63 bases are called correctly: at 1% error that is roughly half of all k-mers, at 8% error well under 1%. The lr-kallisto authors report high accuracy below roughly 10% error and reduced performance above it, so k=63 is the right default for modern high-accuracy basecalls but can collapse on legacy chemistries.
> - **d-list.** Using the genome as a d-list discards reads containing k-mers that are in the genome but not the transcriptome, which is what stops intronic and genomic reads being force-assigned to a transcript. This is far more punitive for long reads than for short ones: a long read only has to contain a single genome-only k-mer to be rejected.
>
> `--kallisto_threshold` will not recover these reads. It bounds the fraction of unmapped k-mers per read, and reads lost here typically have no mapping k-mers at all — raising it from 0.8 to 0.99 changed nothing in the measurements above.
>
> The defaults are deliberately the upstream-recommended, higher-specificity ones. If your mapping rate is far below the rate you get from minimap2 on the same data, try relaxing them one at a time and compare against another quantifier before trusting the counts.

Reads are pseudoaligned as non-strand-specific, because after flexiplex has trimmed the adapter, barcode, UMI and poly-T the remaining cDNA is not in a consistent orientation. `params.stranded` therefore does not apply to this path. If you do want to force a strand, override the `KALLISTO_BUS` arguments with a custom config:

```groovy title="custom.config"
process {
    withName: '.*:QUANTIFY_SCRNA_LRKALLISTO:KALLISTO_BUS' {
        ext.args = { "--threshold ${params.kallisto_threshold} --fr-stranded" }
    }
}
```

> [!NOTE]
> Building a k=63 index is memory-hungry. For a human genome, expect the `KALLISTO_INDEX` process to need substantially more memory than the `process_high` default; see [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources).

The pipeline supports barcode identification and extraction through both `flexiplex` and `blaze` and can be set through `demux_tool_dna` (only works with `flexiplex` for now) and `demux_tool_cdna` parameters. The barcode format can be specified through the `barcode_format` parameter. When working with completely custom barcode structures, you can additionally specify these with `custom_flexiplex_barcode_dna` and `custom_flexiplex_barcode_cdna` parameters. Note: ensure that you are using `flexiplex` as the barcode calling tool. This can be a string formatted as follows `"-x CTACACGACGCTCTTCCGATCT -b ???????????????? -u ?????????? -x TTTCTTATATGGG -f 8 -e 2"`, for more information check the documentation: https://davidsongroup.github.io/flexiplex/

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/scnanoseq -profile <docker/singularity/institute> -params-file params.yaml
```

with

```yaml title="params.yaml"
input: "./samplesheet.csv"
outdir: "./results/"
genome_fasta: "/path/to/genome.fa"
transcript_fasta: "/path/to/transcript.fa"
gtf: "/path/to/file.gtf"
quantifier: "isoquant|oarfish|lrkallisto|isoquant,oarfish"
barcode_format: "10X_3v3"
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/scnanoseq
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/scnanoseq releases page](https://github.com/nf-core/scnanoseq/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://charliecloud.io/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website. We have also provided a pipeline specific example of a custom configuration file in the [Introduction page](https://nf-co.re/scnanoseq/latest/#troubleshooting).

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

## Troubleshooting

If you experience any issues, please make sure to reach out on the [#scnanoseq slack channel](https://nfcore.slack.com/archives/C03TUE2K6NS) or [open an issue on our GitHub repository](https://github.com/nf-core/scnanoseq/issues/new/choose). However, some resolutions for common issues will be noted below:

- Due to the nature of the data this pipeline analyzes, some tools may experience increased runtimes. For some of the custom tools made for this pipeline (`preextract_fastq.py` and `correct_barcodes.py`), we have leveraged the splitting done via the `split_amount` parameter to decrease their overall runtimes. The `split_amount` parameter will split the input FASTQs into a number of FASTQ files, each containing a number of lines based on the value used for this parameter. As a result, it is important not to set this parameter to be too low as doing so would cause the creation of a large number of files the pipeline will be processed. While this value can be highly dependent on the data, a good starting point for an analysis would be to set this value to `500000`. If you find that `PREEXTRACT_FASTQ` and `CORRECT_BARCODES` are still taking long amounts of time to run, it would be worth reducing this parameter to `200000` or `100000`, but keeping the value on the order of hundred of thousands or tens of thousands should help with keeping the total number of processes minimal. An example of setting this parameter to be equal to 500000 is shown below:

```yml title="params.yml"
split_amount: 500000
```

- We have seen a recurrent node failure on slurm clusters that does seem to be related to submission of Nextflow jobs. This issue is not related to this pipeline per se, but rather to Nextflow itself. We are currently working on a resolution. But we have two methods that appear to help overcome should this issue arise:
  1. Provide a custom config that increases the memory request for the job that failed. This may take a couple attempts to find the correct requests, but we have noted that there does appear to be a memory issue occasionally with these errors.
  2. Request an interactive session with a decent amount of time and memory and CPUs in order to run the pipeline on the single node. Note that this will take time as there will be minimal parallelization, but this does seem to resolve the issue.
- We note that umitools dedup can take a large amount of time in order to perform deduplication. One approach we have implemented to assist with speed is to split input files based on chromosome. However for the transcriptome aligned bams, there is some additional work required that involves grouping transcripts into appropriate chromosomes. In order to accomplish this, the pipeline needs to parse the transcript id from the transcriptome FASTA file. The transcript id is often nested in the sequence identifier with additional data and the data is delimited. We have included the delimiters used by reference files obtained from GENCODE, NCBI, and Ensembl. However in case you wish to explicitly control this or if the reference file source uses a different delimiter, you are able to manually set it via the `--fasta_delimiter` parameter.
- We acknowledge that analyzing PromethION data is a common use case for this pipeline. Currently, the pipeline has been developed with defaults to analyze GridION and average sized PromethION data. For cases, where jobs have fail due for larger PromethION datasets, the defaults can be overwritten by a custom configuation file (provided by the `-c` Nextflow option) where resources can be increased (substantially in some cases). Below are some of the overrides we have used, and while these amounts may not work on every dataset, these will hopefully at least note which processes will need to have their resources increased:

```groovy title="custom.config"

process
{
    withName: '.*:.*FASTQC.*'
    {
        cpus = 20
    }
}

process
{
    withName: '.*:BLAZE'
    {
        cpus = 30
    }
}

process
{
    withName: '.*:TAG_BARCODES'
    {
        memory = '60.GB'
    }
}

process
{
    withName: '.*:SAMTOOLS_SORT'
    {
        cpus = 20
    }
}

process
{
    withName: '.*:MINIMAP2_ALIGN'
    {
        cpus = 20
    }
}

process
{
    withName: '.*:ISOQUANT'
    {
        cpus = 30
        memory = '85.GB'
    }
}
```

We further note that while we encourage the use of `split_amount` as discussed above for larger datasets, the pipeline can be executed without enabling this parameter. When doing this, please consider increasing the time limit to `CORRECT_BARCODES` as it can take hours instead of minutes when `split_amount` is disabled:

```groovy title="custom.config"
//NOTE: with split_amount disabled, consider increasing the time limit to CORRECT_BARCODES
process
{
    withName: '.*:CORRECT_BARCODES'
    {
        time = '15.h'
    }
}
```
