/*
 * FastQ QC with NanoPlot, ToulligQC and fastqc
 * subworkflow from nf-core/nanoseq with minor modifications
 * (e.g: addition of ToulligQC)
 * author: @yuukiiwa
 */

params.nanoplot_fastq_options = [:]
params.toulligqc_fastqc_options = [:]
params.fastqc_options         = [:]

include { NANOPLOT     } from '../../modules/nf-core/nanoplot/main'  //addParams( options: params.nanoplot_fastq_options )
include { TOULLIGQC    } from '../../modules/nf-core/toulligqc/main' //addParams( options: params.toulligqc_fastqc_options )                                                                                                                    
include { FASTQC       } from '../../modules/nf-core/fastqc/main'    //addParams( options: params.fastqc_options )

workflow QCFASTQ_NANOPLOT_FASTQC {
    take:
    ch_fastq
    skip_nanoplot
    skip_toulligqc
    skip_fastqc

    main:
    ch_fastq
        .map { ch -> [ ch[0], ch[1] ] }
        .set { ch_fastq }

    /*
     * FastQ QC using NanoPlot
     */
    nanoplot_png     = channel.empty()
    nanoplot_html    = channel.empty()
    nanoplot_txt     = channel.empty()
    nanoplot_log     = channel.empty()
    if (!skip_nanoplot){
        NANOPLOT ( ch_fastq )
        nanoplot_png     = NANOPLOT.out.png
        nanoplot_html    = NANOPLOT.out.html
        nanoplot_txt     = NANOPLOT.out.txt
        nanoplot_log     = NANOPLOT.out.log
    }
    /*
     * FastQ QC using ToulligQC
     */
    toulligqc_report_data   = channel.empty()
    toulligqc_report_html   = channel.empty()
    toulligqc_plots_html    = channel.empty()
    toulligqc_plotly_js     = channel.empty()
    if (!skip_toulligqc){
        TOULLIGQC ( ch_fastq )
        toulligqc_report_data  = TOULLIGQC.out.report_data
        toulligqc_report_html  = TOULLIGQC.out.report_html
        toulligqc_plots_html   = TOULLIGQC.out.plots_html
        toulligqc_plotly_js    = TOULLIGQC.out.plotly_js
    }

    /*
     * FastQ QC using FASTQC
     */
    fastqc_zip     = channel.empty()
    fastqc_html    = channel.empty()
    fastqc_multiqc = channel.empty()
    if (!skip_fastqc){
        FASTQC ( ch_fastq )
        fastqc_zip     = FASTQC.out.zip
        fastqc_html    = FASTQC.out.html
        fastqc_zip
            .map { it -> [ it[1] ] }
            .set { fastqc_zip_only }
        fastqc_html
            .map { it -> [ it[1] ] }
            .set { fastqc_html_only }
        fastqc_multiqc = fastqc_multiqc.mix( fastqc_zip_only, fastqc_html_only )
    }

    emit:
    nanoplot_png
    nanoplot_html
    nanoplot_txt
    nanoplot_log

    toulligqc_report_data
    toulligqc_report_html
    toulligqc_plots_html
    toulligqc_plotly_js

    fastqc_zip
    fastqc_html
    fastqc_multiqc
}
