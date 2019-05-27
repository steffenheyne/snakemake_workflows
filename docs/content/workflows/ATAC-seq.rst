.. _ATAC-seq:

ATAC-seq
========

What it does
------------

The ATAC-seq pipeline takes one or more BAM files and attempts to find accessible regions. If multiple samples and a sample sheet are provided, then CSAW is additionally used to find differentially accessible regions. Prior to finding open/accessible regions, the BAM files are filtered to include only properly paired reads with appropriate fragment sizes (<150 bases by default). These filtered fragments are then used for the remainder of the pipeline.

.. image:: ../images/ATACseq_pipeline.png

.. note:: The **CSAW** step will be skipped if there is no ``sample_info`` tsv file (see :ref:`running_snakePipes`).

Input requirements
------------------

The DNA mapping pipeline generates output that is fully compatible with the ATAC-seq pipeline input requirements!
When running the ATAC-seq pipeline, please specify the output directory of DNA-mapping pipeline as the working directory (``-d``).

* **filtered_bam** directory contains the input BAM files (either filtered or unfiltered, however you prefer).

* **sampleSheet.tsv** (OPTIONAL) is only needed to test for differential binding.

Differential open chromatin analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Similar to differential binding analysis with the ChIP-Seq data. We can perform the differential open chromatin analysis, using the ``--sampleSheet`` option of the ATAC-seq workflow. This requires a sample sheet, which is identical to that required by the ChIP-seq and RNA-seq workflows (see :doc:`ChIP-seq` for details).

An example is below::

    name    condition
    sample1      eworo
    sample2      eworo
    SRR7013047      eworo
    SRR7013048      OreR
    SRR7013049      OreR
    SRR7013050      OreR

.. note:: This sample sheet has the same requirements as the sample sheet in the ChIP-seq workflow, and also uses the same tool (CSAW) with a narrow default window size.

Configuration file
~~~~~~~~~~~~~~~~~~

There is a configuration file in ``snakePipes/workflows/ATACseq/defaults.yaml``::

    ## General/Snakemake parameters, only used/set by wrapper or in Snakemake cmdl, but not in Snakefile
    pipeline: ATAC-seq
    configfile:
    cluster_configfile:
    local: false
    max_jobs: 5
    ## workingdir need to be required DNA-mapping output dir, 'outdir' is set to workingdir internally
    workingdir:
    ## preconfigured target genomes (mm9,mm10,dm3,...) , see /path/to/snakemake_workflows/shared/organisms/
    ## Value can be also path to your own genome config file!
    genome:
    ## The maximum fragment size to retain. This should typically be the size of a nucleosome
    fragmentSize_cutoff: 150
    verbose: false
    # sampleSheet_DB
    sampleSheet:
    # window_size
    window_size: 20
    fragmentCount_cutoff: 1
    #### Flag to control the pipeline entry point
    bam_ext: '.filtered.bam'
    fromBam: 
    ## Bin size of output files in bigWig format
    bw_binsize: 25
    paired: True
    plot_format: png
    ## Median/mean fragment length, only relevant for single-end data (default: 200)
    fragment_length: 
    trim:
    fastqc:
    qval: 0.001

Useful parameters are ``fragmentSize_cutoff`` and ``window_size``, also available from commandline.  

* **window_size**: is the size of windows to test differential binding using CSAW. The default small window size is sufficient for most analysis, since an ATAC-seq peak is sharp.

* **fragmentCount_cutoff**: refers to the minimum number of counts a chromosome must have to be included in the MACS2 analysis. It is introduced to avoid errors in the peak calling step and should only be changed if MACS2 fails.

* **Qval**: a value provided to MACS2 that affects the number and width of the resulting peaks.

Understanding the outputs
---------------------------

Assuming a sample sheet is used, the following will be **added** to the working directory::

    .
    ├── CSAW
    │   ├── CSAW.log
    │   ├── CSAW.session_info.txt
    │   ├── DiffBinding_allregions.bed
    │   ├── DiffBinding_analysis.Rdata
    │   ├── DiffBinding_modelfit.pdf
    │   ├── DiffBinding_scores.txt
    │   ├── DiffBinding_significant.bed
    │   ├── QCplots_first_sample.pdf
    │   ├── QCplots_last_sample.pdf
    │   └── TMM_normalizedCounts.pdf
    ├── deepTools_ATAC
    │   └── plotFingerprint
    │       ├── plotFingerprint.metrics.txt
    │       └── plotFingerprint.png
    ├── MACS2
    │   ├── sample1.filtered.BAM_control_lambda.bdg
    │   ├── sample1.filtered.BAM_peaks.narrowPeak
    │   ├── sample1.filtered.BAM_peaks.xls
    │   ├── sample1.filtered.BAM_summits.bed
    │   ├── sample1.filtered.BAM_treat_pileup.bdg
    │   ├── sample1.short.metrics
    │   ├── sample2.filtered.BAM_control_lambda.bdg
    │   ├── sample2.filtered.BAM_peaks.narrowPeak
    │   ├── sample2.filtered.BAM_peaks.xls
    │   ├── sample2.filtered.BAM_summits.bed
    │   ├── sample2.filtered.BAM_treat_pileup.bdg
    │   └── sample2.short.metrics
    └── MACS2_QC
        ├── sample1.filtered.BAM_peaks.qc.txt
        └── sample2.filtered.BAM_peaks.qc.txt

Currently the ATAC-seq workflow performs detection of open chromatin regions via `MACS2 <https://github.com/taoliu/MACS>`__, and if a sample sheet is provided, the detection of differential open chromatin sites via `CSAW <https://bioconductor.org/packages/release/bioc/html/csaw.html>`__. There are additionally log files in most of the directories. The various outputs are documented in the CSAW and MACS2 documentation.

* **MACS2_QC**: contains a number of QC metrics that we find useful, namely :
    * the number of peaks
    * fraction of reads in peaks (FRiP)
    * percentage of the genome covered by peaks.

* **deepTools_ATAC**: contains the output of `plotFingerPrint <https://deeptools.readthedocs.io/en/develop/content/tools/plotFingerprint.html>`__, which is a useful QC plot to assess signal enrichment between the ATAC-seq samples.


Where to find final bam files and biwgwigs
------------------------------------------

Bam files with the extention filtered.bam are only filtered for PCR duplicates. The final bam files filtered additionally for fragment size and used as direct input to MACS2 are found in the MACS2 folder with the exention ``.short.cleaned.bam``.
Bigwig files calculated from these bam files are found under deepTools_ATAC/bamCompare with the extention ``.filtered.bw``.


Command line options
--------------------

.. argparse::
    :func: parse_args
    :filename: ../snakePipes/workflows/ATAC-seq/ATAC-seq
    :prog: ATAC-seq
    :nodefault:
