rule filterFragments:
    input:
        "filtered_bam/{sample}.filtered.bam"
    output:
        shortBAM = temp(os.path.join(outdir_MACS2, "{sample}.short.bam")),
        metrics = os.path.join(outdir_MACS2, "{sample}.short.metrics")
    params:
        cutoff = fragmentSize_cutoff
    threads: 6
    conda: CONDA_SHARED_ENV
    shell: """
        alignmentSieve --bam {input} \
        --outFile {output.shortBAM} -p {threads} \
        --filterMetrics {output.metrics} \
        --maxFragmentLength {params.cutoff}
        """

# MACS2 BAMPE fails if there is just one fragment mapped
rule filterCoveragePerScaffolds:
    input:
        bam = os.path.join(outdir_MACS2, "{sample}.short.bam")
    output:
        whitelist = os.path.join(outdir_MACS2, "{sample}.chrom.whitelist"),
        shortbai = temp(os.path.join(outdir_MACS2, "{sample}.short.bam.bai")),
        bam = temp(os.path.join(outdir_MACS2, "{sample}.short.cleaned.bam")),
        bai = temp(os.path.join(outdir_MACS2, "{sample}.short.cleaned.bam.bai"))
    params:
        count_cutoff = int(fragmentCount_cutoff) * 2 # must contain more than 2 reads, i.e. 1 fragment
    threads: 6
    conda: CONDA_SHARED_ENV
    shell: """
        samtools index -@ {threads} {input.bam}
        samtools idxstats {input.bam} | awk -v cutoff={params.count_cutoff} \'$3 > cutoff\' | cut -f 1 > {output.whitelist}
        samtools view -@ {threads} -bo {output.bam} {input.bam} $(cat {output.whitelist} | paste -sd\' \')
        samtools index -@ {threads} {output.bam}
        """

# MACS2 BAMPE filter: samtools view -b -f 2 -F 4 -F 8 -F 256 -F 512 -F 2048
rule callOpenChromatin:
    input:
        os.path.join(outdir_MACS2, "{sample}.short.cleaned.bam")
    output:
        peaks = os.path.join(outdir_MACS2, '{sample}.filtered.BAM_peaks.narrowPeak'),
        xls = os.path.join(outdir_MACS2, '{sample}.filtered.BAM_peaks.xls')
    params:
        directory = outdir_MACS2,
        genome_size = int(genome_size),
        name='{sample}',
        qval_cutoff='--qvalue 0.001',
        nomodel='--nomodel',
        write_bdg='--bdg',
        fileformat='--format BAMPE'
    threads: 6
    log:
        out = os.path.join(outdir_MACS2, "logs", "callOpenChromatin", "{sample}_macs2.out"),
        err = os.path.join(outdir_MACS2, "logs", "callOpenChromatin", "{sample}_macs2.err")
    conda: CONDA_ATAC_ENV
    shell: """
        macs2 callpeak --treatment {input} \
            -g {params.genome_size} \
            --name {params.name}.filtered.BAM \
            --outdir {params.directory} \
            {params.fileformat} \
            {params.qval_cutoff} \
            {params.nomodel} \
            {params.write_bdg} > {log.out} 2> {log.err}
        """
