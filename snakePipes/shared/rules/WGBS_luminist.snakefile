import os
import re
from operator import is_not
import tempfile
import pandas


###symlink bams if this is the starting point
if fromBam:
    rule link_bam:
        input:
            indir+"/{sample}"+bam_ext
        output:
            "bams/{sample}.sorted.bam"
        shell:
            "( [ -f {output} ] || ln -s -r {input} {output} ) " #&& touch -h {output}"


if trimReads=='user':
    rule trimReads:
        input:
            R1 = "FASTQ/{sample}"+reads[0]+".fastq.gz",
            R2 = "FASTQ/{sample}"+reads[1]+".fastq.gz"
        output:
            R1cut="FASTQ_Cutadapt/{sample}"+reads[0]+".fastq.gz",
            R2cut="FASTQ_Cutadapt/{sample}"+reads[1]+".fastq.gz"
        log:
            err="FASTQ_Cutadapt/logs/{sample}.trimReads.err",
            out="FASTQ_Cutadapt/logs/{sample}.trimReads.out"
        params:
            adapterSeq=adapterSeq,
            trimThreshold=trimThreshold,
            trimOtherArgs=lambda wildcards: '' if trimOtherArgs is None else str(trimOtherArgs)
        threads: nthreads
        conda: CONDA_SHARED_ENV
        shell: "cutadapt -a {params.adapterSeq} -A {params.adapterSeq} -q {params.trimThreshold} -m 30 -j {threads} {params.trimOtherArgs} -o {output.R1cut} -p {output.R2cut} {input.R1} {input.R2} 1>{log.out} 2>{log.err}"


rule preTrimFastQC:
    input:
        R1="FASTQ/{sample}"+reads[0]+".fastq.gz",
        R2="FASTQ/{sample}"+reads[1]+".fastq.gz"
    output:
        R1fqc="FastQC/{sample}"+reads[0]+"_fastqc.html",
        R2fqc="FastQC/{sample}"+reads[1]+"_fastqc.html"
    log:
        err="FastQC/logs/{sample}.preTrimFastQC.err",
        out="FastQC/logs/{sample}.preTrimFastQC.out"
    params:
        fqcout=os.path.join(outdir,'FastQC')
    threads: nthreads
    conda: CONDA_SHARED_ENV
    shell: "fastqc --outdir {params.fqcout} -t  {threads} {input.R1} {input.R2} 1>{log.out} 2>{log.err}"


rule postTrimFastQC:
    input:
        R1cut="FASTQ_Cutadapt/{sample}"+reads[0]+".fastq.gz",
        R2cut="FASTQ_Cutadapt/{sample}"+reads[1]+".fastq.gz"
    output:
        R1fqc="FastQC_Cutadapt/{sample}"+reads[0]+"_fastqc.html",
        R2fqc="FastQC_Cutadapt/{sample}"+reads[1]+"_fastqc.html"
    log:
        err="FastQC_Cutadapt/logs/{sample}.postTrimFastQC.err",
        out="FastQC_Cutadapt/logs/{sample}.postTrimFastQC.out"
    params:
        fqcout=os.path.join(outdir,'FastQC_Cutadapt')
    threads: nthreads
    conda: CONDA_SHARED_ENV
    shell: "fastqc --outdir {params.fqcout} -t  {threads} {input.R1cut} {input.R2cut} 1>{log.out} 2>{log.err}"


if not fromBam:
    rule map_reads:
        input:
            R1=fastq_dir+"/{sample}"+reads[0]+".fastq.gz",
            R2=fastq_dir+"/{sample}"+reads[1]+".fastq.gz",
            crefG=bwameth_index
        output:
            sbam="bams/{sample}.sorted.bam"
        log:
            err="bams/logs/{sample}.map_reads.err",
            out="bams/logs/{sample}.map_reads.out"
        params:
            #tempdir=tempfile.mkdtemp(suffix='',prefix="{sample}",dir=tempdir),
            sortThreads=min(nthreads,4),
            RG=lambda wildcards: RG_dict[wildcards.sample]
        threads: nthreads
        conda: CONDA_WGBS_ENV
        shell: "tmp_map=$(mktemp -d -p $TMPDIR -t XXXXX.{wildcards.sample});echo $tmp_map; bwameth.py --threads  {threads} --read-group {params.RG} --reference {input.crefG} {input.R1} {input.R2} | samtools sort -T $tmp_map -m 3G -@ {params.sortThreads} -o {output.sbam} 1>{log.out} 2>{log.err}"

rule index_bam:
    input:
        sbam="bams/{sample}.sorted.bam"
    output:
        sbami="bams/{sample}.sorted.bam.bai"
    log:
        err="bams/logs/{sample}.index_bam.err",
        out="bams/logs/{sample}.index_bam.out"
    conda: CONDA_SHARED_ENV
    shell: "samtools index {input.sbam} 1>{log.out} 2>{log.err}"

rule rm_dupes:
    input:
        sbami="bams/{sample}.sorted.bam.bai",
        sbam="bams/{sample}.sorted.bam"
    output:
        rmDupbam="bams/{sample}.PCRrm.bam"
    log:
        err="bams/logs/{sample}.rm_dupes.err",
        out="bams/logs/{sample}.rm_dupes.out"
    params:
        #tempdir=tempfile.mkdtemp(suffix='',prefix='',dir=tempdir)
    threads: nthreads
    conda: CONDA_SAMBAMBA_ENV
    shell: "tmp_dupes=$(mktemp -d -p $TMPDIR -t XXXXX.{wildcards.sample}); echo $tmp_dupes; sambamba markdup --hash-table-size=4194304 --remove-duplicates --tmpdir $tmp_dupes -t {threads} {input.sbam} {output.rmDupbam} 1>{log.out} 2>{log.err}"


rule index_PCRrm_bam:
    input:
        sbam="bams/{sample}.PCRrm.bam"
    output:
        sbami="bams/{sample}.PCRrm.bam.bai"
    params:
    log:
        err="bams/logs/{sample}.index_PCRrm_bam.err",
        out="bams/logs/{sample}.index_PCRrm_bam.out"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell: "samtools index {input.sbam} 1>{log.out} 2>{log.err}"


rule get_ran_CG:
     input:
         refG=genome_fasta
     output:
         pozF="aux_files/genome.poz.gz"
     params:
         methylCtools=os.path.join(workflow_tools,'methylCtools')
     log:
         err="aux_files/logs/get_ran_CG.err"
     threads: nthreads
     conda: CONDA_SHARED_ENV
     shell: """
         set +o pipefail; {params.methylCtools} fapos {input.refG} - | pigz -c -p {threads} > {output.pozF};
         """             


rule calc_Mbias:
    input:
        refG=genome_fasta,
        rmDupBam="bams/{sample}.PCRrm.bam",
        sbami="bams/{sample}.PCRrm.bam.bai"
    output:
        mbiasTXT="QC_metrics/{sample}.Mbias.txt"
    log:
        out="QC_metrics/logs/{sample}.calc_Mbias.out"
    threads: nthreads
    conda: CONDA_WGBS_ENV
    shell: "MethylDackel mbias {input.refG} {input.rmDupBam} {output.mbiasTXT} -@ {threads} 1>{log.out} 2>{output.mbiasTXT}"


if not fromBam:
    rule downsample_reads:
        input:
            R1=fastq_dir+"/{sample}"+reads[0]+".fastq.gz",
            R2=fastq_dir+"/{sample}"+reads[1]+".fastq.gz"
        output:
            R1downsampled="FASTQ_downsampled/{sample}"+reads[0]+".fastq.gz",
            R2downsampled="FASTQ_downsampled/{sample}"+reads[1]+".fastq.gz"
        log:
            err="FASTQ_downsampled/logs/{sample}.downsample_reads.err",
            out="FASTQ_downsampled/logs/{sample}.downsample_reads.out"
        threads: nthreads
        conda: CONDA_SHARED_ENV
        shell: """
                seqtk sample -s 100 {input.R1} 5000000 | pigz -p {threads} -9 > {output.R1downsampled}
                seqtk sample -s 100 {input.R2} 5000000 | pigz -p {threads} -9 > {output.R2downsampled}
                1>{log.out} 2>{log.err}
               """

    rule conv_rate:
        input:
            R1downsampled="FASTQ_downsampled/{sample}"+reads[0]+".fastq.gz",
            R2downsampled="FASTQ_downsampled/{sample}"+reads[1]+".fastq.gz"
        output:
            R12cr="QC_metrics/{sample}.conv.rate.txt"
        params:
            read_root=lambda wildcards: "FASTQ_downsampled/"+wildcards.sample
        log:
            err="QC_metrics/logs/{sample}.conv_rate.err",
            out="QC_metrics/logs/{sample}.conv_rate.out"
        threads: 1
        shell: os.path.join(workflow_tools,'conversionRate_KS.sh ')+ "{params.read_root} {output.R12cr} 1>{log.out} 2>{log.err}"


rule get_flagstat:
    input:
        bam="bams/{sample}.PCRrm.bam"
    output:
        fstat="QC_metrics/{sample}.flagstat"
    log:
        err="QC_metrics/logs/{sample}.get_flagstat.err"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell: "samtools flagstat {input.bam} > {output.fstat} 2>{log.err}"

rule get_flagstat_sorted:
    input:
        bam="bams/{sample}.sorted.bam"
    output:
        fstat="QC_metrics/{sample}.sorted.flagstat"
    log:
        err="QC_metrics/logs/{sample}.sorted.get_flagstat.err"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell: "samtools flagstat {input.bam} > {output.fstat} 2>{log.err}"

rule produce_report:
    input:
        expand("QC_metrics/{sample}.conv.rate.txt",sample=samples) if not fromBam else [],
        mbiasTXT=expand("QC_metrics/{sample}.Mbias.txt",sample=samples),
        fstat=expand("QC_metrics/{sample}.flagstat",sample=samples)
    output:
        QCrep='QC_metrics/QC_report.html'
    params:
        auxdir=os.path.join(outdir,"aux_files")
    log:
        err="QC_metrics/logs/produce_report.err",
        out="QC_metrics/logs/produce_report.out"
    conda: CONDA_RMD_ENV
    threads: 1
    shell: "cp -v " + os.path.join(workflow_rscripts,"WGBS_QC_report_template.Rmd")+ " " + os.path.join("aux_files", "WGBS_QC_report_template.Rmd") + ';Rscript -e "rmarkdown::render(\''+os.path.join(outdir,"aux_files", "WGBS_QC_report_template.Rmd")+'\', params=list(QCdir=\'"' + os.path.join(outdir,"QC_metrics") +'"\' ), output_file =\'"'+ os.path.join(outdir,"QC_metrics",'QC_report.html"\'')+')"' + " 1>{log.out} 2>{log.err}"


if mbias_ignore=="auto":
    rule methyl_extract:
        input:
            rmDupbam="bams/{sample}.PCRrm.bam",
            sbami="bams/{sample}.PCRrm.bam.bai",
            refG=genome_fasta,
            mbiasTXT="QC_metrics/{sample}.Mbias.txt"
        output:
            methTab="methXT/{sample}_CpG.bedGraph"
        params:
            OUTpfx=lambda wildcards,output: os.path.join(outdir,re.sub('_CpG.bedGraph','',output.methTab))
        log:
            err="methXT/logs/{sample}.methyl_extract.err",
            out="methXT/logs/{sample}.methyl_extract.out"
        threads: nthreads
        conda: CONDA_WGBS_ENV
        shell: "mi=$(cat {input.mbiasTXT} | sed 's/Suggested inclusion options: //' );MethylDackel extract  -o {params.OUTpfx} -q 10 -p 20 $mi --minDepth 10 --mergeContext --maxVariantFrac 0.25 --minOppositeDepth 5 -@ {threads} {input.refG} " + os.path.join(outdir,"{input.rmDupbam}") + " 1>{log.out} 2>{log.err}"

else:
    rule methyl_extract:
        input:
            rmDupbam="bams/{sample}.PCRrm.bam",
            sbami="bams/{sample}.PCRrm.bam.bai",
            refG=genome_fasta
        output:
            methTab="methXT/{sample}_CpG.bedGraph"
        params:
            mbias_ignore=mbias_ignore,
            OUTpfx=lambda wildcards,output: os.path.join(outdir,re.sub('_CpG.bedGraph','',output.methTab))
        log:
            err="methXT/logs/{sample}.methyl_extract.err",
            out="methXT/logs/{sample}.methyl_extract.out"
        threads: nthreads
        conda: CONDA_WGBS_ENV
        shell: "MethylDackel extract -o {params.OUTpfx} -q 10 -p 20 {params.mbias_ignore} --minDepth 10 --mergeContext --maxVariantFrac 0.25 --minOppositeDepth 5 -@ {threads} {input.refG} " + os.path.join(outdir,"{input.rmDupbam}") + " 1>{log.out} 2>{log.err}"

if blackList is None:
    rule CpG_filt:
        input:
            methTab="methXT/{sample}_CpG.bedGraph"
        output:
            tabFilt="methXT/{sample}.CpG.filt2.bed"
        params:
            OUTtemp=lambda wildcards,input: os.path.join(outdir,re.sub('_CpG.bedGraph','.CpG.filt.bed',input.methTab))
        log:
            err="methXT/logs/{sample}.CpG_filt.err",
            out="methXT/logs/{sample}.CpG_filt.out"
        threads: 1
        conda: CONDA_WGBS_ENV
        shell: '''awk \'(NR>1)\' {input.methTab} | awk \'{{ print $0, $5+$6, $1\"_\"$2}}\' | tr " " "\t" | sed \'1i chr\tstart\tend\tBeta\tM\tU\tCov\tms\' > {params.OUTtemp};mv -v {params.OUTtemp} {output.tabFilt} 1>{log.out} 2>{log.err}'''

else:
    rule CpG_filt:
        input:
            methTab="methXT/{sample}_CpG.bedGraph",
            blackListF=blackList
        output:
            tabFilt="methXT/{sample}.CpG.filt2.bed"
        params:
            OUTtemp=lambda wildcards,input: os.path.join(outdir,re.sub('_CpG.bedGraph','.CpG.filt.bed',input.methTab))
        log:
            err="methXT/logs/{sample}.CpG_filt.err",
            out="methXT/logs/{sample}.CpG_filt.out"
        threads: 1
        conda: CONDA_WGBS_ENV
        shell: '''awk \'(NR>1)\' {input.methTab} | awk \'{{ print $0, $5+$6, $1\"_\"$2}}\' | tr " " "\t" | sed \'1i chr\tstart\tend\tBeta\tM\tU\tCov\tms\' > {params.OUTtemp};bedtools intersect -v -a {params.OUTtemp} -b {input.blackListF} > {output.tabFilt} 1>{log.out} 2>{log.err}'''


rule make_CG_bed:
    input:
        pozF="aux_files/genome.poz.gz"
    output:
        imdF="aux_files/genome.CpG.bed"
    log:
        err="aux_files/logs/make_CG_bed.err"
    threads: 1
    conda: CONDA_WGBS_ENV
    shell: """
        zcat {input.pozF} | grep "+" | awk '{{OFS="\t"; print $1,$3,$3+1,$4,$6}}' | sort -k 1,1 -k2,2n > {output.imdF}
        """

rule get_CG_per_int:
    input:
        intList=intList,
        refG=genome_fasta,
        imdF="aux_files/genome.CpG.bed"
    output:
        outList=run_int_aggStats(intList)
    log:
        err="aux_files/logs/get_CG_per_int.err"
    params:
        auxshell=lambda wildcards,input,output: ';'.join(["bedtools intersect -wa -a "+ input.imdF + " -b " + bli + ' > ' + oli  for bli,oli in zip(input.intList,output.outList) ])
    threads: 1
    conda: CONDA_WGBS_ENV
    shell: "{params.auxshell} 2>{log.err}"

        
rule on_target_rate:
    input:
        bams=expand("bams/{sample}.PCRrm.bam",sample=samples)
    output:
        tab="custom_stats/on_target_stats.all_reads.txt",
        plot="custom_stats/on_target_stats.all_reads.pdf"
    params:
        targets=intList,
        labels = " ".join(samples)
    log:
        err="custom_stats/logs/on_target_stats.all_reads.err",
        out="custom_stats/logs/on_target_stats.all_reads.out"
    threads: nthreads
    shell:"""
        plotEnrichment -p {threads} \
               -b {input.bams} \
               --plotFile {output.plot}\
               --BED {params.targets} \
               --labels {params.labels} \
               --plotTitle 'Fraction of reads in target regions' \
               --outRawCounts {output.tab} \
               --variableScales > {log.out} 2> {log.err}
        """

rule on_target_rate_mapq:
        input:
            bams=expand("bams/{sample}.PCRrm.bam",sample=samples)
        output:
            tab="custom_stats/on_target_stats.mapq20.txt",
            plot="custom_stats/on_target_stats.mapq20.pdf"
        params:
            targets=intList,
            labels = " ".join(samples)
        log:
            err="custom_stats/logs/on_target_stats.mapq20.err",
            out="custom_stats/logs/on_target_stats.mapq20.out"
        threads: nthreads
        conda: CONDA_SHARED_ENV
        shell:"""
            plotEnrichment -p {threads} \
                   -b {input.bams} \
                   --plotFile {output.plot}\
                   --BED {params.targets} \
                   --labels {params.labels} \
                   --plotTitle 'Fraction of reads in target regions' \
                   --outRawCounts {output.tab} \
                   --variableScales \
                   --minMappingQuality 20 1> {log.out} 2> {log.err}
            """

rule on_target_reads_region:
        input:
            bams=expand("bams/{sample}.PCRrm.bam",sample=samples)
        output:
            tabq20="custom_stats/on_target_stats.per_region.mapq20.tsv",
            tab="custom_stats/on_target_stats.per_region.tsv"
        params:
            targets=intList,
            labels = " ".join(samples)
        log:
            err="custom_stats/logs/on_target_stats.per_region.err",
            out="custom_stats/logs/on_target_stats.per_region.out"
        threads: nthreads
        conda: CONDA_SHARED_ENV
        shell:"""
            multiBamSummary BED-file \
                -b {input.bams} \
                --BED {params.targets} \
                --outRawCounts tmp.q20.tsv \
                --minMappingQuality 20 \
                --labels {params.labels} \
                -p {threads} 1> {log.out} 2> {log.err};
            sort -k1,1 -k2,2n tmp.q20.tsv > {output.tabq20}; rm tmp.q20.tsv; 
            multiBamSummary BED-file \
                -b {input.bams} \
                --BED {params.targets} \
                --outRawCounts tmp.tsv \
                --labels {params.labels} \
                -p {threads} 1>> {log.out} 2>> {log.err};
            sort -k1,1 -k2,2n tmp.tsv > {output.tab}; rm tmp.tsv;
            """

rule methyl_extract_custom:
        input:
            rmDupbam="bams/{sample}.PCRrm.bam",
            sbami="bams/{sample}.PCRrm.bam.bai",
            refG=genome_fasta
        output:
            methTab="custom_stats/{sample}_CpG.bedGraph",
            meanTab="custom_stats/{sample}.mean_methyl_per_region.tsv"
        params:
            mbias_ignore=mbias_ignore,
            targets=intList,
            OUTpfx=lambda wildcards,output: os.path.join(outdir,re.sub('_CpG.bedGraph','',output.methTab))
        log:
            err="custom_stats/logs/{sample}.methyl_extract.err",
            out="custom_stats/logs/{sample}.methyl_extract.out"
        threads: nthreads
        conda: CONDA_WGBS_ENV
        shell: """
            MethylDackel extract -o {params.OUTpfx} -l {params.targets} \
                -q 20 -p 20 --minDepth 1 --mergeContext -@ {threads} \
                {input.refG} {input.rmDupbam} 1>{log.out} 2>{log.err};
            bedtools map -a {params.targets} -b {output.methTab} \
                -c 4 -o mean -prec 4 > {output.meanTab} 2>>{log.err}
            """

rule per_base_cov_custom:
    input:
        bams=expand("bams/{sample}.PCRrm.bam",sample=samples)
    output:
        "custom_stats/coverage_per_base.targets.bed"
    params:
        targets=intList
    conda: CONDA_SHARED_ENV
    log:
        err="custom_stats/logs/coverage_per_base.targets.err",
        out="custom_stats/logs/coverage_per_base.targets.out"
    shell:"""
        cat <(echo -e '#chr\tpos\t'$(echo '{input.bams}' | tr ' ' '\n' | sed 's/.*\///' | sed 's/.PCRrm.bam//g' | tr '\n' '\t')) <(samtools depth -a -q 20 -Q 20 {input.bams} -b <(cat {params.targets} | awk '{{OFS="\t";print $1,$2-1,$3-1}}') ) > {output} 2>{log.err}
        """

rule target_cpgs:
    input:
        "aux_files/genome.CpG.bed"
    params:
        targets=intList
    output:
        "custom_stats/targets.CpG.bed"
    log:
        err="custom_stats/logs/targets.CpG.err"
    conda: CONDA_WGBS_ENV
    shell:"""
        bedtools intersect -a <(cat {params.targets} | awk '{{OFS="\t";$2=$2-1;$3=$3+2; print $0}}') -b {input} -wo | awk '{{OFS="\t"; print $4,$5,$6,$1"_"$2,0,$7,$8}}' > {output} 2>{log.err}
    """

rule target_cpg_coverage:
    input:
        bams=expand("bams/{sample}.PCRrm.bam",sample=samples),
        cpg="custom_stats/targets.CpG.bed"
    output:
        "custom_stats/targets.CpG.coverage.txt"
    log:
        err="custom_stats/logs/targets.CpG.coverage.err"
    conda: CONDA_SHARED_ENV
    shell: """
        cat <(echo -e '#chr\tpos\t'$(echo '{input.bams}' | tr ' ' '\n' | sed 's/.*\///' | sed 's/.PCRrm.bam//g' | tr '\n' '\t')) <(samtools depth -a -q 20 -Q 20 -b {input.cpg} {input.bams}) > {output} 2>{log.err}
    """
    

rule mean_target_coverage:
    input:
        "custom_stats/coverage_per_base.targets.bed"
    output:
        "custom_stats/mean_coverage_per_base.targets.bed"
    params:
        targets=intList
    log:
        err="custom_stats/logs/mean_coverage_per_base.targets.err"
    conda: CONDA_WGBS_ENV
    shell:"""
        cat <(cat {input} | head -n1 | awk '{{OFS="\t";$2="start\tend";print $0}}') \
            <(bedtools map -a {params.targets} -b <(cat {input} | \
              awk '{{OFS="\t";$2=$2-1"\t"$2; print $0}}') \
            -c $(cat {input} | awk '{{}}END{{for (i=4;i<=NF;i++){{printf i","}}; print NF+1}}') \
            -o mean -prec 5) > {output} 2>{log.err}
        """

rule mean_methyl_per_region:
    input:
        tsv=expand("custom_stats/{sample}.mean_methyl_per_region.tsv",sample=samples),
        tab="custom_stats/on_target_stats.per_region.mapq20.tsv",
        meth=expand("custom_stats/{sample}_CpG.bedGraph",sample=samples),
        cpgs="custom_stats/targets.CpG.bed"
        
    output:
        "custom_stats/mean_methyl_per_region.tsv"
    params:
        indir="custom_stats/",
        script = os.path.join(workflow_rscripts,"merge_methyl_data.R")
    log:
        err="custom_stats/logs/mean_methyl_per_region.err",
        out="custom_stats/logs/mean_methyl_per_region.out"
    conda: CONDA_WGBS_ENV
    shell:"""
            Rscript {params.script} {params.indir} {output} 2> {log.err}
        """

#sambamba depth base ../bams/*.PCRrm.bam -L /data_BSpipe/datasets/target_regions/targets.bed -z -m --min-base-quality 20 -F='mapping_quality>20' -t 20
