# Dependencies
import os
import shutil
import pandas as pd

# **********************************************
# Configure Snakemake. Do NOT CHANGE!
localrules: all
configfile: 'Config/config_snake.json'

# Sample sheet & sample list
sampleDF = pd.read_csv('sampleSheet.csv', comment = '#')
sampleList = list(set(list(sampleDF.sampleName)))

# Functions  
def getFastqFiles(wildcards):
    r1 = sampleDF.r1[(sampleDF.sampleName == wildcards.sample)].to_list()
    r2 = sampleDF.r2[(sampleDF.sampleName == wildcards.sample)].to_list()
    return({'r1' : r1[0], 'r2' : r2[0]})
    #return( r1[0], r2[0])

# **********************************************
# Species to be used. EDIT THIS AS NEEDED
species = 'mm10'

# Start of the workflow
rule all:
    input:
        expand('Bam/{sample}_bwa.bam', sample = sampleList),
        #expand('Bam/{sample}_bwa_mapq_{read}.srt.bam.bai', sample = sampleList, read = ['R1', 'R2']),
        expand('Bam/{sample}{step}.flagstat.txt', sample = sampleList, step = ['_bwa', '_bwa_mapq', '_bwa_mapq_R1.srt', '_bwa_mapq_R2.srt', '_bwa_mapq_R1.multi_filt', '_bwa_mapq_R2.multi_filt', '_bwa_mapq_R1.multi_filt.dedup', '_bwa_mapq_R2.multi_filt.dedup']),
        #expand('Matrix/{sample}{filter}.mcool', sample = sampleList, filter = ['.srt', '.multi_filt.dedup']),
        expand('FastQC/{sample}_bwa_{read}_mapStats.tsv', sample = sampleList, read = ['R1', 'R2']),
        expand('FastQC/{sample}_bwa_flagstat_summary.tsv', sample = sampleList),
        expand('Matrix/{sample}_bwa_mapq{filter}.mcool', sample = sampleList, filter = ['.srt', '.multi_filt.dedup']),
        expand('BigWig/{sample}_bwa_mapq_R1{filter}rpgcNorm.bw', sample = sampleList, filter = ['.srt', '.multi_filt.dedup']),
        expand('FastQC/{sample}_{read}_trim_screen.html', sample = sampleList, read = ['R1', 'R2']),
        'FastQC/multiqc_report.html'
        #expand('BigWig/{sample}_bwa_{read}_rpgcNorm.bw', sample = sampleList, read = ['R1', 'R2']),
        #expand('FastQC/{sample}_bwa_{read}_mapStats.tsv', sample = sampleList, read = ['R1', 'R2']),

rule adapter_trim_reads:
    conda: 'trim_galore'
    input:
        unpack(getFastqFiles)
        #r1 = 'Fastq/{sample}_R1.fastq.gz',
        #r2 = 'Fastq/{sample}_R2.fastq.gz'
    output:
        r1 = 'Fastq/{sample}_R1_trim.fastq.gz',
        r2 = 'Fastq/{sample}_R2_trim.fastq.gz',
    threads: 8
    shell:
        """
        trim_galore --basename {wildcards.sample} --cores 4 --phred33 --quality 0 --stringency 10 --length 20 --fastqc -o Fastq/ --paired {input.r1} {input.r2}
        mv Fastq/{wildcards.sample}_val_1.fq.gz {output.r1}
        mv Fastq/{wildcards.sample}_val_2.fq.gz {output.r2}
        mv Fastq/{input.r1}_trimming_report.txt FastQC/
        mv Fastq/{input.r2}_trimming_report.txt FastQC/
        """

        #trim_galore --basename {wildcards.sample} --cores 4 --phred33 --quality 0 --stringency 10 --length 20 --fastqc -o Fastq/ --paired {input}
#        mv Fastq/{wildcards.sample}_R1.fastq.gz_trimming_report.txt FastQC/
#        mv Fastq/{wildcards.sample}_R2.fastq.gz_trimming_report.txt FastQC/

rule align_reads_bwa:
    conda: 'bwa_0_7_18'
    input:
        #getFastqFiles
        r1 = 'Fastq/{sample}_R1_trim.fastq.gz',
        r2 = 'Fastq/{sample}_R2_trim.fastq.gz',
    output:
        temp('Bam/{sample}_bwa.sam')
    params:
        index = config['index'][species]['bwa']
    threads: 12
    shell:
        """
        bwa mem -SP -t 12 {params.index} {input} 2>>FastQC/{wildcards.sample}_bwa_log.txt > {output}
        """

#bwa mem -SP -t 12 ./version0.6.0/genome.fa mESC_8_tn5_1hr_HICAR_sample_REP1_T1_1.trim.fastq.gz mESC_8_tn5_1hr_HICAR_sample_REP1_T1_2.trim.fastq.gz
#bwa mem -A1 -B4 -E50 -L0 {params.index} {input} 2>>mate_{wildcards.read}.log > {output}

rule convert_sam_filter_and_sort_bam:
    conda: 'samtools'
    input:
        'Bam/{sample}_bwa.sam'
    output:
        bam = 'Bam/{sample}_bwa.bam',
        mapq =     temp('Bam/{sample}_bwa_mapq.bam'),
        r1_coord = temp('Bam/{sample}_bwa_mapq_R1.srt.bam'),
        r2_coord = temp('Bam/{sample}_bwa_mapq_R2.srt.bam'),
        r1_name =  temp('Bam/{sample}_bwa_mapq_R1.nSrt.bam'),
        r2_name =  temp('Bam/{sample}_bwa_mapq_R2.nSrt.bam'),
        bamIndex =      'Bam/{sample}_bwa.bam.bai',
        mapqIndex =     'Bam/{sample}_bwa_mapq.bam.bai',
        r1_coordIndex = 'Bam/{sample}_bwa_mapq_R1.srt.bam.bai',
        r2_coordIndex = 'Bam/{sample}_bwa_mapq_R2.srt.bam.bai',
    threads: 8
    params: 
        mapq = 5
    shell:
        """
        samtools view -@ {threads} -hb {input} | samtools sort -@ {threads} -o {output.bam}
        samtools view -@ {threads} -hb -q {params.mapq} {input} | samtools sort - > {output.mapq}
        samtools view -@ {threads} -hb -f 64 {output.mapq} | samtools sort -@ {threads} - > {output.r1_coord}
        samtools view -@ {threads} -hb -f 128 {output.mapq} | samtools sort -@ {threads} - > {output.r2_coord}
        samtools sort -@ {threads} -n {output.r1_coord} -o {output.r1_name}
        samtools sort -@ {threads} -n {output.r2_coord} -o {output.r2_name}

        samtools index {output.bam}
        samtools index {output.mapq}
        samtools index {output.r1_coord}
        samtools index {output.r2_coord}
        """

rule collect_multi_map_stats:
    conda: 'base'
    input:
        'Bam/{sample}_bwa_mapq_{read}.nSrt.bam'
    output:
        'FastQC/{sample}_bwa_{read}_mapStats.tsv'
    params:
        source_dir = 'Src/'
    shell:
        """
        python3.12 {params.source_dir}/filter_and_check_multimappers_v3.py -i {input} -o {output}
        """

rule filter_multimappers:
    conda: 'base'
    input:
        r1 = 'Bam/{sample}_bwa_mapq_R1.nSrt.bam',
        r2 = 'Bam/{sample}_bwa_mapq_R2.nSrt.bam'
    output:
        r1_filter = temp('Bam/{sample}_bwa_mapq_R1.nSrt.filtered.bam'),
        r2_filter = temp('Bam/{sample}_bwa_mapq_R2.nSrt.filtered.bam'),
        r1_sync = temp('Bam/{sample}_bwa_mapq_R1.nSrt.sync.bam'),
        r2_sync = temp('Bam/{sample}_bwa_mapq_R2.nSrt.sync.bam'),
    params:
        source_dir = 'Src/'
    shell:
        """
        python3.12 {params.source_dir}/filter_and_output_multi_mappers.py -i {input.r1} -o {output.r1_filter}
        python3.12 {params.source_dir}/filter_and_output_multi_mappers.py -i {input.r2} -o {output.r2_filter}
        python3.12 {params.source_dir}/sync_bams.py -r1 {output.r1_filter} -r2 {output.r2_filter} -o1 {output.r1_sync} -o2 {output.r2_sync}
        """

rule sort_filtered_reads:
    conda: 'samtools'
    input:
        'Bam/{sample}_bwa_mapq_{read}.nSrt.filtered.bam'
    output:
        bam = temp('Bam/{sample}_bwa_mapq_{read}.multi_filt.bam'),
        bam_ind = 'Bam/{sample}_bwa_mapq_{read}.multi_filt.bam.bai'
    threads: 8
    shell:
        """
        samtools sort -@ {threads} {input} -o {output.bam}
        samtools index {output.bam}
        """

rule remove_duplicates:
    conda: 'picard'
    input:
        bam = 'Bam/{sample}_bwa_mapq_{read}.multi_filt.bam',
        #r2_bam = 'Bam/{sample}_bwa_mapq_R2.multi_filt.bam',
        bam_ind = 'Bam/{sample}_bwa_mapq_{read}.multi_filt.bam.bai'
    output:
        bam = 'Bam/{sample}_bwa_mapq_{read}.multi_filt.dedup.bam',
        #bam_ind = 'Bam/{sample}_bwa_mapq_{read}.multi_filt.dedup.bam.bai',
        metrics = 'Bam/{sample}_bwa_mapq_{read}.multi_filt.dedup.metrics'
    params:
        #picard = 'picard',
        tmp_dir = 'Tmp/'
    shell:
        """
        picard MarkDuplicates I={input.bam} O={output.bam} M={output.metrics} REMOVE_DUPLICATES=true TMP_DIR={params.tmp_dir} ASSUME_SORTED=true
        """
        #samtools index {output.bam}

rule index_bam:
    conda: 'samtools'
    input:
        'Bam/{sample}_bwa_mapq_{read}.multi_filt.dedup.bam'
    output:
        'Bam/{sample}_bwa_mapq_{read}.multi_filt.dedup.bam.bai'
    shell:
        """
        samtools index {input} 
        """

rule run_flagstat:
    conda: 'samtools'
    input:
        'Bam/{sample}{step}.bam'
        #expand('Bam/{{sample}}{step}.bam', step = [
        #    '', '_mapq', '_mapq_R1.srt', '_mapq_R2.srt',
        #     '_mapq_R1.multi_filt', '_mapq_R2.multi_filt',
        #      '_mapq_R1.multi_filt.dedup', '_mapq_R2.multi_filt.dedup']),
    output:
        'Bam/{sample}{step}.flagstat.txt'
    shell:
        """
        samtools flagstat {input} > {output}
        """

#        for input_bam in {input}; do
#            step=$(basename ${{input_bam}} .bam | sed 's/{wildcards.sample}_bwa//')
#            samtools flagstat ${{input_bam}} > Bam/{wildcards.sample}_bwa${{step}}.flagstat.txt
#        done

#rule run_flagstat:
#    conda: 'samtools'
#    input:
#        'Bam/{sample}_bwa{step}.bam'
#    output:
#        'Bam/{sample}_bwa{step}.flagstat.txt'
#    shell:
#        """
#        samtools flagstat {input} > {output}
#        """

rule collect_flagstat:
    conda: 'base'
    input:
        align =       'Bam/{sample}_bwa.flagstat.txt',
        mapq =        'Bam/{sample}_bwa_mapq.flagstat.txt',
        r1_map =      'Bam/{sample}_bwa_mapq_R1.srt.flagstat.txt',
        r2_map =      'Bam/{sample}_bwa_mapq_R2.srt.flagstat.txt',
        r1_filtered = 'Bam/{sample}_bwa_mapq_R1.multi_filt.flagstat.txt',
        r2_filtered = 'Bam/{sample}_bwa_mapq_R2.multi_filt.flagstat.txt',
        r1_dupsRem =  'Bam/{sample}_bwa_mapq_R1.multi_filt.dedup.flagstat.txt',
        r2_dupsRem =  'Bam/{sample}_bwa_mapq_R2.multi_filt.dedup.flagstat.txt',
    output:
        'FastQC/{sample}_bwa_flagstat_summary.tsv'
    shell:
        """
        echo "step\tdescription\tcount1\tcount1_percent\tcount2\tcount2_percent\ttotal > {output}
        python3.12 Src/restruct_flagstat_v2.py -i {input.align} -o FastQC/{wildcards.sample}_align_flagstat_summary.tsv -s all_align
        python3.12 Src/restruct_flagstat_v2.py -i {input.mapq} -o FastQC/{wildcards.sample}_mapq_flagstat_summary.tsv -s all_mapq
        python3.12 Src/restruct_flagstat_v2.py -i {input.r1_map} -o FastQC/{wildcards.sample}_R1_mapq_flagstat_summary.tsv -s r1_mapq
        python3.12 Src/restruct_flagstat_v2.py -i {input.r2_map} -o FastQC/{wildcards.sample}_R2_mapq_flagstat_summary.tsv -s r2_mapq
        python3.12 Src/restruct_flagstat_v2.py -i {input.r1_filtered} -o FastQC/{wildcards.sample}_R1_multi_filter_flagstat_summary.tsv -s r1_multi_filter
        python3.12 Src/restruct_flagstat_v2.py -i {input.r2_filtered} -o FastQC/{wildcards.sample}_R2_multi_filter_flagstat_summary.tsv -s r2_multi_filter
        python3.12 Src/restruct_flagstat_v2.py -i {input.r1_dupsRem} -o FastQC/{wildcards.sample}_R1_dupsRem_flagstat_summary.tsv -s r1_dupsRem
        python3.12 Src/restruct_flagstat_v2.py -i {input.r2_dupsRem} -o FastQC/{wildcards.sample}_R2_dupsRem_flagstat_summary.tsv -s r2_dupsRem

        for f in FastQC/{wildcards.sample}_*_flagstat_summary.tsv; do
            tail -n +2 $f >> {output}
        done
        """

rule build_matrix:
    conda: 'hicexplorer_v2_2'
    input:
        readOne = 'Bam/{sample}_bwa_mapq_R1{filter}.bam',
        readTwo = 'Bam/{sample}_bwa_mapq_R2{filter}.bam',
    output:
        #'Matrix/{sample}_{aligner}.h5'
        'Matrix/{sample}_bwa_mapq{filter}.mcool'
    params:
        mse_seq = 'TTAA',
        nla_seq = 'CATG',
        mse_dangling = 'TA',
        nla_dangling = 'CATG',
        mse_cuts = config['index'][species]['mse_cuts'],
        nla_cuts = config['index'][species]['nla_cuts'],
        #mse_cuts = '/athena/apostoloulab/scratch/collab/Genome/Mus_Musculus/mm10_msei_rest_site_positions.bed',
        #nla_cuts = '/athena/apostoloulab/scratch/collab/Genome/Mus_Musculus/mm10_nlaIII_rest_site_positions.bed',
    threads: 8
    shell:
        """
        hicBuildMatrix --samFiles {input.readOne} {input.readTwo} \\
            --restrictionSequence {params.mse_seq} {params.nla_seq} \\
            --danglingSequence {params.mse_dangling} {params.nla_dangling} \\
            --restrictionCutFile {params.mse_cuts} {params.nla_cuts} \\
            --threads {threads} \\
            --inputBufferSize 400000 \\
            --minMappingQuality 1 \\
            -o {output} \\
            --QCfolder ./hicQC_{wildcards.sample}{wildcards.filter} \\
            --binSize 5000 10000 20000 50000 100000
        """

rule make_bigwig:
    conda: "deeptools_3_5_1",
    input:
        bam = 'Bam/{sample}_bwa_mapq_R1{filter}.bam',
        bam_ind = 'Bam/{sample}_bwa_mapq_R1{filter}.bam.bai',
    output:
        'BigWig/{sample}_bwa_mapq_R1{filter}rpgcNorm.bw',
    params:
        genomeSize = config['index'][species]['genomeSize'],
        #blacklist = config['index']['mm10']['blacklist'],
    shell:
        """
        bamCoverage --bam {input.bam} -o {output} -p 10 \\
          --binSize 10 --ignoreForNormalization chrX chrM --normalizeUsing RPGC \\
          --effectiveGenomeSize {params.genomeSize} --extendReads 200 --ignoreDuplicates
        """

rule fastqc:
    conda: 'fastqc'
    input:
        #fq_raw = unpack(getFastqFiles),
        unpack(getFastqFiles),
        r1_trim = 'Fastq/{sample}_R1_trim.fastq.gz',
        r2_trim = 'Fastq/{sample}_R2_trim.fastq.gz'
    output:
        'FastQC/{sample}_R1_trim_fastqc.html',
        'FastQC/{sample}_R2_trim_fastqc.html',
        'FastQC/{sample}_R1_trim_fastqc.zip',
        'FastQC/{sample}_R2_trim_fastqc.zip'
    threads: 8
    shell:
        """
        fastqc -t {threads} -o FastQC/ {input.r1} {input.r2}
        fastqc -t {threads} -o FastQC/ {input.r1_trim} {input.r2_trim}
        """

rule fastq_screen:
    conda: 'fastq_screen'
    input:
        'Fastq/{sample}_{read}_trim.fastq.gz',
    output:
        'FastQC/{sample}_{read}_trim_screen.html',
    params:
        conf = 'Config/fastq_screen.conf',
        fastqc_dir = 'FastQC/'
    threads: 8
    shell:
        """
        fastq_screen --conf {params.conf} --threads {threads} --outdir {params.fastqc_dir} {input}
        """

rule multiqc:
    conda: 'multiqc'
    input:
        expand('FastQC/{sample}_{read}_trim_fastqc.zip', sample = sampleList, read = ['R1', 'R2']),
        expand('FastQC/{sample}_{read}_trim_screen.html', sample = sampleList, read = ['R1', 'R2']),
    output:
        'FastQC/multiqc_report.html'
    params:
        fastqc_dir = 'FastQC/'
    shell:
        """
        multiqc {params.fastqc_dir} -o {params.fastqc_dir}
        """


#rule merge_bam:
#    conda: 'samtools'
#    input:
#        readOne = 'Bam/{sample}_{aligner}_R1.bam',
#        readTwo = 'Bam/{sample}_{aligner}_R2.bam',
#    output:
#        unSortBam = temp('Bam/{sample}_{aligner}_{merge}_unsorted.bam'),
#        nameSortBam = 'Bam/{sample}_{aligner}_{merge}_nameSort.bam',
#        posSortBam = 'Bam/{sample}_{aligner}_{merge}_posSort.bam'
#    threads: 8
#    shell:
#        """
#        samtools merge -@ {threads} {output.unSortBam} {input.readOne} {input.readTwo} 
#        samtools sort -@ {threads} -n {output.unSortBam} -o {output.posSortBam}
#        samtools sort -@ {threads} {output.unSortBam} -o {output.nameSortBam}
#        """
#
#rule convert_to_bed:
#    conda: 'bedtools'
#    input: 
#        'Bam/{sample}_{aligner}_merge_posSort.bam'
#    output:
#        'Bed/{sample}_{aligner}_merge.bedpe'
#    shell:
#        """
#        bedtools bamtobed -bedpe -i {input} > {output}
#        """

# hicFindRestSite

#rule make_mcool:
#    conda: 'hicexplorer_v2_2'
#    input:
#        'Matrix/{sample}_{aligner}.h5'
#    output:
#        'Matrix/{sample}_{aligner}.mcool'
#    params:
#        cs = '/home/fs01/cmu4001/chromSizes_hiGlass.txt'
#    shell:
#        """
#        hicConvertFormat -m {input} --inputFormat h5 --outputFormat mcool -o {output} -cs {params.cs} --resolutions 5000 10000 20000 50000 100000
#        """
#

#rule make_diagnostic_plot:
#    input:
#    params:
#    shell:
#        """
#        hicCorrectMatrix diagnostic_plot -m {input} -o {output}
#        """


#        hicBuildMatrix --samFiles {input.readOne} {input.readTwo} \


#                 --binSize 10000 \
#                 --restrictionSequence GATC \
#                 --danglingSequence GATC \
#                 --restrictionCutFile {params.cutSites} \
#                 --threads 4 \
#                 --inputBufferSize 100000 \
#                 --outBam hic.bam \
#                 -o hic_matrix.h5 \
#                 --QCfolder ./hicQC
#

#rule name_sort_reads:
#    input:
#
#    output:
#
#    shell:
#        """
#        source rep_aware_pipeline/condaInit.sh
#        conda activate {params.samtools}
#
#        samtools sort -@ {threads} -n {input} -o Bam/{wildcards.sample}_nameSort.final.bam
#
#        """
#"""
#STAR --runMode alignReads --genomeLoad NoSharedMemory --genomeDir {params.genomeDir} \\
#--readFilesIn {input} --readFilesCommand zcat --outFileNamePrefix Bam/{wildcards.sample}_trim_starLocal. \\
#--runThreadN {threads} --outSAMtype BAM SortedByCoordinate \\
#--outFilterMultimapNmax 5000 --outSAMmultNmax 1 \\
#--outFilterMismatchNmax 3 --outMultimapperOrder Random --winAnchorMultimapNmax 5000 \\
#--alignEndsType Local --alignIntronMax 1 --alignMatesGapMax 350 \\
#--seedSearchStartLmax 30 --alignTranscriptsPerReadNmax 30000 \\
#--alignWindowsPerReadNmax 30000 --alignTranscriptsPerWindowNmax 300 \\
#--seedPerReadNmax 3000 --seedPerWindowNmax 300 --seedNoneLociPerWindow 1000 \\
#--alignSJoverhangMin 999 --alignSJDBoverhangMin 999
#"""
#rule makeBigWig:
#    #### MAKE BIGWIGS FROM BAM FILES
#    # Requires: Deeptools
#    input:
#        'Bam/{sample}.final.bam',
#    output:
#        'BigWig/{sample}_1xCovNorm.bw'
#    params:
#        deeptools = config['software']['deeptools'],
#        genomeSize = config['index']['mm10']['genomeSize'],
#        blacklist = config['index']['mm10']['blacklist'],
#    shell:
#        """
#        source rep_aware_pipeline/condaInit.sh
#        conda activate {params.deeptools}
#
#        bamCoverage --bam {input} -o {output} -p 10 \\
#          --binSize 10 --ignoreForNormalization chrX chrM --normalizeUsing RPGC \\
#          --effectiveGenomeSize {params.genomeSize} --extendReads 200 --ignoreDuplicates \\
#          --blackListFileName {params.blacklist} 
#
#        """
#rule align_reads_star:
#    conda: "star_2_7_10"
#    input:
#        readOne = 'Fastq/{sample}_{read}_trim.fastq.gz',
#        #readTwo = 'Fastq/{sample}_R2_trim.fastq.gz',
#    output:
#        #star = 'Bam/{sample}_star_{read}.Aligned.out.bam',
#        newname = 'Bam/{sample}_star_{read}.bam',
#    params:
#        #genomeDir = config['index']['mm10']['starIndex'], 
#        star = config['software']['star'],
#        genomeDir = config['index']['mm10']['starRepAware']
#    threads: 1
#    shell:
#        """
#        STAR --runMode alignReads --genomeLoad NoSharedMemory --genomeDir {params.genomeDir} \\
#        --readFilesIn {input} \\
#        --readFilesCommand zcat --outFileNamePrefix Bam/{wildcards.sample}_star_{wildcards.read}. \\
#        --runThreadN 1 --outSAMtype BAM Unsorted \\
#        --outFilterMultimapNmax 5000 --outSAMmultNmax 1 \\
#        --outFilterMismatchNmax 3  --outMultimapperOrder Random --winAnchorMultimapNmax 5000 \\
#        --alignEndsType Local --alignIntronMax 1 --alignMatesGapMax 350 \\
#        --seedSearchStartLmax 30 --alignTranscriptsPerReadNmax 30000 \\
#        --outSAMunmapped Within \\
#        --alignWindowsPerReadNmax 30000 --alignTranscriptsPerWindowNmax 300 \\
#        --seedPerReadNmax 3000 --seedPerWindowNmax 300 --seedNoneLociPerWindow 1000 \\
#        --alignSJoverhangMin 999 --alignSJDBoverhangMin 999
#
#        mv Bam/{wildcards.sample}_star_{wildcards.read}.Aligned.out.bam {output.newname}
#rule resort_bam:
#    conda: 'samtools'
#    input:
#        'Bam/{sample}_{aligner}_{read}.bam',
#    output:
#        bam = temp('Bam/{sample}_{aligner}_{read}_coordSort.bam'),
#        bam_ind = temp('Bam/{sample}_{aligner}_{read}_coordSort.bam.bai')
#    threads: 4
#    shell:
#        """
#        samtools sort -@ {threads} {input} > {output.bam}
#        samtools index {output.bam}
#        """
        #expand('Matrix/{sample}_{aligner}.mcool', sample = sampleList, read = ['R1', 'R2'], aligner = ['bwa', 'star']),
        #expand('BigWig/{sample}_{aligner}_{read}_rpgcNorm.bw', sample = sampleList, read = ['R1', 'R2'], aligner = ['bwa', 'star']),
        #expand('Matrix_Mapq_Test/{sample}_bwa_{mapq}.mcool', sample = sampleList, mapq = ['1', '5', '10', '15']),
        #expand('Matrix_Mapq_Test/{sample}_star_{mapq}.mcool', sample = sampleList, mapq = ['1', '2', '3', '4', '5', '10', '15']),
        #expand('Bed/{sample}_{aligner}_merge.bedpe', sample = sampleList, aligner = ['bwa', 'star'], bamType = ['R1', 'R2', 'merge']),
        #expand('Bed/{sample}_{aligner}_{bamType}.bedpe', sample = sampleList, aligner = ['bwa', 'star'], bamType = ['R1', 'R2', 'merge']),
        #expand('Fastq/{sample}_{read}.fastq.gz', sample = sampleList, read = ['R1', 'R2']),
        #expand('Bam/{sample}_{aligner}_{read}.bam', sample = sampleList, read = ['R1', 'R2'], aligner = ['bwa', 'star']),
        #expand('Matrix/{sample}_{aligner}.h5', sample = sampleList, read = ['R1', 'R2'], aligner = ['bwa', 'star']),    
        #input_arr=({input})
        #for ((i=0; i<${{#input_arr[@]}}; i++)); do
        #    in_file=${{input_arr[$i]}}
        #    out_summary=${{in_file/Bam/FastQC}}
        #    out_summary=${{out_summary/.flagstat.txt/_mapStats.tsv}}
        #    step_name=$(basename $in_file .flagstat.txt | sed 's/{wildcards.sample}_bwa_//')
        #    python3.12 Src/restruct_flagstat.py -i $in_file -s $step_name -o $out_summary
        #    tail -n +2 $out_summary >> {output}
        #done
        #input_arr=({input})
        #for ((i=0; i<${{#input_arr[@]}}; i++)); do
        #    in_file=${{input_arr[$i]}}
        #    step_name = $(basename $in_file .flagstat.txt | sed 's/{wildcards.sample}_bwa_//')
        #done
        #python3.12 Src/restruct_flagstat.py -i {input.align} {input.mapq} {input.r1_map} {input.r2_map} {input.r1_filtered} {input.r2_filtered} {input.r1_dupsRem} {input.r2_dupsRem} -o {output}

        #samtools flagstat {input.align} > FastQC/{wildcards.sample}_align_flagstat.txt
        #samtools flagstat {input.mapq} > FastQC/{wildcards.sample}_mapq_flagstat.txt
        #samtools flagstat {input.r1_map} > FastQC/{wildcards.sample}_r1_mapq_flagstat.txt
        #samtools flagstat {input.r2_map} > FastQC/{wildcards.sample}_r2_map2_flagstat.txt
        #samtools flagstat {input.r1_filtered} > FastQC/{wildcards.sample}_flagstat.txt
        #samtools flagstat {input.r2_filtered} > FastQC/{wildcards.sample}_flagstat.txt

        #align_index =       'Bam/{sample}_bwa.bam.bai',
        #mapq_index =        'Bam/{sample}_bwa_mapq.bam.bai',
        #r1_index =          'Bam/{sample}_bwa_mapq_R1.srt.bam.bai',
        #r2_index =          'Bam/{sample}_bwa_mapq_R2.srt.bam.bai',
        #r1_filtered_index = 'Bam/{sample}_bwa_R1.filtered.bam.bai',
        #r2_filtered_index = 'Bam/{sample}_bwa_R2.filtered.bam.bai',
        #r1_dupsRem_index =  'Bam/{sample}_bwa_mapq_R1.filter.dedup.bam.bai',
        #r2_dupsRem_index =  'Bam/{sample}_bwa_mapq_R2.filter.dedup.bam.bai',#        """
#rule build_matrix_mapq_test:
#    conda: 'hicexplorer_v2_2'
#    input:
#        readOne = 'Bam/{sample}_{aligner}_R1.bam',
#        readTwo = 'Bam/{sample}_{aligner}_R2.bam',
#    output:
#        #'Matrix/{sample}_{aligner}.h5'
#        'Matrix_Mapq_Test/{sample}_{aligner}_{mapq}.mcool'
#    params:
#        rd_seq = 'TTAA',
#        dangling_seq = 'TA',
#        cutSites = '/athena/apostoloulab/scratch/collab/Genome/Mus_Musculus/mm10_msei_rest_site_positions.bed',
#        nThreads = 8,
#    shell:
#        """
#        if [[ ! -d ./hicQC/{wildcards.sample}_{wildcards.aligner}_{wildcards.mapq} ]]; then
#            mkdir ./hicQC/{wildcards.sample}_{wildcards.aligner}_{wildcards.mapq}
#        fi
#
#        hicBuildMatrix --samFiles {input.readOne} {input.readTwo} \\
#            --restrictionSequence {params.rd_seq} \\
#            --danglingSequence {params.dangling_seq} \\
#            --restrictionCutFile {params.cutSites} \\
#            --threads {params.nThreads} \\
#            --inputBufferSize 400000 \\
#            --minMappingQuality {wildcards.mapq} \\
#            -o {output} \\
#            --QCfolder ./hicQC/{wildcards.sample}_{wildcards.aligner}_{wildcards.mapq} \\
#            --binSize 5000 10000 20000 50000 100000
#
#        """
#rule run_flagstat_initial:
#    conda: 'samtools'
#    input:
#        'Bam/{sample}_bwa.bam'
#    output:
#        'Bam/{sample}_bwa.flagstat.txt'
#    shell:
#        """
#        samtools flagstat {input} > {output}
#        """
#