#!/usr/bin/env python3

import os
import re
import sys
import tempfile
import pandas as pd
from snake_setup import set_config, load_samples

BASE = workflow.basedir

# Define path to conda environment specifications
ENVS = f'{BASE}/workflow/envs'
# Defne path to custom scripts directory
SCRIPTS = f'{BASE}/workflow/scripts'

if not config:
    configfile: f'{BASE}/config/config.yaml'

# Defaults configuration file - use empty string to represent no default value.
default_config = {
    'workdir':           workflow.basedir,
    'tmpdir':            tempfile.gettempdir(),
    'threads':           workflow.cores,
    'data':              ''          ,
    'paired':            ''          ,
    'computeMatrixScale':
        {'exon':         True        ,},
    'genome':
        {'build':        'genome'    ,
         'index':        ''          ,
         'annotation':   ''          ,
         'sequence':     ''          ,
         'genes':        ''          ,
         'blacklist':    None        ,},
    'cutadapt':
        {'forwardAdapter': 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA',
         'reverseAdapter': 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT',
         'overlap':         3                                 ,
         'errorRate':       0.1                               ,
         'minimumLength':   0                                 ,
         'qualityCutoff':  '0,0'                              ,
         'GCcontent':       50                                ,},
    'fastq_screen':         None                              ,
}

config = set_config(config, default_config)

workdir: config['workdir']
BUILD = config['genome']['build']

# Read path to samples in pandas
samples = load_samples(config['data'])

# Extract sample names
SAMPLES = list(samples['sample'].unique())

wildcard_constraints:
    single = '|'.join(samples['single']),
    sample = '|'.join(samples['sample']),
    rep = '|'.join(samples['rep'].unique())

rule all:
    input:
        ['qc/multiqc', 'counts/countMatrix-merged.txt']

rule fastqc:
    input:
        lambda wc: samples.xs(wc.single, level = 2)['path']
    output:
        html = 'qc/fastqc/{single}.raw_fastqc.html',
        zip = 'qc/fastqc/unmod/{single}.raw.fastqc.zip'
    group:
        'processFASTQ'
    log:
        'logs/fastqc/{single}.log'
    wrapper:
        '0.49.0/bio/fastqc'


rule modifyFastQC:
    input:
        'qc/fastqc/unmod/{single}.raw.fastqc.zip'
    output:
        'qc/fastqc/{single}.raw_fastqc.zip'
    params:
        name = lambda wc: f'{wc.single}'
    group:
        'processFASTQ'
    log:
        'logs/modifyFastQC/{single}.raw.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/modifyFastQC.py {input} {output} {params.name} &> {log}'


def cutadaptOutput():
    if config['paired']:
        return ['fastq/trimmed/{sample}-R1.trim.fastq.gz',
                'fastq/trimmed/{sample}-R2.trim.fastq.gz']
    else:
        return ['fastq/trimmed/{sample}-R1.trim.fastq.gz']


def cutadaptCmd():
    if config['paired']:
        cmd = ('cutadapt -a {params.forwardAdapter} -A {params.reverseAdapter} '
            '-o {output.trimmed[0]} -p {output.trimmed[1]} ')
    else:
        cmd = 'cutadapt -a {params.forwardAdapter} -o {output.trimmed[0]} '

    cmd += ('--overlap {params.overlap} --error-rate {params.errorRate} '
        '--minimum-length {params.minimumLength} '
        '--quality-cutoff {params.qualityCutoff} '
        '--gc-content {params.GCcontent} '
        '--cores {threads} {input} > {output.qc} 2> {log}')
    return cmd


rule cutadapt:
    input:
        lambda wc: samples.xs(wc.sample, level=1)['path']
    output:
        trimmed = cutadaptOutput(),
        qc = 'qc/cutadapt/unmod/{sample}.cutadapt.txt'
    params:
        forwardAdapter = config['cutadapt']['forwardAdapter'],
        reverseAdapter = config['cutadapt']['reverseAdapter'],
        overlap = config['cutadapt']['overlap'],
        errorRate = config['cutadapt']['errorRate'],
        minimumLength = config['cutadapt']['minimumLength'],
        qualityCutoff = config['cutadapt']['qualityCutoff'],
        GCcontent = config['cutadapt']['GCcontent']
    group:
        'processFASTQ'
    log:
        'logs/cutadapt/{sample}.log'
    conda:
        f'{ENVS}/cutadapt.yaml'
    threads:
        config['threads']
    shell:
        cutadaptCmd()


rule modifyCutadapt:
    input:
        rules.cutadapt.output.qc
    output:
        'qc/cutadapt/{sample}.cutadapt.txt'
    group:
        'processFASTQ'
    log:
        'logs/modifyCutadapt/{sample}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/modifyCutadapt.py {wildcards.sample} {input} '
        '> {output} 2> {log}'

if config['fastq_screen']:
    rule fastqScreen:
        input:
            'fastq/trimmed/{single}.trim.fastq.gz'
        output:
            txt = 'qc/fastq_screen/{single}.fastq_screen.txt',
            png = 'qc/fastq_screen/{single}.fastq_screen.png'
        params:
            fastq_screen_config = config['fastq_screen'],
            subset = 100000,
            aligner = 'bowtie2'
        group:
            'processFASTQ'
        log:
            'logs/fastq_screen/{single}.log'
        threads:
            config['threads']
        wrapper:
            '0.49.0/bio/fastq_screen'


rule fastQCtrimmed:
    input:
        'fastq/trimmed/{single}.trim.fastq.gz'
    output:
        html = 'qc/fastqc/{single}.trim_fastqc.html',
        zip = 'qc/fastqc/{single}.trim_fastqc.zip'
    group:
        'processFASTQ'
    log:
        'logs/fastqc_trimmed/{single}.log'
    wrapper:
        '0.49.0/bio/fastqc'


def hisat2Cmd():
    if config['paired']:
        return ('hisat2 -x {params.index} -1 {input.reads[0]} '
            '-2 {input.reads[1]} --threads {threads} '
            '--summary-file {output.qc} > {output.sam} 2> {log}')
    else:
        return ('hisat2 -x {params.index} -p {threads} '
            '-U {input.reads[0]} --summary-file {output.qc} '
            '> {output.sam} 2> {log}')


rule hisat2:
    input:
        reads = rules.cutadapt.output.trimmed
    output:
        sam = pipe('mapped/{sample}.sam'),
        qc = 'qc/hisat2/{sample}.hisat2.txt'
    params:
        index = config['genome']['index']
    log:
        'logs/hisat2/{sample}.log'
    conda:
        f'{ENVS}/hisat2.yaml'
    threads:
        max(1, (config['threads'] / 2) - 1)
    shell:
        hisat2Cmd()


rule fixBAM:
    input:
        rules.hisat2.output.sam
    output:
        pipe('mapped/{sample}.fixmate.bam')
    log:
        'logs/fixBAM/{sample}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools fixmate -O bam,level=0 -m {input} {output} &> {log}'


rule sortBAM:
    input:
        rules.fixBAM.output
    output:
        'mapped/{sample}.sort.bam'
    log:
        'logs/sortBAM/{sample}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    threads:
        max(1, (config['threads'] / 2) - 1)
    shell:
        'samtools sort -@ {threads} {input} > {output} 2> {log}'


rule markdupBAM:
    input:
        rules.sortBAM.output
    output:
        bam = 'mapped/{sample}.markdup.bam',
        qc = 'qc/deduplicate/{sample}.txt'
    log:
        'logs/markdupBAM/{sample}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    threads:
        config['threads']
    shell:
        'samtools markdup -@ {threads} '
        '-sf {output.qc} {input} {output.bam} &> {log}'


rule indexBAM:
    input:
        'mapped/{sample}.{stage}.bam'
    output:
        'mapped/{sample}.{stage}.bam.bai'
    log:
        'logs/indexBAM/{sample}-{stage}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    threads:
        min(1, (config['threads'] / 2) - 1,)
    shell:
        'samtools index -@ {threads} {input} &> {log}'


rule samtoolsStats:
    input:
        rules.markdupBAM.output.bam
    output:
        'qc/samtools/stats/{sample}.stats.txt'
    group:
        'samQC'
    log:
        'logs/samtools_stats/{sample}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools stats {input} > {output} 2> {log}'


rule samtoolsIdxstats:
    input:
        bam = 'mapped/{sample}.markdup.bam',
        index = 'mapped/{sample}.markdup.bam.bai'
    output:
        'qc/samtools/idxstats/{sample}.idxstats.txt'
    group:
        'samQC'
    log:
        'logs/samtools_idxstats/{sample}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools idxstats {input.bam} > {output} 2> {log}'


rule samtoolsFlagstat:
    input:
        rules.markdupBAM.output.bam
    output:
        'qc/samtools/flagstat/{sample}.flagstat.txt'
    group:
        'samQC'
    log:
        'logs/samtoolsFlagstat/{sample}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools flagstat {input} > {output} 2> {log}'


def setBlacklistCommand():
    if config['genome']['blacklist']:
        cmd = ('bedtools merge -d {params.distance} -i {input} '
               '> {output} 2> {log}')
    else:
        cmd = 'touch {output} &> {log}'
    return cmd


rule processBlacklist:
    input:
        config['genome']['blacklist'] if config['genome']['blacklist'] else []
    output:
        f'genome/{BUILD}-blacklist.bed'
    params:
        # Max distance between features allowed for features to be merged.
        distance = 1000
    log:
        f'logs/processBlacklist/{BUILD}.log'
    conda:
        f'{ENVS}/bedtools.yaml'
    shell:
        setBlacklistCommand()


rule estimateReadFiltering:
    input:
        bam = rules.markdupBAM.output.bam,
        index = 'mapped/{sample}.markdup.bam.bai',
        blacklist = rules.processBlacklist.output
    output:
        'qc/deeptools/estimateReadFiltering/{sample}.txt'
    params:
        minMapQ = 15,
        binSize = 10000,
        distanceBetweenBins = 0,
        properPair = '--samFlagInclude 2' if config['paired'] else '',
    log:
        'logs/estimateReadFiltering/{sample}.log'
    conda:
        f'{ENVS}/deeptools.yaml'
    threads:
        config['threads']
    shell:
        'estimateReadFiltering --bamfiles {input.bam} --outFile {output} '
        '--binSize {params.binSize} --blackListFileName {input.blacklist} '
        '--distanceBetweenBins {params.distanceBetweenBins} '
        '--minMappingQuality {params.minMapQ} --ignoreDuplicates '
        '--samFlagExclude 260 {params.properPair} '
        '--numberOfProcessors {threads} &> {log}'


rule alignmentSieve:
    input:
        bam = rules.markdupBAM.output.bam,
        index = 'mapped/{sample}.markdup.bam.bai',
        blacklist = rules.processBlacklist.output
    output:
        bam = 'mapped/{sample}.filtered.bam',
        qc = 'qc/deeptools/{sample}-filter-metrics.txt'
    params:
        minMapQ = 15,
        maxFragmentLength = 2000,
        properPair = '--samFlagInclude 2' if config['paired'] else '',
    log:
        'logs/alignmentSieve/{sample}.log'
    conda:
        f'{ENVS}/deeptools.yaml'
    threads:
        config['threads']
    shell:
        'alignmentSieve --bam {input.bam} --outFile {output.bam} '
        '--minMappingQuality {params.minMapQ} --ignoreDuplicates '
        '--samFlagExclude 260 {params.properPair} '
        '--maxFragmentLength {params.maxFragmentLength} '
        '--blackListFileName {input.blacklist} '
        '--numberOfProcessors {threads} --filterMetrics {output.qc} &> {log}'


rule multiBamSummary:
    input:
        bams = expand('mapped/{sample}.filtered.bam', sample=SAMPLES),
        indexes = expand('mapped/{sample}.filtered.bam.bai', sample=SAMPLES)
    output:
        'qc/deeptools/multiBamSummary.npz'
    params:
        binSize = 10000,
        distanceBetweenBins = 0,
        labels = ' '.join(SAMPLES),
        extendReads = 150
    log:
        'logs/multiBamSummary.log'
    conda:
        f'{ENVS}/deeptools.yaml'
    threads:
        config['threads']
    shell:
        'multiBamSummary bins --bamfiles {input.bams} --outFileName {output} '
        '--binSize {params.binSize} --labels {params.labels} '
        '--distanceBetweenBins {params.distanceBetweenBins} '
        '--extendReads {params.extendReads} --numberOfProcessors {threads} &> {log}'


rule plotCorrelation:
    input:
        rules.multiBamSummary.output
    output:
        plot = 'qc/deeptools/plotCorrelation.png',
        matrix = 'qc/deeptools/plotCorrelation.tsv'
    params:
        corMethod = 'pearson',
        colorMap = 'viridis',
        labels = ' '.join(SAMPLES)
    log:
        'logs/plotCorrelation.log'
    conda:
        f'{ENVS}/deeptools.yaml'
    shell:
        'plotCorrelation --corData {input} --corMethod {params.corMethod} '
        '--whatToPlot heatmap --labels {params.labels} --skipZeros '
        '--colorMap {params.colorMap} --plotNumbers '
        '--plotFile {output.plot} --outFileCorMatrix {output.matrix} &> {log}'


def setColours(samples):
    """ Find group and assign to specific colour. """
    colours = ''
    colourPool = ['#E69F00', '#56B4E9', '#009E73', '#F0E442',
                  '#0072B2', '#D55E00', '#CC79A7']
    # Add double quotes for compatibility with shell
    colourPool = [f'"{colour}"' for colour in colourPool]
    usedColours = {}
    for sample in samples:
        group = sample.split('-')[0]
        if group not in usedColours:
            usedColours[group] = colourPool[0]
            # Remove from pool
            colourPool = colourPool[1:]
        colours += f'{usedColours[group]} '
    return f'{colours}'


rule plotPCA:
    input:
        rules.multiBamSummary.output
    output:
        plot = 'qc/deeptools/plotPCA.png',
        data = 'qc/deeptools/plotPCA.tab'
    params:
        labels = ' '.join(SAMPLES),
        colours = setColours(SAMPLES)
    log:
        'logs/plotPCA.log'
    conda:
        f'{ENVS}/deeptools.yaml'
    shell:
        'plotPCA --corData {input} --colors {params.colours} '
        '--labels {params.labels} --plotFile {output.plot} '
        '--outFileNameData {output.data} &> {log}'


rule plotCoverage:
    input:
        bams = expand('mapped/{sample}.filtered.bam', sample=SAMPLES),
        indexes = expand('mapped/{sample}.filtered.bam.bai', sample=SAMPLES)
    output:
        plot = 'qc/deeptools/plotCoverage.png',
        data = 'qc/deeptools/plotCoverage.tab',
        info = 'qc/deeptools/plotCoverage.info'
    params:
        nSamples = 1000000,
        labels = ' '.join(SAMPLES)
    log:
        'logs/plotCoverage.log'
    conda:
        f'{ENVS}/deeptools.yaml'
    threads:
        config['threads']
    shell:
        'plotCoverage --bamfiles {input.bams} --labels {params.labels} '
        '--plotFile {output.plot} --outRawCounts {output.data} '
        '--numberOfSamples {params.nSamples} --numberOfProcessors {threads} '
        '> {output.info} 2> {log}'


rule bamCoverage:
    input:
        bam  = rules.alignmentSieve.output.bam,
        index  = f'{rules.alignmentSieve.output.bam}.bai'
    output:
        'bigwig/{sample}.filtered.bigwig'
    params:
        binSize = 50,
    log:
        'logs/bamCoverage/{sample}.log'
    conda:
        f'{ENVS}/deeptools.yaml'
    threads:
        config['threads']
    shell:
        'bamCoverage --bam {input.bam} --outFileName {output} '
        '--binSize {params.binSize} --numberOfProcessors {threads} &> {log}'


rule computeMatrixScaledGroups:
    input:
        expand('bigwig/{sample}.filtered.bigwig', sample=SAMPLES)
    output:
        scaledGZ = 'deeptools/computeMatrixGroups/matrix-scaled.gz',
        scaled = 'deeptools/computeMatrixGroups/matrix-scaled.tab',
        sortedRegions = 'deeptools/computeMatrixGroups/genes-scaled.bed'
    params:
        binSize = 50,
        regionBodyLength = 5000,
        upstream = 2000,
        downstream = 2000,
        metagene = '--metagene' if config['computeMatrixScale']['exon'] else '',
        samplesLabel = ' '.join(SAMPLES),
        genes = config['genome']['genes'],
        averageType = 'mean'
    log:
        'logs/computeMatrixScaledGroups.log'
    conda:
        f'{ENVS}/deeptools.yaml'
    threads:
        config['threads']
    shell:
        'computeMatrix scale-regions --scoreFileName {input} '
        '--regionsFileName {params.genes} --outFileName {output.scaledGZ} '
        '--outFileNameMatrix {output.scaled} --skipZeros '
        '--regionBodyLength {params.regionBodyLength} '
        '--samplesLabel {params.samplesLabel} --binSize {params.binSize} '
        '--averageTypeBins {params.averageType} {params.metagene} '
        '--upstream {params.upstream} --downstream {params.downstream} '
        '--outFileSortedRegions {output.sortedRegions} '
        '--numberOfProcessors {threads} &> {log}'


rule computeMatrixReferenceGroups:
    input:
        expand('bigwig/{sample}.filtered.bigwig', sample=SAMPLES)
    output:
        referenceGZ = 'deeptools/computeMatrixGroups/matrix-reference.gz',
        reference = 'deeptools/computeMatrixGroups/matrix-reference.tab',
        sortedRegions = 'deeptools/computeMatrixGroups/genes-reference.bed'
    params:
        binSize = 50,
        upstream = 5000,
        downstream = 5000,
        samplesLabel = ' '.join(SAMPLES),
        genes = config['genome']['genes'],
        averageType = 'mean',
        referencePoint = 'TSS'
    log:
        'logs/computeMatrixReferenceGroups.log'
    conda:
        f'{ENVS}/deeptools.yaml'
    threads:
        config['threads']
    shell:
        'computeMatrix reference-point --scoreFileName {input} '
        '--regionsFileName {params.genes} --outFileName {output.referenceGZ} '
        '--outFileNameMatrix {output.reference} --skipZeros '
        '--samplesLabel {params.samplesLabel} --binSize {params.binSize} '
        '--averageTypeBins {params.averageType} '
        '--referencePoint {params.referencePoint} '
        '--upstream {params.upstream} --downstream {params.downstream} '
        '--outFileSortedRegions {output.sortedRegions} '
        '--numberOfProcessors {threads} &> {log}'


rule plotProfileGroups:
    input:
        'deeptools/computeMatrixGroups/matrix-{mode}.gz'
    output:
        plot = 'qc/deeptools/compareGroup/plotProfileGroup-{mode}.png',
        data = 'qc/deeptools/compareGroup/plotProfileGroup-data-{mode}.tab',
        bed = 'qc/deeptools/compareGroup/plotProfileGroup-regions-{mode}.bed'
    params:
        dpi = 300,
        plotsPerRow = 2,
        averageType = 'mean',
        referencePoint = 'TSS'
    log:
        'logs/plotProfileGroups-{mode}.log'
    conda:
        f'{ENVS}/deeptools.yaml'
    shell:
        'plotProfile --matrixFile {input} --outFileName {output.plot} '
        '--outFileSortedRegions {output.bed} --outFileNameData {output.data} '
        '--dpi {params.dpi} --averageType {params.averageType} '
        '--refPointLabel {params.referencePoint} '
        '--numPlotsPerRow {params.plotsPerRow} &> {log}'


rule plotHeatmapGroups:
    input:
        'deeptools/computeMatrixGroups/matrix-{mode}.gz'
    output:
        plot = 'qc/deeptools/compareGroup/plotHeatmapGroup-{mode}.png',
        data = 'qc/deeptools/compareGroup/plotHeatmapGroup-data-{mode}.tab',
        bed = 'qc/deeptools/compareGroup/plotHeatmapGroup-regions-{mode}.bed'
    params:
        dpi = 300,
        zMax = 3,
        zMin = -3,
        kmeans = 3,
        width = max(4, len(SAMPLES) * 2),
        colorMap = 'RdBu_r',
        averageType = 'mean',
        referencePoint = 'TSS',
        interpolationMethod = 'auto'
    log:
        'logs/plotHeatmapGroups-{mode}.log'
    conda:
        f'{ENVS}/deeptools.yaml'
    threads:
        config['threads']
    shell:
        'plotHeatmap --matrixFile {input} --outFileName {output.plot} '
        '--outFileSortedRegions {output.bed} --outFileNameMatrix {output.data} '
        '--interpolationMethod {params.interpolationMethod} '
        '--dpi {params.dpi} --kmeans {params.kmeans} '
        '--zMin {params.zMin} --zMax {params.zMax} '
        '--colorMap {params.colorMap} --refPointLabel {params.referencePoint} '
        '--averageTypeSummaryPlot {params.averageType} '
        '--heatmapWidth {params.width} '
        '--refPointLabel {params.referencePoint} &> {log}'


rule featureCounts:
    input:
        bam = rules.sortBAM.output,
        gtf = config['genome']['annotation']
    output:
        counts = 'counts/{sample}.featureCounts.txt',
        qc = 'qc/featurecounts/{sample}.featureCounts.txt.summary'
    log:
        'logs/featureCounts/{sample}.log'
    conda:
        f'{ENVS}/subread.yaml'
    shell:
        '(featureCounts --primary -C -t exon -g gene_id '
        '-a {input.gtf} -o {output.counts} {input.bam} '
        '&& mv {output.counts}.summary {output.qc}) &> {log}'


rule mergeFeatureCounts:
    input:
        expand('counts/{sample}.featureCounts.txt', sample=SAMPLES)
    output:
        'counts/countMatrix-merged.txt'
    params:
        samples = SAMPLES
    log:
        'logs/mergeFeatureCounts.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/mergeCounts.py {input} --samples {params.samples} '
        '> {output} 2> {log}'


rule multiQC:
    input:
        expand('qc/fastqc/{single}.raw_fastqc.zip',
            single=samples['single']),
        expand('qc/fastq_screen/{single}.fastq_screen.txt',
            single=samples['single']) if config['fastq_screen'] else [],
        expand('qc/cutadapt/{sample}.cutadapt.txt',
            sample=SAMPLES),
        expand('qc/fastqc/{single}.trim_fastqc.zip',
            single=samples['single']),
        expand('qc/hisat2/{sample}.hisat2.txt',
            sample=SAMPLES),
        #expand('qc/deeptools/estimateReadFiltering/{sample}.txt',
        #    sample=SAMPLES),
        'qc/deeptools/plotCorrelation.tsv',
        'qc/deeptools/plotPCA.tab',
        'qc/deeptools/plotCoverage.info',
        'qc/deeptools/plotCoverage.tab',
        'qc/deeptools/compareGroup/plotProfileGroup-data-scaled.tab',
        expand('qc/samtools/stats/{sample}.stats.txt',
            sample=SAMPLES),
        expand('qc/samtools/idxstats/{sample}.idxstats.txt',
            sample=SAMPLES),
        expand('qc/samtools/flagstat/{sample}.flagstat.txt',
            sample=SAMPLES),
        expand('qc/featurecounts/{sample}.featureCounts.txt.summary',
            sample=SAMPLES)
    output:
        directory('qc/multiqc')
    log:
        'logs/multiqc/multiqc.log'
    conda:
        f'{ENVS}/multiqc.yaml'
    shell:
        'multiqc --outdir {output} '
        '--force --config {BASE}/config/multiqc_config.yaml {input} '
        '&> {log}'
