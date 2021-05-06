#!/usr/bin/env python3

import os
import re
import sys
import tempfile
import pandas as pd
from snake_setup import set_config, Samples

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
    'strand':            ''          ,
    'sample':            0           ,
    'QConly':            False       ,
    'genome':
        {'build':         'genome'    ,
         'transcriptome': ''          ,
         'gff3':          ''          ,
         'sequence':      ''          ,
         'index':         None        ,},
    'cutadapt':
        {'forwardAdapter': 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA',
         'reverseAdapter': 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT',
         'overlap':         3                                 ,
         'errorRate':       0.1                               ,
         'minimumLength':   1                                 ,
         'qualityCutoff':  '0,0'                              ,
         'GCcontent':       50                                ,},
    'fastqScreen':         None                              ,
}

config = set_config(config, default_config)

workdir: config['workdir']
BUILD = config['genome']['build']

RNASeq = Samples(config['data'], config['paired'])

# Extract sample names
SAMPLES = RNASeq.samples()

wildcard_constraints:
    sample = '|'.join(RNASeq.samples())

# Create tmpdir to ensure it is set
os.makedirs(config['tmpdir'], exist_ok=True)

# Minimum length must be atleast 1 (if 0 this break seqkt)
if config['cutadapt']['minimumLength'] < 1:
    config['cutadapt']['minimumLength'] = 1

if config['paired']:
    assert config['strand'] in ['RF', 'FR', 'unstranded']
    reads = ['R1', 'R2']
else:
    assert config['strand'] in ['R', 'F', 'unstranded']
    reads = ['R1']

rule all:
    input:
        'qc/multiqc'


rule bgzipGenome:
    input:
        config['genome']['sequence']
    output:
        f'genome/{config["genome"]["build"]}.fa.gz'
    group:
        'prepareGenome'
    log:
        'logs/bgzipGenome.log'
    conda:
        f'{ENVS}/tabix.yaml'
    shell:
        '(zcat -f {input} | bgzip > {output}) 2> {log}'


rule indexGenome:
    input:
        rules.bgzipGenome.output
    output:
        f'{rules.bgzipGenome.output}.fai'
    group:
        'prepareGenome'
    log:
        'logs/indexGenome.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools faidx {input} &> {log}'


rule getChromSizes:
    input:
        rules.indexGenome.output
    output:
        f'genome/{config["genome"]["build"]}.chrom.sizes'
    group:
        'prepareGenome'
    log:
        'logs/getChromSizes.log'
    shell:
        'cut -f 1,2 {input} > {output} 2> {log}'


rule fastqc:
    input:
        lambda wc: RNASeq.path(wc.sample, [wc.read])
    output:
        html = 'qc/fastqc/{sample}-{read}.raw_fastqc.html',
        zip = 'qc/fastqc/unmod/{sample}-{read}.raw.fastqc.zip'
    group:
        'processFASTQ'
    log:
        'logs/fastqc/{sample}-{read}.log'
    wrapper:
        '0.49.0/bio/fastqc'


rule modifyFastQC:
    input:
        'qc/fastqc/unmod/{sample}-{read}.raw.fastqc.zip'
    output:
        'qc/fastqc/{sample}-{read}.raw_fastqc.zip'
    group:
        'processFASTQ'
    log:
        'logs/modifyFastQC/{sample}-{read}.raw.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/modifyFastQC.py {input} {output} '
        '{wildcards.sample}-{wildcards.read} &> {log}'


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
        lambda wc: RNASeq.path(wc.sample, reads)
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


if config['fastqScreen']:
    rule fastqScreen:
        input:
            'fastq/trimmed/{sample}-{read}.trim.fastq.gz'
        output:
            txt = 'qc/fastq_screen/{sample}-{read}.fastq_screen.txt',
            png = 'qc/fastq_screen/{sample}-{read}.fastq_screen.png'
        params:
            fastq_screen_config = config['fastqScreen'],
            subset = 100000,
            aligner = 'bowtie2'
        group:
            'processFASTQ'
        log:
            'logs/fastqScreen/{sample}-{read}.log'
        threads:
            config['threads']
        wrapper:
            '0.49.0/bio/fastq_screen'


rule fastQCtrimmed:
    input:
        'fastq/trimmed/{sample}-{read}.trim.fastq.gz'
    output:
        html = 'qc/fastqc/{sample}-{read}.trim_fastqc.html',
        zip = 'qc/fastqc/{sample}-{read}.trim_fastqc.zip'
    group:
        'processFASTQ'
    log:
        'logs/fastqc_trimmed/{sample}-{read}.log'
    wrapper:
        '0.49.0/bio/fastqc'


rule buildKallistoIndex:
    input:
        config['genome']['transcriptome']
    output:
        f'kallisto/{config["genome"]["build"]}-transcripts.idx'
    log:
        'logs/buildKallistoIndex.log'
    conda:
        f'{ENVS}/kallisto.yaml'
    shell:
        'kallisto index -i {output} {input} &> {log}'


def kallistoCmd():
    if config['paired']:
        return ('kallisto quant -i {input.index} -o kallisto/{wildcards.sample}/ '
                '-b {params.bootstraps} {params.strand} {input.reads} &> {log} '
                '&& cp {log} {output.qc}')
    else:
        return ('kallisto quant -i {input.index} -o kallisto/{wildcards.sample}/ '
                '-b {params.bootstraps} {params.strand} '
                '--single -l {params.fragmentLength} -s {params.fragmentSD} '
                '--threads {threads} {input.reads} &> {log} '
                '&& cp {log} {output.qc}')


def getKallistoStrand():
    if config['strand'] in ['F', 'FR']:
        return '--fr-stranded'
    elif config['strand'] in ['R', 'RF']:
        return '--rf-stranded'
    else:
        return ''


rule kallistoQuant:
    input:
        index = rules.buildKallistoIndex.output,
        reads = rules.cutadapt.output.trimmed
    output:
        qc = 'qc/kallisto/{sample}.stdout',
        abundancesH5 = 'kallisto/{sample}/abundance.h5',
        abundancesTSV = 'kallisto/{sample}/abundance.tsv',
        runInfo = 'kallisto/{sample}/run_info.json'
    params:
        bootstraps = 100,
        fragmentLength = 200,
        fragmentSD = 2,
        strand = getKallistoStrand()
    log:
        'logs/kallistoQuant/{sample}.log'
    conda:
        f'{ENVS}/kallisto.yaml'
    threads:
        config['threads']
    shell:
        kallistoCmd()


def setBlacklistCommand():
    if config['filtering']['blackListFileName']:
        cmd = ('bedtools merge -d {params.distance} -i {input} '
               '> {output} 2> {log}')
    else:
        cmd = 'touch {output} &> {log}'
    return cmd


rule sampleReads:
    input:
        'fastq/trimmed/{sample}-{read}.trim.fastq.gz'
    output:
        'fastq/sampled/{sample}-{read}.trim.fastq.gz'
    params:
        seed = 42,
        nReads = config['sample']
    log:
        'logs/sampleReads/{sample}-{read}.log'
    conda:
        f'{ENVS}/seqtk.yaml'
    shell:
        '(seqtk sample -s {params.seed} {input} {params.nReads} '
        '| gzip > {output}) 2> {log}'


def hisat2Cmd():
    if config['paired']:
        command = 'hisat2 -1 {input[0]} -2 {input[1]} '
    else:
        command = 'hisat2 -U {input[0]} '
    if config['strand'] != 'unstranded':
        command += '--rna-strandness {params.strand} '
    command += ('-x {params.index} --threads {threads} '
                '--summary-file {output.qc} > {output.sam} 2> {log}')
    return command


def hisat2Input(wc):
    if config['sample'] > 0:
        return expand(
            'fastq/sampled/{sample}-{read}.trim.fastq.gz',
            sample=wc.sample, read=reads)
    else:
        return expand(
            'fastq/trimmed/{sample}-{read}.trim.fastq.gz',
            sample=wc.sample, read=reads)


rule hisat2:
    input:
        hisat2Input
    output:
        sam = pipe('aligned/{sample}.mapped.sam'),
        qc = 'qc/hisat2/{sample}.hisat2.txt'
    params:
        index = config['genome']['index'],
        strand = config['strand']
    log:
        'logs/hisat2/{sample}.log'
    conda:
        f'{ENVS}/hisat2.yaml'
    threads:
        #max(1, (config['threads'] / 2) - 1)
        max(1, config['threads'] - 1)
    shell:
        hisat2Cmd()


rule fixBAM:
    input:
        rules.hisat2.output.sam
    output:
        temp('aligned/{sample}.fixmate.bam')
    log:
        'logs/fixBAM/{sample}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        #'samtools fixmate -O bam,level=0 -m {input} {output} &> {log}'
        'samtools fixmate -O bam -m {input} {output} &> {log}'


rule sortBAM:
    input:
        rules.fixBAM.output
    output:
        temp('aligned/{sample}.sort.bam')
    log:
        'logs/sortBAM/{sample}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    threads:
        config['threads']
    shell:
        'samtools sort -@ {threads} {input} > {output} 2> {log}'


rule markdupBAM:
    input:
        rules.sortBAM.output
    output:
        bam = 'aligned/{sample}.markdup.bam',
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
        'aligned/{sample}.{stage}.bam'
    output:
        'aligned/{sample}.{stage}.bam.bai'
    log:
        'logs/indexBAM/{sample}-{stage}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    threads:
        max(1, (config['threads'] / 2) - 1,)
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
        bam = 'aligned/{sample}.markdup.bam',
        index = 'aligned/{sample}.markdup.bam.bai'
    output:
        'qc/samtools/idxstats/{sample}.idxstats.txt'
    group:
        'samQC'
    log:
        'logs/samtoolsIdxstats/{sample}.log'
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


rule multiBamSummary:
    input:
        bams = expand('aligned/{sample}.markdup.bam',
            sample=SAMPLES),
        indexes = expand('aligned/{sample}.markdup.bam.bai',
            sample=SAMPLES)
    output:
        'qc/deeptools/multiBamSummary.npz'
    params:
        binSize = 10000,
        distanceBetweenBins = 0,
        labels = ' '.join(SAMPLES),
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
        '--numberOfProcessors {threads} &> {log}'


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


def setColours(groups):
    """ Find group and assign to specific colour. """
    colours = ''
    colourPool = ['#E69F00', '#56B4E9', '#009E73', '#F0E442',
                  '#0072B2', '#D55E00', '#CC79A7']
    # Add double quotes for compatibility with shell
    colourPool = [f'"{colour}"' for colour in colourPool]
    usedColours = {}
    for group in groups:
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
        colours = setColours(RNASeq.groups())
    log:
        'logs/plotPCA.log'
    conda:
        f'{ENVS}/deeptools.yaml'
    shell:
        'plotPCA --corData {input} --colors {params.colours} '
        '--labels {params.labels} --plotFile {output.plot} '
        '--outFileNameData {output.data} &> {log}'


def getStrand():
    if config['strand'] in ['F', 'FR']:
        return 1
    elif config['strand'] in ['R', 'RF']:
        return 2
    else:
        return 0


rule gff3ToGTF:
    input:
        config['genome']['gff3']
    output:
        f'annotation/{config["genome"]["build"]}.gtf'
    log:
        'logs/gff3ToGTF.log'
    conda:
        f'{ENVS}/gffread.yaml'
    shell:
        'gffread <(zcat -f {input}) -T -o {output} &> {log}'


rule featureCounts:
    input:
        bam = rules.markdupBAM.output.bam,
        gtf = rules.gff3ToGTF.output
    output:
        counts = 'aligned/{sample}.featureCounts.txt',
        qc = 'qc/featurecounts/{sample}.featureCounts.txt.summary'
    params:
        paired = '-pC' if config['paired'] else '',
        strand = getStrand()
    log:
        'logs/featureCounts/{sample}.log'
    conda:
        f'{ENVS}/subread.yaml'
    threads:
        1
    shell:
        '(featureCounts {params.paired} -a {input.gtf} -o {output.counts} '
        '{input.bam} -t exon -g gene_id -T {threads} -s {params.strand} '
        '&& mv {output.counts}.summary {output.qc}) &> {log}'


rule gff3ToGenePred:
    input:
        config['genome']['gff3']
    output:
        pipe(f'annotation/{config["genome"]["build"]}.unsorted.genePred')
    log:
        'logs/gff3ToGenePred.log'
    conda:
        f'{ENVS}/UCSCtools.yaml'
    shell:
        'gff3ToGenePred -geneNameAttr=gene_name {input} {output} &> {log}'


rule sortGenePred:
    input:
        rules.gff3ToGenePred.output
    output:
        f'annotation/{config["genome"]["build"]}.genePred'
    log:
        'logs/sortGenePred.log'
    conda:
        f'{ENVS}/UCSCtools.yaml'
    shell:
        'sort -k2,2 -k4n,4n {input} > {output} 2> {log}'


rule genePredToBed:
    input:
        rules.sortGenePred.output
    output:
        f'annotation/{config["genome"]["build"]}.bed12'
    log:
        'logs/genePredToBed.log'
    conda:
        f'{ENVS}/UCSCtools.yaml'
    shell:
        'genePredToBed {input} {output} &> {log}'


def getBams():
    bams = expand('subsampled/{sample}.markdup.bam',
        sample=SAMPLES)
    return ','.join(bams)


rule readDistribution:
    input:
        bed12 = rules.genePredToBed.output,
        bam = 'aligned/{sample}.markdup.bam',
        index = 'aligned/{sample}.markdup.bam.bai'
    output:
        'qc/rseqc/readDistribution/{sample}.txt'
    log:
        'logs/readDistribution/{sample}.log'
    conda:
        f'{ENVS}/rseqc.yaml'
    shell:
        'read_distribution.py -r {input.bed12} -i {input.bam} '
        '> {output} 2> {log} '


rule geneBodyCoverage:
    input:
        bed12 = rules.genePredToBed.output,
        bams = expand('aligned/{sample}.markdup.bam',
            sample=SAMPLES),
        indexes = expand('aligned/{sample}.markdup.bam.bai',
            sample=SAMPLES)
    output:
        'qc/rseqc/geneBodyCoverage/{sample}.txt',
        'qc/rseqc/geneBodyCoverage/{sample}.pdf'
    params:
        bams = getBams(),
        prefix = lambda wc: f'qc/rseqc/geneBodyCoverage/{wc.sample}'
    log:
        'logs/geneBodyCoverage/{sample}.log'
    conda:
        f'{ENVS}/rseqc.yaml'
    shell:
        'geneBody_coverage.py -r {input.bed12} -i {params.bams} '
        '-o {params.prefix} &> {log} '


rule innerDistance:
    input:
        bed12 = rules.genePredToBed.output,
        bam = 'aligned/{sample}.markdup.bam',
        index = 'aligned/{sample}.markdup.bam.bai'
    output:
        'qc/rseqc/innerDistance/{sample}.inner_distance_freq.txt',
        'qc/rseqc/innerDistance/{sample}.inner_distance_plot.pdf',
        'qc/rseqc/innerDistance/{sample}.inner_distance_plot.r',
        'qc/rseqc/innerDistance/{sample}.inner_distance.txt'
    params:
        lowerBound = -250,
        upperBound = 250,
        sampleSize = 1000000,
        stepSize = 5,
        mapQual = 30,
        prefix = lambda wc: f'qc/rseqc/innerDistance/{wc.sample}'
    log:
        'logs/innerDistance/{sample}.log'
    conda:
        f'{ENVS}/rseqc.yaml'
    shell:
        'inner_distance.py -r {input.bed12} -i {input.bam} '
        '-l {params.lowerBound} -u {params.upperBound} '
        '-k {params.sampleSize} -s {params.stepSize} '
        '-q {params.mapQual} -o {params.prefix} &> {log} '


rule readGC:
    input:
        bam = 'aligned/{sample}.markdup.bam',
        index = 'aligned/{sample}.markdup.bam.bai'
    output:
        'qc/rseqc/readGC/{sample}.GC.xls',
        'qc/rseqc/readGC/{sample}.GC_plot.r',
        'qc/rseqc/readGC/{sample}.GC_plot.pdf'
    params:
        mapQual = 30,
        prefix = lambda wc: f'qc/rseqc/readGC/{wc.sample}'
    log:
        'logs/readGC/{sample}.log'
    conda:
        f'{ENVS}/rseqc.yaml'
    shell:
        'read_GC.py -i {input.bam} '
        '-q {params.mapQual} -o {params.prefix} &> {log} '


rule readDuplication:
    input:
        bam = 'aligned/{sample}.markdup.bam',
        index = 'aligned/{sample}.markdup.bam.bai'
    output:
        'qc/rseqc/readDuplication/{sample}.pos.DupRate.xls',
        'qc/rseqc/readDuplication/{sample}.seq.DupRate.xls',
        'qc/rseqc/readDuplication/{sample}.DupRate_plot.r',
        'qc/rseqc/readDuplication/{sample}.DupRate_plot.pdf'
    params:
        mapQual = 30,
        prefix = lambda wc: f'qc/rseqc/readDuplication/{wc.sample}'
    log:
        'logs/readDuplication/{sample}.log'
    conda:
        f'{ENVS}/rseqc.yaml'
    shell:
        'read_duplication.py -i {input.bam} '
        '-q {params.mapQual} -o {params.prefix} &> {log}'


rule junctionAnnotation:
    input:
        bed12 = rules.genePredToBed.output,
        bam = 'aligned/{sample}.markdup.bam',
        index = 'aligned/{sample}.markdup.bam.bai'
    output:
        expand('qc/rseqc/junctionAnnotation/{{sample}}.{ext}',
            ext=['junction.bed', 'junction_plot.r',
                 'splice_events.pdf', 'junction.Interact.bed',
                 'junction.xls', 'splice_junction.pdf']),
        txt = 'qc/rseqc/junctionAnnotation/{sample}.txt'
    params:
        minIntron = 50,
        mapQual = 30,
        prefix = lambda wc: f'qc/rseqc/junctionAnnotation/{wc.sample}'
    log:
        'logs/junctionAnnotation/{sample}.log'
    conda:
        f'{ENVS}/rseqc.yaml'
    shell:
        'junction_annotation.py -r {input.bed12} -i {input.bam} '
        '-q {params.mapQual} -m {params.minIntron} -o {params.prefix} '
        '&> {log} && cp {log} {output.txt}'


rule junctionSaturation:
    input:
        bed12 = rules.genePredToBed.output,
        bam = 'aligned/{sample}.markdup.bam',
        index = 'aligned/{sample}.markdup.bam.bai'
    output:
        expand('qc/rseqc/junctionSaturation/{{sample}}.{ext}',
            ext=['junctionSaturation_plot.r',
                 'junctionSaturation_plot.pdf'])
    params:
        lowerBound = 5,
        upperBound = 100,
        percentStep = 5,
        minIntron = 50,
        minSpliceRead = 1,
        mapQual = 30,
        prefix = lambda wc: f'qc/rseqc/junctionSaturation/{wc.sample}'
    log:
        'logs/junctionSaturation/{sample}.log'
    conda:
        f'{ENVS}/rseqc.yaml'
    shell:
        'junction_saturation.py -r {input.bed12} -i {input.bam} '
        '-l {params.lowerBound} -u {params.upperBound} '
        '-s {params.percentStep} -m {params.minIntron} '
        '-q {params.mapQual} -v {params.minSpliceRead} '
        '-o {params.prefix} &> {log}'


rule inferExperiment:
    input:
        bed12 = rules.genePredToBed.output,
        bam = 'aligned/{sample}.markdup.bam',
        index = 'aligned/{sample}.markdup.bam.bai'
    output:
        'qc/rseqc/inferExperiment/{sample}.txt'
    log:
        'logs/inferExperiment/{sample}.log'
    conda:
        f'{ENVS}/rseqc.yaml'
    shell:
        'infer_experiment.py -r {input.bed12} -i {input.bam} '
        '> {output} 2> {log} '


rule bamStat:
    input:
        bam = 'aligned/{sample}.markdup.bam',
        index = 'aligned/{sample}.markdup.bam.bai'
    output:
        'qc/rseqc/bamStat/{sample}.txt'
    params:
        mapQual = 30,
    log:
        'logs/readDuplication/{sample}.log'
    conda:
        f'{ENVS}/rseqc.yaml'
    shell:
        'bam_stat.py -i {input.bam} -q {params.mapQual} > {output} 2> {log}'


rule multiQC:
    input:
        expand('qc/fastqc/{single}.raw_fastqc.zip',
            single=RNASeq.singles()),
        expand('qc/fastq_screen/{single}.fastq_screen.txt',
            single=RNASeq.singles()) if config['fastqScreen'] else [],
        expand('qc/cutadapt/{sample}.cutadapt.txt', sample=SAMPLES),
        expand('qc/fastqc/{single}.trim_fastqc.zip',
            single=RNASeq.singles()),
        expand('qc/kallisto/{sample}.stdout', sample=SAMPLES),
        expand('qc/hisat2/{sample}.hisat2.txt', sample=SAMPLES),
        expand('qc/samtools/stats/{sample}.stats.txt', sample=SAMPLES),
        expand('qc/samtools/idxstats/{sample}.idxstats.txt',
            sample=SAMPLES),
        expand('qc/samtools/flagstat/{sample}.flagstat.txt',
            sample=SAMPLES),
        'qc/deeptools/plotCorrelation.tsv',
        'qc/deeptools/plotPCA.tab',
        expand('qc/rseqc/{tool}/{sample}.txt', sample=SAMPLES,
            tool=['inferExperiment', 'readDistribution',
                  'junctionAnnotation', 'bamStat']),#, 'geneBodyCoverage']),
        expand('qc/rseqc/junctionSaturation/{sample}.junctionSaturation_plot.r',
            sample=SAMPLES),
        expand('qc/rseqc/readGC/{sample}.GC.xls', sample=SAMPLES),
        expand('qc/rseqc/readDuplication/{sample}.pos.DupRate.xls',
            sample=SAMPLES),
        expand('qc/rseqc/innerDistance/{sample}.inner_distance_freq.txt',
            sample=SAMPLES) if config['paired'] else [],
        expand('qc/featurecounts/{sample}.featureCounts.txt.summary',
            sample=SAMPLES) if not config['QConly'] else [],
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
