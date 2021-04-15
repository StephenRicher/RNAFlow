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
    'QConly':            False       ,
    'genome':
        {'build':         'genome'    ,
         'transcriptome': ''          ,
         'annotation':    ''          ,
         'sequence':      ''          ,
         'genes':         ''          ,},
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

RNASeq = Samples(config['data'], config['paired'])

# Extract sample names
SAMPLES = RNASeq.samples()

wildcard_constraints:
    sample = '|'.join(RNASeq.samples())

# Create tmpdir to ensure it is set
os.makedirs(config['tmpdir'], exist_ok=True)

rule all:
    input:
        'qc/multiqc',
        #'counts/countMatrix-merged.txt' if not config['QConly'] else []


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
        lambda wc: RNASeq.path(wc.sample, ['R1', 'R2'])
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
            'fastq/trimmed/{sample}-{read}.trim.fastq.gz'
        output:
            txt = 'qc/fastq_screen/{sample}-{read}.fastq_screen.txt',
            png = 'qc/fastq_screen/{sample}-{read}.fastq_screen.png'
        params:
            fastq_screen_config = config['fastq_screen'],
            subset = 100000,
            aligner = 'bowtie2'
        group:
            'processFASTQ'
        log:
            'logs/fastq_screen/{sample}-{read}.log'
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
                '--genomebam --gtf {input.gtf} --chromosomes {input.chromSizes} '
                '-b {params.bootstraps} {input.reads} &> {log} '
                '&& cp {log} {output.qc}')
    else:
        return ('kallisto quant -i {input.index} -o kallisto/{wildcards.sample}/ '
                '--genomebam --gtf {input.gtf} --chromosomes {input.chromSizes} '
                '-b {params.bootstraps} '
                '--single -l {params.fragmentLength} -s {params.fragmentSD} '
                '--threads {threads} {input.reads} &> {log} '
                '&& cp {log} {output.qc}')


rule kallistoQuant:
    input:
        index = rules.buildKallistoIndex.output,
        reads = rules.cutadapt.output.trimmed,
        gtf = config['genome']['annotation'],
        chromSizes = rules.getChromSizes.output
    output:
        qc = 'qc/kallisto/{sample}.stdout',
        abundancesH5 = 'kallisto/{sample}/abundance.h5',
        abundancesTSV = 'kallisto/{sample}/abundance.tsv',
        runInfo = 'kallisto/{sample}/run_info.json',
        pseudoAlign = 'kallisto/{sample}/pseudoalignments.bam',
        pseudoAlignIdx = 'kallisto/{sample}/pseudoalignments.bam.bai'
    params:
        bootstraps = 30,
        fragmentLength = 200,
        fragmentSD = 2
    log:
        'logs/kallistoQuant/{sample}.log'
    conda:
        f'{ENVS}/kallisto.yaml'
    threads:
        config['threads']
    shell:
        kallistoCmd()


rule samtoolsStats:
    input:
        rules.kallistoQuant.output.pseudoAlign
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
        bam = rules.kallistoQuant.output.pseudoAlign,
        index = rules.kallistoQuant.output.pseudoAlignIdx
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
        rules.kallistoQuant.output.pseudoAlign
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
        bams = expand('kallisto/{sample}/pseudoalignments.bam',
            sample=SAMPLES),
        indexes = expand('kallisto/{sample}/pseudoalignments.bam.bai',
            sample=SAMPLES)
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


def setBlacklistCommand():
    if config['filtering']['blackListFileName']:
        cmd = ('bedtools merge -d {params.distance} -i {input} '
               '> {output} 2> {log}')
    else:
        cmd = 'touch {output} &> {log}'
    return cmd

if False:
    rule processBlacklist:
        input:
            config['filtering']['blackListFileName'] if config['filtering']['blackListFileName'] else []
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


rule multiQC:
    input:
        expand('qc/fastqc/{single}.raw_fastqc.zip',
            single=RNASeq.singles()),
        expand('qc/fastq_screen/{single}.fastq_screen.txt',
            single=RNASeq.singles()) if config['fastq_screen'] else [],
        expand('qc/cutadapt/{sample}.cutadapt.txt',
            sample=SAMPLES),
        expand('qc/fastqc/{single}.trim_fastqc.zip',
            single=RNASeq.singles()),
        expand('qc/kallisto/{sample}.stdout', sample=SAMPLES),
        expand('qc/samtools/stats/{sample}.stats.txt', sample=SAMPLES),
        expand('qc/samtools/idxstats/{sample}.idxstats.txt',
            sample=SAMPLES),
        expand('qc/samtools/flagstat/{sample}.flagstat.txt',
            sample=SAMPLES),
        'qc/deeptools/plotCorrelation.tsv',
        'qc/deeptools/plotPCA.tab'
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
