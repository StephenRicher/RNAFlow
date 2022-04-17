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
    'ignoreChrom':       None        ,
    'genome':
        {'build':         'genome'    ,
         'transcriptome': ''          ,
         'gff3':          ''          ,
         'index':         ''          ,},
    'cutadapt':
        {'forwardAdapter': 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA',
         'reverseAdapter': 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT',
         'overlap':         3                                 ,
         'errorRate':       0.1                               ,
         'minimumLength':   1                                 ,
         'qualityCutoff':  '0,0'                              ,
         'GCcontent':       50                                ,},
    'misc':
        {'inferExperimentSampleSize': 200000                 ,
         'sample'                   : 0                      ,
         'geneBodyCoverage'         : True                   ,
         'readDuplication'          : True                   ,},
    'fastqScreen':         None                              ,
    'multiQCconfig':       None                              ,
}

config = set_config(config, default_config)

workdir: config['workdir']
BUILD = config['genome']['build']

RNASeq = Samples(config['data'], config['paired'], config['strand'])

# Extract sample names
SAMPLES = RNASeq.samples()

strandedness = RNASeq.getStrandedness()

wildcard_constraints:
    sample = '|'.join(RNASeq.samples())

# Create tmpdir to ensure it is set
os.makedirs(config['tmpdir'], exist_ok=True)

# Minimum length must be atleast 1 (if 0 this break seqkt)
if config['cutadapt']['minimumLength'] < 1:
    config['cutadapt']['minimumLength'] = 1

for strand in strandedness.values():
    if config['paired']:
        assert strand in ['RF', 'FR', 'unstranded']
        reads = ['R1', 'R2']
    else:
        assert strand in ['R', 'F', 'unstranded']
        reads = ['R1']

rule all:
    input:
        'qc/multiQC',
        expand('bigWig/{sample}.bigWig', sample=SAMPLES)


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
        return [temp('fastq/trimmed/{sample}-R1.trim.fastq.gz'),
                temp('fastq/trimmed/{sample}-R2.trim.fastq.gz')]
    else:
        return [temp('fastq/trimmed/{sample}-R1.trim.fastq.gz')]


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
            txt = 'qc/fastqScreen/{sample}-{read}.fastq_screen.txt',
            png = 'qc/fastqScreen/{sample}-{read}.fastq_screen.png'
        params:
            config = config['fastqScreen'],
            subset = 100000,
        group:
            'processFASTQ'
        log:
            'logs/fastqScreen/{sample}-{read}.log'
        conda:
            f'{ENVS}/fastqScreen.yaml'
        threads:
            config['threads']
        shell:
            '{SCRIPTS}/fastqScreen.py {input} {params.config} '
            '--subset {params.subset} --threads {threads} '
            '--plotOut {output.png} --dataOut {output.txt} &> {log}'


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
    cmd = ('kallisto quant -i {input.index} -o kallisto/{wildcards.sample}/ '
           '-b {params.bootstraps} --seed {params.seed} {params.strand} '
           '--threads {threads} ')
    if not config['paired']:
        cmd += '--single -l {params.fragmentLength} -s {params.fragmentSD} '
    cmd += '{input.reads} &> {log} && cp {log} {output.qc}'
    return cmd


def getKallistoStrand(wc):
    strand = strandedness[wc.sample]
    if strand in ['F', 'FR']:
        return '--fr-stranded'
    elif strand in ['R', 'RF']:
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
        seed = 42,
        fragmentSD = 2,
        bootstraps = 100,
        fragmentLength = 200,
        strand = getKallistoStrand
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
        nReads = config['misc']['sample']
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
    command += ('{params.strand} -x {params.index} --threads {threads} '
                '--summary-file {output.qc} > {output.sam} 2> {log}')
    return command


def hisat2Input(wc):
    if config['misc']['sample'] > 0:
        return expand(
            'fastq/sampled/{sample}-{read}.trim.fastq.gz',
            sample=wc.sample, read=reads)
    else:
        return expand(
            'fastq/trimmed/{sample}-{read}.trim.fastq.gz',
            sample=wc.sample, read=reads)


def getHisat2Strand(wc):
    strand = strandedness[wc.sample]
    if strand != 'unstranded':
        return f'--rna-strandness {strand} '
    else:
        return ''


rule hisat2:
    input:
        hisat2Input
    output:
        sam = pipe('aligned/{sample}.mapped.sam'),
        qc = 'qc/hisat2/{sample}.hisat2.txt'
    params:
        index = config['genome']['index'],
        strand = getHisat2Strand
    log:
        'logs/hisat2/{sample}.log'
    conda:
        f'{ENVS}/hisat2.yaml'
    threads:
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
        'aligned/{sample}.sort.bam'
    output:
        bam = 'aligned/{sample}.markdup.bam',
        qc = 'qc/markDuplicates/{sample}.markDuplicates.log'
    params:
        tmp = config['tmpdir']
    log:
        'logs/markDuplicates/{sample}.log'
    conda:
        f'{ENVS}/picard.yaml'
    shell:
        'picard MarkDuplicates I={input} O={output.bam} M={output.qc} '
        'VALIDATION_STRINGENCY=LENIENT TMP_DIR={params.tmp} &> {log}'


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
        config['threads']
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
        'logs/samtoolsStats/{sample}.log'
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

def getStrand(wc):
    strand = strandedness[wc.sample]
    if strand in ['F', 'FR']:
        return 1
    elif strand in ['R', 'RF']:
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
        counts = 'featureCounts/{sample}.featureCounts.txt',
        qc = 'featureCounts/{sample}.featureCounts.txt.summary'
    params:
        paired = '-pC' if config['paired'] else '',
        strand = getStrand
    log:
        'logs/featureCounts/{sample}.log'
    conda:
        f'{ENVS}/subread.yaml'
    threads:
        config['threads']
    shell:
        'featureCounts {params.paired} -a {input.gtf} -o {output.counts} '
        '{input.bam} -t exon -g gene_id -T {threads} -s {params.strand} '
        '&> {log}'


rule generateScaleFactor:
    input:
        expand('featureCounts/{sample}.featureCounts.txt', sample=SAMPLES)
    output:
        'featureCounts/allScaleFactors.tsv'
    log:
        'logs/generateScaleFactor.log'
    shell:
        'Rscript {SCRIPTS}/getScaleFactors.R {input} > {output} 2> {log} '
        

rule multiBamSummary:
    input:
        bams = expand('aligned/{sample}.markdup.bam',
            sample=SAMPLES),
        indexes = expand('aligned/{sample}.markdup.bam.bai',
            sample=SAMPLES)
    output:
        'qc/deepTools/multiBamSummary.npz'
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

def setStrandFilter(wc):
    strand = getStrand(wc)
    if strand == 0:
        return ''
    elif strand == 1:
        return '--filterRNAstrand reverse'
    else:
        return '--filterRNAstrand forward'

def ignoreChrom(wc):
    if config['ignoreChrom'] is None:
        return ''
    else:
        chroms = []
        with open(config['ignoreChrom']) as fh:
            for chrom in fh:
                chroms.append(chrom.strip())
        if len(chroms) == 0:
            return ''
        else:
            return f'--ignoreForNormalization {" ".join(chroms)}'

def getScaleFactor(wc):
    with open('featureCounts/allScaleFactors.tsv') as fh:
        for line in fh:
            sample, factor = line.strip().split()
            if sample == wc.sample:
                return factor

rule bamCoverage:
    input:
        bam = 'aligned/{sample}.markdup.bam',
        index = 'aligned/{sample}.markdup.bam.bai',
        scaleFactors = rules.generateScaleFactor.output
    output:
        'bigWig/{sample}.bigWig'
    params:
        binSize = 50,
        scaleFactor = getScaleFactor,
        filterRNAstrand = setStrandFilter,
        ignoreChromForNormalisation = ignoreChrom
    log:
        'logs/bamCoverage/{sample}.log'
    conda:
        f'{ENVS}/deeptools.yaml'
    threads:
        config['threads']
    shell:
        'bamCoverage --bam {input.bam} --outFileName {output} '
        '{params.filterRNAstrand} '
        '--binSize {params.binSize} --scaleFactor {params.scaleFactor} '
        '--skipNonCoveredRegions --numberOfProcessors {threads} '
        '--verbose &> {log}'


rule plotCorrelation:
    input:
        rules.multiBamSummary.output
    output:
        plot = 'qc/deepTools/plotCorrelation.png',
        matrix = 'qc/deepTools/plotCorrelation.tsv'
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
    colourPool = ['#a6cee3', '#1f78b4', '#b2df8a', '#33a02c',
                  '#fb9a99', '#e31a1c', '#fdbf6f', '#ff7f00',
                  '#cab2d6', '#6a3d9a', '#ffff99', '#b15928',
                  '#000000', '#ff0000']
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
        plot = 'qc/deepTools/plotPCA.png',
        data = 'qc/deepTools/plotPCA.tab'
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


rule readDistribution:
    input:
        bed12 = rules.genePredToBed.output,
        bam = 'aligned/{sample}.markdup.bam',
        index = 'aligned/{sample}.markdup.bam.bai'
    output:
        'qc/RSeQC/readDistribution/{sample}.txt'
    log:
        'logs/readDistribution/{sample}.log'
    conda:
        f'{ENVS}/rseqc.yaml'
    shell:
        'read_distribution.py -r {input.bed12} -i {input.bam} '
        '> {output} 2> {log} '


def getBams():
    bams = expand('aligned/{sample}.markdup.bam', sample=SAMPLES)
    return ','.join(bams)


rule geneBodyCoverage:
    input:
        bed12 = rules.genePredToBed.output,
        bams = expand('aligned/{sample}.markdup.bam',
            sample=SAMPLES),
        indexes = expand('aligned/{sample}.markdup.bam.bai',
            sample=SAMPLES)
    output:
        'qc/RSeQC/geneBodyCoverage/data.geneBodyCoverage.txt',
        'qc/RSeQC/geneBodyCoverage/data.geneBodyCoverage.r',
        'qc/RSeQC/geneBodyCoverage/data.geneBodyCoverage.curves.pdf'
    params:
        bams = getBams(),
        prefix = lambda wc: f'qc/RSeQC/geneBodyCoverage/data'
    log:
        'logs/geneBodyCoverage.log'
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
        'qc/RSeQC/innerDistance/{sample}.inner_distance_freq.txt',
        'qc/RSeQC/innerDistance/{sample}.inner_distance_plot.pdf',
        'qc/RSeQC/innerDistance/{sample}.inner_distance_plot.r',
        'qc/RSeQC/innerDistance/{sample}.inner_distance.txt'
    params:
        lowerBound = -250,
        upperBound = 250,
        sampleSize = 1000000,
        stepSize = 5,
        mapQual = 30,
        prefix = lambda wc: f'qc/RSeQC/innerDistance/{wc.sample}'
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
        'qc/RSeQC/readGC/{sample}.GC.xls',
        'qc/RSeQC/readGC/{sample}.GC_plot.r',
        'qc/RSeQC/readGC/{sample}.GC_plot.pdf'
    params:
        mapQual = 30,
        prefix = lambda wc: f'qc/RSeQC/readGC/{wc.sample}'
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
        'qc/RSeQC/readDuplication/{sample}.pos.DupRate.xls',
        'qc/RSeQC/readDuplication/{sample}.seq.DupRate.xls',
        'qc/RSeQC/readDuplication/{sample}.DupRate_plot.r',
        'qc/RSeQC/readDuplication/{sample}.DupRate_plot.pdf'
    params:
        mapQual = 30,
        prefix = lambda wc: f'qc/RSeQC/readDuplication/{wc.sample}'
    log:
        'logs/readDuplication/{sample}.log'
    conda:
        f'{ENVS}/rseqc.yaml'
    threads:
        config['threads']
    shell:
        'read_duplication.py -i {input.bam} '
        '-q {params.mapQual} -o {params.prefix} &> {log}'


rule junctionAnnotation:
    input:
        bed12 = rules.genePredToBed.output,
        bam = 'aligned/{sample}.markdup.bam',
        index = 'aligned/{sample}.markdup.bam.bai'
    output:
        expand('qc/RSeQC/junctionAnnotation/{{sample}}.{ext}',
            ext=['junction.bed', 'junction_plot.r',
                 'splice_events.pdf', 'junction.Interact.bed',
                 'junction.xls', 'splice_junction.pdf']),
        txt = 'qc/RSeQC/junctionAnnotation/{sample}.txt'
    params:
        minIntron = 50,
        mapQual = 30,
        prefix = lambda wc: f'qc/RSeQC/junctionAnnotation/{wc.sample}'
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
        expand('qc/RSeQC/junctionSaturation/{{sample}}.{ext}',
            ext=['junctionSaturation_plot.r',
                 'junctionSaturation_plot.pdf'])
    params:
        lowerBound = 5,
        upperBound = 100,
        percentStep = 5,
        minIntron = 50,
        minSpliceRead = 1,
        mapQual = 30,
        prefix = lambda wc: f'qc/RSeQC/junctionSaturation/{wc.sample}'
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
        'qc/RSeQC/inferExperiment/{sample}.txt'
    params:
        sampleSize = config['misc']['inferExperimentSampleSize']
    log:
        'logs/inferExperiment/{sample}.log'
    conda:
        f'{ENVS}/rseqc.yaml'
    shell:
        'infer_experiment.py -r {input.bed12} -i {input.bam} -s 75000 '
        '> {output} 2> {log} '


rule bamStat:
    input:
        bam = 'aligned/{sample}.markdup.bam',
        index = 'aligned/{sample}.markdup.bam.bai'
    output:
        'qc/RSeQC/bamStat/{sample}.txt'
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
        expand('qc/fastqScreen/{single}.fastq_screen.txt',
            single=RNASeq.singles()) if config['fastqScreen'] else [],
        expand('qc/cutadapt/{sample}.cutadapt.txt', sample=SAMPLES),
        expand('qc/fastqc/{single}.trim_fastqc.zip',
            single=RNASeq.singles()),
        expand('qc/kallisto/{sample}.stdout', sample=SAMPLES),
        expand('qc/hisat2/{sample}.hisat2.txt', sample=SAMPLES),
        expand('qc/samtools/{tool}/{sample}.{tool}.txt',
            sample=SAMPLES, tool=['stats', 'idxstats', 'flagstat']),
        expand('qc/markDuplicates/{sample}.markDuplicates.log', sample=SAMPLES),
        #'qc/deepTools/plotCorrelation.tsv',
        #'qc/deepTools/plotPCA.tab',
        expand('qc/RSeQC/{tool}/{sample}.txt', sample=SAMPLES,
            tool=['readDistribution', 'junctionAnnotation',
                  'bamStat', 'inferExperiment']),
        ('qc/RSeQC/geneBodyCoverage/data.geneBodyCoverage.txt'
            if config['misc']['geneBodyCoverage'] else []),
        expand('qc/RSeQC/junctionSaturation/{sample}.junctionSaturation_plot.r',
            sample=SAMPLES),
        expand('qc/RSeQC/readGC/{sample}.GC.xls', sample=SAMPLES),
        expand('qc/RSeQC/readDuplication/{sample}.pos.DupRate.xls',
            sample=SAMPLES) if config['misc']['readDuplication'] else [],
        expand('qc/RSeQC/innerDistance/{sample}.inner_distance_freq.txt',
            sample=SAMPLES) if config['paired'] else [],
        expand('featureCounts/{sample}.featureCounts.txt.summary',
            sample=SAMPLES),
    output:
        directory('qc/multiQC')
    params:
        config = f'--config {config["multiQCconfig"]}' if config["multiQCconfig"] else ''
    log:
        'logs/multiQC/multiQC.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        'multiqc --outdir {output} --force {params.config} {input} &> {log}'
