import pandas as pd

wildcard_constraints:
    group = '[^\/-]+',
    sample = '[^\-]+-\d+',
    single = '[^\-]+-\d+-R[12]',
    rep = '\d+',
    read = 'R[12]'

BASE = workflow.basedir

ENVS = f'{BASE}/workflow/envs'
SCRIPTS = f'{BASE}/workflow/scripts'

configfile : f'{BASE}/config/config.yaml'
GENOME = config['genome']
BUILD = config['build']
ANNOTATION = config['annotation']

try:
    OVERHANG = config['read_length'] - 1
except KeyError:
    OVERHANG = 100 - 1
try:
    workdir : config['outdir']
except KeyError:
    pass

DATA_TABLE = config['data']
DATA = pd.read_table(config['data'], sep = ',', dtype = {'rep' : str})

# Validate read file input with wildcard definitions
if not DATA['group'].str.match(r'[^-\/]+').all():
    sys.exit(f'Invalid group definition in {DATA_TABLE}.')
if not DATA['rep'].str.match(r'\d+').all():
    sys.exit(f'Invalid replicate definition in {DATA_TABLE}.')
if not DATA['read'].str.match(r'R[12]').all():
    sys.exit(f'Invalid read definition in {DATA_TABLE}.')

# Define single and sample names from definitions.
DATA['single'] = (DATA[['group', 'rep', 'read']]
    .apply(lambda x: '-'.join(x), axis = 1))
DATA['sample'] = (DATA[['group', 'rep']]
    .apply(lambda x: '-'.join(x), axis = 1))
DATA = DATA.set_index(['group', 'sample', 'single'], drop = False)

# Extract groups and replicates.
GROUPS = {}
for group in DATA['group']:
    GROUPS[group] = list(DATA.loc[group]['rep'].unique())
# Extract sample names
SAMPLES = list(DATA['sample'].unique())
# Define READS (R1 and R2)
READS = ['R1', 'R2']

THREADS = 1

rule all:
    input:
        ['multiqc_report.html', f'genome/{BUILD}.fa.gz.fai']

rule bgzip_genome:
    input:
        GENOME
    output:
        f'genome/{BUILD}.fa.gz'
    log:
        f'logs/bgzip_genome/{BUILD}.log'
    conda:
        f'{ENVS}/tabix.yaml'
    shell:
        '(zcat -f {input} | bgzip --stdout > {output}) 2> {log}'

rule faidx:
    input:
        rules.bgzip_genome.output
    output:
         multiext(f'{rules.bgzip_genome.output}', '.fai', '.gzi')
    log:
        'logs/index_genome/index_genome.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools faidx {input} &> {log}'

rule STAR_index:
    input:
        genome = rules.bgzip_genome.output,
    output:
        directory('genome/idx')
    log:
        'logs/build_STAR_index/STAR_index.log'
    conda:
        f'{ENVS}/STAR.yaml'
    threads:
        THREADS
    shell:
        'STAR '
            '--runMode genomeGenerate --runThreadN {threads} '
            '--genomeDir {output} --genomeFastaFiles {input.genome} '
            '--sjdbOverhang {OVERHANG} --sjdbGTFfile {ANNOTATION} '
        '&> {log}'

rule fastqc:
    input:
        lambda wc: DATA.xs(wc.single, level = 2)['path']
    output:
        html = 'qc/{single}-fastqc.html',
        zip = 'qc/{single}-fastqc.zip'
    log:
        'logs/fastqc/{single}.log'
    wrapper:
        '0.49.0/bio/fastqc'

#rule trimmomatic:
#    input:
#        lambda wc: DATA.xs(wc.sample, level = 1)['path']
#    output:
#        r1_paired = 'fastq/{sample}-R1-trim.fq.gz',
#        r1_unpaired = 'fastq/{sample}-R1-trim-unpaired.fq.gz',
#        r2_paired = 'fastq/{sample}-R2-trim.fq.gz',
#        r2_unpaired = 'fastq/{sample}-R2-trim-unpaired.fq.gz',
#        qc = 'qc/{sample}-trimmomatic.txt'
#    log:
#        'logs/trimmomatic/{sample}.log'
#    conda:
#        f'{ENVS}/trimmomatic.yaml'
#    threads:
#        THREADS
#    shell:
#        'trimmomatic PE '
#            '{input} -threads {threads} -trimlog {output.qc} '
#            '{output.r1_paired} {output.r1_unpaired} '
#            '{output.r2_paired} {output.r2_unpaired} '
#            'ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 '
#            'LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20 '
#        '&> {log}'

rule cutadapt:
    input:
        lambda wc: DATA.xs(wc.sample, level = 1)['path']
    output:
        r1 = 'fastq/trimmed/{sample}-R1-trim.fq.gz',
        r2 = 'fastq/trimmed/{sample}-R2-trim.fq.gz',
        qc = 'qc/{sample}-cutadapt.txt'
    params:
        adapters = '-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA '
                   '-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT',
        others = '--minimum-length 20 --quality-cutoff 20 '
                 '--gc-content 46 --overlap 6 --error-rate 0.1'
    log:
        'logs/cutadapt/{sample}.log'
    conda:
        f'{ENVS}/cutadapt.yaml'
    threads:
        THREADS
    shell:
        'cutadapt '
            '{params.adapters} {params.others} --cores {THREADS} '
            '-o {output.r1} -p {output.r2} '
            '{input} > {output.qc} '
        '2> {log}'

rule fastqc_trimmed:
    input:
        'fastq/trimmed/{sample}-{read}-trim.fq.gz'
    output:
        html = 'qc/{sample}-{read}-trim-fastqc.html',
        zip = 'qc/{sample}-{read}-trim-fastqc.zip'
    log:
        'logs/fastqc_trimmed/{sample}-{read}.log'
    wrapper:
        '0.49.0/bio/fastqc'

rule STAR_map:
    input:
        fastq = [rules.cutadapt.output.r1, rules.cutadapt.output.r2],
        idx_dir = rules.STAR_index.output
    output:
        bam = 'mapped/{sample}.out.bam',
        qc = 'qc/{sample}-STAR.txt'
    log:
        'logs/STAR_map/{sample}.log'
    conda:
        f'{ENVS}/STAR.yaml'
    threads:
        THREADS
    shell:
        'STAR '
            '--runMode alignReads --runThreadN {threads} '
            '--genomeDir {input.idx_dir} '
            '--readFilesIn {input.fastq} --readFilesCommand zcat '
            '--outFileNamePrefix mapped/{wildcards.sample} '
            '--outSAMtype BAM SortedByCoordinate '
        '&> {log}; cp {log} {output.qc}'

rule hisat2_map:
    input:
        fastq_r1 = rules.cutadapt.output.r1,
        fastq_r2 = rules.cutadapt.output.r2
    output:
        sam = 'mapped/{sample}.sam',
        qc = 'qc/{sample}-hisat2.txt'
    params:
        index = '/media/stephen/Data/genomes/GRCm38/index/GRCm38'
    log:
        'logs/hisat2_map/{sample}.log'
    conda:
        f'{ENVS}/hisat2.yaml'
    threads:
        12
    shell:
        'hisat2 '
            '-x {params.index} -p 4 '
            '-1 {input.fastq_r1} -2 {input.fastq_r2} '
            '--summary-file {output.qc} '
        '> {output.sam} '
        '2> {log}'

rule samtools_sort:
    input:
        rules.hisat2_map.output.sam
    output:
        'mapped/{sample}.bam',
    log:
        'logs/samtools_sort/{sample}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    threads:
        12
    shell:
        'samtools sort -@ 2 {input} > {output} 2> {log}'


rule index_bam:
    input:
        rules.samtools_sort.output
        #rules.STAR_map.output.bam
    output:
        f'{rules.samtools_sort.output}.bai'
        #f'{rules.STAR_map.output.bam}.bai'
    log:
        'logs/index_bam/{sample}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    threads:
        12
    shell:
        'samtools index -@ 2 {input} &> {log}'

rule samtools_stats:
    input:
        rules.samtools_sort.output
        #rules.STAR_map.output.bam
    output:
        'qc/{sample}-stats.txt'
    group:
        'SAM_QC'
    log:
        'logs/samtools_stats/{sample}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools stats {input} > {output} 2> {log}'

rule samtools_idxstats:
    input:
        rules.samtools_sort.output,
        #rules.STAR_map.output.bam,
        rules.index_bam.output
    output:
        'qc/{sample}-idxstats.txt'
    group:
        'SAM_QC'
    log:
        'logs/samtools_idxstats/{sample}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools stats {input} > {output} 2> {log}'

rule samtools_flagstat:
    input:
        rules.samtools_sort.output
        #rules.STAR_map.output.bam
    output:
        'qc/{sample}-flagstat.txt'
    group:
        'SAM_QC'
    log:
        'logs/samtools_flagstat/{sample}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools flagstat {input} > {output} 2> {log}'

rule bamqc:
    input:
        rules.samtools_sort.output
        #rules.STAR_map.output.bam
    output:
        file = 'qc/{sample}-bamqc.html',
        dir = directory('qc/bamqc/{sample}/')
    log:
        'logs/bamqc/{sample}.log'
    conda:
        f'{ENVS}/qualimap.yaml'
    threads:
        12
    shell:
        'qualimap bamqc '
            '-bam {input} -gff {ANNOTATION} '
            '-outdir {output.dir} -outfile {output.file} '
            '-nt 2 --outside-stats '
        '&> {log}'

rule rnaqc:
    input:
        rules.samtools_sort.output
        #rules.STAR_map.output.bam
    output:
        file = 'qc/{sample}-rnaqc.html',
        dir = directory('qc/rnaqc/{sample}')
    log:
        'logs/rnaqc/{sample}.log'
    conda:
        f'{ENVS}/qualimap.yaml'
    threads:
        12
    shell:
        'qualimap rnaqc '
            '-bam {input} -gff {ANNOTATION} '
            '-outdir {output.dir} -outfile {output.file} '
            '-nt 2 --outside-stats '
        '&> {log}'

rule multiqc:
    input:
        expand('qc/{sample}-{read}-fastqc.html',
            sample = SAMPLES, read = READS),
        expand('qc/{sample}-cutadapt.txt', sample = SAMPLES),
        expand('qc/{sample}-{read}-trim-fastqc.html',
            sample = SAMPLES, read = READS),
        #expand('qc/{sample}-STAR.txt', sample = SAMPLES),
        expand('qc/{sample}-stats.txt', sample = SAMPLES),
        expand('qc/{sample}-idxstats.txt', sample = SAMPLES),
        expand('qc/{sample}-flagstat.txt', sample = SAMPLES),
        expand('qc/{sample}-bamqc.html', sample = SAMPLES),
        expand('qc/{sample}-rnaqc.html', sample = SAMPLES)
    output:
        'multiqc_report.html'
    log:
        'logs/multiqc/multiqc.log'
    conda:
        f'{ENVS}/multiqc.yaml'
    shell:
        'multiqc qc/ &> {log}'
