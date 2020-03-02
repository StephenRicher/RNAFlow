import pandas as pd

wildcard_constraints:
    group = '[^\/\s.-]+',
    sample = '[^\/\s.-]+-\d+',
    single = '[^\/\s.-]+-\d+-R[12]',
    rep = '\d+',
    read = 'R[12]'

BASE = workflow.basedir

ENVS = f'{BASE}/workflow/envs'
SCRIPTS = f'{BASE}/workflow/scripts'

configfile : f'{BASE}/config/config.yaml'
GENOME = config['genome']
BUILD = config['build']
FASTQ_SCREEN_CONFIG = config['fastq_screen_config']
ANNOTATION = config['annotation'] # Must NOT be gzipped

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
if not DATA['group'].str.match(r'[^\/\s.-]+').all():
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
        ['multiqc_report.html', 'qc/multibamqc', f'genome/{BUILD}.fa.gz.fai']

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
        html = 'qc/{single}.fastqc.html',
        zip = 'qc/{single}.fastqc.zip'
    log:
        'logs/fastqc/{single}.log'
    wrapper:
        '0.49.0/bio/fastqc'

rule cutadapt:
    input:
        lambda wc: DATA.xs(wc.sample, level = 1)['path']
    output:
        r1 = 'fastq/trimmed/{sample}-R1.trim.fastq.gz',
        r2 = 'fastq/trimmed/{sample}-R2.trim.fastq.gz',
        qc = 'qc/{sample}.cutadapt.txt'
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

rule fastq_screen:
    input:
        'fastq/trimmed/{sample}-{read}.trim.fastq.gz'
    output:
        txt = 'qc/{sample}-{read}.fastq_screen.txt',
        png = 'qc/{sample}-{read}.fastq_screen.png'
    params:
        fastq_screen_config = FASTQ_SCREEN_CONFIG,
        subset = 100000,
        aligner = 'bowtie2'
    log:
        'logs/fastq_screen/{sample}-{read}.log'
    threads:
        8
    wrapper:
        "0.49.0/bio/fastq_screen"

rule fastqc_trimmed:
    input:
        'fastq/trimmed/{sample}-{read}.trim.fastq.gz'
    output:
        html = 'qc/{sample}-{read}.trim.fastqc.html',
        zip = 'qc/{sample}-{read}.trim.fastqc.zip'
    log:
        'logs/fastqc_trimmed/{sample}-{read}.log'
    wrapper:
        '0.49.0/bio/fastqc'

rule STAR_map:
    input:
        fastq = [rules.cutadapt.output.r1, rules.cutadapt.output.r2],
        idx_dir = rules.STAR_index.output
    output:
        bam = 'mapped/{sample}.sorted.bam',
        qc = 'qc/{sample}.STAR.txt'
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
        bam = 'mapped/{sample}.bam',
        qc = 'qc/{sample}.hisat2.txt'
    params:
        index = '/media/stephen/Data/genomes/GRCm38/index/GRCm38'
    log:
        'logs/hisat2_map/{sample}.log'
    conda:
        f'{ENVS}/hisat2.yaml'
    threads:
        12
    shell:
        '(hisat2 '
            '-x {params.index} -p 6 '
            '-1 {input.fastq_r1} -2 {input.fastq_r2} '
            '--summary-file {output.qc} '
        '| samtools view -b > {output.bam}) '
        '2> {log}'

rule samtools_sort:
    input:
        rules.hisat2_map.output.bam
    output:
        'mapped/{sample}.sort.bam',
    log:
        'logs/samtools_sort/{sample}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    threads:
        12
    shell:
        'samtools sort -@ 6 {input} > {output} 2> {log}'

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
        'samtools index -@ 6 {input} &> {log}'

rule samtools_stats:
    input:
        rules.samtools_sort.output
        #rules.STAR_map.output.bam
    output:
        'qc/{sample}.stats.txt'
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
        'qc/{sample}.idxstats.txt'
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
        'qc/{sample}.flagstat.txt'
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
        directory('qc/bamqc/{sample}')
    resources:
        mem_mb = 3000
    log:
        'logs/bamqc/{sample}.log'
    conda:
        f'{ENVS}/qualimap.yaml'
    threads:
        12
    shell:
        'qualimap bamqc '
            '-bam {input} -gff {ANNOTATION} -outdir {output} '
            '-nt 6 --outside-stats --java-mem-size={resources.mem_mb}M '
            '--paint-chromosome-limits '
        '&> {log}'

rule rnaqc:
    input:
        rules.samtools_sort.output
        #rules.STAR_map.output.bam
    output:
        directory('qc/rnaqc/{sample}')
    resources:
        mem_mb = 3000
    log:
        'logs/rnaqc/{sample}.log'
    conda:
        f'{ENVS}/qualimap.yaml'
    shell:
        'qualimap rnaseq '
            '-bam {input} -gtf {ANNOTATION} -outdir {output} '
            '--java-mem-size={resources.mem_mb}M '
        '&> {log}'

rule multibamqc_config:
    input:
        expand('qc/bamqc/{sample}', sample = SAMPLES)
    output:
        'qc/bamqc/multibamqc.config'
    group:
        'multiBAM_QC'
    log:
        'logs/multibamqc_config/multibamqc_config.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/multibamqc_config.py {input} > {output} 2> {log}'

rule multibamqc:
    input:
        rules.multibamqc_config.output
    output:
        directory('qc/multibamqc')
    group:
        'multiBAM_QC'
    log:
        'logs/multibamqc/multibamqc.log'
    conda:
        f'{ENVS}/qualimap.yaml'
    shell:
        'qualimap multi-bamqc '
            '--data {input} -outdir {output} '
        '&> {log}'

rule multiqc:
    input:
        expand('qc/{sample}-{read}.fastqc.zip',
            sample = SAMPLES, read = READS),
        expand('qc/{sample}-{read}.fastq_screen.txt',
            sample = SAMPLES, read = READS),
        expand('qc/{sample}.cutadapt.txt', sample = SAMPLES),
        expand('qc/{sample}-{read}.trim.fastqc.zip',
            sample = SAMPLES, read = READS),
        #expand('qc/{sample}-STAR.txt', sample = SAMPLES),
        expand('qc/{sample}.stats.txt', sample = SAMPLES),
        expand('qc/{sample}.idxstats.txt', sample = SAMPLES),
        expand('qc/{sample}.flagstat.txt', sample = SAMPLES),
        expand('qc/bamqc/{sample}', sample = SAMPLES),
        expand('qc/rnaqc/{sample}', sample = SAMPLES),
    output:
        'multiqc_report.html'
    log:
        'logs/multiqc/multiqc.log'
    conda:
        f'{ENVS}/multiqc.yaml'
    shell:
        'multiqc qc/ &> {log}'
