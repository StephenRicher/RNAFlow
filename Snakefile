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
        ['qc/multiqc', 'qc/multibamqc', f'genome/{BUILD}.fa.gz.fai']

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
        gtf = ANNOTATION
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
            '--sjdbOverhang {OVERHANG} --sjdbGTFfile {input.gtf} '
        '&> {log}'

rule fastqc:
    input:
        lambda wc: DATA.xs(wc.single, level = 2)['path']
    output:
        html = 'qc/fastqc/{single}.raw_fastqc.html',
        zip = 'qc/fastqc/unmod/{single}.raw.fastqc.zip'
    log:
        'logs/fastqc/{single}.log'
    wrapper:
        '0.49.0/bio/fastqc'

# Modify FASTQ filename to match {sample}-{read} for multiQC
rule modify_fastqc:
    input:
        'qc/fastqc/unmod/{sample}-{read}.raw.fastqc.zip'
    output:
        'qc/fastqc/{sample}-{read}.raw_fastqc.zip'
    params:
        name = lambda wc: f'{wc.sample}-{wc.read}'
    log:
        'logs/modify_fastqc/{sample}-{read}.raw.log'
    conda:
        f'{ENVS}/coreutils.yaml'
    shell:
        '{SCRIPTS}/modify_fastqc.sh {input} {output} {params.name} &> {log}'

rule cutadapt:
    input:
        lambda wc: DATA.xs(wc.sample, level = 1)['path']
    output:
        r1 = 'fastq/trimmed/{sample}-R1.trim.fastq.gz',
        r2 = 'fastq/trimmed/{sample}-R2.trim.fastq.gz',
        qc = 'qc/cutadapt/unmod/{sample}.cutadapt.txt'
    group:
        'cutadapt'
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
            '-o {output.r1} -p {output.r2} {input} > {output.qc} '
        '2> {log}'

rule modify_cutadapt:
    input:
        rules.cutadapt.output.qc
    output:
        'qc/cutadapt/{sample}.cutadapt.txt'
    group:
        'cutadapt'
    log:
        'logs/modify_cutadapt/{sample}.log'
    conda:
        f'{ENVS}/coreutils.yaml'
    shell:
        'awk -v sample={wildcards.sample} '
            '-f {SCRIPTS}/modify_cutadapt.awk > {output} 2> {log}'

rule fastq_screen:
    input:
        'fastq/trimmed/{sample}-{read}.trim.fastq.gz'
    output:
        txt = 'qc/fastq_screen/{sample}-{read}.fastq_screen.txt',
        png = 'qc/fastq_screen/{sample}-{read}.fastq_screen.png'
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
        html = 'qc/fastqc/{sample}-{read}.trim_fastqc.html',
        zip = 'qc/fastqc/{sample}-{read}.trim_fastqc.zip'
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
        'qc/samtools/stats/{sample}.stats.txt'
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
        'qc/samtools/idxstats/{sample}.idxstats.txt'
    group:
        'SAM_QC'
    log:
        'logs/samtools_idxstats/{sample}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools idxstats {input} > {output} 2> {log}'

rule samtools_flagstat:
    input:
        rules.samtools_sort.output
        #rules.STAR_map.output.bam
    output:
        'qc/samtools/flagstat/{sample}.flagstat.txt'
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
        bam = rules.samtools_sort.output,
        #bam = rules.STAR_map.output.bam,
        gtf = ANNOTATION
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
            '-bam {input.bam} -gff {input.gtf} -outdir {output} '
            '-nt 6 --outside-stats --java-mem-size={resources.mem_mb}M '
            '--paint-chromosome-limits '
        '&> {log}'

rule rnaqc:
    input:
        bam = rules.samtools_sort.output,
        #rules.STAR_map.output.bam,
        gtf = ANNOTATION
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
            '-bam {input.bam} -gtf {input.gtf} -outdir {output} '
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

rule featureCounts:
    input:
        rules.samtools_sort.output,
        #rules.STAR_map.output.bam,
        gtf = ANNOTATION
    output:
        'counts/{sample}.featureCounts.txt'
    log:
        'logs/featureCounts/{sample}.log'
    conda:
        f'{ENVS}/subread.yaml'
    shell:
        'featureCounts '
            '--primary -C -t exon -g gene_id '
            '-a {input.gtf} -o {output} {input} '
        '&> {log}'

rule multiqc:
    input:
        expand('qc/fastqc/{sample}-{read}.raw_fastqc.zip',
            sample = SAMPLES, read = READS),
        expand('qc/fastq_screen/{sample}-{read}.fastq_screen.txt',
            sample = SAMPLES, read = READS),
        expand('qc/cutadapt/{sample}.cutadapt.txt', sample = SAMPLES),
        expand('qc/fastqc/{sample}-{read}.trim_fastqc.zip',
            sample = SAMPLES, read = READS),
        #expand('qc/{sample}-STAR.txt', sample = SAMPLES),
        expand('qc/samtools/stats/{sample}.stats.txt', sample = SAMPLES),
        expand('qc/samtools/idxstats/{sample}.idxstats.txt', sample = SAMPLES),
        expand('qc/samtools/flagstat/{sample}.flagstat.txt', sample = SAMPLES),
        expand('qc/bamqc/{sample}', sample = SAMPLES),
        expand('qc/rnaqc/{sample}', sample = SAMPLES),
        expand('qc/rnaqc/{sample}', sample = SAMPLES),
        expand('counts/{sample}.featureCounts.txt', sample = SAMPLES)
    output:
        directory('qc/multiqc')
    log:
        'logs/multiqc/multiqc.log'
    conda:
        f'{ENVS}/multiqc.yaml'
    shell:
        'multiqc -outdir {output}'
            '--force --config {BASE}/config/multiqc_config.yaml {input} '
        '&> {log}'
