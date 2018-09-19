
from os.path import dirname
from os.path import join as pjoin

# parameters
INPUT_DIR = 'fast5'
#ONT_KIT = 'FLO-MIN106'
ALBACORE_CONFIG_FILE = 'r94_450bps_linear.cfg'
THREADS = 2
LOWER_BOUND_LENGTH = '300'
KEEP_PERCENT = '90'
CONDA = 'source activate humansv && '
CONDA3 = 'source activate humansvp3 && '
CUSTOM_CONDA3 = 'source activate /home/dkoppstein/envs/py3 && '
HUMAN_SV_SOURCE_FILE = '/mnt/humanSV/.humansv'
CHROM_20_SIZE = '63m'

# rule basecall:
#     output: '1_basecalls/.sentinel'
#     threads: THREADS
#     shell:
#         '{CONDA3} read_fast5_basecaller.py '
#         '-t {threads} '
#         '-i {INPUT_DIR} '
#         '-c {ALBACORE_CONFIG_FILE} '
#         '-s 1_basecalled '
#         '-o fastq'

# rule poretools_extract:
#     '{CONDA3} poretools fastq fast5


# rule concat:
#     input: rules.basecall.output
#     output: '2_concat/combined.fastq'
#     shell:
#         'cat 1_basecalled/workspace/fail/*.fastq ' # if using both passed and fail
#         '1_basecalled/workspace/pass/*.fastq > {output}'

# rule filter:
#     input: rules.filter.output
#     output: '2_filtered/something.txt'
#     shell:
#         '{CONDA3} NanoFilt '
#        ''

rule porechop:
    input: 'data/chr20.fastq'
    output: '1_porechop/output.fastq.gz'
    shell:
        '{CUSTOM_CONDA3} porechop '
        '-i {input} '
        '| gzip -9 > {output}'

rule plot_quals:
    input: '{somefile}.fastq'
    output: '{somefile}/NanoPlot-report.html'
    threads: THREADS
    run:
        outdir = pjoin(dirname(str())
        shell('{CUSTOM_CONDA3} NanoPlot '
              '--fastq {input} '
              '-t {threads} '
              '--loglength '
              '-o {wildcards.somefile} '
              '--maxlength 100000 '
              '--percentqual '
              '--plots hex dot'

rule filter:
    input:
        fastq=rules.porechop.output,
        reference='data/chr20.fa.gz'
    output: '2_filtered/output.fastq.gz'
    shell:
        '{CUSTOM_CONDA3} filtlong '
        '-a {input.reference} '
        '--min_length {LOWER_BOUND_LENGTH} '
        '--keep_percent {KEEP_PERCENT} '
        '{input.fastq} '
        '| gzip -9 > {output}

rule canu:
    input: rules.filter.output
    output: '3_canu/output.fastq'
    shell:
        '{CONDA} canu -d 3_canu '
        '-p output '
        'genomeSize={CHROM_20_SIZE} '
        '-nanopore-raw '
        ' {input}'
