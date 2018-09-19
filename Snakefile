
from os.path import dirname

# parameters
INPUT_DIR = 'fast5'
#ONT_KIT = 'FLO-MIN106'
ALBACORE_CONFIG_FILE = 'r94_450bps_linear.cfg'
THREADS = 2
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
    input: 'input/chr20.fastq'
    output: '1_porechop/output.fastq'
    shell:
        '{CUSTOM_CONDA3} porechop '
        '-i {input} '
        '-o 1_porechop/output.fastq'

rule canu:
    input: rules.porechop.output
    output: '2_canu/output.fastq'
    shell:
        '{CONDA} canu -d 2_canu '
        '-p output '
        '-nanopore-raw '
        'genomeSize={CHROM_20_SIZE} {input}'
