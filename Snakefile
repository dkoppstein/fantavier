
# parameters
INPUT_DIR = 'fast5'
#ONT_KIT = 'FLO-MIN106'
ALBACORE_CONFIG_FILE = 'r94_450bps_linear.cfg'
THREADS = 2
CONDA = 'source activate humansv && '
CONDA3 = 'source activate humansvp3 && '

rule basecall:
    output: '1_basecalls/'
    threads: THREADS
    shell:
        '{CONDA3} read_fast5_basecaller.py '
        '-t {threads} '
        '-i {INPUT_DIR} '
        '-c {ALBACORE_CONFIG_FILE} '
        '-s 1_basecalled/'

rule filter:
    input: rules.basecall.output
    output: '2_filtered/something.txt'
    shell:
        '{CONDA3} NanoFilt '
#        ''
