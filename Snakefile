
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
QUAST = 'source activate /home/cgross/envs/quast && '
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

rule preprocess:
    input: rules.porechop.output
    output: '2_preprocessed/output.fastq.gz'
    shell:
        '{CONDA3} gunzip -c {input} | '
        'python3 scripts/deduplicate.py '
        '| gzip -9 > {output}'

rule plot_quals:
    input: '{somefile}.fastq'
    output: '{somefile}/NanoPlot-report.html'
    threads: THREADS
    shell:
        '{CUSTOM_CONDA3} NanoPlot '
        '--fastq {input} '
        '-t {threads} '
        '--loglength '
        '-o {wildcards.somefile} '
        '--maxlength 100000 '
        '--percentqual '
        '--plots hex dot'

rule filter:
    input:
        fastq=rules.preprocess.output,
        reference='db/chr20_GRCh38.fa.gz'
    output: '2_filtered/output.fastq.gz'
    shell:
        '{CUSTOM_CONDA3} filtlong '
        '-a {input.reference} '
        '--min_length {LOWER_BOUND_LENGTH} '
        '--keep_percent {KEEP_PERCENT} '
        '{input.fastq} '
        '| gzip -9 > {output}'

rule canu:
    input: rules.filter.output
    output: '3_canu/output.fasta'
    shell:
        '{CONDA} canu -d 3_canu '
        '-p output '
        'genomeSize={CHROM_20_SIZE} '
        '-nanopore-raw '
        ' {input}'

# from https://github.com/caspargross/hybridAssembly/blob/master/main.nf
rule miniasm:
    input: rules.filter.output
    output:
        paf='4_miniasm/overlap.paf',
        fasta='4_miniasm/assembled.fasta',
        gfa='4_miniasm/miniasm_graph.gfa'
    shell:
        "set +u; {CONDA} minimap2 -x ava-ont {input} {input} > {output.paf}; "
        "{CONDA} miniasm -f {input} {output.paf} > {output.gfa}; "
        "bash scripts/fold_fasta.sh {output.gfa} > {output.fasta}"

rule racon:
    input:
        assembly=rules.miniasm.output.fasta,
        filtered=rules.filter.output
    output:
        paf='5_racon/assembly_map.paf',
        consensus='5_racon/assembly_consensus.fasta'
    threads: THREADS
    shell:
        'set +u; {CONDA} minimap2 -x map-ont '
        '-t {threads} {input.assembly} {input.filtered} '
        '> {output.paf}; '
        '{CUSTOM_CONDA3} racon -t {threads} '
        '{input.filtered} {output.paf} {input.assembly} {output.consensus}'

# rule nucmer:
#     input:
#     output:
#     shell:
#         'nucmer -maxmatch -l 100 '
#         '-c 500 -p my_fav_alignment '
#         'data/Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa data/pilon_yeast_Aug2018.fasta'
