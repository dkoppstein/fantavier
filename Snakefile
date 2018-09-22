
from os.path import dirname
from os.path import join as pjoin

shell.prefix("set +u; ")

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
CONDA_QUAST = 'source activate /home/cgross/envs/quast && '
HUMAN_SV_SOURCE_FILE = '/mnt/humanSV/.humansv'
CHROM_20_SIZE = '63m'
REFERENCE_GENOME_UNZIPPED = 'db/chr20.fa'
REFERENCE_GENOME = REFERENCE_GENOME_UNZIPPED + '.gz'

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

NANOPLOT_CMD = ('{CUSTOM_CONDA3} NanoPlot '
                '--fastq {input} '
                '-t {threads} '
                '--loglength '
                '-o {params.outdir} '
                '--maxlength 100000 '
                '--percentqual '
                '--plots hex dot')

rule raw_plot_quals:
    input: 'data/chr20.fastq'
    output: 'plots/raw_quals/NanoPlot-report.html'
    params:
        outdir='plots/raw_quals'
    threads: THREADS
    shell: NANOPLOT_CMD

rule porechop_plot_quals:
    input: rules.porechop.output
    output: 'plots/porechop_quals/NanoPlot-report.html'
    params:
        outdir='plots/porechop_quals'
    threads: THREADS
    shell: NANOPLOT_CMD

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
        "{CONDA} minimap2 -x ava-ont {input} {input} > {output.paf}; "
        "{CONDA} miniasm -f {input} {output.paf} > {output.gfa}; "
        "bash scripts/fold_fasta.sh {output.gfa} > {output.fasta}"

RACON_COMMAND = ('{CONDA} minimap2 -x map-ont '
                 '-t {threads} {input.assembly} {input.filtered} '
                 '> {output.paf}; '
                 '{CUSTOM_CONDA3} racon -t {threads} '
                 '{input.filtered} {output.paf} {input.assembly} > {output.consensus}')

rule racon:
    input:
        assembly=rules.miniasm.output.fasta,
        filtered=rules.filter.output
    output:
        paf='5_racon/assembly_map.paf',
        consensus='5_racon/assembly_consensus.fasta'
    threads: THREADS
    shell: RACON_COMMAND

rule racon_round2:
    input:
        assembly=rules.racon.output.consensus,
        filtered=rules.filter.output
    output:
        paf='6_racon_round2/assembly_map.paf',
        consensus='6_racon_round2/assembly_consensus.fasta'
    threads: THREADS
    shell: RACON_COMMAND

rule gunzip:
    input: REFERENCE_GENOME
    output: temp(REFERENCE_GENOME_UNZIPPED)
    shell: 'gunzip -c {input} > {output}'

rule nucmer:
    input:
        reference=REFERENCE_GENOME_UNZIPPED,
        consensus=rules.racon_round2.output.consensus
    output: '7_nucmer/nucmer_alignment.delta'
    shell:
        '{CUSTOM_CONDA3} nucmer -maxmatch -l 100 '
        '-c 500 -p 7_nucmer/nucmer_alignment '
        '{input.reference} {input.consensus}'

rule assemblytics:
    input: rules.nucmer.output
    output: '8_assemblytics/assemblytics.sentinel'
    shell:
        'export PATH=/mnt/humanSV/software/Assemblytics:$PATH; '
        '{CONDA} Assemblytics {input} 8_assemblytics/assemblytics 10000 /mnt/humanSV/software/Assemblytics; '
        'touch {output}'

rule quast:
    input:
        reference=REFERENCE_GENOME_UNZIPPED,
        filtered=rules.filter.output,
        assembly=rules.racon.output.consensus
    output:
        report='9_quast_results/latest/report.txt'
    threads: THREADS
    shell:
        '{CONDA_QUAST} quast {input.filtered} -R {input.reference} --large -o 9_quast '

## alignment branch
rule ngmlr:
    input:
        reference='db/chr20.fa',
        reads=rules.filter.output
    output: 'alignment/ngmlr/output.sorted.bam'
    threads: THREADS
    shell:
        'rm -rf alignment/ngmlr/*.bam; {CONDA} gunzip -c {input.reads} > alignment/ngmlr/unzipped_reads.tmp.fastq; '
        'ngmlr -t {threads} -r {input.reference} -x ont -q alignment/ngmlr/unzipped_reads.tmp.fastq -o alignment/ngmlr/output.tmp.sam; '
        'rm -f alignment/ngmlr/unzipped_reads.tmp.fastq; '
        'samtools sort -@ {threads} -o output.sorted.bam alignment/ngmlr/output.tmp.sam; '
        'rm -f alignment/ngmlr/output.tmp.sam'

rule samtools_index:
    input: rules.ngmlr.output
    output: 'alignment/ngmlr/output.sorted.bam.bai'
    shell:
        '{CONDA} samtools index {input}'

rule sniffles:
    input:
        bam=rules.ngmlr.output,
        idx=rules.samtools_index.output
    output: 'alignment/sniffles/output.vcf'
    shell:
        '{CONDA} source {HUMAN_SV_SOURCE_FILE} && sniffles -m {input.bam} -v {output}'

rule survivor:
    input: rules.assemblytics.output
    output: '10_survivor/assemblytics_sv.vcf'
    params:
        bedfile=str(rules.assemblytics.output).replace('sentinel', 'Assemblytics_structural_variants.bed')
    shell:
        '{CONDA_QUAST} SURVIVOR convertAssemblytics '
        '{params.bedfile} 0 {output}'


rule vcf2tab:
    input: '{somefile}.vcf'
    output: '{somefile}.tab'
    shell:
        '{CONDA} python scripts/vcf2tab.py {input} > {output}'

rule svpv:
    input:
        bam=rules.ngmlr.output,
        vcf=rules.sniffles.output
    output: '11_svpv/sentinel.txt'
    shell:
        '{CONDA_QUAST} SVPV -vcf {input.vcf} -aln {input.bam} -o 11_svpv -samples output; '
        'touch {output}'


rule all:
    input:
        assemblytics=str(rules.survivor.output).replace('vcf', 'tab'),
        sniffles=str(rules.sniffles.output).replace('vcf', 'tab'),
        plot_quals=[rules.porechop_plot_quals.output, rules.raw_plot_quals.output]
