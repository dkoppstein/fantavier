
# parameters
ONT_KIT = 'FLO-MIN106'

rule basecall:
    input: 'rel3-fast5-chr20.part05.tar'
    output: ''
    threads: 2
    shell:
        'read_fast5_basecaller.py '
        '-t {threads} '
        '-i rel3 '
        '-c r94_450bps_linear.cfg'


read_fast5_basecaller.py -t 4 -s ./BaseCalling/ -i fast5 -c r94_450bps_linear.cfg
