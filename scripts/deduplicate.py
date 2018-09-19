
from Bio import SeqIO
import argparse
import sys
from collections import defaultdict

def main(argv):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--infile', type=argparse.FileType('r'), default=sys.stdin)
    parser.add_argument('-o', '--outfile', type=argparse.FileType('w'), default=sys.stdout)
    opts = parser.parse_args(argv)
    d = defaultdict(int)
    for read in SeqIO.parse(opts.infile, 'fastq'):
        d[read.description] += 1
        read.description = str(read.description) + '-' + str(d[read.description])
        print(read.format('fastq'), file=opts.outfile)

if __name__ == '__main__':
    main(sys.argv[1:])
