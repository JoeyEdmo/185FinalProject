import argparse

'''
mostly lets me define the arg parsing in it's own file so it doesn't clog up my main function

'''


def handle(parser):
    # Input
    parser.add_argument('fq', help="fastq file", \
        type=str)

    # Output
    parser.add_argument("-d", "--directory", help="directory to write output to " \
        "Default: stdout", metavar="DIRECTORY", type=str, required=False)

    # Other options
    parser.add_argument("-f", "--fasta-ref", \
        help="faidx indexed reference sequence file", \
        type=str, metavar="FILE", required=False)
    parser.add_argument("-r", "--region", help="region in which pileup is " \
        "generated. Format chr:start-end", \
        type=str, metavar="REG", required=False)
    
    return parser
    