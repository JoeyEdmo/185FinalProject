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

    
    return parser
    