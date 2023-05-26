#!/usr/bin/env python3

"""
Command-line script to perform  pileup of BAM files

Similar to samtools mpileup
"""

import argparse
import os
import Bio
import sys
import parseArgs

def main():
	parser = argparse.ArgumentParser(
		prog="ffqc",
		description="Command-line script to visualize quality of fastqc file"
	)
	#push extra code into it's own file.
	parser = parseArgs.handle(parser)
	args = parser.parse_args()
	print(args)

if __name__ == "__main__":
    main()
