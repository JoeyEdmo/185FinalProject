#!/usr/bin/env python

"""
Command-line script to perform  pileup of BAM files

Similar to samtools mpileup
"""
import numpy as np
import pandas as pd
import argparse
import os
import Bio
import sys
from . import parseArgs
from . import quality

def main():
	parser = argparse.ArgumentParser(
		prog="ffqc",
		description="Command-line script to visualize quality of fastqc file"
	)
	#push extra code into it's own file.
	parser = parseArgs.handle(parser)
	args = parser.parse_args()
	df = quality.dfScores(args.fq) #get df of all reads
	quality.mainPlot(df, args)#use df to plot main plot

if __name__ == "__main__":
    main()
