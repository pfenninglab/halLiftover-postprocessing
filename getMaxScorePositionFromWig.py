import sys
import gzip
import argparse
from argparse import Namespace
import os
import pybedtools as bt
import shutil
from getMaxScorePositionFromBedgraph import getMaxScorePositionFromBedgraph

def parseArgument():
        # Parse the input
	parser=argparse.ArgumentParser(description=\
		"Get the position with the max. score from a wig file for each region; if there is a tie, choose the center-most score")
	parser.add_argument("--bedFileName", required=True,\
		help='Bed file with the regions; should be sorted by chromosome, start, end')
	parser.add_argument("--wigFileName", required=True,\
		help='wig file name with the scores')
	parser.add_argument("--chromSizesFileName", required=True,\
		help='chrom sizes file name for the query genome')
	parser.add_argument("--bigwigFileName", required=False,\
		help='bigwig file with the scores that will be outputted as intermediate')
	parser.add_argument("--bedgraphFileName", required=False,\
		help='bedgraph file with the scores that will be outputted as intermediate; will be sorted by chromosome, start, end')
	parser.add_argument("--gz", action="store_true", required=False,\
		help='The input bed file is gzipped and the intermediate sorted bedgraph file will be gzipped')
	parser.add_argument("--highestScoreLocationFileName", required=True,\
		help='bed file where position with the highest score that is closest to the center will be written')
	options = parser.parse_args()
	return options

def getMaxScorePositionFromWig(options):
	# Get the position with the max. score from a wig file for each region; if there is a tie, choose the center-most score
	bigwigFileName = options.bigwigFileName
	wigFileNameElements = options.wigFileName.strip().split(".")
	wigFileNamePrefix = ".".join(wigFileNameElements[0:-1])
	if not bigwigFileName:
		# Create the name of the intermediate bigwig file
		bigwigFileName = wigFileNamePrefix + ".bw"
	os.system(" ".join(["wigToBigWig", options.wigFileName, options.chromSizesFileName, bigwigFileName]))
	unsortedBedgraphFileName = wigFileNamePrefix + ".bedgraph"
	os.system(" ".join(["bigWigToBedGraph", bigwigFileName, unsortedBedgraphFileName]))
	bedgraphFileName = options.bedgraphFileName
	if not bedgraphFileName:
		# Create the name of the intermediate bedgraph file
		bedgraphFileName = wigFileNamePrefix + "_sorted.bedgraph"
	elif options.gz:
		# Remove the gz from the end of the bedgraph file name
		bedgraphFileName = bedgraphFileName[0:-3]
	bt.BedTool(unsortedBedgraphFileName).sort().saveas(bedgraphFileName)
	os.remove(unsortedBedgraphFileName)
	if options.gz:
		# gzip the bedgraph file
		gzipBedgraphFileName = options.bedgraphFileName
		if not gzipBedgraphFileName:
			# Create the name of the intermediate gzipped bedgraph file
			gzipBedgraphFileName = bedgraphFileName + ".gz"
		with open(bedgraphFileName, "rb") as f_in, gzip.open(gzipBedgraphFileName, "wb") as f_out:
			# Compress the bedgraph file
			shutil.copyfileobj(f_in, f_out)
		os.remove(bedgraphFileName)
		bedgraphFileName = gzipBedgraphFileName
	maxScorePositionOptions = Namespace(bedFileName=options.bedFileName, bedgraphFileName=bedgraphFileName, gz=options.gz, \
                highestScoreLocationFileName=options.highestScoreLocationFileName)
	getMaxScorePositionFromBedgraph(maxScorePositionOptions)

if __name__=="__main__":
	bt.helpers.set_tempdir("/scratch")
	options = parseArgument()
	getMaxScorePositionFromWig(options)
	bt.helpers.cleanup()
