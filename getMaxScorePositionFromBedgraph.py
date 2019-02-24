import sys
import gzip
import argparse

def parseArgument():
        # Parse the input
        parser=argparse.ArgumentParser(description=\
                        "Get the position with the max. score from a bedgraph file for each region; if there is a tie, choose the center-most score")
        parser.add_argument("--bedFileName", required=True,\
                        help='Bed file with the regions; should be sorted by chromosome, start, end')
        parser.add_argument("--bedgraphFileName", required=True,\
                        help='bedgraph file with the scores; should be sorted by chromosome, start, end')
        parser.add_argument("--gz", action="store_true", required=False,\
                        help='The input bed and bedgraph files are gzipped')
        parser.add_argument("--highestScoreLocationFileName", required=True,\
                        help='bed file where position with the highest score that is closest to the center will be written')
        options = parser.parse_args()
        return options

def getMaxScorePositionFromBedgraph(options):
	# Get the position with the max. score from a bedgraph file for each region; if there is a tie, choose the center-most score
	bedFile = None
	bedgraphFile = None
	if options.gz == 1:
		# Open the input files using the gzip library
		bedFile = gzip.open(options.bedFileName)
                bedgraphFile = gzip.open(options.bedgraphFileName)
	else:
		# Open the input files normally
		bedFile = open(options.bedFileName)
		bedgraphFile = open(options.bedgraphFileName)
	highestScoreLocationFile = open(options.highestScoreLocationFileName, 'w+')
	bedgraphChrom = "chr0"
	bedgraphStart = 0
	bedgraphEnd = 0
	bedgraphMid = 0
	bgs = 0.0
	bedgraphLine = ""
	lastChrom = "chr0"
	for line in bedFile:
		# Iterate through the regions and find the position with the maximum score for each
		lineElements = line.split("\t")
		chrom = lineElements[0]
		if chrom != lastChrom:
			lastChrom = chrom
		start = int(lineElements[1])
		end = int(lineElements[2])
		peakName = lineElements[3]
		mid = int(round(float(end + start)/2.0))
		bedgraphScores = []
		highestbedgraphScorePosition = (chrom, mid, mid + 1, peakName, 0.0)
		highestbedgraphScore = 0.0
		stopReached = False
		while chrom > bedgraphChrom:
			# At the wrong chromosome, so read through bedgraph file until the next chromosome is reached
			bedgraphLine = bedgraphFile.readline()
			bedgraphLineElements = bedgraphLine.split("\t")
			if len(bedgraphLineElements) < 4:
				# Record the center because there is no score information
				highestScoreLocationFile.write("\t".join([highestbedgraphScorePosition[0], str(highestbedgraphScorePosition[1]), \
					str(highestbedgraphScorePosition[2]), highestbedgraphScorePosition[3], str(highestbedgraphScorePosition[4])])\
					+ "\n")
				stopReached = True
				break
			bedgraphChrom = bedgraphLineElements[0]  
			bedgraphStart = int(bedgraphLineElements[1])
			bedgraphEnd = int(bedgraphLineElements[2])
			bedgraphMid = int(round(float(bedgraphStart + bedgraphEnd)/2.0))
		while (not stopReached) and ((bedgraphMid < start) and (bedgraphChrom == chrom)):
			# Go through bases until the start is reached
			bedgraphLine = bedgraphFile.readline()
			bedgraphLineElements = bedgraphLine.split("\t")
			if len(bedgraphLineElements) < 4:
                                # Record the center because there is no score information
                                highestScoreLocationFile.write("\t".join([highestbedgraphScorePosition[0], str(highestbedgraphScorePosition[1]), \
					str(highestbedgraphScorePosition[2]), highestbedgraphScorePosition[3], str(highestbedgraphScorePosition[4])])\
					+ "\n")
				stopReached = True
				break
			bedgraphChrom = bedgraphLineElements[0]
			bedgraphStart = int(bedgraphLineElements[1])
			bedgraphEnd = int(bedgraphLineElements[2])
                        bedgraphMid = int(round(float(bedgraphStart + bedgraphEnd)/2.0))
			bgs = float(bedgraphLineElements[3])
		if (not stopReached) and ((bedgraphChrom > chrom) or (bedgraphMid >= end)):
                        # Record the center because there is no score information
                        highestScoreLocationFile.write("\t".join([highestbedgraphScorePosition[0], str(highestbedgraphScorePosition[1]), \
				str(highestbedgraphScorePosition[2]), highestbedgraphScorePosition[3], str(highestbedgraphScorePosition[4])]) + "\n")
                        continue
		while (not stopReached) and ((bedgraphMid < end) and (bedgraphChrom == chrom)):
			# ASSUMES THAT THERE ARE NO IDENTICAL REGIONS
			# Gets scores of every base in the region
			assert(bedgraphMid >= start)
			bedgraphScores.append(bgs)
			if bgs > highestbedgraphScore:
                        	# The current score is the highest bedgraph score so far
                                highestbedgraphScore = bgs
                                highestbedgraphScorePosition = (chrom, bedgraphMid, bedgraphMid + 1, peakName, bgs)
			if bgs == highestbedgraphScore:
				# The current score is equal to previous highest bedgraph score, so check if position is closer to center
				if abs(bedgraphMid - mid) < abs(highestbedgraphScorePosition[1] - mid):
					# The current position is closer to the center, so change highest bedgraph score position to current poistion
					highestbedgraphScorePosition = (chrom, bedgraphMid, bedgraphMid + 1, peakName, bgs)
			bedgraphLine = bedgraphFile.readline()
			bedgraphLineElements = bedgraphLine.split("\t")
			if len(bedgraphLineElements) < 4:
                                # Record the center because score information is not available for the entire region
                                highestScoreLocationFile.write("\t".join([highestbedgraphScorePosition[0], str(highestbedgraphScorePosition[1]), \
                                        str(highestbedgraphScorePosition[2]), highestbedgraphScorePosition[3], str(highestbedgraphScorePosition[4])])\
                                        + "\n")
				stopReached = True
				break
			bedgraphChrom = bedgraphLineElements[0]
			bedgraphStart = int(bedgraphLineElements[1])
			bedgraphEnd = int(bedgraphLineElements[2])
                        bedgraphMid = int(round(float(bedgraphStart + bedgraphEnd)/2.0))
			bgs = float(bedgraphLineElements[3])
		if (len(bedgraphScores) > 0) and (not stopReached):
			# Scores are available for the entire region, so record the position with the highest score
			highestScoreLocationFile.write("\t".join([highestbedgraphScorePosition[0], str(highestbedgraphScorePosition[1]), \
				str(highestbedgraphScorePosition[2]), highestbedgraphScorePosition[3], str(highestbedgraphScorePosition[4])]) + "\n")
	bedFile.close()
	bedgraphFile.close()
	highestScoreLocationFile.close()

if __name__=="__main__":
	options = parseArgument()
	getMaxScorePositionFromBedgraph(options)
