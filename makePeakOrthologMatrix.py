import sys
import argparse

def parseArgument():
        # Parse the input
        parser=argparse.ArgumentParser(description=\
                        "Make a matrix that is species by peaks, where a 1 indicates that the peak ortholog is present in the species")
        parser.add_argument("--bedFileName", required=True,\
                        help='bed file with list of original peaks')
        parser.add_argument("--orthologsFileNameListFileName", required=True,\
                        help='File name with list of files with orthologs')
        parser.add_argument("--outputFileName", required=True,\
                        help='File where ortholog matrix will be recorded')
        options = parser.parse_args()
        return options

def makePeakOrthologMatrix(options):
	# Make a matrix that is species by peaks, where a 1 indicates that the peak ortholog is present in the species and 0 indicates not present
	bedFile = open(options.bedFileName)
	peakList = [line.strip().split("\t")[3] for line in bedFile]
	bedFile.close()
	outputFile = open(options.outputFileName, 'w+')
	outputFile.write("Species File Name")
	for peak in peakList:
		# Iterate through the peaks and add the names to the output file
		outputFile.write("\t" + peak)
	orthologsFileNameListFile = open(options.orthologsFileNameListFileName)
	for line in orthologsFileNameListFile:
		# Iterate through the species and make a line in the table for each
		outputFile.write("\n")
		orthologsFileName = line.strip()
		outputFile.write(orthologsFileName)
		orthologsFile = open(orthologsFileName)
		orthologsPeakList = [line.strip().split("\t")[3] for line in orthologsFile]
		for peak in peakList:
			# Iterate through the peaks and fill out the matrix based on which peaks have orthologs in the species
			if peak in orthologsPeakList:
				# The current peak has an ortholog in the species
				outputFile.write("\t" + "1")
			else:
				# The current peak does not have an ortholog in the species
				outputFile.write("\t" + "0")
		orthologsFile.close()
	orthologsFileNameListFile.close()
	outputFile.close()

if __name__ == "__main__":
        options = parseArgument()
        makePeakOrthologMatrix(options)
