import sys
import argparse
import os

def parseArgument():
        # Parse the input
        parser=argparse.ArgumentParser(description=\
                        "Make a matrix for each bed file that is species by peaks, where a 1 indicates that the peak ortholog is present in the species")
        parser.add_argument("--bedFileNameToOrthologsFileNameFileName", required=True,\
                        help='2-column file, where the 1st is query bed file names and the 2nd is corresponding target species ortholog bed file names')
        parser.add_argument("--outputFileNamePrefix", required=True,\
                        help='Prefix of file where ortholog matrix will be recorded, should not end in _')
	parser.add_argument("--outputFileNameSuffixDictFileName", required=True,\
                        help='2-column file, where the 1st is the query file name and the 2nd is the suffix of file where ortholog matrix will be recorded')
        options = parser.parse_args()
        return options

def makePeakOrthologMatrixList(options):
	# Make a matrix for each bed file that is species by peaks, where a 1 indicates that the peak ortholog is present in the species
	bedFileNameToOrthologsFileNameFile = open(options.bedFileNameToOrthologsFileNameFileName)
	bedFileNameToOrthologsFileNameDict = {}
	for line in bedFileNameToOrthologsFileNameFile:
		# Iterate through the query bed and orthologs files and make a dictionary mapping the bed files to the orthologs files
		lineElements = line.strip().split()
		bedFileName = lineElements[0]
		if bedFileName in bedFileNameToOrthologsFileNameDict:
			# The query bed file has already been added to the dictonary, so add the ortholog file to its value list
			bedFileNameToOrthologsFileNameDict[bedFileName].append(lineElements[1])
		else:
			# Make a new entry for the query bed file whose value is a list consisting of the ortholog file
			bedFileNameToOrthologsFileNameDict[bedFileName] = [lineElements[1]]
	bedFileNameToOrthologsFileNameFile.close()
	bedFileNameToOutputFileNameDict = {}
	outputFileNameSuffixDictFile = open(options.outputFileNameSuffixDictFileName)
	for line in outputFileNameSuffixDictFile:
		# Iterate through the file mapping query bed files to output file name suffixes and add each mapping to the dictionary
		lineElements = line.strip().split("\t")
		outputFileName = options.outputFileNamePrefix + lineElements[1]
		bedFileNameToOutputFileNameDict[lineElements[0]] = outputFileName
	for bedFileName in bedFileNameToOrthologsFileNameDict:
		# Iterate through the query bed files and make a matrix for each
		bedFile = open(bedFileName)
		peakList = [line.strip().split("\t")[3] for line in bedFile]
		bedFile.close()
		outputFile = open(bedFileNameToOutputFileNameDict[bedFileName], 'w+')
		outputFile.write("Species File Name")
		for peak in peakList:
			# Iterate through the peaks and add the names to the output file
			outputFile.write("\t" + peak)
		for orthologsFileName in bedFileNameToOrthologsFileNameDict[bedFileName]:
			# Iterate through the species and make a line in the table for each
			outputFile.write("\n")
			outputFile.write(orthologsFileName)
			orthologsPeakList = []
			if os.path.isfile(orthologsFileName):
				# At least 1 peak has an ortholog
				orthologsFile = open(orthologsFileName)
				orthologsPeakList = [line.strip().split("\t")[3] for line in orthologsFile]
				orthologsFile.close()
			else:
				print(orthologsFileName + " does not exist")
			for peak in peakList:
				# Iterate through the peaks and fill out the matrix based on which peaks have orthologs in the species
				if peak in orthologsPeakList:
					# The current peak has an ortholog in the species
					outputFile.write("\t" + "1")
				else:
					# The current peak does not have an ortholog in the species
					outputFile.write("\t" + "0")
		outputFile.close()

if __name__ == "__main__":
        options = parseArgument()
        makePeakOrthologMatrixList(options)
