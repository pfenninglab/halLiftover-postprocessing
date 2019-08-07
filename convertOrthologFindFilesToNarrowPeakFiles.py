import sys
import argparse

def parseArgument():
        # Parse the input
	parser=argparse.ArgumentParser(description=\
                "Make a script that will conver the outputs of orthologFind.py to narrowPeak format")
	parser.add_argument('-oFileListFileName', help='list of bed files that were outputs of orthologFind.py', \
                required=True)
	parser.add_argument('-convertedFileNameSuffix', help='suffix that will be used for narrowPeak output files, should not begin with .', \
                required=False, default="narrowPeak.gz")
	parser.add_argument("-scriptFileName", required=True,\
                help='Name of the file where the script will be recorded')
	options = parser.parse_args()
	return options

def convertOrthologFindFilesToNarrowPeakFiles(options):
	# Make a script that will conver the outputs of orthologFind.py to narrowPeak format
	oFileListFile = open(options.oFileListFileName)
	scriptFile = open(options.scriptFileName, 'w+')
	for line in oFileListFile:
		# Iterate through the files from orthologFind.py and add a line to the script to convert each
		oFile = line.strip()
		oFileElements = oFile.split(".")
		convertedFileName = ".".join(oFileElements[0:-1]) + "." + options.convertedFileNameSuffix
		scriptFile.write(" ".join(["awk \'BEGIN{OFS=\"\\t\"} {print $1, $2, $3, $5, \".\", \".\", \".\", \".\", \".\", $8}\'", oFile, "| gzip >", \
			convertedFileName]) + "\n")
	oFileListFile.close()
	scriptFile.close()

if __name__=="__main__":
	options = parseArgument()
	convertOrthologFindFilesToNarrowPeakFiles(options)
