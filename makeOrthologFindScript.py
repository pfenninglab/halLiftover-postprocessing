import sys
import argparse

def parseArgument():
        # Parse the input
	parser=argparse.ArgumentParser(description=\
		"Make a script that will run orthologFind.py on a list of target, summit, query file combinations")
	parser.add_argument('-max_len', required=False, \
		help='maximum number of base pairs of the ortholog')
	parser.add_argument('-max_frac', required=False, \
		help='maximum percentage of original peak of the ortholog')
	parser.add_argument('-protect_dist', required=False, help='summit protection distance', \
		default=50)
	parser.add_argument('-min_len', required=False, \
		help='minimum number of base pairs of the ortholog')
	parser.add_argument('-min_frac', required=False, \
		help='minimum percentage of original peak of the ortholog')
	parser.add_argument('-qFileListFileName', help='list of input bed files, file can have other columns if the bed file name is in the first column', \
		required=True)
	parser.add_argument('-tFileListFileName', help='list of input mapped bed files, lines should correspond to lines of qFileListFile', \
		required=True)
	parser.add_argument('-sFileListFileName', help='list of input mapped-summit bed files, lines should correspond to lines of qFileListFile', \
		required=True)
	parser.add_argument('-oFileNameSuffix', help='suffix to add to target file names to create the output file names, should end in .bed', \
		required=True)
	parser.add_argument('-mult_keepone', action="store_true", \
		help='if a region\'s summit maps to multiple positions in the target species, use the first position in file specified in -sFile', \
                required=False)
	parser.add_argument('-narrowPeak', action="store_true", help='output files in narrowPeak format', \
                required=False)
	parser.add_argument("-codePath", required=True,\
		help='Path to orthologFind.py')
	parser.add_argument("-scriptFileName", required=True,\
		help='Name of the file where the script will be recorded')
	options = parser.parse_args()
	return options

def makeOrthologFindScript(options):
	# Make a script that will run orthologFind.py on a list of target, summit, query file combinations
	qFileListFile = open(options.qFileListFileName)
	tFileListFile = open(options.tFileListFileName)
	sFileListFile = open(options.sFileListFileName)
	scriptFile = open(options.scriptFileName, 'w+')
	for qFileStr, tFileStr, sFileStr in zip(qFileListFile, tFileListFile, sFileListFile):
		qFile = qFileStr.strip().split("\t")[0]
		tFile = tFileStr.strip()
		sFile = sFileStr.strip()
		tFileNameElements = tFile.split(".")
		oFile = ".".join(tFileNameElements[0:-1]) + "_" + options.oFileNameSuffix
		cmd = options.codePath + "/orthologFind.py"
		scriptFile.write(" ".join(["python", cmd, "-qFile", qFile, "-tFile", tFile, "-sFile", sFile, "-oFile", oFile]))
		if options.max_len != None:
			# Add the max_len option
			scriptFile.write(" -max_len " + options.max_len)
		if options.max_frac != None:
			# Add the max_frac option
			scriptFile.write(" -max_frac " + options.max_frac)
		if options.protect_dist != None:
			# Add the protect_dist option
			scriptFile.write(" -protect_dist " + options.protect_dist)
		if options.min_len != None:
			# Add the min_len option
			scriptFile.write(" -min_len " + options.min_len)
		if options.min_frac != None:
			# Add the min_frac option
			scriptFile.write(" -min_frac " + options.min_frac)
		if options.mult_keepone:
                        # Add the min_frac option
                        scriptFile.write(" -mult_keepone ")
		if options.narrowPeak:
			# Add the narrowPeak option
			scriptFile.write(" -narrowPeak")
		scriptFile.write("\n")
	qFileListFile.close()
	tFileListFile.close()
	sFileListFile.close()
	scriptFile.close()

if __name__=="__main__":
        options = parseArgument()
        makeOrthologFindScript(options)
