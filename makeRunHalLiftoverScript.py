import sys
import argparse

def parseArgument():
        # Parse the input
        parser=argparse.ArgumentParser(description=\
                        "Make a script that will run halLiftover on a list of files and map each file to a list of species")
        parser.add_argument("--bedToSpeciesFileName", required=True,\
                        help='List of bed files and their corresponding species')
        parser.add_argument("--speciesToLiftFileName", required=True,\
                        help='File with target species, add a second column to specifiy a different output suffix than the species name')
        parser.add_argument("--CactusFileName", required=True,\
                        help='Name of the Cactus file')
        parser.add_argument("--halLiftoverPath", required=False, default="/home/ikaplow/RegulatoryElementEvolutionProject/src/hal/bin",\
                        help='Path to hal-Liftover executable')
	parser.add_argument("--gz", action="store_true", required=False,\
                        help='The input file is gzipped')
        parser.add_argument("--scriptFileName", required=True,\
                        help='Name of the file where the script will be recorded')
        options = parser.parse_args()
        return options

def makeRunHalLiftoverScript(options):
	# Make a script that will run halLiftover on a list of files and map each file to a list of species
	bedToSpeciesFile = open(options.bedToSpeciesFileName)
	speciesToLiftFile = open(options.speciesToLiftFileName)
	speciesToLiftDict = {}
	for line in speciesToLiftFile:
		# Iterate through the species to lift and add each and its corresponding output suffix to the dictionary
		lineElements = line.strip().split("\t")
		if len(lineElements) > 1:
			# The output suffix for the species is different from the species name
			speciesToLiftDict[lineElements[0]] = lineElements[1]
		else:
			speciesToLiftDict[lineElements[0]] = lineElements[0]
	speciesToLiftFile.close()
	halLiftoverCmd = options.halLiftoverPath + "/halLiftover"
	scriptFile = open(options.scriptFileName, 'w+')
	for line in bedToSpeciesFile:
		# Iterate through the query files and write a line in the script for lifting over each to each target species
		lineElements = line.strip().split("\t")
		bedFileNameElements = lineElements[0].split(".")
		outputFileNamePrefix = ".".join(bedFileNameElements[0:-2])
		for species in speciesToLiftDict:
			# Iterate through the species and write a line in the script for lifting over the query species to each target species
			if species == lineElements[1]:
				# The current species is the same as the species of the bed file, so skip it
				continue
			outputFileName = outputFileNamePrefix + "_" + speciesToLiftDict[species] + ".bed"
			if options.gz:
				# The input bed file is gzipped
				scriptFile.write(" ".join(["zcat", lineElements[0], "|", halLiftoverCmd, "--inBedVersion 4", options.CactusFileName, \
					lineElements[1], "stdin", species, outputFileName]) + "\n")
			else:
				# The input bed file is not gzipped
				scriptFile.write(" ".join([halLiftoverCmd, "--inBedVersion 4", options.CactusFileName, lineElements[1], \
					lineElements[0], species, outputFileName]) + "\n")
	bedToSpeciesFile.close()
	scriptFile.close()

if __name__=="__main__":
        options = parseArgument()
	makeRunHalLiftoverScript(options)
