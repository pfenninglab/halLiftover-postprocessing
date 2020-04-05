import sys
import argparse

def parseArgument():
        # Parse the input
        parser=argparse.ArgumentParser(description=\
                        "Make a script that will run halLiftover on a single file and map it to a list of species")
        parser.add_argument("--bedFileName", required=True,\
                        help='bed file')
        parser.add_argument("--querySpecies", required=True,\
                        help='Name of species that will be lifted')
        parser.add_argument("--speciesToLiftFileName", required=True,\
                        help='File with target species, add a second column to specifiy a different output suffix than the species name')
        parser.add_argument("--CactusFileName", required=True,\
                        help='Name of the Cactus file')
        parser.add_argument("--halLiftoverPath", required=False, default="/home/ikaplow/RegulatoryElementEvolutionProject/src/hal/bin",\
                        help='Path to hal-Liftover executable')
        parser.add_argument("--numInputFilePartsToRemoveForOutput", type=int, required=False, default=2,\
                        help='Number of parts of the input file name to remove from the end when creating the output file name')
        parser.add_argument("--gz", action="store_true", required=False,\
                        help='The input file is gzipped')
        parser.add_argument("--scriptFileName", required=True,\
                        help='Name of the file where the script will be recorded')
        options = parser.parse_args()
        return options

def makeRunHalLiftoverSingleBedScript(options):
	# Make a script that will run halLiftover on a list of files and map each file to a list of species
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
	bedFileNameElements = options.bedFileName.split(".")
	outputFileNamePrefix = ".".join(bedFileNameElements[0:0-options.numInputFilePartsToRemoveForOutput])
	for species in speciesToLiftDict:
		# Iterate through the species and write a line in the script for lifting over the query species to each target species
		if species == options.querySpecies:
			# The current species is the same as the species of the bed file, so skip it
			continue
		outputFileName = outputFileNamePrefix + "_" + speciesToLiftDict[species] + ".bed"
		if options.gz:
			# The input bed file is gzipped
			scriptFile.write(" ".join(["zcat", options.bedFileName, "|", halLiftoverCmd, "--bedType 4", options.CactusFileName, \
				options.querySpecies, "stdin", species, outputFileName]) + "\n")
		else:
			# The input bed file is not gzipped
			scriptFile.write(" ".join([halLiftoverCmd, "--bedType 4", options.CactusFileName, options.querySpecies, \
				options.bedFileName, species, outputFileName]) + "\n")
	scriptFile.close()

if __name__=="__main__":
        options = parseArgument()
        makeRunHalLiftoverSingleBedScript(options)
