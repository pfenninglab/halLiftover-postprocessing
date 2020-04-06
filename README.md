# HALPER (***hal***Liftover ***P***ostprocessing for the ***E***volution of ***R***egulatory Elements)


## Running HALPER
* `python orthologFind.py` using python3 (see "Example Run of HALPER" below for an example)


## Introduction
HALPER is designed for constructing contiguous orthologs from the outputs of halLiftover.  While it was originally designed for contructing orthologs of transcription factor ChIP-seq and open chromatin peaks, it can be applied to any genomic regions of interest. Since HALPER relies on halLiftover, the assembly of the query and target genomic regions must be in a Cactus alginment hal file.


## Dependencies
* Python version 3.7 (https://www.python.org/downloads/release/python-371/)
* Python libraries `matplotlib` and `numpy`
	* numpy (http://www.numpy.org/)
		* HALPER has been tested using numpy versions 1.14.3, 1.16.0, 1.16.4, 1.16.6, and 1.18.2
	* matplotlib (https://matplotlib.org/downloads.html)
		* HALPER has been tested using matplotlib versions 1.5.1, 2.2.3, 3.2.1
* HALPER has been tested on Linux (CentOS 6, CentOS 7, and Ubuntu 18.04.4), Windows (Windows 10), and Mac


## Tips for Installing the HAL Format API
* To install, follow the instructions in this website: https://github.com/ComparativeGenomicsToolkit/hal
	* For detailed installation tips, follow the instructions in https://github.com/pfenninglab/halLiftover-postprocessing/blob/master/halliftoverInstallationSpecifics.md
* HALPER has been tested using outputs from HAL Format API Release 2.1


## Program Parameters 
* -qFile: bed file with query regions (used as input to halLiftover) containing (at least) the following information: chromosome_name, start, end, region name
	* The 1st 4 columns **MUST** be in standard bed format
	* The names in column 4 must be unique -- these names will be used in HALPER
	
* -tFile: bed file of the file specified in -qFile mapped to the target species using halLiftover (no modifications to the output from halLiftover are necessary) 
	* Line format must be: ` chr_name    peak_start    peak_end    peak_name ` (the output file from halLiftover should conform to this format)
	* Examples:
```
		chr8	55610267	55610335	peak0
		chr8	55610240	55610267	peak0
		chr8	55610220	55610240	peak0
		chr8	55610191	55610220	peak0
		chr8	55610183	55610190	peak0 
```

* -sFile: bed file of the peak summits file specified in -qFile mapped to the target species using halLiftover
	* Line format must be: ` chr_name    peak_start    peak_end    peak_name` (the output file from halLiftover should conform to this format)
	* To obtain peak summits from a narrowPeak file, add the second column to the tenth column
	* See "Preparing Histone Modification Data for HALPER" below for instructions for how to create the file for -sFile when using this program with histone modification ChIP-seq peaks or regions without peak summits
	* Examples:
```
		chr8	55609835	55609836	peak0
		chr8	55609437	55609438	peak1
		chr8	55591653	55591654	peak2
		chr8	55592205	55592206	peak4
		chr8	55536703	55536704	peak6
		chr8	55499203	55499204	peak8
		chr8	55473539	55473540	peak9 
```
	
* -narrowPeak: output files in narrowPeak format (optional argument)

* -oFile: output file name (the chromosome name and all positions in the file specified in -oFile are from the target species)
	* Line format (from left to right, if -narrowPeak option is not used): 
```
		chr_name 
		target_ortholog_start 
		target_ortholog_end 
		summit_position 
		peakname 
		target_ortholog_length 
		query_ortholog_length 
		target_summit_to_target_ortholog_start_length
		target_summit_to_target_ortholog_end_length
```

* -mult_keepone: if a region's summit maps to multiple positions in the target species, use the first position in file specified in -sFile; otherwise, such a region is discarded
* -max_len: ortholog length must be less or equal to max_len
* -max_frac: ortholog length must be less or equal to max_frac * region length 
	* provide either max_len or max_frac
* -min_len: ortholog length must be greater or equal to min_len
* -min_frac: ortholog length must be greater or equal to min_frac * region length 
	* provide either min_len or min_frac
* -protect_dist: the ortholog length in each direction from the ortholog of the summit must be at least proct_dist 
			![alt text](https://github.com/pfenninglab/multiple_alignment-python/blob/master/min_proct_dist.png)


## Example Run of HALPER
Running these examples requires the files in the examples directory and 10plusway-master.hal, a Cactus alignment with 12 mammals that can be obtained from the authors of the paper describing Cactus (see "Relevant Publications" below).  One can compare the outputs of each step to the files with the corresponding names in the examples directory.  To run only HALPER (not halLiftover), go directly to #4.
1.  Run halLiftover on the file from the query species (example is in narrowPeak format, so columns not in standard bed format are first removed) to obtain the regions' orthologs in the target species:
```
	[directory with hal]/hal/bin/halLiftover --bedType 4 [directory with Cactus alignment]/10plusway-master.hal Human [directory with halLiftover-postprocessing]/halLiftover-postprocessing/examples/hg38Peaks.bed Mouse hg38Peaks_halLiftovermm10.bed
```
2.  Get the peak summits (example is for a narrowPeak file, see "Preparing Histone Modification Data for HALPER" below for how to do this for histone modification ChIP-seq peaks or genomic regions without summits):
```
	awk 'BEGIN{OFS="\t"}{print $1, $2+$10, $2+$10+1, $4}' [directory with halLiftover-postprocessing]/halLiftover-postprocessing/examples/hg38Peaks.bed > hg38Peaks_summits.bed
```
3.  Run halLiftover on the peak summits to obtain their orthologs in the target species:
```
	[directory with hal]/hal/bin/halLiftover [directory with Cactus alignment]/10plusway-master.hal Human hg38Peaks_summits.bed Mouse hg38Peaks_summits_halLiftovermm10.bed
```
4.  Run HALPER (note that there is only one '-' for the parameter names):
```
	python [directory with halLiftover-postprocessing]/orthologFind.py -max_len 1000 -min_len 50 -protect_dist 5 -qFile [directory with halLiftover-postprocessing]/halLiftover-postprocessing/examples/hg38Peaks.bed -tFile hg38Peaks_halLiftovermm10.bed -sFile  hg38Peaks_summits_halLiftovermm10.bed -oFile hg38Peaks_halLiftovermm10_summitExtendedMin50Max1000Protect5.bed -mult_keepone
```
* Examples of output without -narrowPeak option:
```
		chr8	55609305	55610335	55609835	peak0	1031	1019	530	500
		chr8	55609305	55610335	55609437	peak1	1031	1019	132	898
```
* Examples of output with -narrowPeak option (columns 5-9 do not have meaningful values):
```
		chr8	55609305	55610335	peak0	-1	.	-1	-1	-1	530
		chr8	55609305	55610335	peak1	-1	.	-1	-1	-1	132
```


## Output Files Produced by HALPER
* File with coherent target species orthologs (name specified in -oFile)
* File with target species orthologs that did not meet all of the criteria specified by the user (name is the name specified in -oFile + ".failed")
* File with histogram of lengths of coherent target species orthologs (name is the name specified in -oFile + ".png")
* File with histogram of lengths of query species regions with coherent orthologs in the target species (name is the name specified in -oFile + "-peak.png", not created when the -narrowPeak option is used)
	* Note: To obtain ortholog length histograms when running orthologFind.py on a cluster, submit the job (or open the interactive session in which the program will be run) using x11 forwarding (option for slurm is --x11).


## Preparing Histone Modification Data for HALPER
Starting to construct target species orthologs with the target species orthologs of peak summits is sub-optimal for histone modification ChIP-seq data because, in this data, TFs are thought to bind, not where there are large numbers of reads, but in the valleys between the parts of regions with large numbers of reads. A reasonable place to start with histone modification data, therefore, is the location within the region that has the largest number of species in the alignment, as this is likely to be an important part of the region. If there are multiple such locations, which often happens, then choosing the one that is closest to the center makes sense because the centers of the histone modification regions tend to be more important than their edges. This same approach can be used for other genomic regions that do not have summits.

Here are the dependencies required for making an -sFile using this process:
* HAL Format API (https://github.com/ComparativeGenomicsToolkit/hal)
* wigToBigWig and bigWigToBedGraph (http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/)
* pybedtools (https://daler.github.io/pybedtools/main.html)
	* This has been tested using pybedtools version 0.8.0

Here is how to make an -sFile using this process:

1.  Get the alignment depth for the query species:
```
[directory with hal]/hal/bin/halAlignmentDepth --outWiggle [alignmentDepthFileName] [cactusFileName] [speciesName]
```
This can require up to 8 gigabytes for a hal file with 35 species.  Running this on 35 species can take over a week, and the output files can be at least a few gigabytes.  For a larger hal file, one can run halAlignmentDepth on each genomic region instead of on the entire genome.

2.  Convert the alignment depth file from a wig file to a bigwigh file:
```
wigToBigWig [alignmentDepthFileName] [chromSizesFileName] [alignmentDepthBigwigFileName]
```
This can require up to 64 gigabytes for an alignment depth file produced from a full genome.  Note that the chromosome naming conventions used in the alignment depth file and the chrom sizes file need to be the same, which might require converting the chromosome names in the chrom sizes file.

3.  Convert the alignment depth bigwig file to a bedgraph file:
```
bigWigToBedGraph [alignmentDepthBigwigFileName] [alignmentDepthBedgraphFileName]
```

4.  Sort the bedgraph file by chromosome, start, end:
```
sort -k1,1 -k2,2n -k3,3n [alignmentDepthBedgraphFileName] > [sortedAlignmentDepthBedgraphFileName]
```
The bedgraph files can be gzipped so that they take up less space.

5.  Get the file that will be used for starting the ortholog extension for each region using the scores in the bedgraph file:
```
python [directory with halLiftover-postprocessing]/getMaxScorePositionFromBedgraph.py --bedFileName [file with regions you will be getting scores for, will be -qFile for next step] --bedgraphFileName [sortedAlignmentDepthBedgraphFileName] --highestScoreLocationFileName [where the positions with the highest scores will be recored, you can map this with hal-liftover to create -sFile for the next step] --gz
```
This program requires the bed file and the bedgraph file to be sorted and not contain duplicated entires.  Leave out --gz if the file with the regions and the alignment depth bedgraph file are not gzipped.  Note that this program is compatible with both python version 2 and python version 3 while orthologFind.py is compatible with only python verison 3.

Alternatively, steps 2-5 can be replaced with the following script that combines them:
```
python [directory with halLiftover-postprocessing]/getMaxScorePositionFromWig.py --bedFileName [file with regions you will be getting scores for, will be -qFile for next step] --wigFileName [alignmentDepthFileName] --chromSizesFileName [chromSizesFileName] --highestScoreLocationFileName [where the positions with the highest scores will be recored, you can map this with hal-liftover to create -sFile for the next step] --gz
```
This program requires the bed file to be sorted and not contain duplicated rows.  Leave out --gz if the bed file is not gzipped.  This program is compatible with both python version 2 and python version 3.  Note that this script runs UCSC tools internally that sometimes fail silently; therefore, check the sorted bedgraph file when it finishes and re-run it with more memory alloted if that file is not large.

6.  Use halLiftover to map the positions where the highest scores are recorded to the target species.  This will create your -sFile for orthologFind.py.


## Additional Utilities
* makeRunHalLiftoverSingleBedScript.py: Makes a script that will run halLiftover on a single file and map the regions in it to a list of species
* makeRunHalLiftoverScript.py: Makes a script that will run halLiftover on a list of files and map the regions in each file to a list of species
* makeOrthologFindSingleBedScript.py: Makes a script that will run orthologFind.py on a list of target, summit file combinations for a single query file
* makeOrthologFindScript.py: Makes a script that will run orthologFind.py on a list of target, summit, query file combinations
* makePeakOrthologMatrix.py: Uses a list of outputs of orthologFind.py run on the same query file to make a species x peak matrix, where the entry at i,j is a 1 if species i has an ortholog of peak j and a 0 otherwise


## Relevant Publications
* Manuscript describing the Cactus alignment method: Benedict Paten, Dent Earl, Ngan Nguyen, Mark Diekhans, Daniel Zerbino and David Haussler. Cactus: Algorithms for genome multiple sequence alignment. Genome Research, Volume 21, Issue 9, 10 June 2011, Pages 1512-1528.
* Manuscript describing method for creating Cactus alignments for hundreds of species: Joel Armstrong, Glenn Hickey, Mark Diekhans, Alden Deran, Qi Fang, Duo Xie, et al. Progressive alignment with Cactus: a multiple-genome aligner for the thousand-genome era. bioRxiv, 9 August 2019.
* Manuscript describing HAL Format API: Glenn Hickey, Benedict Paten, Dent Earl, Daniel Zerbino, and David Haussler. HAL: A Hierarchical Format for Storing and Analyzing Multiple Genome Alignments. Bioinformatics, Volume 29, Issue 10, 15 May 2013, Pages 1341–1342.


## Contributors
* Erin Zhang (xiaoyuz1@andrew.cmu.edu)
* Irene Kaplow (ikaplow@cs.cmu.edu)
* Morgan Wirthlin (mwirthlin@cmu.edu)
* Andreas Pfenning (apfenning@cmu.edu)
