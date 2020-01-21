## HALPER


# Running HALPER
* `python orthologFind.py` using python3


# Introduction
HALPER is designed for constructing coherent orthologs from the outputs of halLiftover.  While it was originally designed for contructing orthologs of transcription factor ChIP-seq and open chromatin peaks, it can be applied to any genomic regions of interest. Since HALPER relies on halLiftover, the assembly of the query and target genomic regions must be in a Cactus alginment hal file.  (If the assemblies are not in a Cactus alignment, liftOver can be used to map regions to the closest available assembly in a Cactus alignment.)


# Dependencies
* Python version 3.7 (https://www.python.org/downloads/release/python-371/)
* Python libraries `matplotlib` and `numpy`
	* matplotlib (https://matplotlib.org/downloads.html)
	* numpy (http://www.numpy.org/)


# Tips for Installing hal toolkit
* To install, follow the instructions in this website: https://github.com/ComparativeGenomicsToolkit/hal
	* For detailed installation tips, follow the instructions in https://github.com/pfenninglab/halLiftover-postprocessing/blob/master/halliftoverInstallationSpecifics.txt
* hal toolkit requires gcc >= 4.9 


# Program Parameters 
* -qFile: the query bed file (used as input to halLiftover) containing information on (at least): chromosome_name, start, end, region name
	* The 1st 4 columns **MUST** be in standard bed format
	* The names in column 4 must be unique -- these names will be used in HALPER
	
* -tFile: bed file of the file specified in -qFile mapped to the target species using halLiftover 
	* Line format must be: ` chr_name    peak_start    peak_end    peak_name ` (halLiftover should output file conforming to this format)
	* Examples:
		```
		chr8	55610267	55610335	peak0
		chr8	55610240	55610267	peak0
		chr8	55610220	55610240	peak0
		chr8	55610191	55610220	peak0
		chr8	55610183	55610190	peak0 
		```

* -sFile: bed file of the peak summits file specified in -qFile mapped to the target species using halLiftover
	* Line format must be: ` chr_name    peak_start    peak_end    peak_name`. halLiftover should output file conforming to this format. 
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
	* See Preparing Histone Modification Data for HALPER below for instructions for how to create the sFile if you are using this program with histone modification ChIP-seq peaks or regions without peak summits

* -oFile: output file name
	* Line format (from left to right): 
		```
		chr_name 
		ortholog_start 
		ortholog_end 
		summit_position 
		peakname 
		ortholog_length 
		original_peak_length 
		summit_to_ortholog_start_length
		summit_to_ortholog_end_length
		```
	The chromosome name and all positions in oFile are from the target species.
	* Examples:
		```
		chr8	55609305	55610335	55609835	peak0	1031	1019	530	500
		chr8	55609305	55610335	55609437	peak1	1031	1019	132	898
		```


* -max_len: ortholog length must be less or equal to max_len
* -max_frac: ortholog length must be less or euqal to max_frac * peak length 
	* provide either max_len or max_frac
* -min_len: ortholog length must be greater or equal to min_len
* -min_frac: ortholog length must be greater or equal to min_frac * peak length 
	* provide either min_len or min_frac
* -protect_dist: the ortholog length in each direction from the ortholog of the summit must be at least proct_dist 
			![alt text](https://github.com/pfenninglab/multiple_alignment-python/blob/master/min_proct_dist.png)


# Example Run of HALPER
```
python orthologFind.py -max_len 1000 -min_len 50 -protect_dist 5 -qFile hg38Peaks.bed -tFile hg38Peaks_halLiftovermm10.bed -sFile  hg38Peaks_summits_halLiftovermm10.bed -oFile hg38Peaks_halLiftovermm10_summitExtendedMin50Max1000Protect5.bed
```
Note how there is only one '-' (dash) for the parameter name. 


# Output Files Produced by HALPER
* File with coherent orhtologs (name specified in -oFile)
* File with orthologs that did not meet all of the criteria specified by the user (name is the name specified in -oFile + ".failed")
* File with histogram of ortholog lengths of all orthologs, including those that did not meet the criteria specified by the user (name is the name specified in -oFile + ".png")
* File with historgram of ortholog lengths of coherent orthologs (name is the name specified in -oFile + "-peak.png")
	* Note: To obtain ortholog length histograms when running orthologFind.py on a cluster, submit the job (or open the interactive session in which the program will be run) using the --x11 option.


# Preparing Histone Modification Data for HALPER
There are many reasons that starting with the summits is sub-optimal for histone modifcation data.  Unlike for TF ChIP-seq and open chromatin data, where for which the motifs are known to be clustered around motif summits, TFs are thought not to bind where there are large numbers of reads in histone modification datas but in the valleys between the regions with large numbers of reads.  In addition, the summit locations produced by MACS2, a commonly used peak caller for histone modification data, are thought to be unreliable.  A reasonable place to start with histone modification data, therefore, is the location within the region that has the largest number of species in the alignment, as this is likely to be an important part of the region.  If there are multiple such locations, which often happens, then choosing the one closest to the center makes sense because the centers of the histone modification regions tend to be more important than their edges.  Here is how to make an -sFile that contains these locations:

1.  Get the alignment depth for your species of interest:
```
halAlignmentDepth --outWiggle [alignmentDepthFileName] [cactusFileName] [speciesName]
```
This can require up to 8 gigabytes.  This can take a long time (days), and these files take up a few gigabytes, so, if you are in the Pfenning Lab, find out if someone else in the lab has already done this for your species of interest before doing this.

2.  Convert the alignment depth file from a wig file to a bigwigh file:
```
wigToBigWig [alignmentDepthFileName] [chromSizesFileName] [alignmentDepthBigwigFileName]
```
This can require up to 64 gigabytes.  Note that the chromosome naming conventions might be different from those in the chrom sizes file name.  Some people in the Pfenning Lab have experience converting chromosome names, so ask around for advice if you find that this is an issue.

3.  Convert the alignment depth bigwig file to a bedgraph file:
```
bigWigToBedGraph [alignmentDepthBigwigFileName] [alignmentDepthBedgraphFileName]
```

4.  Sort the bedgraph file by chromosome, start, end:
```
sort -k1,1 -k2,2n -k3,3n [alignmentDepthBedgraphFileName] > [sortedAlignmentDepthBedgraphFileName]
```
You can gzip the bedgraph files so that they do not take up too much space.

5.  Get the file that will be used for starting the ortholog extension for each region using the scores in the bedgraph file:
```
python getMaxScorePositionFromBedgraph.py --bedFileName [file with regions you will be getting scores for, will be -qFile for next step] --bedgraphFileName [sortedAlignmentDepthBedgraphFileName] --highestScoreLocationFileName [where the positions with the highest scores will be recored, you can map this with hal-liftover to create -sFile for the next step] --gz
```
This program requires the bed file and the bedgraph file to be sorted and not contain duplicated entires.  You should leave out --gz if the file with the regions and the alignment depth bedgraph file are not gzipped.  Note that this program is compatible with both python version 2 and python version 3 while orthologFind.py is compatible with only python verison 3.

Alternatively, you can replace steps 2-5 with the following script that combines them:
```
python getMaxScorePositionFromWig.py --bedFileName [file with regions you will be getting scores for, will be -qFile for next step] --wigFileName [alignmentDepthFileName] --chromSizesFileName [chromSizesFileName] --highestScoreLocationFileName [where the positions with the highest scores will be recored, you can map this with hal-liftover to create -sFile for the next step] --gz
```
This program requires the bed file to be sorted and not contain duplicated entires.  You should leave out --gz if the file with the regions is not gzipped.  This program is compatible with both python version 2 and python version 3.  Note that this script runs UCSC tools internally that sometimes fail silently; if you are running it, you should check the sorted bedgraph file when it finishes and re-run it with more memory alloted if that file is not large.

6.  Use hal-liftover to map the positions where the highest scores are recorded to the target species.  This will create your -sFile for the program.

