## halLiftover-postprocessing

# Running Program
* `python orthologFind.py` using python3


# Introduction
This tool is helpful if you have a bed file of open data, and you want to find orthologs of these peaks. Since this tool relies on halLiftover, the assembly of your peak data must be one of the available ones in the cactus file that you are using for halLiftover.

These are the species in the 12-way mammalian cactus file, which is the most common cactus file used by the Pfenning Lab:
* Tree_shrew tupChi1
* Kangaroo_rat dipOrd1
* Human hg38
* Chimp panTro6
* Rhesus rheMac8
* Mouse mm10
* Rat rn6
* Dog canFam3
* Cat felCat8
* Pig susScr11
* Cow bosTau8
* Horse equCab3

If your assembly is not in the cactus file that you are using for halLiftover, use UCSC-Liftover to lift the peaks over to the correct assembly.

# Dependencies
* Python version 3.7 (https://www.python.org/downloads/release/python-371/)
* Python libraries `matplotlib` and `numpy`
	* matplotlib (https://matplotlib.org/downloads.html)
	* numpy (http://www.numpy.org/)
# Tips for installing hal toolkits
* To install, follow the instructions in this website: https://github.com/ComparativeGenomicsToolkit/hal
	* For details about installation on the Lane cluster, follow the instructions in https://github.com/pfenninglab/halLiftover-postprocessing/blob/master/halliftoverInstallationSpecifics.txt
* On Lane cluster, load gcc 4.9 module 
* Have a correct ~/.bashrc is very helpful, here is an example:
```
export PATH=/home/xiaoyuz1/multalign/hal/hdf5-1.10.1/hdf5/bin:${PATH}
export h5prefix=-prefix=/home/xiaoyuz1/multalign/hal/hdf5
export sonLibRootPath=/home/xiaoyuz1/multalign/hal/sonLib
export PATH=/home/xiaoyuz1/multalign/hal/bin:${PATH}
export PYTHONPATH=/home/xiaoyuz1/multalign:${PYTHONPATH}
export PYTHONPATH=/home/xiaoyuz1/multalign/hal:${PYTHONPATH}

```

# Running halLiftover 
* Run halLiftover on tFile 
	* To halLiftover, an example sbatch job would be:
	```
	#!/bin/bash
	#SBATCH -e mappingHal.err
	#SBATCH -o mappingHal.out
	#SBATCH -t 48:00:00
	#SBATCH -n 1
	#SBATCH -p pfen1
	#SBATCH -c 2
	#SBATCH --mem=8G
	halLiftover /projects/pfenninggroup/machineLearningForComputationalBiology/alignCactus/10plusway-master.hal \
	Mouse \
	"/projects/pfenninggroup/machineLearningForComputationalBiology/alignCactus/Cortex_AgeB_ATAC/Cortex_AgeB_ATAC_out_ppr.IDR0.1.filt.narrowPeak.summit" \
	Human \
	"/projects/pfenninggroup/machineLearningForComputationalBiology/alignCactus/Cortex_AgeB_ATAC/Cortex_AgeB_ATAC_out_ppr.IDR0.1.filt.toHuman.narrowPeak.summit"
	```

 

# Program Parameters 
* -max_len: ortholog length must be less or equal to max_len
* -max_frac: ortholog length must be less or euqal to max_frac * peak length 
	* provide either max_len or max_frac
* -min_len: ortholog length must be greater or equal to min_len
* -min_frac: ortholog length must be greater or equal to min_frac * peak length 
	* provide either min_len or min_frac
* -protect_dist: the ortholog length in each direction from the ortholog of the summit must be at least proct_dist 
			![alt text](https://github.com/pfenninglab/multiple_alignment-python/blob/master/min_proct_dist.png)

* -qFile: the original bed file containing information on (at least): chromosome_name, start, end, peak name -- The 1st 4 columns **MUST** be in standard bed format. 
	

* -tFile: bed file of the halLiftover-ed result for each peak 
	* Line format must be: ` chr_name    peak_start    peak_end    peak_name `. halLiftover should output file conforming to this format. 
	* Last column must be of the format "peak[number]", for example, "peak0"
	* Examples:
		```
		chr8	55610267	55610335	peak0
		chr8	55610240	55610267	peak0
		chr8	55610220	55610240	peak0
		chr8	55610191	55610220	peak0
		chr8	55610183	55610190	peak0 
		```

* -sFile: bed file of the halLiftover-ed result for each peak summit
	* Line format must be: ` chr_name    peak_start    peak_end    peak_name`. halLiftover should output file conforming to this format. 
	* Last column must be of the format "peak[number]", for example, "peak0"
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
	* Examples:
		```
		chr8	55609305	55610335	55609835	peak0	1031	1019	530	500
		chr8	55609305	55610335	55609437	peak1	1031	1019	132	898
		```
	* Along with the output file, another file is generated, which has name oFile.failed
		* Line format is the same as oFile 
		* oFile.failed would contain all the orthologs that are not valid when judged against the parameters users have entered
		* Examples:
		```
		chr4	76114521	76116294	76116110	peak13	1774	544	1589	184
		chr4	76477814	76478909	76477950	peak49	1096	459	136	959
		chr4	76809080	76816154	76816028	peak66	7075	288	6948	126
		chr6	53443438	53447196	53443500	peak69	3759	750	62	3696
		```



# Example run of the program 
```
python orthologFind.py -max_len 1000 -min_len 50 -protect_dist 5 -tFile HumanMiddleFrontalGyrusDNase_ppr.IDR0.1.filt_brainUpEnhancerStrictShort_withPeakOffetsAndNames.bed -qFile HumanMiddleFrontalGyrusDNase_ppr.IDR0.1.filt_brainUpEnhancerStrictShort_mm10Hal.bed -sFile  HumanMiddleFrontalGyrusDNase_ppr.IDR0.1.filt_brainUpEnhancerStrictShort_summits_mm10Hal.bed -oFile HumanMiddleFrontalGyrusDNase_ppr.IDR0.1.filt_brainUpEnhancerStrictShort_mm10Hal_summitExtendedMin50Max1000Protect5.bed
```
Note how there is only one '-' (dash) for the parameter name. 

Note: When running the program on the Lane cluster, submit the job (or open the interactive session in which the program will be run) using the --x11 option.

