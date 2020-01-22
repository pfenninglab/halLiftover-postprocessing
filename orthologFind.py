import argparse
import sys
import matplotlib.pyplot as plt
import numpy as np
import subprocess
import os

from orthologFindHelper import *
from quickSort import quicksort
from tupleMergeSort import *


'''
organizing mapped peak file into dictionary
key: peakname
value: [(peak_start,peak_end,chr_name),..., (peak_start,peak_end,chr_name)]
	sorted list
'''

#chr8	55610267	55610365	Peak_12452	1000	.	55610267
def create_tFile_dict(tFileH): 
	tFileH.seek(0)
	tFile_segDict={} 
	for line in tFileH:
		strList=line.split("\t")
		t_chrName=strList[0]
		t_segStart=int(strList[1])
		t_segEnd=int(strList[2])
		t_segName=strList[3].strip()
		t_segName_list = tFile_segDict.get(t_segName,[])
		t_segName_list.append((t_segStart,t_segEnd,t_chrName))
		tFile_segDict[t_segName]=t_segName_list
	for key, value in tFile_segDict.items():
		merge_sort(value, cmp_tuple)
		if(not sortedSeg(value)):
			print("Fatal Error: list not sorted for "+key)
			return {}
	return tFile_segDict



def num_segments_hist(dict_segtFile):
	numFragmentsDict={}
	plt.figure(1)
	hist_len = []
	for key, value in dict_segtFile.items():
		length = len(value)
		numFragmentsDict[key] = length
		hist_len.append(length)
	binwidth=np.linspace(0, 175, num=20)
	n, bins, patches =plt.hist(hist_len, edgecolor='black',bins=binwidth)
	title="Number of Segments of Mapped Peaks"
	plt.title(title)
	plt.xlabel('Number of Fragments')
	plt.ylabel('Count')
	plt.savefig("num_frags_mapped_peaks.png")
	plt.close()
	return n
	


''' if a summit maps multiple places, see if all segments (peak_start,peak_end,chr_name)
are adjacent  '''
def adj_pos(arr):
    n = len(arr)
    quicksort(arr, 0, n - 1)
    for i in range(1, n):
        if(i == 0):
            continue
        else:
            if(not(arr[i][0] == arr[i - 1][1] or arr[i][0] == arr[i - 1][1] + 1)):
                return False


'''
As we go through line by line
2 dicts
	- unique summit mapping
		key: peak_name
		value: (mapped_s,mapped_e,chr_name)
	- multiple summit mapping
		key: peak_name
		value: [(mapped_s,mapped_e,chr_name)...]
'''
def create_SFile_dict(FileH,mult_keepone):
	FileH.seek(0)
	peak_summit = {}
	multpeak_dict = {}
	# accumulating values
	multpeak_pos_list = []
	num_multpeak = 0
	num_multpeak_nonad = 0
	# specially dealing with first line
	first_ln_list = (FileH.readline()).split("\t")
	last_peak_name = first_ln_list[3].strip()
	last_chrstart = int(first_ln_list[1])
	last_chrend = int(first_ln_list[2])
	last_chrname = first_ln_list[0]
	firstline=True
	#
	FileH.seek(0)
	for line in FileH:
		strList = line.split("\t")
		chr_name = strList[0]
		mapped_s = int(strList[1])
		mapped_e = int(strList[2])
		peak_name = strList[3].strip()
		if peak_name != last_peak_name:
			if(multpeak_pos_list != []):
				if(not adj_pos(multpeak_pos_list)):
				    num_multpeak_nonad += 1
				multpeak_dict[last_peak_name] = multpeak_pos_list
				multpeak_pos_list = []
				# Delete the peak whose summit maps to multiple places
				if not mult_keepone: 
					peak_summit.pop(last_peak_name, None)
			# At a new peak, so add its mapped summit to the list
			peak_summit[peak_name] = (mapped_s, mapped_e, chr_name)
		else:
			if firstline:
				# Not at the occurrence of a repeated peak because at the first line
                                peak_summit[last_peak_name] = (last_chrstart, last_chrend, last_chrname)
                                firstline = False
			else:
				if(multpeak_pos_list == []):
					num_multpeak += 1
					multpeak_pos_list.append((last_chrstart, last_chrend, last_chrname))
				multpeak_pos_list.append((mapped_s, mapped_e, chr_name))
		last_peak_name = peak_name
		last_chrstart = mapped_s
		last_chrend = mapped_e
		last_chrname = chr_name
	# Repeat it one more time in case last line is a peak whose summit maps to multiple places
	if(multpeak_pos_list != []):
		if(not adj_pos(multpeak_pos_list)):
			num_multpeak_nonad += 1
		multpeak_dict[last_peak_name] = multpeak_pos_list
		# Delete the peak whose summit maps to multiple places
		if not mult_keepone: 
			peak_summit.pop(last_peak_name, None)
	return (peak_summit, multpeak_dict)


'''
The sorted list of (peak_start,peak_end,chr_name) for each peak
might not be contiguous, so we fill in the gaps between adjacent
segments on the same chromosome
'''
def process_search_seg(L):
	last_seg_s = L[0][0]
	last_seg_e = L[0][1]
	last_chrname = L[0][2]
	res = []
	#
	res.append(L[0])
	for seg in L[1:]:
		seg_s = seg[0]
		seg_e = seg[1]
		seg_chrname = seg[2]
		if(str_cmp(last_chrname,seg_chrname)==0):
			res.append((last_seg_e,seg_s,seg_chrname))
		res.append((seg_s,seg_e,seg_chrname))
		last_seg_s = seg_s
		last_seg_e = seg_e
		last_chrname = seg_chrname
	return res

'''
locate where the mapped-summit is in the sorted
list of (peak_start,peak_end,chr_name) of a given peak,
and extend left and right to incldue all segments on the
same chromosome
'''
def extend_summit(q_peak_list,summit_seg):
	q_peak_list_proc = process_search_seg(q_peak_list)
	n=len(q_peak_list_proc)
	#find in this arr of (summit_start,summit_end,summit_chrname) corresponding
	s_index = binsearch_summitseg(q_peak_list_proc,summit_seg,0,n-1)	
	if(s_index==-1): return()	
	############################separately deal with the (s,e,chrname) that includes the summit_seg
	summit_ortho_s=q_peak_list_proc[s_index][0]
	summit_ortho_e=q_peak_list_proc[s_index][1]
	summit_chrname=summit_seg[2]
	summit_s = q_peak_list_proc[s_index][0]
	summit_e = q_peak_list_proc[s_index][1]
	#
	l_index=s_index-1
	r_index=s_index+1
	#
	summit_q_pos = summit_seg[0] + (summit_seg[1] - summit_seg[0])//2
	l_len=0
	r_len= 0
	sum_len = summit_ortho_e - summit_ortho_s +1
	l_deadend=l_index<0
	r_deadend=r_index>=n
	###################
	while(not(l_deadend and r_deadend )): #as long as you can still extend to one side of the list
		if(not l_deadend):	
			l_seg_s=q_peak_list_proc[l_index][0]
			l_seg_chrname=q_peak_list_proc[l_index][2]
			if(l_seg_chrname == summit_chrname):
				l_len = summit_s - l_seg_s 
				summit_ortho_s = l_seg_s 	
		if(not r_deadend):
			r_seg_e=q_peak_list_proc[r_index][1]
			r_seg_chrname=q_peak_list_proc[r_index][2]
			if(r_seg_chrname == summit_chrname):
				r_len = r_seg_e - summit_e 
				summit_ortho_e = r_seg_e
		l_index-=1
		r_index+=1
		if(l_index<0):
			l_deadend=True
		if(r_index>=n):
			r_deadend=True
	sum_len += l_len + r_len 
	l_len += summit_q_pos - summit_s 
	r_len += summit_e - summit_q_pos 
	return(summit_ortho_s,summit_q_pos,summit_ortho_e,sum_len,l_len,r_len)

'''
test if a ortholog is valid against the user parameters
'''
def validOrtholog(summit_ortho_info,max_len,min_len,proct_dist, peak_name):
	#summit_ortho_info:
	## summit_ortho_s,summit_q_pos, summit_ortho_e,sum_len,l_len,r_len
	sum_len = summit_ortho_info[3]
	l_len = summit_ortho_info[4]
	r_len = summit_ortho_info[5]
	if(sum_len > max_len):
		# print("max_len is"+str(max_len))
		# print("peak "+str(peak_name)+" sum len is "+str(sum_len))
		return False
	if(sum_len < min_len):
		# print("%.4f min_len" % (this_min_len))
		# print("peak "+str(peak_name)+" sum len is "+str(sum_len))
		return False
	if(not(l_len >= proct_dist and r_len>=proct_dist)):
		# print("peak "+str(peak_name)+" l_len is "+str(l_len)+" r_len is "+str(r_len))
		return False
	return True

'''
make a histogram of all the valid orthologs
x label: ortholog length
y label: number of orthologs
'''
def make_hist(oFile,outname,bin_max,narrowPeak=False):
	oFileH = open(oFile,"r")
	plt.figure(1)
	hist_len = []
	peaks_len = []
	for line in oFileH:
		strList=line.split("\t")
		ortholog_len = int(strList[2]) - int(strList[1])
		hist_len.append(ortholog_len)
		if not narrowPeak:
			# The output file is not in narrowPeak format, so get the query species region length
			peak_len = int(strList[6])
			peaks_len.append(peak_len)
	binwidth=np.linspace(0, bin_max, num=20)
	fig=plt.hist(hist_len, edgecolor='black',bins=binwidth)
	title="Orthologs"
	plt.title(title)
	plt.xlabel('Length')
	plt.ylabel('Count')
	plt.savefig(outname+".png")
	plt.close()
	if not narrowPeak:
		# The output file is not in narrowPeak format, so plot the query species region lengths
		plt.figure(2)
		binwidth=np.linspace(0, bin_max, num=20)
		fig=plt.hist(peaks_len, edgecolor='black',bins=binwidth)
		title="Peaks"
		plt.title(title)
		plt.xlabel('Length')
		plt.ylabel('Count')
		plt.savefig(outname+"-peak.png")
		plt.close()
	oFileH.close()

def make_hist_peaks(oFile,outname,bin_max):
	oFileH = open(oFile,"r")
	plt.figure(1)
	hist_len = []
	for line in oFileH:
		strList=line.split("\t")
		peak_len = int(strList[3])
		hist_len.append(peak_len)
	binwidth=np.linspace(0, bin_max, num=20)
	fig=plt.hist(hist_len, edgecolor='black',bins=binwidth)
	title="All Peaks"
	plt.title(title)
	plt.xlabel('Length')
	plt.ylabel('Count')
	plt.savefig(outname+"-all-peaks.png")
	plt.close()
	oFileH.close()

'''
finding valid orthologs and then plot the histogram
'''
def ortholog_find(file_H,max_len,alen,min_len,blen,proct_dist,mult_keepone=False,narrowPeak=False):
	qFileH = open(file_H[0],"r+")
	tFileH = open(file_H[1],"r+")
	sFileH = open(file_H[2],"r+")
	oFileH = open(file_H[3],"w+")
	qFileH.seek(0) 
	qFile_fix_name=file_H[0]+".fixed"
	qFile_failed_name = file_H[3]+".failed"
	#
	qFile_FH = open(qFile_failed_name, "w+")
	dict_ortholog={}
	# 
	dict_segtFile = create_tFile_dict(tFileH)
	if(dict_segtFile=={}):
		print("Fatal Error")
		return 1
	dict_summit = create_SFile_dict(sFileH,mult_keepone)[0]
	#
	for line in qFileH: # qFileH's first four columns are chromosome name, region start, region end, region name
		# Iterate through the query peaks and construct coherent orthologs of their corresponding target peaks if possible
		strList=line.strip().split("\t")
		chr_name=strList[0]
		peak_s=int(strList[1])
		peak_e=int(strList[2])
		peak_len=peak_e-peak_s
		peak_name=strList[3]
		#if given fraction, calculate max_len 
		if(not alen):
			this_max_len = max_len*(peak_e-peak_s+1)
		else:
			this_max_len = max_len
		if(not blen):
			this_min_len = min_len * (float(peak_e-peak_s+1))
		else:
			this_min_len = min_len
		#key:peak_name, value:list of (s,e,chr_name) sorted wrt s
		q_peak_list = dict_segtFile.get(peak_name,[]) #q_segStart,q_segEnd,q_chrName
		summit_seg = dict_summit.get(peak_name,()) #mapped_summit_start, end, chr_name
		if(q_peak_list==[] or summit_seg==()):
			continue
		#
		q_extent=extend_summit(q_peak_list,summit_seg)
		if(q_extent == ()):
			continue
		# summit_ortho_s,summit_q_pos,summit_ortho_e,sum_len,l_len,r_len
		ortho_s=q_extent[0]
		ortho_e=q_extent[2]
		ortho_len = q_extent[3]
		summit_q_pos = q_extent[1]
		newLineList = []
		if not narrowPeak:
			# Format the output line in the original format for orthologFind.py
			newLineList = [summit_seg[2],str(ortho_s),str(ortho_e),str(summit_q_pos),peak_name,str(ortho_len)]
			newLineList.append(str(peak_len))
			newLineList.append(str(q_extent[-2]))
			newLineList.append(str(q_extent[-1]))
		else:
			# Format the output line in narrowPeak format
			newLineList = [summit_seg[2],str(ortho_s),str(ortho_e),peak_name,"-1",".","-1","-1","-1",str(q_extent[-2])]
		newLine = fromStringListToStr(newLineList)
		if(validOrtholog(q_extent,this_max_len,this_min_len,proct_dist,peak_name)):
			oFileH.write(newLine)
		else:
			qFile_FH.write(newLine)
	tFileH.close()
	qFileH.close()
	sFileH.close()
	oFileH.close()
	qFile_FH.close()
	make_hist(file_H[3],file_H[3],2500,narrowPeak=narrowPeak)
	return 0



def main(argv):
	parser = argparse.ArgumentParser(description='Ortholog Find')
	parser.add_argument('-max_len',
		help='maximum number of base pairs of the ortholog')
	
	parser.add_argument('-max_frac',
		help='maximum percentage of original peak of the ortholog')
	
	parser.add_argument('-protect_dist',help='summit protection distance',
	default=50)
	
	parser.add_argument('-min_len',
		help='minimum number of base pairs of the ortholog')
	
	parser.add_argument('-min_frac',
		help='minimum percentage of original peak of the ortholog')
	
	parser.add_argument('-qFile', help='input bed file', 
		required=True)
	
	parser.add_argument('-tFile', help='input mapped bed file',
		required=True)
	
	parser.add_argument('-sFile', help='input mapped-summit bed file',
		required=True)

	parser.add_argument('-oFile', help='out bed file name',
		required=True)
	parser.add_argument('-mult_keepone', action="store_true", \
		help='keep one position to use for a peak whose summit maps to multiple places; otherwise, such a peak is discarded', \
		required=False)
	parser.add_argument('-narrowPeak', action="store_true", \
		help='output file in narrowPeak format, string columns other than 1-4 and 10 will be ., number columns other than 1-4 and 10 will be -1',
                required=False)
	args = parser.parse_args()

	if(args.max_len is None and args.max_frac is None):
		print("Error: Must supply max_len or max_frac")
		exit(1)
	alen=True
	if(args.max_len is None):
		max_len=float(args.max_frac)
		alen=False
	else:
		max_len=int(args.max_len)
	#
	if(args.min_len is None and args.min_frac is None):
		print("Error: Must supply min_len or min_frac")
		exit(1)
	blen=True
	if(args.min_len is None):
		min_len=float(args.min_frac)
		blen=False
	else:
		min_len=int(args.min_len)
	file_H=[]
	file_H.append(args.qFile)
	file_H.append(args.tFile)
	file_H.append(args.sFile)
	file_H.append(args.oFile)
	if(not check_valid_files(args.tFile)):
		print("Error: tFile is empty")
		exit(1)
	if(not check_valid_files(args.qFile)):
		print("Error: qFile is empty")
		exit(1)
	if(not check_valid_files(args.sFile)):
		print("Error: sFile is empty")
		exit(1)
	ortholog_find(file_H,max_len,alen,min_len,blen,int(args.protect_dist),mult_keepone=args.mult_keepone,narrowPeak=args.narrowPeak);
		

	

if __name__ == "__main__":
   main(sys.argv[1:])
