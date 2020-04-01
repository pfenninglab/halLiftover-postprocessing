import subprocess
import os

def check_valid_files(file):
	if (os.path.getsize(file) <=0):
		return False 
	return True

''' comparison functions for 2 strings '''
def str_cmp(s1,s2):
	s1_len=len(s1)
	s2_len=len(s2)
	if(s1_len > s2_len):
		return 1
	else:
		if(s1_len < s2_len):
			return -1
		else:
			if(s1 > s2):
				return 1
			else:
				if(s1 < s2):
					return -1
				else:
					return 0

''' comparison functions for 2 tuples (chr_start,chr_end,chr_name) 
first sort by chromosome names and then by peak_start
'''
def cmp_tuple(t1,t2): 
	peak_start = t1[0]
	chr_name = t1[2]
	peak_start2 = t2[0]
	chr_name2 = t2[2]
	if(str_cmp(chr_name,chr_name2)==0):
		if peak_start < peak_start2:
			return -1
		elif peak_start > peak_start2:
			return 1
		else:
			return 0
	else:
		return str_cmp(chr_name,chr_name2)


''' Concatenate a list of strings into one string '''
def fromStringListToStr(strL):
    re = ""
    num = 0
    for s in strL:
        if num == 0:
            re = re + s
        else:
            re = re + "\t" + s
        num = 1
    return re + "\n"


''' t1 is the summit tuple, t2 is one of the other tuples 
	for a peak
	you want to see if (s_start,s_end) is in between
	the range (seg_start,seg_end) and both are on the 
	same chromosome
'''
def cmp_tuple_summit(t1,t2): #(seg_start,seg_end,chr_name)
	s_start = t1[0]
	s_end = t1[1]
	s_chr_name = t1[2]
	#
	seg_start = t2[0]
	seg_end = t2[1]
	chr_name = t2[2]
	#
	if str_cmp(s_chr_name,chr_name)==0:
		if (s_start >= seg_start) and (s_end <= seg_end):
			return 0
		elif s_start < seg_start:
			return -1
		else:
			return 1
	else:
		return str_cmp(s_chr_name,chr_name)

'''
search for the segment (peak_start,peak_end,chr_name)
in a list for a given peak that contains
the (mapped_summit_start,mapped_summit_end,mapped_summit_chr_name)
'''
def binsearch_summitseg(L,summit_seg,low,high):
	while(low<=high):
		mid=low+(high-low)//2
		t2 = L[mid]
		cmp_res = cmp_tuple_summit(summit_seg,t2)
		# print("low is"+str(low)+" high is "+str(high))
		# print("comparing "+str(summit_seg)+" and "+str(t2))
		# print("\t result is:"+str(cmp_res))
		if (cmp_res==0):
			return mid
		else:
			if (cmp_res==-1): #summit_seg<t2
				high=mid-1
			else:
				low=mid+1
	return -1

'''
collect a list of peaksnames from a file
'''
def find_all_peaknames(fileH): #assume last columns are peak names 
	fileH.seek(0)
	name_l = {}
	for line in fileH:
		strList=line.split("\t")
		peakname=strList[-1].strip()
		peak_times = name_l.get(peakname, 0)
		peak_times +=1
		name_l[peakname]  = peak_times
	peaknames = list(name_l.keys())
	merge_sort(peaknames, str_cmp)
	return peaknames

'''
find peaks that are not mapped by hal-liftover
'''
def find_notmapped_peaks_h(qFileName, total_peaks):
	fileH = open(qFileName,"r+")
	name_l = find_all_peaknames(fileH)
	out_text = (subprocess.check_output(['wc','-l',qFileName])).decode('utf-8')
	total_peaks = int(out_test.split()[0])
	not_mapped = []
	last_peaknum=int(name_l[0][4:])
	for peakname in name_l:
		peaknum = int(peakname[4:])
		if peaknum > last_peaknum+1:
			for not_mapped_num in range(last_peaknum+1,peaknum):
				not_mapped_name = "peak"+str(not_mapped_num)
				not_mapped.append(not_mapped_name)
		last_peaknum = peaknum
	for peak_num in range(last_peaknum+1,total_peaks):
		not_mapped_name = "peak"+str(peak_num)
		not_mapped.append(not_mapped_name)
	return not_mapped