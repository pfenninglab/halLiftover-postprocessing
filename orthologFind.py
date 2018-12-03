import argparse
import sys
import matplotlib.pyplot as plt
import numpy as np
import subprocess

'''
def create_halLiftover_sh(slurm_params)
'''

def partition(arr,low,high):
    resultarr=arr
    i = ( low-1 )         # index of smaller element
    pivot = resultarr[high][0]     # pivot
 
    for j in range(low , high):
 
        # If current element is smaller than or
        # equal to pivot
        if resultarr[j][0] <= pivot:
         
            # increment index of smaller element
            i = i+1
            resultarr[i],resultarr[j] = resultarr[j],resultarr[i]
 
    resultarr[i+1],resultarr[high] = resultarr[high],resultarr[i+1]
    return (resultarr, i+1 )
 
# The main function that implements QuickSort
# arr[] --> Array to be sorted,
# low  --> Starting index,
# high  --> Ending index
 
# Function to do Quick sort
def quicksort(arr,low,high):
    if low < high:
 
        # pi is partitioning index, arr[p] is now
        # at right place
        pi = partition(arr,low,high)
 
        # Separately sort elements before
        # partition and after partition
        quicksort(pi[0], low, pi[1]-1)
        quicksort(pi[0], pi[1]+1, high)

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

def columnChrNameStartEnd(f, f2):
    f.seek(0)
    f2.seek(0)
    index = 0
    peakName = "peak0"
    for line in f:
        peakName = "peak" + str(index)
        strList = line.split("\t")
        newStrList = strList[0:3]
        newStrList.append(peakName)
        newLine = fromStringListToStr(newStrList)
        f2.write(newLine)
        index = index + 1
    f2.close()
    f.close()

def summitPlusMinusLength(f, f2, slen, summit):
    f.seek(0)
    f2.seek(0)
    slen = int(slen)
    index = 0
    peakName = "peak0"
    summit_offset = 0
    if(summit):
        summit_offset = 1
    for line in f:
        peakName = "peak" + str(index)
        strList = line.split("\t")
        peakStart = int(strList[1])
        summitDisFromStart = int(strList[9])
        summitStart = peakStart + summitDisFromStart - slen + summit_offset
        summitEnd = peakStart + summitDisFromStart + slen

        newLineList = [strList[0], str(summitStart), str(summitEnd), peakName]
        newLine = fromStringListToStr(newLineList)
        f2.write(newLine)
        index = index + 1
    f2.close()
    f.close()

''' preposessing of revtFile if needed
peakName0
peakName0
peakName0
--> 
peakName0_1
peakName0_2
peakName0_3
'''
def assignPeakNameSuffix(fname, fname2):
    f = open(fname, "r+")
    f2 = open(fname2, "x+")
    curPeakName = ""
    acc = 0
    for line in f:
        strList = line.split("\t")  # CAUTION: assume always be delimited by tab
        peakNamePrefix = strList[-1]
        if(peakNamePrefix != curPeakName):
            curPeakName = peakNamePrefix
            acc = 0
        strList[-1] = peakNamePrefix[:-1] + "_" + str(acc)
        acc+=1
        newLine = fromStringListToStr(strList)
        f2.write(newLine)
    f2.close()
    f.close()


def preprocess_tFile(tFileH, outf):
	tFileH.seek(0)
	outH=open(outf,"w+")
	index = 0
	peakName = "peak0"
	for line in tFileH:
		peakName = "peak" + str(index)
		strList=line.split("\t")
		peak_s=int(strList[1])
		peak_e=int(strList[2])
		#
		newLineList = strList[0:3]
		newLineList.append(str(peak_e-peak_s+1))
		newLineList.append(peakName)
		newLine = fromStringListToStr(newLineList)
		outH.write(newLine)
		index = index + 1
	outH.close()

'''
assume that is already preprocessed 
segName is unique e.g. chrName0_0 only appears once  
look-up dictionary of revtFile
	key: chrName
	value: a list: (segStart, segEnd, segName), sorted according to segStart
		**** ASSUME that there is no overlapping segment
'''
def cmp_tuple_b4(t1,t2): #(chr_start,chr_end,chr_name)
	seg_start = t1[0]
	chr_name = t1[2]
	seg_start2 = t2[0]
	chr_name2 = t2[2]
	if chr_name < chr_name2:
		return -1 #t1<t2
	elif chr_name > chr_name2:
		return 1 #t1>t2
	else:
		if seg_start < seg_start2:
			return -1
		elif seg_start > seg_start2:
			return 1
		else:
			return 0

def cmp_tuple(t1,t2): #(chr_start,chr_end,chr_name)
	seg_start = t1[0]
	chr_name = t1[2]
	seg_start2 = t2[0]
	chr_name2 = t2[2]
	if(str_cmp(chr_name,chr_name2)==0):
		if seg_start < seg_start2:
			return -1
		elif seg_start > seg_start2:
			return 1
		else:
			return 0
	else:
		return str_cmp(chr_name,chr_name2)



def merge_sort(arr, cmp_func): 
    if len(arr) >1: 
        mid = len(arr)//2 #Finding the mid of the array 
        L = arr[:mid] # Dividing the array elements  
        R = arr[mid:] # into 2 halves 
  
        merge_sort(L, cmp_func) # Sorting the first half 
        merge_sort(R, cmp_func) # Sorting the second half 
  
        i = j = k = 0
          
        # Copy data to temp arrays L[] and R[] 
        while i < len(L) and j < len(R): 
            # if L[i] < R[j]:
            if (cmp_func(L[i],R[j]) ==-1):
                arr[k] = L[i] 
                i+=1
            else: 
                arr[k] = R[j] 
                j+=1
            k+=1
          
        # Checking if any element was left 
        while i < len(L): 
            arr[k] = L[i] 
            i+=1
            k+=1
          
        while j < len(R): 
            arr[k] = R[j] 
            j+=1
            k+=1

def sortedSeg(L):
	last_s = L[0][0]
	last_chr_name = L[0][2]
	valid=True
	for (start,end,chr_name) in L:
		if((str_cmp(chr_name,last_chr_name)==-1)):
			print("last chr_name is greater than chr_name")
			print("\t Chr name is "+chr_name)
			print("\t last_chr_name "+last_chr_name)
			valid=False
		else:
			if(start < last_s and (str_cmp(chr_name,last_chr_name)==0)):
				print("Start is less than last_s")
				print("\t Chr name is "+chr_name)
				print("\t last_chr_name "+last_chr_name)
				print("\t Start is "+str(start))
				print("\t last_s is "+str(last_s))
				valid = False
		last_s = start
		last_chr_name = chr_name
	return valid

def check_qFile_sorted(qdict):
	not_valid=[]
	for key,value in qdict.items():
		if(not sortedSeg(value)):
			not_valid.append(key)
	return not_valid

'''
organizing query file into dictionary
'''
def create_qFile_dict(qFileH): 
	qFileH.seek(0)
	qFile_segDict={} #key: seg_name, value: seg_start,seg_end,chr_name
	for line in qFileH:
		strList=line.split("\t")
		q_chrName=strList[0]
		q_segStart=int(strList[1])
		q_segEnd=int(strList[2])
		q_segName=strList[3][:-1]
		# 
		q_segName_list = qFile_segDict.get(q_segName,[])
		q_segName_list.append((q_segStart,q_segEnd,q_chrName))
		qFile_segDict[q_segName]=q_segName_list
	for key, value in qFile_segDict.items():
		merge_sort(value, cmp_tuple)
		if(not sortedSeg(value)):
			print("Fatal Error: list not sorted for "+key)
			return {}
	return qFile_segDict


def find_all_peaknames(fileH): #assume last columns are peak names 
	fileH.seek(0)
	name_l = {}
	for line in fileH:
		strList=line.split("\t")
		peakname=strList[-1][:-1]
		peak_times = name_l.get(peakname, 0)
		peak_times +=1
		name_l[peakname]  = peak_times
	peaknames = list(name_l.keys())
	merge_sort(peaknames, str_cmp)
	return peaknames

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

def num_segments_hist(dict_segqFile):
	numFragmentsDict={}
	plt.figure(1)
	hist_len = []
	for key, value in dict_segqFile.items():
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
	



def adj_pos(arr):
    n = len(arr)
    quicksort(arr, 0, n - 1)
    #for i in range(n):
    #    print ('{0}'.format(arr[i]))
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
def create_SFile_dict(FileH):
	FileH.seek(0)
	peak_summit = {}
	multpeak_dict = {}
	# accumulating values
	multpeak_pos_list = []
	num_multpeak = 0
	num_multpeak_nonad = 0
	# specially dealing with first line
	first_ln_list = (FileH.readline()).split("\t")
	last_peak_name = first_ln_list[3][0:-1]
	last_chrstart = int(first_ln_list[1])
	last_chrend = int(first_ln_list[2])
	last_chrname = first_ln_list[0]
	firstline=True
	#
	for line in FileH:
		strList = line.split("\t")
		chr_name = strList[0]
		mapped_s = int(strList[1])
		mapped_e = int(strList[2])
		peak_name = strList[3][0:-1]
		if peak_name != last_peak_name:
			if(multpeak_pos_list != []):
				if(not adj_pos(multpeak_pos_list)):
				    num_multpeak_nonad += 1
				multpeak_dict[last_peak_name] = multpeak_pos_list
				multpeak_pos_list = []
			else:
				if firstline:
					peak_summit[last_peak_name] = (last_chrstart, last_chrend, last_chrname)
					firstline = False;
				peak_summit[peak_name] = (mapped_s, mapped_e, chr_name)
		else:
			if(multpeak_pos_list == []):
				num_multpeak += 1
				multpeak_pos_list.append((last_chrstart, last_chrend, last_chrname))
			multpeak_pos_list.append((mapped_s, mapped_e, chr_name))
		last_peak_name = peak_name
		last_chrstart = mapped_s
		last_chrend = mapped_e
		last_chrname = chr_name
	return (peak_summit, multpeak_dict)



# t1 is the summit tuple, t2 is the other tuples for the peak
def cmp_tuple_summit(t1,t2): #(seg_start,seg_end,chr_name)
	s_start = t1[0]
	s_end = t1[1]
	s_chr_name = t1[2]
	#
	seg_start2 = t2[0]
	seg_end2 = t2[1]
	chr_name2 = t2[2]
	#
	if str_cmp(s_chr_name,chr_name2)==0:
		if (s_start >= seg_start2) and (s_end <= seg_end2):
			return 0
		elif s_start < seg_start2:
			return -1
		else:
			return 1
	else:
		return str_cmp(s_chr_name,chr_name2)


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

def make_hist(oFile,outname,bin_max):
	oFileH = open(oFile,"r")
	plt.figure(1)
	hist_len = []
	peaks_len = []
	for line in oFileH:
		strList=line.split("\t")
		ortholog_len = int(strList[5])
		peak_len = int(strList[6])
		hist_len.append(ortholog_len)
		peaks_len.append(peak_len)
	binwidth=np.linspace(0, bin_max, num=20)
	fig=plt.hist(hist_len, edgecolor='black',bins=binwidth)
	title="Orthologs"
	plt.title(title)
	plt.xlabel('Length')
	plt.ylabel('Count')
	plt.savefig(outname+".png")
	plt.close()
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
def double_valid(oFile,fixedFile,max_len,alen,min_len,blen,proct_dist):
	oFileH = open(oFile,"r")
	fixedH = open(fixedFile,"r")
	for line in oFileH:
		strList=line.split("\t")
'''

def ortholog_find(file_H,max_len,alen,min_len,blen,proct_dist):
	tFileH = open(file_H[0],"r+")
	qFileH = open(file_H[1],"r+")
	sFileH = open(file_H[2],"r+")
	oFileH = open(file_H[3],"w+")
	tFileH.seek(0) #tFileH has 5 fields: chr_name, peak_s, peak_e, peak_summit_d, peak_name
	tFile_fix_name=file_H[0]+".fixed"
	tFile_failed_name = file_H[3]+".failed"
	#
	tFile_FH = open(tFile_failed_name, "w+")
	# chrname, start, end, length, peakname 
	preprocess_tFile(tFileH,tFile_fix_name)
	tFileH.close()
	tFileH = open(tFile_fix_name,"r+")
	dict_ortholog={}
	# 
	dict_segqFile = create_qFile_dict(qFileH)
	if(dict_segqFile=={}):
		print("Fatal Error")
		return 1
	dict_summit = create_SFile_dict(sFileH)[0]
	#
	# test_trial=1000
	for line in tFileH:
		# if(test_trial == 0):
		# 	break
		strList=line.split("\t")
		chr_name=strList[0]
		peak_s=int(strList[1])
		peak_e=int(strList[2])
		peak_len=int(strList[-2])
		peak_name = strList[-1][0:-1]
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
		q_peak_list = dict_segqFile.get(peak_name,[]) #q_segStart,q_segEnd,q_chrName
		summit_seg = dict_summit.get(peak_name,()) #mapped_summit_start, end, chr_name
		if(q_peak_list==[] or summit_seg==()):
			# test_trial = test_trial-1
			continue
		#
		q_extent=extend_summit(q_peak_list,summit_seg)
		if(q_extent == ()):
			# test_trial = test_trial-1
			continue
		# summit_ortho_s,summit_q_pos,summit_ortho_e,sum_len,l_len,r_len
		ortho_s=q_extent[0]
		ortho_e=q_extent[2]
		ortho_len = q_extent[3]
		summit_q_pos = q_extent[1]
		newLineList = [summit_seg[2],str(ortho_s),str(ortho_e),str(summit_q_pos),peak_name,str(ortho_len)]
		newLineList.append(str(peak_len))
		newLineList.append(str(q_extent[-2]))
		newLineList.append(str(q_extent[-1]))
		newLine = fromStringListToStr(newLineList)
		if(validOrtholog(q_extent,this_max_len,this_min_len,proct_dist,peak_name)):
			oFileH.write(newLine)
		else:
			tFile_FH.write(newLine)
		# test_trial = test_trial-1
	tFileH.close()
	qFileH.close()
	sFileH.close()
	oFileH.close()
	tFile_FH.close()
	make_hist(file_H[3],file_H[3],2500)
	return 0
	# return dict_ortholog



def main(argv):
	parser = argparse.ArgumentParser(description='Ortholog Find')
	parser.add_argument('--max_len',
		help='maximum number of base pairs of the ortholog')
	
	parser.add_argument('--max_frac',
		help='maximum percentage of original peak of the ortholog')
	
	parser.add_argument('--protect_dist',help='summit protection distance',
	default=50)
	
	parser.add_argument('--min_len',
		help='minimum number of base pairs of the ortholog')
	
	parser.add_argument('--min_frac',
		help='minimum percentage of original peak of the ortholog')
	
	parser.add_argument('-tFile', help='input bed file', 
		required=True)
	
	parser.add_argument('-qFile', help='input mapped bed file',
		required=True)
	
	parser.add_argument('-sFile', help='input mapped-summit bed file',
		required=True)

	parser.add_argument('-oFile', help='out bed file name',
		required=True)
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
	file_H.append(args.tFile)
	file_H.append(args.qFile)
	file_H.append(args.sFile)
	file_H.append(args.oFile)
	ortholog_find(file_H,max_len,alen,min_len,blen,int(args.protect_dist))

	'''if max_len is passed in as percentage of original peak, 
	then first calculate that number and then find ortholog'''

	

if __name__ == "__main__":
   main(sys.argv[1:])
