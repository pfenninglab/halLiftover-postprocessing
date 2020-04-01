from orthologFindHelper import *

'''
mergesort a list of (peak_start,peak_end,chr_name)
for a given peak
first by chr_name, and then by peak_start 
'''
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
'''
Check if all the segments (peak_start,peak_end,chr_name)
in a list for a peak are sorted
'''
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
'''
collecting a list of peaks that are not sorted
'''
def check_qFile_sorted(qdict):
	not_valid=[]
	for key,value in qdict.items():
		if(not sortedSeg(value)):
			not_valid.append(key)
	return not_valid