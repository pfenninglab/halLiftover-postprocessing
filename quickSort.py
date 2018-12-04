import argparse
import sys
import orthologFindHelper
import subprocess

'''
partition the list based on pivot, the first element
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