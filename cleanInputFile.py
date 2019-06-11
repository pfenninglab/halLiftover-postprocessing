from orthologFindHelper import *


'''
* f is a file handler of the input file
* f2 is a file handler of the output file
For each line in f, take the first 3 columns, 
peakName on the 4th column is assigned by the user.
Each row in the output file will have 5 fields.
chr_name peak_start peak_end peak_length peak_name 
'''
def preprocess_qFile(qFileH, outf):
    qFileH.seek(0)
    outH=open(outf,"w+")
    for line in qFileH:
        strList=line.strip().split("\t")
        peakName = strList[3]
        peak_s=int(strList[1])
        peak_e=int(strList[2])
        #
        newLineList = strList[0:3]
        newLineList.append(str(peak_e-peak_s+1))
        newLineList.append(peakName)
        newLine = fromStringListToStr(newLineList)
        outH.write(newLine)
    outH.close()

'''
* f is a file handler of the input file
* f2 is a file handler of the output file
For each line in f, take the first 3 columns, assign
a peakname in the form peakx, where x is some integer
enumerated from 0.
'''
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

'''
* f is a file handler of the input file
* f2 is a file handler of the output file
* slen, int, is how many bases to the left and right 
of the summit do you want to take to be in the output
    for summit, slen=1, but peak_offset is 0
    we need this because for cactus, if we have
    coordinate range (x,y), y is not included.
* summit, bool, indicates whether we are taking the
summit of each peak or not
'''
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

''' 
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

