#!/usr/bin/env python

import argparse
from Bio import SeqIO
import re
import bisect

def findIDs(fwdPath, fwdBarcode):
    
    id = []
    handle = open(fwdPath, 'rU')
    
    for num, record in enumerate(SeqIO.parse(handle, 'fastq')):
        seq = str(record.seq)
        barcode = seq[0:5]
        
        if barcode == fwdBarcode:
            print num, record.id
            id.append(record.id[0:-2])  # Last two characters are /[1 or 2]

    return id

def writeSeqs(ids, inPath, outPath):
    
    header = ''
    data = ''
    
    inFile = open(inPath, 'r')
    outFile = open(outPath, 'a')
    
    # details
    ids = sorted(ids)
    
    i = 0
    
    for line in inFile:
        if line[0] == "@":
            if header != '':
                i+= 1
                
                print i, header
                barcode = header
                if isInSet(ids, barcode):
                    outFile.write(header + "\n" + data)
    
            header = line.strip()
            data = ''
        else:
            data += line
    
    barcode = header
    if barcode in ids:
        outFile.write(header + "\n" + data)
        
def isInSet(sortedList, seq):
    a = bisect.bisect_left(sortedList, seq)
    
    return (a < len(sortedList) and seq == sortedList[a])

argp = argparse.ArgumentParser( prog = "demultiplex_by_fwd_barcode.py")
argp.add_argument("-f", dest = "fwdPath", metavar = "fwdPath", 
                  type = str, help = "Path to forward reads")
argp.add_argument("-r", dest = "revPath", metavar = "revPath", 
                  type = str, help = "Path to reverse reads")
argp.add_argument("-fo", dest = "foutPath", metavar = "foutPath", 
                  type = str, help = "Output path for demultiplexed forward reads")
argp.add_argument("-ro", dest = "routPath", metavar = "routPath", 
                  type = str, help = "Output path for demultiplexed reverse reads")
argp.add_argument("-b", dest = "fwdBarcode", metavar = "fwdBarcode", 
                  type = str, help = "Forward barcode to look for")

def _main( ):
    args = argp.parse_args( )
    [fwdPath, revPath, foutPath, routPath, fwdBarcode] = [args.fwdPath, args.revPath, args.foutPath, args.routPath, args.fwdBarcode]
    
    allIDs = findIDs(fwdPath, fwdBarcode)
            
    fwdIDs = ['@' + item + '/1' for item in allIDs]
    revIDs = ['@' + item + '/2' for item in allIDs]
 
    writeSeqs(fwdIDs, fwdPath, foutPath)
    writeSeqs(revIDs, revPath, routPath)
    
if __name__ == "__main__":
    _main( )