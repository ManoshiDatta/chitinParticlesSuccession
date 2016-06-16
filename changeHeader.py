#!/usr/bin/env python

import sys
import pandas as pd
import argparse
from Bio import SeqIO
import re

def readBarcodes(filePath):
    
    barcodeDF = pd.io.parsers.read_csv(filePath, delimiter='\t')
    sampleBarcodes = barcodeDF['Barcode'].tolist()
    sampleIDs = barcodeDF['SampleID'].tolist()
    
    sampleDict = dict(zip(sampleIDs, sampleBarcodes))
    
    return sampleDict

def writeFastaNewHeader(inPath, outPath, barcodeDict):
    
    seqDict = SeqIO.to_dict(SeqIO.parse(inPath, 'fasta'))
    
    with open(outPath, 'a') as f:
        for sampleID, barcode in barcodeDict.iteritems():
            print 'SampleID is ' + str(sampleID) + ', barcode is ' + barcode
            records = [r for r in seqDict.keys() if r.rfind(barcode) != -1]
             
            for ind, r in enumerate(records):
                seq = str(seqDict[r].seq)
             
                header = '>' + str(sampleID) + "_" +str(ind) + ';barcodelabel=' + barcode + '\n'
                str2write = header + seq + '\n'
                
                f.write(header + seq + '\n')
                   
    f.close()
        

argp = argparse.ArgumentParser( prog = "changeHeader.py",
    description = "Changes sequence headers from Illumina header to SeqID + barcodelabel.")
argp.add_argument( "-b",        dest = "barcodePath",
    type = str, help = "Path to tab-delimited file containing barcodes (Column 1: Sample ID, Column 2: Barcode)" )
argp.add_argument( "-i",        dest = "inputPath",
    type = str, help = "Path to input .fastq file" )
argp.add_argument( "-o",        dest = "outputPath",
    type = str, help = "Path to output .fastq file" )

def _main():
    
    args = argp.parse_args( )
    barcodePath, inPath, outPath = [args.barcodePath, args.inputPath, args.outputPath]

    barcodeDict = readBarcodes(barcodePath)
    
    writeFastaNewHeader(inPath, outPath, barcodeDict)

if __name__ == "__main__":
    
    _main()