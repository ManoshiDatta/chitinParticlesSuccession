#!/usr/bin/env python

import sys
import pandas as pd
import argparse
import csv

def readOTUtable(filePath):
    
    with open(filePath) as f:
        reader = csv.reader(f, delimiter='\t')
        d = list(reader)

    return d

def makeBarcodeDict(barcodePath):
    
    df = pd.io.parsers.read_csv(barcodePath, delimiter='\t', index_col=1, header=False)
    dfDict = df.to_dict()
    
    return dfDict['SampleID']
    
def changeHeader(headerList, barcodeDict):
    
    newHeader = []
    
    newHeader.append(headerList[0])
    headerList.pop(0)
    
    for elem in headerList:
        
        newHeader.append(barcodeDict[elem])
    
    return newHeader

def writeNewOTUTable(OTUList, outputPath):
    
    with open(outputPath, 'wb') as f:
        
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(OTUList)
        
    f.close()

argp = argparse.ArgumentParser( prog = "change_OTU_table_headers.py",
    description = "Changes OTU table headers from barcodes to sample identifiers")
argp.add_argument( "-b",        dest = "barPath",
    type = str, help = "Path to tab-delimited file with barcodes and sample identifiers")
argp.add_argument( "-i",        dest = "inputPath",
    type = str, help = "Path to tab-delimited OTU table" )
argp.add_argument( "-o",        dest = "outputPath",
    type = str, help = "Output file path" )

def _main():
    
    args = argp.parse_args( )
    barcodePath, inputPath, outputPath = [args.barPath, args.inputPath, args.outputPath]
    
    OTUtableList = readOTUtable(inputPath)
    barcodeDict = makeBarcodeDict(barcodePath)
    
    newHeader = changeHeader(OTUtableList[0], barcodeDict)
    
    OTUtableList[0] = newHeader
    
    writeNewOTUTable(OTUtableList, outputPath)

if __name__ == "__main__":
    
    _main()