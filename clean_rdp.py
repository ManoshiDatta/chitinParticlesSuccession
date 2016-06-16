#! /usr/bin/env python

import sys
import csv
import argparse

argp = argparse.ArgumentParser( prog = "clean_rdp.py",
    description = "Removes RDP classifications with confidence < 80%.")
argp.add_argument( "-i",        dest = "inputPath",
    type = str, help = "Path to input .fastq file" )
argp.add_argument( "-t",        dest = "threshold",
    type = float, help = "Threshold for confidence of classification")

def _main():
    
    args = argp.parse_args( )
    inPath, thres = [args.inputPath, args.threshold]
    
    fileIn = open(inPath, 'rU')
    
    for astrLine in csv.reader(fileIn, delimiter='\t'):
        astrTaxon = []
        adTaxon = []
         
        for i in range( 2, len( astrLine ), 3 ):
            astrTaxon.append( astrLine[i] )
            
            try:
                adTaxon.append( float(astrLine[i + 2]) )
            except:
                continue
         
        fUnclassified = False
        
        while adTaxon[-1] < args.threshold:
             
            fUnclassified = True
            adTaxon.pop( )
            astrTaxon.pop( )
        
        if fUnclassified:
            astrTaxon.append( "unclassified" )
        
        print( "\t".join( [astrLine[0], "|".join( astrTaxon[1:len(adTaxon)+1] )] ) )

if __name__ == "__main__":
    
    _main()
