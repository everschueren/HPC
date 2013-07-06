#!/usr/bin/env python
# encoding: utf-8
"""
ParsePprospector.py

Created by erik verschueren on 2013-03-13.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

import sys
import getopt

help_message = '''
The help message goes here.
'''

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


# split the description in uniprot_id+description
def parseProspector(fileName, outputFileName):
    with open(fileName) as f:
        outputFile = open(outputFileName,'w')
        header = next(f)
        header = header.split('\t')
        header[len(header)-1] = "uniprot_id"
        header.append("description\n") 
        header = "\t".join(header)
        outputFile.write(header)
        for line in f:
            lineSplit = line.split('\t')
            lastEntry = lineSplit.pop(len(lineSplit)-1)
            lastEntrySplit = lastEntry.split()
            uniprot_id = lastEntrySplit.pop(0)
            description = " ".join(lastEntrySplit)
            outputFile.write("\t".join(lineSplit)+"\t"+uniprot_id+"\t"+description+"\n")
        outputFile.close()
    f.close()

def main(argv=None):
    
    verbose=False
    output = ""
    fileName = ""
    
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "hofs:v", ["help", "output=", "file="])
        except getopt.error, msg:
            raise Usage(msg)
    
        # option processing
        for option, value in opts:
            if option == "-v":
                verbose = True
            if option in ("-h", "--help"):
                raise Usage(help_message)
            if option in ("-o", "--output"):
                output = value
            if option in ("-f", "--file"):
                fileName = value
         
    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        print >> sys.stderr, "\t for help use --help"
        return 2
    
    ## call functions here
    parseProspector(fileName, output)
    
if __name__ == "__main__":
    sys.exit(main())
