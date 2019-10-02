#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 14:45:04 2019

@author: nkleinbo
"""

import getopt, sys
import os

fullCmdArguments = sys.argv
argumentList = fullCmdArguments[1:]

unixOptions = "hf:t:o:b:a:l:"
gnuOptions = ["help", "fastq=", "tdna=", "output=", "blast=", "assembly=", "line="]

try:
    arguments, values = getopt.getopt(argumentList, unixOptions, gnuOptions)
except getopt.error as err:
    print (str(err))
    sys.exit(2)
    
for currentArgument, currentValue in arguments:
    if currentArgument in ("-h", "--help"):
        print ("displaying help")
    elif currentArgument in ("-f", "--fastq"):
        print (("Location of fastq file: (%s)") % (currentValue))
        fastqfile = currentValue
    elif currentArgument in ("-t", "--tdna"):
        print (("Location of tdnafile: (%s)") % (currentValue))
        tdnafile = currentValue
    elif currentArgument in ("-b", "--blast"):
        print (("Location of blastfile: (%s)") % (currentValue))
        blastfile = currentValue
    elif currentArgument in ("-o", "--output"):
        print (("Results will be stored in: (%s)") % (currentValue))
        outfile = currentValue
    elif currentArgument in ("-a", "--assembly"):
        print (("Assembly results will be stored in: (%s)") % (currentValue))
        assemblyprefix = currentValue
    elif currentArgument in ("-l", "--line"):
        print (("Lineid is: (%s)") % (currentValue))
        lineid = currentValue
        
     

if not (os.path.isfile(blastfile)):
    print "BLASTING..."
    blastcommand = 'gzip -dc '+fastqfile+' | seqtk seq -A | blastn -db '+tdnafile+ ' -out '+blastfile+' -outfmt 6 -num_threads 4 -perc_identity 90'
    os.system(blastcommand)
    print "Done."
else:
    print "Skipping BLAST.";

print "Creating idfile..."
idfile = blastfile+".ids"
os.system("cut -f1 "+blastfile+" | sort | uniq > "+idfile)
print "Done."

if not (os.path.isfile(outfile)):
    print "Grepping sequences..."
    os.system("zgrep -A 3 --no-group-separator -f "+idfile+" "+fastqfile+" > "+outfile)
    print "Done."
else:
    print "Skipping grepping sequences."

#print "Starting asssembly..."
#assemblydir=assemblyprefix+"_canu"
#canu_command="canu -d "+assemblydir+" -p "+lineid+" -nanopore-raw "+outfile+ " 'genomeSize=100k' 'useGrid=1' 'gridEngineThreadsOption=-pe multislot THREADS' 'gridEngineMemoryOption=-l m_mem_total=MEMORY' 'gridEngineSubmitCommand=qsub -P denbi -l idle=1'"
#os.system(canu_command)
#print "Done."



#qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore