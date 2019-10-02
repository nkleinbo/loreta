#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 14:45:04 2019

@author: nkleinbo
"""

import getopt, sys
import os

#move to seperate module later:
#wether to be verbose or not:
BE_VERBOSE = True;
#no of CPUs to use at max:
NO_CPUS = 64;
#BLAST parameters:
BLAST_PARAMS = " -perc_identity 90"

#runs BLAST and returns an ID file with all hits
def run_blast(fastqfile, blastfile, tdnafile):
    #check if makeblastdb has been run
    if not (os.path.isfile(tdnafile+".nin")):
        os.system("makeblastdb -dbtype nucl -in "+tdnafile)
    #run blast if results do not already exist
    if not (os.path.isfile(blastfile)):
        if BE_VERBOSE: print ("BLASTING "+fastqfile+" ...")
        #blastcommand = "gzip -dc "+fastqfile+" | seqtk seq -A | blastn -db "+tdnafile+" -out "+blastfile+" -outfmt 6 -num_threads "+str(NO_CPUS)+" "+BLAST_PARAMS
        blastcommand = "blastn -query "+fastqfile+" -db "+tdnafile+" -out "+blastfile+" -outfmt 6 -num_threads "+str(NO_CPUS)+" "+BLAST_PARAMS
        print(blastcommand)
        os.system(blastcommand)
        if BE_VERBOSE:  "Done BLASTING "+fastqfile+"."
    else:
        if BE_VERBOSE: print ("Skipping BLAST for "+fastqfile+".");
    if BE_VERBOSE: print ("Creating idfile for "+blastfile+"...")
    idfile = blastfile+".ids"
    os.system("cut -f1 "+blastfile+" | sort | uniq > "+idfile)
    if BE_VERBOSE:  "Done IDFILE "+fastqfile+"."
    return idfile

def filter_sequences(idfile, fastqfile, filtered_fastq):
    if not (os.path.isfile(filtered_fastq)):
        if BE_VERBOSE: print ("Grepping sequences from "+fastqfile+"...")
        os.system("grep -A 3 --no-group-separator -f "+idfile+" "+fastqfile+" > "+filtered_fastq)
        if BE_VERBOSE: print ("Done grepping sequences from "+fastqfile+".")
    else:
        if BE_VERBOSE: print ("Skipping grepping sequences from "+fastqfile+".")

def summarise_results(outdir):
    print ("Summary "+outdir)


def analyse_insertion (line, fastqfile, tdnafile, outdir):
    #run BLAST
    blastfile = os.path.join(outdir, line+".bls");
    idfile = run_blast(fastqfile, blastfile, tdnafile)
    #filter FASTQ
    filtered_fastq = os.path.join(outdir, line+".fastq")
    filter_sequences(idfile, fastqfile, filtered_fastq)
     




#print "Starting asssembly..."
#assemblydir=assemblyprefix+"_canu"
#canu_command="canu -d "+assemblydir+" -p "+lineid+" -nanopore-raw "+outfile+ " 'genomeSize=100k' 'useGrid=1' 'gridEngineThreadsOption=-pe multislot THREADS' 'gridEngineMemoryOption=-l m_mem_total=MEMORY' 'gridEngineSubmitCommand=qsub -P denbi -l idle=1'"
#os.system(canu_command)
#print "Done."



#qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore