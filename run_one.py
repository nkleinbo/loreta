#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 10:17:35 2019

@author: nkleinbo
"""

import getopt, sys
import os
import analyse_insertion as ai
import re

fullCmdArguments = sys.argv
argumentList = fullCmdArguments[1:]

unixOptions = "hf:t:o:w:a:g:"
gnuOptions = ["help", "fastqpath=", "tdnafile=", "outputpath=", "webpath=", "allfasta="]

#this file is not working yet!

try:
    arguments, values = getopt.getopt(argumentList, unixOptions, gnuOptions)
except getopt.error as err:
    print (str(err))
    sys.exit(2)

genomepath = None
    
for currentArgument, currentValue in arguments:
    if currentArgument in ("-h", "--help"):
        print ("displaying help")
    elif currentArgument in ("-f", "--fastqpath"):
        print (("Location of fastq file: (%s)") % (currentValue))
        fastqpath = currentValue
    elif currentArgument in ("-t", "--tdna"):
        print (("Location of tdnafile: (%s)") % (currentValue))
        tdnafile = currentValue
    elif currentArgument in ("-o", "--outputpath"):
        print (("Results will be stored in: (%s)") % (currentValue))
        outputpath = currentValue
    elif currentArgument in ("-w", "--webpath"):
        print (("Webpages will be stored in: (%s)") % (currentValue))
        webpath = currentValue
    elif currentArgument in ("-a", "--allfasta"):
        print (("Location of fasta containing all sequences of interest: (%s)") % (currentValue))
        allfasta = currentValue
    elif currentArgument in ("-g", "--genome_path"):
        print (("Location of fasta containing assemblies (instead of fastq reads): (%s)") % (currentValue))
        genomepath = currentValue
        


assembly_extension = ""
if (genomepath is not None):
    input_dict = {}
    assembly_extension = "_assembly"
    for f in os.listdir(genomepath):
        prefix, ext = os.path.splitext(f)
        if not ext == ".fasta":
                continue
        input_dict[prefix] = os.path.join(genomepath, f)
else:
    input_dict = {}
    for f in os.listdir(fastqpath):
        prefix, ext = os.path.splitext(f)
        input_dict[prefix] = os.path.join(fastqpath, f)
        
fastq_dict = {}
for f in os.listdir(fastqpath):
    prefix, ext = os.path.splitext(f)
    fastq_dict[prefix] = os.path.join(fastqpath, f)
    
out_dict = {}
for line in input_dict:
    outdir = os.path.join(outputpath, line+assembly_extension)
    out_dict[line] = outdir

web_dict = {}
for line in input_dict:
    webdir = os.path.join(webpath, line+assembly_extension)
    web_dict[line] = webdir
    
#print (fastq_dict)
#print(out_dict)

#run analysis for each line:
for line in input_dict:
    #_thread.start_new_thread(ai.analyse_insertion(line, fastq_dict[line], tdnafile, out_dict[line]))
    if(genomepath is None):
        (contigfile, references_file, references_file_no_tdna, idfile) = ai.analyse_insertion(line, fastq_dict[line], tdnafile, allfasta, out_dict[line], web_dict[line])
        if(os.path.isfile(contigfile)):
            (contigfile2, references_file2, references_file_no_tdna2, idfile2) = ai.analyse_insertion(line, input_dict[line], contigfile, allfasta, out_dict[line]+"_with_flanking", web_dict[line]+"_with_flanking")
            (contigfile3, references_file3, references_file_no_tdna3, idfile3) = ai.analyse_insertion(line, input_dict[line], references_file, allfasta, out_dict[line]+"_with_flanking_from_ref", web_dict[line]+"_with_flanking_from_ref")
            (contigfile4, references_file4, references_file_no_tdna4, idfile4) = ai.analyse_insertion(line, input_dict[line], references_file_no_tdna, allfasta, out_dict[line]+"_with_flanking_no_tdna", web_dict[line]+"_with_flanking_no_tdna", filter_reads=idfile)
    else:
        (contigfile, references_file, references_file_no_tdna) = ai.analyse_insertion_assembly(line, fastq_dict[line], tdnafile, allfasta, out_dict[line], web_dict[line], assembly=input_dict[line])
    #break
    #contigfile2 = ai.analyse_insertion(line, fastq_dict[line], contigfile, allfasta, out_dict[line]+"_2", web_dict[line]+"_2")
    

#for line in fastq_dict:
#    ai.summarise_results(out_dict[line])
