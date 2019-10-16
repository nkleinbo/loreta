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
        print (("Location of fastq files: (%s)") % (currentValue))
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
        
if(genomepath is None):
    input_dict = {}
    for f in os.listdir(fastqpath):
        prefix, ext = os.path.splitext(f)
        input_dict[prefix] = os.path.join(fastqpath, f)
else:
    input_dict = {}
    for f in os.listdir(genomepath):
        #if not (re.match(".*fasta", f)): continue
        prefix, ext = os.path.splitext(f)
        input_dict[prefix+"_assembly"] = os.path.join(genomepath, f)
    
out_dict = {}
for line in input_dict:
    outdir = os.path.join(outputpath, line)
    out_dict[line] = outdir

web_dict = {}
for line in input_dict:
    webdir = os.path.join(webpath, line)
    web_dict[line] = webdir
    
#print (fastq_dict)
#print(out_dict)

#run analysis for each line:
for line in input_dict:
    #_thread.start_new_thread(ai.analyse_insertion(line, fastq_dict[line], tdnafile, out_dict[line]))
    if(genomepath is None):
        (contigfile, references_file) = ai.analyse_insertion(line, input_dict[line], tdnafile, allfasta, out_dict[line], web_dict[line])
        (contigfile2, references_file2) = ai.analyse_insertion(line, input_dict[line], references_file, allfasta, out_dict[line]+"_with_flanking", web_dict[line]+"_with_flanking")
        
    else:
        (contigfile, references_file) = ai.analyse_insertion(line, input_dict[line], tdnafile, allfasta, out_dict[line], web_dict[line], assembly_input=True)
    break
    #contigfile2 = ai.analyse_insertion(line, fastq_dict[line], contigfile, allfasta, out_dict[line]+"_2", web_dict[line]+"_2")
    

#for line in fastq_dict:
#    ai.summarise_results(out_dict[line])
