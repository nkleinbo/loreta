#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 14:45:04 2019

@author: nkleinbo
"""

import os
import re


#move to seperate module later:
#wether to be verbose or not:
BE_VERBOSE = True;
#no of CPUs to use at max:
NO_CPUS = 64;
#BLAST parameters:
BLAST_PARAMS = " -perc_identity 85"
#expected genome size:
GENOME_SIZE=145000000


#runs BLAST and returns an ID file with all hits
def run_blast(fastqfile, blastfile, dbfile):
    #check if makeblastdb has been run
    if not (os.path.isfile(dbfile+".nin")):
        os.system("makeblastdb -dbtype nucl -in "+dbfile)
    #run blast if results do not already exist
    if not (os.path.isfile(blastfile)):
        if BE_VERBOSE: print ("BLASTING "+fastqfile+" ...")
        #blastcommand = "gzip -dc "+fastqfile+" | seqtk seq -A | blastn -db "+tdnafile+" -out "+blastfile+" -outfmt 6 -num_threads "+str(NO_CPUS)+" "+BLAST_PARAMS
        (prefix, extension) = os.path.splitext(fastqfile)
        blastcommand = "seqtk seq -A "+fastqfile+" |  parallel -S : --block 10M --recstart '>' --pipe blastn -query - -db "+dbfile+" -out "+blastfile+" -outfmt 6 "+BLAST_PARAMS
        if(extension == "fasta"):
            blastcommand = "cat "+fastqfile+" |  parallel -S : --block 10M --recstart '>' --pipe blastn -query - -db "+dbfile+" -out "+blastfile+" -outfmt 6 "+BLAST_PARAMS
        print(blastcommand)
        os.system(blastcommand)
        if BE_VERBOSE:  "Done BLASTING "+fastqfile+"."
    else:
        if BE_VERBOSE: print ("Skipping BLAST for "+fastqfile+".");
    if BE_VERBOSE: print ("Creating idfile for "+blastfile+"...")
    idfile = blastfile+".ids"
    os.system("cut -f1 "+blastfile+" | sort | uniq > "+idfile)
    if BE_VERBOSE: print("Done IDFILE "+fastqfile+".")
    return idfile

def run_canu(filtered_fastq, statistics, outdir):
    assemblydir = os.path.join(outdir, "assembly")
    projectname = "assembly"
    contigfile = os.path.join(assemblydir, projectname+".contigs.fasta")
    longest_read = statistics["filtered_fastq_stats"]["longest_read"]
    coverage_for_longest_read = statistics["filtered_fastq_stats"]["coverage_for_longest_read"]
    genome_size = 0
    if (coverage_for_longest_read > 10):
        genome_size=longest_read
    else:
        genome_size=int((longest_read * coverage_for_longest_read) / 11)
    if not (os.path.isfile(contigfile)):
        if BE_VERBOSE: print ("Starting asssembly for "+filtered_fastq+"...")
        canu_command="canu -d "+assemblydir+" -p "+projectname+" -nanopore-raw "+filtered_fastq+ " 'stopOnLowCoverage=5' 'genomeSize="+str(genome_size)+"' 'useGrid=0' 'corMhapFilterThreshold=0.0000000002' 'ovlMerThreshold=500' 'corMhapOptions=--threshold 0.80 --num-hashes 512 --num-min-matches 3 --ordered-sketch-size 1000 --ordered-kmer-size 14 --min-olap-length 2000 --repeat-idf-scale 50' 'correctedErrorRate=0.17'"
        os.system(canu_command)
        if BE_VERBOSE: print ("Done asssembly for "+filtered_fastq+".")
    else:
        if BE_VERBOSE: print ("Skipping asssembly for "+filtered_fastq+"...")
    return contigfile
#os.system(canu_command)
#print "Done."


def filter_sequences(idfile, fastqfile, filtered_fastq):
    if not (os.path.isfile(filtered_fastq)):
        if BE_VERBOSE: print ("Grepping sequences from "+fastqfile+"...")
        os.system("grep -A 3 --no-group-separator -f "+idfile+" "+fastqfile+" > "+filtered_fastq)
        if BE_VERBOSE: print ("Done grepping sequences from "+fastqfile+".")
    else:
        if BE_VERBOSE: print ("Skipping grepping sequences from "+fastqfile+".")
    os.system("wc -l "+fastqfile)

def get_fastq_stats(fastqfile):
    stats = {}
    number_reads = 0
    longest_read = 0
    total_bases = 0
    with open(fastqfile, 'r') as fh:
        i = 0
        for line in fh:
            i += 1
            if(i % 4 == 0):
                number_reads += 1
            elif(i % 2 == 0):
                bases = len(line)
                total_bases += bases
                if(bases > longest_read):
                    longest_read = bases
    mean_length = total_bases / number_reads
    coverage_for_longest_read = total_bases / longest_read
    coverage_for_genome = total_bases / GENOME_SIZE
    stats["number_reads"] = number_reads
    stats["longest_read"] = longest_read
    stats["total_bases"] = total_bases
    stats["mean_length"] = mean_length
    stats["coverage_for_longest_read"] = coverage_for_longest_read
    stats["coverage_for_genome"] = coverage_for_genome
    return stats
    

def summarise_results(outdir):
    print ("Summary "+outdir)

def html_summary(line, statistics, webdir):
    if not (os.path.isdir(webdir)):
        os.mkdir(webdir)    
    index_html = os.path.join(webdir, "index.html")
    fh = open(index_html, "w+")
    fh.write("<h1>Results for "+line+"</h1>")
    for head in statistics:
        fh.write("<h2>"+head+"</h2>")
        stats = statistics[head]
        fh.write("<table>")
        for subhead in stats:
            fh.write("<tr>")
            fh.write("<td><b>"+subhead+"</b></td>")
            fh.write("<td>"+str(stats[subhead])+"</b></td>")
            fh.write("</tr>")
        fh.write("</table>")
    fh.close()
    
    return index_html
    

def get_annotation_from_blast_result(contigfile, blastfile):
    contigs = {}
    with open(contigfile, 'r') as fh:
        for line in fh:
            test = re.match("^>(tig\d+)\slen=(\d+).*", line)
            if test is not None: 
                test2 = test.groups()
                test3 = test2[0]
                test4 = test2[1]
                print(test3+"/"+test4)
                


def analyse_insertion (line, fastqfile, tdnafile, allfasta, outdir, webdir):
    if not (os.path.isdir(outdir)):
        os.mkdir(outdir)
    statistics= {}
    #get statistics for all reads
    #fastq_stats = get_fastq_stats(fastqfile)
    #statistics["fastq_stats"] = fastq_stats
    #run BLAST
    blastfile = os.path.join(outdir, line+".bls");
    idfile = run_blast(fastqfile, blastfile, tdnafile)
    #filter FASTQ
    filtered_fastq = os.path.join(outdir, line+".fastq")
    filter_sequences(idfile, fastqfile, filtered_fastq)
    filtered_fastq_stats = get_fastq_stats(filtered_fastq)
    statistics["filtered_fastq_stats"] = filtered_fastq_stats
    #run canu:
    contigfile = run_canu(filtered_fastq, statistics, outdir)
    
    blast_vs_allfasta = os.path.join(outdir, line+"_assembly_vs_allfasta.bls");
    run_blast(contigfile, blast_vs_allfasta, allfasta)

    get_annotation_from_blast_result(contigfile, blast_vs_allfasta)
    
        
    html_summary(line, statistics, webdir)
     

#qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore