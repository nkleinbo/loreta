#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 14:45:04 2019

@author: nkleinbo
"""

import os
import re
import visualisation as vis

#move to seperate module later:
#wether to be verbose or not:
BE_VERBOSE = True;
#no of CPUs to use at max:
NO_CPUS = 64;
#BLAST parameters:
BLAST_PARAMS = " -perc_identity 85"
#expected genome size:
GENOME_SIZE=145000000
#write fastq files WITHOUT T-DNA:
WRITE_FASTQ_WITHOUT_TDNA = False;
#allowed overlap of BLAST results for contig annotation:
ALLOWED_OVERLAP = 10
#allowed identity of BLAST results for contig annotation:
ALLOWED_IDENTITY = 90
#use corrected reads as input:
CORRECTED_READS = True


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
        if(extension == ".fasta"):
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
        if CORRECTED_READS:
            canu_command="canu trim-assemble -d "+assemblydir+" -p "+projectname+" -nanopore-corrected "+filtered_fastq+ " 'stopOnLowCoverage=5' 'genomeSize="+str(genome_size)+"' 'useGrid=0' 'corMhapFilterThreshold=0.0000000002' 'ovlMerThreshold=500' 'corMhapOptions=--threshold 0.80 --num-hashes 512 --num-min-matches 3 --ordered-sketch-size 1000 --ordered-kmer-size 14 --min-olap-length 2000 --repeat-idf-scale 50' 'correctedErrorRate=0.17'"
        else:
            canu_command="canu -d "+assemblydir+" -p "+projectname+" -nanopore-raw "+filtered_fastq+ " 'stopOnLowCoverage=5' 'genomeSize="+str(genome_size)+"' 'useGrid=0' 'corMhapFilterThreshold=0.0000000002' 'ovlMerThreshold=500' 'corMhapOptions=--threshold 0.80 --num-hashes 512 --num-min-matches 3 --ordered-sketch-size 1000 --ordered-kmer-size 14 --min-olap-length 2000 --repeat-idf-scale 50' 'correctedErrorRate=0.17'"
        os.system(canu_command)
        if BE_VERBOSE: print ("Done asssembly for "+filtered_fastq+".")
    else:
        if BE_VERBOSE: print ("Skipping asssembly for "+filtered_fastq+"...")
    return contigfile
#os.system(canu_command)
#print "Done."


def filter_sequences(idfile, fastqfile, filtered_fastq, write_remaining=False):
    if not (os.path.isfile(filtered_fastq)):
        if BE_VERBOSE: print ("Grepping sequences from "+fastqfile+"...")
        os.system("grep -A 3 --no-group-separator -f "+idfile+" "+fastqfile+" > "+filtered_fastq)
        if BE_VERBOSE: print ("Done grepping sequences from "+fastqfile+".")
    else:
        if BE_VERBOSE: print ("Skipping grepping sequences from "+fastqfile+".")
    if(write_remaining):
        prefix, ext = os.path.splitext(fastqfile)
        remaining_fastqfile = prefix+".notdna.fastq"
        remaining_idfile = prefix+".notdna.ids"
        if not (os.path.isfile(remaining_fastqfile)):
            if BE_VERBOSE: print ("Grepping non-T-DNA sequences from "+fastqfile+"...")
            os.system("grep '^@.*runid' "+fastqfile+" | grep -v -f "+idfile+" > "+remaining_idfile)
            os.system("grep -A 3 --no-group-separator -f "+remaining_idfile+" "+fastqfile+" > "+remaining_fastqfile) 
            if BE_VERBOSE: print ("Done grepping non-T-DNA sequences from "+fastqfile+".")
        else:
            if BE_VERBOSE: print ("Skipping grepping non-T-DNA sequences from "+fastqfile+".")
        
    
    


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

def html_summary(lineid, statistics, webdir):
    if not (os.path.isdir(webdir)):
        os.mkdir(webdir)    
    index_html = os.path.join(webdir, "index.html")
    fh = open(index_html, "w+")
    fh.write("<h1>Results for "+lineid+"</h1>")
    for head in statistics:
        fh.write("<h2>"+head+"</h2>")
        stats = statistics[head]
        if(head == "assembly_statistics"):
            fh.write("<table>")
            header_written = False
            for l in stats:
                fh.write("<tr>")
                for c in l:
                    if not header_written:
                        fh.write("<th>"+str(c)+"</th>")
                    else:
                        fh.write("<td>"+str(c)+"</td>")
                header_written = True
                fh.write("</tr>")
            fh.write("</table>")
        elif(head == "contigs_and_blast_results"):
            for contig in stats:
                fh.write("<h3>"+contig+", Length: "+stats[contig]['length']+"</h3>")
                fh.write("<table>")
                first_hit = stats[contig]["hits"][0]
                fh.write("<tr>")
                for desc in first_hit:
                    fh.write("<th>"+str(desc)+"</th>")
                fh.write("</tr>")
                for hit in stats[contig]["hits"]:
                    fh.write("<tr>")
                    for d in hit:
                        fh.write("<td>"+str(hit[d])+"</td>")
                    fh.write("</tr>")
                fh.write("</table>")
                fh.write("<h4>Visualisation for "+contig+"</h4>")
                fh.write("<img src='"+stats[contig]["image"]+"'/>")
        else:
            fh.write("<table>")
            for subhead in stats:
                fh.write("<tr>")
                fh.write("<td><b>"+subhead+"</b></td>")
                fh.write("<td>"+str(stats[subhead])+"</b></td>")
                fh.write("</tr>")
            fh.write("</table>")
    fh.close()
    
    return index_html
    

def hit_overlaps_other_hit(hit, hits):
    s = int(hit["qstart"])
    e = int(hit["qend"])
    for h in hits:
        overlap = 0
        contained = False
        s2 = int(h["qstart"])
        e2 = int(h["qend"])
        if s < s2:
            overlap = e - s2
        else:
            if e < e2:
                overlap = e - s
                contained = True
            else:
                overlap = e2 -s
        if overlap > ALLOWED_OVERLAP:
            #print ("Overlap: "+str(s)+"/"+str(e)+", "+str(s2)+"/"+str(e2))
            return True
    return False
        
    

def get_annotation_from_blast_result(contigfile, blastfile):
    contigs = {}
    #sort blastfile by length/identity:
    sortedblastfile = blastfile+".sorted"
    os.system("sort -nr -k 4 -k 3 "+blastfile+" > "+sortedblastfile)
    with open(contigfile, 'r') as fh:
        for line in fh:
            match_line = re.match("^>(tig\d+)\slen=(\d+).*", line)
            if match_line is not None: 
                tig_len = match_line.groups()
                tig = tig_len[0]
                length = tig_len[1]
                length_dic = {"length": length}
                contigs[tig] = length_dic
    for contig in contigs:
        hits = []
        borders = []
        with open(sortedblastfile, 'r') as fh:
            for line in fh:
                if(re.match("^"+contig, line)):
                    hit = {}
                    (hit["qaccver"], hit["saccver"], hit["pident"], hit["length"], hit["mismatch"], hit["gapopen"], hit["qstart"], hit["qend"], hit["sstart"], hit["send"], hit["evalue"], hit["bitscore"])= line.split()
                    if(float(hit["pident"]) < ALLOWED_IDENTITY):
                        #print ("Skipping due to identity hit "+hit["qaccver"]+" on "+hit["qaccver"]+" on "+hit["saccver"]+" length "+hit["length"])
                        continue
                    #print ("Examining hit "+hit["qaccver"]+" on "+hit["saccver"]+" length "+hit["length"])
                    if(hit["qaccver"] in ["LB","RB"]):
                        borders.append(hit)
                    elif not hit_overlaps_other_hit(hit, hits):
                        #print ("Hit "+hit["qaccver"]+" on "+hit["saccver"]+" length "+hit["length"]+" is ok!")
                        hits.append(hit)
                    #else:
                        #print ("Hit "+hit["qaccver"]+" on "+hit["saccver"]+" length "+hit["length"]+" overlaps!")
        fh.close()
        for b in borders:
            hits.append(b)
        contigs[contig]["hits"] = hits
    print (contigs)
    return contigs

     

def analyse_insertion (lineid, fastqfile, tdnafile, allfasta, outdir, webdir):
    if not (os.path.isdir(outdir)):
        os.mkdir(outdir)
    statistics= {}
    #get statistics for all reads
    fastq_stats = get_fastq_stats(fastqfile)
    statistics["fastq_stats"] = fastq_stats
    #run BLAST
    blastfile = os.path.join(outdir, lineid+".bls");
    idfile = run_blast(fastqfile, blastfile, tdnafile)
    #filter FASTQ
    filtered_fastq = os.path.join(outdir, lineid+".fastq")
    filter_sequences(idfile, fastqfile, filtered_fastq, WRITE_FASTQ_WITHOUT_TDNA)
    filtered_fastq_stats = get_fastq_stats(filtered_fastq)
    statistics["filtered_fastq_stats"] = filtered_fastq_stats
    #run canu:
    contigfile = run_canu(filtered_fastq, statistics, outdir)

    #get assembly statistics:
    (p1,s1) = os.path.splitext(contigfile)
    assemblytigfile = p1+".layout.tigInfo"
    assembly_statistics = []
    with open (assemblytigfile, "r") as fh:
        for line in fh:
            assembly_statistics_line = line.split()
            assembly_statistics.append(assembly_statistics_line)
    fh.close()
    statistics["assembly_statistics"] = assembly_statistics
    
    
    blast_vs_allfasta = os.path.join(outdir, lineid+"_assembly_vs_allfasta.bls");
    run_blast(contigfile, blast_vs_allfasta, allfasta)

    contigs = get_annotation_from_blast_result(contigfile, blast_vs_allfasta)
    
    statistics["contigs_and_blast_results"] = contigs
    
    for c in contigs:        
        imagename = os.path.join(webdir, c+".png")
        if BE_VERBOSE: print ("Creating image "+imagename+"for "+lineid+".")
        vis.draw_insertion(imagename, contigs[c])
        if BE_VERBOSE: print ("Done creating image "+imagename+"for "+lineid+".")
        statistics["contigs_and_blast_results"][c]["image"] = (os.path.split(imagename))[1]
    html_summary(lineid, statistics, webdir)
    return contigfile
     

#qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore