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
BLAST_PARAMS = " -perc_identity 85 -evalue 1e-150 -min_raw_gapped_score 1000"
#expected genome size:
GENOME_SIZE=145000000
#write fastq files WITHOUT T-DNA:
WRITE_FASTQ_WITHOUT_TDNA = False;
#allowed overlap of BLAST results for contig annotation:
ALLOWED_OVERLAP = 10
#allowed identity of BLAST results for contig annotation:
ALLOWED_IDENTITY = 90
#use corrected reads as input:
CORRECTED_READS = False
#bed quality value cutoff for mapping
MAPPING_QUAL_CUTOFF = 60
#mapping length cutoff
MAPPING_LENGTH_CUTOFF = 200
#path to stylesheet:
STYLESHEET = "style/style.css"

def get_id_file(blastfile, filter_reads = None):
    if BE_VERBOSE: print ("Creating idfile for "+blastfile+"...")
    idfile = blastfile+".ids"
    idfile_tmp = blastfile+".tmp.ids"
    os.system("cut -f1 "+blastfile+" | sort | uniq > "+idfile_tmp)
    if(filter_reads is not None):
        os.system("comm -23 "+idfile_tmp+" "+filter_reads+" > "+idfile)
    else:
        os.system("cat "+idfile_tmp+" > "+idfile)
    if BE_VERBOSE: print("Done IDFILE "+blastfile+".")
    return idfile 

#runs BLAST and returns an ID file with all hits
def run_blast(fastqfile, blastfile, dbfile, nofilters=False):
    #check if makeblastdb has been run
    if not (os.path.isfile(dbfile+".nin")):
        os.system("makeblastdb -dbtype nucl -in "+dbfile)
    #run blast if results do not already exist
    if not (os.path.isfile(blastfile)):
        if BE_VERBOSE: print ("BLASTING "+fastqfile+" ...")
        #blastcommand = "gzip -dc "+fastqfile+" | seqtk seq -A | blastn -db "+tdnafile+" -out "+blastfile+" -outfmt 6 -num_threads "+str(NO_CPUS)+" "+BLAST_PARAMS
        (prefix, extension) = os.path.splitext(fastqfile)
        if(nofilters):
            blastparams = ""
        else:
            blastparams = BLAST_PARAMS
        blastcommand = "seqtk seq -A "+fastqfile+" |  parallel -S : --block 10M --recstart '>' --pipe blastn -query - -db "+dbfile+" -outfmt 6 "+blastparams+" > "+blastfile
        if(extension == ".fasta"):
            blastcommand = "cat "+fastqfile+" |  parallel -S : --block 10M --recstart '>' --pipe blastn -query - -db "+dbfile+" -outfmt 6 "+blastparams+" > "+blastfile
        #print(blastcommand)
        os.system(blastcommand)
        if BE_VERBOSE:  "Done BLASTING "+fastqfile+"."
    else:
        if BE_VERBOSE: print ("Skipping BLAST for "+fastqfile+".");
    return blastfile

def run_minimap2(fastqfile, contigfile, bedfile):
    if not (os.path.isfile(bedfile)):
        os.system("minimap2 "+contigfile+" "+fastqfile+" -a -t "+str(NO_CPUS)+" -x map-ont | samtools view -bS - | samtools sort - | bedtools bamtobed -cigar > "+bedfile)
    return bedfile
    

def run_canu(filtered_fastq, statistics, outdir):
    assemblydir = os.path.join(outdir, "assembly")
    projectname = "assembly"
    contigfile = os.path.join(assemblydir, projectname+".contigs.fasta")
    unassembledfile = os.path.join(assemblydir, projectname+".unassembled.fasta")
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
            canu_command="canu -assemble -d "+assemblydir+" -p "+projectname+" -nanopore-corrected "+filtered_fastq+ " 'stopOnLowCoverage=5' 'genomeSize="+str(genome_size)+"' 'useGrid=0' 'corMhapFilterThreshold=0.0000000002' 'ovlMerThreshold=500' 'corMhapOptions=--threshold 0.80 --num-hashes 512 --num-min-matches 3 --ordered-sketch-size 1000 --ordered-kmer-size 14 --min-olap-length 2000 --repeat-idf-scale 50' 'correctedErrorRate=0.17'"
        else:
            canu_command="canu -d "+assemblydir+" -p "+projectname+" -nanopore-raw "+filtered_fastq+ " 'stopOnLowCoverage=5' 'genomeSize="+str(genome_size)+"' 'useGrid=0' 'corMhapFilterThreshold=0.0000000002' 'ovlMerThreshold=500' 'corMhapOptions=--threshold 0.80 --num-hashes 512 --num-min-matches 3 --ordered-sketch-size 1000 --ordered-kmer-size 14 --min-olap-length 2000 --repeat-idf-scale 50' 'correctedErrorRate=0.17'"
        os.system(canu_command)
        if BE_VERBOSE: print ("Done asssembly for "+filtered_fastq+".")
    else:
        if BE_VERBOSE: print ("Skipping asssembly for "+filtered_fastq+"...")
    return (contigfile, unassembledfile)
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
    if number_reads > 0:
        mean_length = total_bases / number_reads
        coverage_for_longest_read = total_bases / longest_read
        coverage_for_genome = total_bases / GENOME_SIZE
    else: 
        mean_length = 0
        coverage_for_longest_read = 0
        coverage_for_genome = 0
    stats["number_reads"] = number_reads
    stats["longest_read"] = longest_read
    stats["total_bases"] = total_bases
    stats["mean_length"] = mean_length
    stats["coverage_for_longest_read"] = coverage_for_longest_read
    stats["coverage_for_genome"] = coverage_for_genome
    return stats
    

def html_summary(lineid, statistics, webdir):
    if not (os.path.isdir(webdir)):
        os.mkdir(webdir)    
    index_html = os.path.join(webdir, "index.html")
    style = ""
    with open(STYLESHEET, 'r') as fh:
        for line in fh:
            style += line
    fh = open(index_html, "w+")
    fh.write("<!doctype html>\n")
    fh.write("<html>\n")
    fh.write("<head>\n")
    fh.write(" <style>\n")
    fh.write(style)
    fh.write("</style>\n")
    fh.write("</head>\n")

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
                if (len(stats[contig]["hits"]) == 0):
                    continue
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
        elif(head == "unassembled_and_blast_results"):
            for contig in stats:
                if (len(stats[contig]["hits"]) == 0):
                    continue
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
    fh.write("</html>\n")
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
 
       
def get_length_for_reads(fastqfile, fasta=False):
    reads = {}
    with open (fastqfile, "r") as fh:
        linecount = 0;
        read = ""
        for line in fh:
            linecount += 1
            if linecount == 1:
                match_line = re.match("^@(\S+)\srunid=.*", line)
                if match_line is not None: 
                    read_match = match_line.groups()
                    read = read_match[0]
            elif linecount == 2:
                    reads[read] = {"length": len(line)}
            elif linecount == 4:
                    linecount = 0
    return reads

def get_mappings_from_blast_results(fastqfile, contigfile, blastfile, fasta=False):
    #get length for reads:
    
    reads = get_length_for_reads(fastqfile, fasta=fasta)

    mappings = get_non_overlapping_hits_from_blast_result(reads, blastfile)        
    
    contigs = get_contigs_with_length(contigfile)
    mappings_by_contig = {}
    for c in contigs:
        mappings_by_contig[c] = {}
    for read in mappings:
        hits = mappings[read]["hits"]
        length = mappings[read]["length"]
        for hit in hits:
            contig_hit = hit["saccver"]
            if not read in mappings_by_contig[contig_hit].keys():
                mappings_by_contig[contig_hit][read] = {}
                mappings_by_contig[contig_hit][read]["hits"] = []
            mappings_by_contig[contig_hit][read]["length"] = length
            mappings_by_contig[contig_hit][read]["hits"].append(hit)
                
    return mappings_by_contig

def get_mappings_from_bedfile(fastqfile, contigfile, mapping_bed_file, fasta=False):
    reads = get_length_for_reads(fastqfile, fasta=fasta)
    contigs = get_contigs_with_length(contigfile)
    mappings_by_contig = {}
    
    for c in contigs:
        mappings_by_contig[c] = {}
#    
#    print(fastqfile)
#    print(reads)
#    print(contigs)
#    
    with open (mapping_bed_file, "r") as fh:
        for line in fh:
            (tig, start, stop, read, qual, strand, cigar) = line.split()
            if int(qual) < MAPPING_QUAL_CUTOFF:
                continue
            read_length = reads[read]["length"]
            max_mapping_length = abs(int(stop) - int(start))
            if(max_mapping_length < MAPPING_LENGTH_CUTOFF):
                continue
            if not read in mappings_by_contig[tig].keys():
                mappings_by_contig[tig][read] = {}
                mappings_by_contig[tig][read]["hits"] = []
                mappings_by_contig[tig][read]["length"] = read_length
                mappings_by_contig[tig][read]["max_mapping_length"] = max_mapping_length
            if(strand == "-"):
                tmp = start
                start = stop
                stop = tmp
                
            qstart = 0
            qend = read_length

            #check cigar string:
            #print(cigar)
            match_cigar_front = re.match("^(\d+[HS])", cigar)
            if match_cigar_front is not None:
                groups = match_cigar_front.groups()
                #print("CIGAR GROUPS FRONT: "+str(groups))
                if(groups[0] != ""):
                    qstart = int(groups[0][:-1])
            match_cigar_end = re.search("(\d+[HS])$", cigar)
            if match_cigar_end is not None:
                groups = match_cigar_end.groups()
                #print("CIGAR GROUPS END: "+str(groups))
                if(groups[0] != ""):
                    qend = int(read_length) - int(groups[0][:-1])

            hit = {"qaccver": read,
                   "saccver": tig,
                   "sstart": start,
                   "send": stop,
                   "qstart": qstart,
                   "qend": qend
                    }
            mappings_by_contig[tig][read]["hits"].append(hit)
            if max_mapping_length >  mappings_by_contig[tig][read]["max_mapping_length"]:
                mappings_by_contig[tig][read]["max_mapping_length"] = max_mapping_length
    #print(mappings_by_contig)
    return mappings_by_contig
    

def get_contigs_with_length(contigfile):
    contigs = {}
    with open(contigfile, 'r') as fh:
        for line in fh:
            match_line = re.match("^>(tig\d+)\slen=(\d+).*", line)
            if match_line is not None: 
                tig_len = match_line.groups()
                tig = tig_len[0]
                length = tig_len[1]
                length_dic = {"length": length}
                contigs[tig] = length_dic
            else:
                match_line = re.match("^>(tig\d+):(\d+)-(\d+)", line)
                if match_line is not None: 
                    tig_len = match_line.groups()
                    tig = tig_len[0]+":"+tig_len[1]+"-"+tig_len[2]
                    length = int(tig_len[2])-int(tig_len[1])
                    length_dic = {"length": str(length)}
                    contigs[tig] = length_dic
    return contigs

def get_non_overlapping_hits_from_blast_result(contigs, blastfile):    
    #sortedblastfile = blastfile+".sorted"
    #os.system("sort -nr -k 12 -k 3 "+blastfile+" > "+sortedblastfile)
    for contig in contigs:
        hits = []
        #borders = []
        with open(blastfile, 'r') as fh:
            for line in fh:
                if(re.match("^"+contig, line)):
                    hit = {}
                    (hit["qaccver"], hit["saccver"], hit["pident"], hit["length"], hit["mismatch"], hit["gapopen"], hit["qstart"], hit["qend"], hit["sstart"], hit["send"], hit["evalue"], hit["bitscore"])= line.split()
                    if(float(hit["pident"]) < ALLOWED_IDENTITY):
                        #print ("Skipping due to identity hit "+hit["qaccver"]+" on "+hit["qaccver"]+" on "+hit["saccver"]+" length "+hit["length"])
                        continue
                    #print ("Examining hit "+hit["qaccver"]+" on "+hit["saccver"]+" length "+hit["length"])
                   # if(hit["qaccver"] in ["LB","RB"]):
                        #borders.append(hit)
                    elif not hit_overlaps_other_hit(hit, hits):
                        #print ("Hit "+hit["qaccver"]+" on "+hit["saccver"]+" length "+hit["length"]+" is ok!")
                        hits.append(hit)
                    #else:
                        #print ("Hit "+hit["qaccver"]+" on "+hit["saccver"]+" length "+hit["length"]+" overlaps!")
        fh.close()
#        for b in borders:
#            hits.append(b)
        contigs[contig]["hits"] = hits
    return contigs
    

def get_annotation_from_blast_result(contigfile, blastfile, allfasta, references_base):
    contigs = get_contigs_with_length(contigfile)
    #sort blastfile by length/identity:
    contigs = get_non_overlapping_hits_from_blast_result(contigs, blastfile) 
    
    #write fasta with references:
    bedfile = blastfile+".bed"
    references_file = references_base+".fasta"
    fh = open(bedfile, "w")
    for c in contigs:
        for hit in contigs[c]["hits"]:
            (subject, start, end) = (hit["saccver"], hit["sstart"], hit["send"])
            if int(start) < int(end):
                fh.write(subject+"\t"+str(start)+"\t"+str(end)+"\n")
            else:
                fh.write(subject+"\t"+str(end)+"\t"+str(start)+"\n")
    fh.close()
    
    bedfile_no_tdna = blastfile+"_no_tdna.bed"
    references_file_no_tdna = references_base+"_no_tdna.fasta"
    fh = open(bedfile_no_tdna, "w")
    for c in contigs:
        for hit in contigs[c]["hits"]:
            (subject, start, end) = (hit["saccver"], hit["sstart"], hit["send"])
            if(vis.BLAST_FEATURE_MAPPING[subject] == "tdna"):
                continue
            if int(start) < int(end):
                fh.write(subject+"\t"+str(start)+"\t"+str(end)+"\n")
            else:
                fh.write(subject+"\t"+str(end)+"\t"+str(start)+"\n")
    fh.close()   
    os.system("bedtools getfasta -fi "+allfasta+" -bed "+bedfile+" -fo "+references_file)
    os.system("bedtools getfasta -fi "+allfasta+" -bed "+bedfile_no_tdna+" -fo "+references_file_no_tdna)
    
    
    return (contigs, references_file, references_file_no_tdna)


def create_contigfile_from_assembly(assemblyfile, blastfile, outdir, lineid, extension=50000):
    hits = []
    with open(blastfile, 'r') as fh:
        for line in fh:
            hit = {}
            split = line.split()
            if len(split) < 12: continue
            if not (re.match("tig", line)): continue
            (hit["qaccver"], hit["saccver"], hit["pident"], hit["length"], hit["mismatch"], hit["gapopen"], hit["qstart"], hit["qend"], hit["sstart"], hit["send"], hit["evalue"], hit["bitscore"])= split
            if(float(hit["pident"]) > 95.0 and int(hit["length"]) > 50):
                hits.append(hit)
    tdna_ranges = []
    for hit in hits:
        (subjectname, qstart, qend) = (hit["qaccver"], int(hit["qstart"]), int(hit["qend"]))
        #check if already annotated:
        already_there = False
        for tdna in tdna_ranges:
            (tdna_subject,tdna_start, tdna_end) = tdna
            if(tdna_subject == subjectname):
                contained = False
                if qstart < tdna_start:
                    overlap = qend - tdna_start
                else:
                    if qend < tdna_end:
                        overlap = qend - qstart
                        contained = True
                    else:
                        overlap = tdna_end - qstart
                if(overlap > -10000 and not contained):
                    already_there = True
                    if(qstart < tdna_start):
                        tdna[1] = qstart
                    if(qend > tdna_end):
                        tdna[2] = qend
        if not already_there:
            new_tdna = [subjectname, qstart, qend]
            tdna_ranges.append(new_tdna)
    #merge nearby:
    unique_tdnas = tdna_ranges
    finished = False
    while not finished:
        finished = True
        for tdna1 in unique_tdnas:
            (tdna_subject1,tdna_start1, tdna_end1) = tdna1
            for tdna2 in tdna_ranges:
                (tdna_subject2,tdna_start2, tdna_end2) = tdna2
                if(tdna_subject1 == tdna_subject2 and tdna_start1 == tdna_start2 and tdna_end1 == tdna_end2):
                    continue
                if(tdna_subject1 == tdna_subject2):
                    contained = False
                    if tdna_start1 < tdna_start2:
                        overlap = tdna_end1 - tdna_start2
                    else:
                        if tdna_end1 < tdna_end2:
                            overlap = tdna_end1 - tdna_start1
                            contained = True
                        else:
                            overlap = tdna_end2 - tdna_start1
                    if(overlap > -10000 and not contained):
                        finished = False
                        if(tdna_start2 < tdna_start1):
                            tdna1[1] = tdna_start2
                        if(tdna_end1 > tdna_end2):
                            tdna1[2] = tdna_end2
        #remove duplicates:
        new_unique = []
        for tdna1 in unique_tdnas:
            (tdna_subject1,tdna_start1, tdna_end1) = tdna1
            its_in = False;
            for tdna2 in new_unique:
                (tdna_subject2,tdna_start2, tdna_end2) = tdna2
                if(tdna_subject1 == tdna_subject2 and tdna_start1 == tdna_start2 and tdna_end1 == tdna_end2):
                    its_in = True
            if not its_in:
                new_unique.append(tdna1)
        unique_tdnas = new_unique
        tdna_ranges = new_unique
    
    os.system("samtools faidx "+assemblyfile)
    contigs_length = {}
    with open(assemblyfile+".fai", 'r') as fh:
        for line in fh:
            (contig, c_length, a, b, c) = line.split()
            contigs_length[contig] = int(c_length)
  
    bedfile = blastfile+".bed"
    fh = open(bedfile, "w")
    for tdna in unique_tdnas:
        (tdna_subject,tdna_start, tdna_end) = tdna
        start = tdna_start-extension
        if start < 0: start = 0
        end = tdna_end+extension
        if end > contigs_length[tdna_subject]: end = contigs_length[tdna_subject]
        fh.write(tdna_subject+"\t"+str(start)+"\t"+str(end)+"\n")
    fh.close()
    contigfile = os.path.join(outdir, lineid+"_contigs.fasta")
    os.system("bedtools getfasta -fi "+assemblyfile+" -bed "+bedfile+" -fo "+contigfile)
    return contigfile        
    
   

     
def analyse_insertion_assembly(lineid, fastqfile, tdnafile, allfasta, outdir, webdir, assembly):
    if not (os.path.isdir(outdir)):
        os.mkdir(outdir)
    statistics= {}
    fastq_stats = get_fastq_stats(fastqfile)
    statistics["fastq_stats"] = fastq_stats
    #run BLAST
    blastfile = os.path.join(outdir, lineid+".bls");

    #create contigfile from assembly, +- 50kb from insertion
    run_blast(assembly, blastfile, tdnafile)
    contigfile = create_contigfile_from_assembly(assembly, blastfile, outdir, lineid)
    if(os.path.isfile(contigfile)):
        #blast contigs vs all possible targets:
        blast_vs_allfasta = os.path.join(outdir, lineid+"_assembly_vs_allfasta.bls");
        run_blast(contigfile, blast_vs_allfasta, allfasta)
    
        #get annotated contigs and a file containing all matched sequences from possible targets
        references_base = os.path.join(outdir, lineid+"_with_flanking")
        (contigs, references_file, references_file_no_tdna) = get_annotation_from_blast_result(contigfile, blast_vs_allfasta, allfasta, references_base)
        statistics["contigs_and_blast_results"] = contigs
        mapping_bed_file = os.path.join(outdir, lineid+"_reads_vs_assembly.bed")

        #map reads back to assembly:
        mapping_bed_file = os.path.join(outdir, lineid+"_reads_vs_assembly.bed")
        run_minimap2(fastqfile, contigfile, mapping_bed_file)
        mappings_by_contig = get_mappings_from_bedfile(fastqfile, contigfile, mapping_bed_file)
            

  
        for c in contigs: 
            if not (os.path.isdir(webdir)):
                os.mkdir(webdir)    
            c_name = c.replace(":", "_")
            imagename = os.path.join(webdir, c_name+".png")
            if BE_VERBOSE: print ("Creating image "+imagename+" for "+lineid+".")
            if(mappings_by_contig is not None):
                vis.draw_insertion(imagename, contigs[c], mappings_by_contig[c])
            else:
                vis.draw_insertion(imagename, contigs[c], None)
            if BE_VERBOSE: print ("Done creating image "+imagename+" for "+lineid+".")
            statistics["contigs_and_blast_results"][c]["image"] = (os.path.split(imagename))[1]
    html_summary(lineid, statistics, webdir)
    return (contigfile, references_file, references_file_no_tdna)      


def analyse_insertion (lineid, fastqfile, tdnafile, allfasta, outdir, webdir, filter_reads = None):
    if not (os.path.isdir(outdir)):
        os.mkdir(outdir)
    statistics= {}
    fastq_stats = get_fastq_stats(fastqfile)
    statistics["fastq_stats"] = fastq_stats
    #run BLAST
    blastfile = os.path.join(outdir, lineid+".bls");

    run_blast(fastqfile, blastfile, tdnafile)
    idfile = get_id_file(blastfile, filter_reads)
    #filter FASTQ
    filtered_fastq = os.path.join(outdir, lineid+".fastq")
    filter_sequences(idfile, fastqfile, filtered_fastq, WRITE_FASTQ_WITHOUT_TDNA)
    filtered_fastq_stats = get_fastq_stats(filtered_fastq)
    statistics["filtered_fastq_stats"] = filtered_fastq_stats
    #run canu:
    (contigfile, unassembledfile) = run_canu(filtered_fastq, statistics, outdir)
    #get assembly statistics:
    (p1,s1) = os.path.splitext(contigfile)
    assemblytigfile = p1+".layout.tigInfo"
    assembly_statistics = []
    if(os.path.isfile(assemblytigfile)):
        with open (assemblytigfile, "r") as fh:
            for line in fh:
                assembly_statistics_line = line.split()
                assembly_statistics.append(assembly_statistics_line)
        fh.close()
        statistics["assembly_statistics"] = assembly_statistics
    (references_file, references_file_no_tdna) = (None, None)

    
    #annotate contigs
    if(os.path.isfile(contigfile)):
        #blast contigs vs all possible targets:
        blast_vs_allfasta = os.path.join(outdir, lineid+"_assembly_vs_allfasta.bls");
        run_blast(contigfile, blast_vs_allfasta, allfasta)
    
        #get annotated contigs and a file containing all matched sequences from possible targets
        references_base = os.path.join(outdir, lineid+"_with_flanking")
        (contigs, references_file, references_file_no_tdna) = get_annotation_from_blast_result(contigfile, blast_vs_allfasta, allfasta, references_base)
        statistics["contigs_and_blast_results"] = contigs
        
        #map reads back to assembly
        mapping_bed_file = os.path.join(outdir, lineid+"_reads_vs_assembly.bed")
        run_minimap2(filtered_fastq, contigfile, mapping_bed_file)
        mappings_by_contig = get_mappings_from_bedfile(filtered_fastq, contigfile, mapping_bed_file)

        for c in contigs: 
            if not (os.path.isdir(webdir)):
                os.mkdir(webdir)    
            c_name = c.replace(":", "_")
            imagename = os.path.join(webdir, c_name+".png")
            if BE_VERBOSE: print ("Creating image "+imagename+" for "+lineid+".")
            if(mappings_by_contig is not None):
                vis.draw_insertion(imagename, contigs[c], mappings_by_contig[c])
            else:
                vis.draw_insertion(imagename, contigs[c], None)
            if BE_VERBOSE: print ("Done creating image "+imagename+" for "+lineid+".")
            statistics["contigs_and_blast_results"][c]["image"] = (os.path.split(imagename))[1]
    
    #annotate unassembled
    if(os.path.isfile(unassembledfile)):
        #blast contigs vs all possible targets:
        blast_vs_allfasta_ua = os.path.join(outdir, lineid+"_unassembled_vs_allfasta.bls");
        run_blast(unassembledfile, blast_vs_allfasta_ua, allfasta)
    
        #get annotated contigs and a file containing all matched sequences from possible targets
        references_base_ua = os.path.join(outdir, lineid+"unassembled_with_flanking")
        (contigs_ua, references_file_ua, references_file_no_tdna_ua) = get_annotation_from_blast_result(unassembledfile, blast_vs_allfasta_ua, allfasta, references_base_ua)
        statistics["unassembled_and_blast_results"] = contigs_ua
        
        #map reads back to assembly
        mapping_bed_file_ua = os.path.join(outdir, lineid+"_reads_vs_unassembled.bed")
        run_minimap2(filtered_fastq, unassembledfile, mapping_bed_file_ua)
        mappings_by_contig_ua = get_mappings_from_bedfile(filtered_fastq, unassembledfile, mapping_bed_file_ua)

        for c in contigs_ua: 
            if not (os.path.isdir(webdir)):
                os.mkdir(webdir)    
            c_name = c.replace(":", "_")
            imagename = os.path.join(webdir, c_name+".png")
            if BE_VERBOSE: print ("Creating image "+imagename+" for "+lineid+".")
            if(mappings_by_contig_ua is not None):
                vis.draw_insertion(imagename, contigs_ua[c], mappings_by_contig_ua[c])
            else:
                vis.draw_insertion(imagename, contigs_ua[c], None)
            if BE_VERBOSE: print ("Done creating image "+imagename+" for "+lineid+".")
            statistics["unassembled_and_blast_results"][c]["image"] = (os.path.split(imagename))[1]


            
    
    html_summary(lineid, statistics, webdir)
    return (contigfile, references_file, references_file_no_tdna, idfile)
     

#qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore