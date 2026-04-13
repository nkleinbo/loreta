#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import re
import subprocess
import visualisation as vis
from Config import (BE_VERBOSE, NO_CPUS, BLAST_PARAMS, GENOME_SIZE, 
                    WRITE_FASTQ_WITHOUT_TDNA, ALLOWED_OVERLAP, 
                    ALLOWED_IDENTITY, CORRECTED_READS, MAPPING_QUAL_CUTOFF, 
                    MAPPING_LENGTH_CUTOFF, STYLESHEET, RERUN_BLAST, 
                    RERUN_ASSEMBLY, RERUN_MAPPING, RECREATE_FASTQ, RERUN_ALL)

# --- HELPER FUNCTIONS ---

def run_cmd(cmd):
    """Wrapper for subprocess to replace os.system"""
    subprocess.run(cmd, shell=True)

def get_id_file(blastfile, filter_reads=None):
    if BE_VERBOSE: print(f"Creating idfile for {blastfile}...")
    idfile = f"{blastfile}.ids"
    idfile_tmp = f"{blastfile}.tmp.ids"
    
    run_cmd(f"cut -f1 {blastfile} | sort | uniq > {idfile_tmp}")
    
    if filter_reads is not None:
        run_cmd(f"comm -23 {idfile_tmp} {filter_reads} > {idfile}")
    else:
        run_cmd(f"cat {idfile_tmp} > {idfile}")
        
    if BE_VERBOSE: print(f"Done IDFILE {blastfile}.")
    return idfile 

def run_blast(fastqfile, blastfile, dbfile, nofilters=False):
    if not os.path.isfile(f"{dbfile}.nin"):
        run_cmd(f"makeblastdb -dbtype nucl -in {dbfile}")
        
    if not os.path.isfile(blastfile) or RERUN_BLAST or RERUN_ALL:
        if BE_VERBOSE: print(f"BLASTING {fastqfile} ...")
        prefix, extension = os.path.splitext(fastqfile)
        blastparams = "" if nofilters else BLAST_PARAMS
        
        if extension == ".fasta":
            blastcommand = f"cat {fastqfile} | parallel --block 10M --recstart '>' --pipe blastn -query - -db {dbfile} -outfmt 6 {blastparams} > {blastfile}"
        else:
            blastcommand = f"seqtk seq -A {fastqfile} | parallel --block 10M --recstart '>' --pipe blastn -query - -db {dbfile} -outfmt 6 {blastparams} > {blastfile}"
            
        run_cmd(blastcommand)
        if BE_VERBOSE: print(f"Done BLASTING {fastqfile}.")
    else:
        if BE_VERBOSE: print(f"Skipping BLAST for {fastqfile}.")
    return blastfile

def run_minimap2(fastqfile, contigfile, bedfile):
    if not os.path.isfile(bedfile) or RERUN_MAPPING or RERUN_ALL:
        run_cmd(f"minimap2 {contigfile} {fastqfile} -a -t {NO_CPUS} -x map-ont | samtools view -bS - | samtools sort - | bedtools bamtobed -cigar > {bedfile}")
    return bedfile

def create_fasta_from_fastq(fastqfile):
    prefix, ext = os.path.splitext(fastqfile)
    fastafile = f"{prefix}.fasta"
    run_cmd(f"seqtk seq -A {fastqfile} > {fastafile}")
    return fastafile

def run_canu(filtered_fastq, statistics, outdir):
    assemblydir = os.path.join(outdir, "assembly")
    projectname = "assembly"
    contigfile = os.path.join(assemblydir, f"{projectname}.contigs.fasta")
    unassembledfile = os.path.join(assemblydir, f"{projectname}.unassembled.fasta")
    
    longest_read = statistics["filtered_fastq_stats"]["longest_read"]
    coverage = statistics["filtered_fastq_stats"]["coverage_for_longest_read"]
    genome_size = longest_read if coverage > 10 else int((longest_read * coverage) / 11)
    
    if os.path.isfile(contigfile) and (RERUN_ASSEMBLY or RERUN_ALL):
        run_cmd(f"rm -rf {assemblydir}")
        
    if not os.path.isfile(contigfile):
        if BE_VERBOSE: print(f"Starting asssembly for {filtered_fastq}...")
        
        base_opts = f"-d {assemblydir} -p {projectname} 'stopOnLowCoverage=5' 'useGrid=0' 'corMhapFilterThreshold=0.0000000002' 'ovlMerThreshold=500' 'corMhapOptions=--threshold 0.80 --num-hashes 512 --num-min-matches 3 --ordered-sketch-size 1000 --ordered-kmer-size 14 --min-olap-length 2000 --repeat-idf-scale 50' 'correctedErrorRate=0.17'"
        
        if CORRECTED_READS:
            canu_command = f"canu -assemble {base_opts} -nanopore-corrected {filtered_fastq} 'genomeSize={genome_size}'"
        else:
            canu_command = f"canu {base_opts} -nanopore-raw {filtered_fastq} 'genomeSize=10k' -fast 'corOutCoverage=200'"
            
        run_cmd(canu_command)
        if BE_VERBOSE: print(f"Done asssembly for {filtered_fastq}.")
    else:
        if BE_VERBOSE: print(f"Skipping asssembly for {filtered_fastq}...")
        
    return contigfile, unassembledfile

def filter_sequences(idfile, fastqfile, filtered_fastq, write_remaining=False):
    if os.path.isfile(filtered_fastq) and (RECREATE_FASTQ or RERUN_ALL):
        run_cmd(f"rm {filtered_fastq}")
        
    if not os.path.isfile(filtered_fastq):
        if BE_VERBOSE: print(f"Grepping sequences from {fastqfile}...")
        run_cmd(f"grep -A 3 --no-group-separator -f {idfile} {fastqfile} > {filtered_fastq}")
        if BE_VERBOSE: print(f"Done grepping sequences from {fastqfile}.")
    else:
        if BE_VERBOSE: print(f"Skipping grepping sequences from {fastqfile}.")
        
    if write_remaining:
        prefix, ext = os.path.splitext(fastqfile)
        remaining_fastqfile = f"{prefix}.notdna.fastq"
        remaining_idfile = f"{prefix}.notdna.ids"
        
        if os.path.isfile(remaining_fastqfile) and (RECREATE_FASTQ or RERUN_ALL):
            run_cmd(f"rm {remaining_fastqfile}")
            
        if not os.path.isfile(remaining_fastqfile):
            if BE_VERBOSE: print(f"Grepping non-T-DNA sequences from {fastqfile}...")
            run_cmd(f"grep '^@.*runid' {fastqfile} | grep -v -f {idfile} > {remaining_idfile}")
            run_cmd(f"grep -A 3 --no-group-separator -f {remaining_idfile} {fastqfile} > {remaining_fastqfile}") 
            if BE_VERBOSE: print(f"Done grepping non-T-DNA sequences from {fastqfile}.")
        else:
            if BE_VERBOSE: print(f"Skipping grepping non-T-DNA sequences from {fastqfile}.")

def get_fastq_stats(fastqfile):
    number_reads, longest_read, total_bases = 0, 0, 0
    with open(fastqfile, 'r') as fh:
        for i, line in enumerate(fh, 1):
            if i % 4 == 0:
                number_reads += 1
            elif i % 2 == 0:
                bases = len(line.strip())
                total_bases += bases
                if bases > longest_read:
                    longest_read = bases
                    
    mean_length = total_bases / number_reads if number_reads > 0 else 0
    cov_longest = total_bases / longest_read if longest_read > 0 else 0
    cov_genome = total_bases / GENOME_SIZE if GENOME_SIZE > 0 else 0
    
    return {
        "number_reads": number_reads,
        "longest_read": longest_read,
        "total_bases": total_bases,
        "mean_length": mean_length,
        "coverage_for_longest_read": cov_longest,
        "coverage_for_genome": cov_genome
    }

def html_summary(lineid, statistics, webdir):
    os.makedirs(webdir, exist_ok=True)
    index_html = os.path.join(webdir, "index.html")
    
    with open(STYLESHEET, 'r') as fh:
        style = fh.read()
        
    html = [
        "<!doctype html>", "<html>", "<head>", f"<style>\n{style}</style>", "</head>",
        f"<body>\n<h1>Results for {lineid}</h1>"
    ]
    
    for head, stats in statistics.items():
        html.append(f"<h2>{head}</h2>")
        if head == "assembly_statistics":
            html.append("<table>")
            for i, row in enumerate(stats):
                tag = "th" if i == 0 else "td"
                row_html = "".join(f"<{tag}>{c}</{tag}>" for c in row)
                html.append(f"<tr>{row_html}</tr>")
            html.append("</table>")
            
        elif head in ["contigs_and_blast_results", "unassembled_and_blast_results", "reads_and_blast_results"]:
            for contig, data in stats.items():
                if not data.get("hits"): continue
                html.append(f"<h3>{contig}, Length: {data['length']}</h3>")
                html.append("<table>")
                
                first_hit = data["hits"][0]
                html.append("<tr>" + "".join(f"<th>{desc}</th>" for desc in first_hit.keys()) + "</tr>")
                
                for hit in data["hits"]:
                    html.append("<tr>" + "".join(f"<td>{val}</td>" for val in hit.values()) + "</tr>")
                    
                html.append("</table>")
                if "image" in data:
                    html.append(f"<h4>Visualisation for {contig}</h4><img src='{data['image']}'/>")
        else:
            html.append("<table>")
            for subhead, val in stats.items():
                html.append(f"<tr><td><b>{subhead}</b></td><td>{val}</td></tr>")
            html.append("</table>")
            
    html.append("</body></html>\n")
    
    with open(index_html, "w") as fh:
        fh.write("\n".join(html))
    return index_html

def hit_overlaps_other_hit(hit, hits):
    s, e = int(hit["qstart"]), int(hit["qend"])
    for h in hits:
        s2, e2 = int(h["qstart"]), int(h["qend"])
        overlap = e - s2 if s < s2 else (e - s if e < e2 else e2 - s)
        if overlap > ALLOWED_OVERLAP:
            return True
    return False

def get_length_for_reads(fastqfile, fasta=False):
    reads = {}
    read_pattern = re.compile(r"^@(\S+)\srunid=.*")
    with open(fastqfile, "r") as fh:
        read_id = ""
        for i, line in enumerate(fh, 1):
            if i % 4 == 1:
                match = read_pattern.match(line)
                if match: read_id = match.group(1)
            elif i % 4 == 2 and read_id:
                reads[read_id] = {"length": len(line.strip())}
    return reads

def get_mappings_from_bedfile(fastqfile, contigfile, mapping_bed_file, fasta=False):
    reads = get_length_for_reads(fastqfile, fasta=fasta)
    contigs = get_contigs_with_length(contigfile)
    mappings_by_contig = {c: {} for c in contigs}
    
    # Vorab kompilierte Regex für massiven Speedup
    cigar_front_pattern = re.compile(r"^(\d+[HS])")
    cigar_end_pattern = re.compile(r"(\d+[HS])$")
    
    with open(mapping_bed_file, "r") as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < 7: continue
            tig, start, stop, read, qual, strand, cigar = parts[:7]
            
            if int(qual) < MAPPING_QUAL_CUTOFF: continue
            
            read_length = reads.get(read, {"length": 0})["length"]
            if read_length == 0: continue
            
            start, stop = int(start), int(stop)
            max_mapping_length = abs(stop - start)
            
            if max_mapping_length < MAPPING_LENGTH_CUTOFF: continue
            
            if read not in mappings_by_contig.get(tig, {}):
                if tig not in mappings_by_contig: mappings_by_contig[tig] = {}
                mappings_by_contig[tig][read] = {"hits": [], "length": read_length, "max_mapping_length": 0}
                
            if strand == "-": start, stop = stop, start
                
            qstart, qend = 0, read_length
            
            m_front = cigar_front_pattern.match(cigar)
            if m_front and m_front.group(1): qstart = int(m_front.group(1)[:-1])
                
            m_end = cigar_end_pattern.search(cigar)
            if m_end and m_end.group(1): qend = read_length - int(m_end.group(1)[:-1])

            hit = {"qaccver": read, "saccver": tig, "sstart": start, "send": stop, "qstart": qstart, "qend": qend}
            mappings_by_contig[tig][read]["hits"].append(hit)
            mappings_by_contig[tig][read]["max_mapping_length"] = max(mappings_by_contig[tig][read]["max_mapping_length"], max_mapping_length)
            
    return mappings_by_contig

def get_contigs_with_length(contigfile):
    contigs = {}
    lengthfile = f"{contigfile}.length"
    run_cmd(f"infoseq -noheading -only -name -length {contigfile} > {lengthfile}")
    
    if os.path.exists(lengthfile):
        with open(lengthfile, 'r') as fh:
            for line in fh:
                parts = line.split()
                if len(parts) == 2:
                    contigs[parts[0]] = {"length": parts[1]}
    return contigs

def get_non_overlapping_hits_from_blast_result(contigs, blastfile):    
    sortedblastfile = f"{blastfile}.sorted"
    run_cmd(f"sort -gr -k 12 {blastfile} > {sortedblastfile}")
    
    for contig in contigs:
        hits = []
        with open(sortedblastfile, 'r') as fh:
            for line in fh:
                if line.startswith(contig):
                    parts = line.split()
                    if len(parts) < 12: continue
                    hit = {k: v for k, v in zip(["qaccver", "saccver", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"], parts)}
                    
                    if float(hit["pident"]) >= ALLOWED_IDENTITY and not hit_overlaps_other_hit(hit, hits):
                        hits.append(hit)
        contigs[contig]["hits"] = hits
    return contigs

def get_annotation_from_blast_result(contigfile, blastfile, allfasta, references_base):
    contigs = get_contigs_with_length(contigfile)
    contigs = get_non_overlapping_hits_from_blast_result(contigs, blastfile) 

    bedfile = f"{blastfile}.bed"
    bedfile_no_tdna = f"{blastfile}_no_tdna.bed"
    
    with open(bedfile, "w") as fh, open(bedfile_no_tdna, "w") as fh_no:
        for c, data in contigs.items():
            for hit in data.get("hits", []):
                subject, start, end = hit["saccver"], int(hit["sstart"]), int(hit["send"])
                s, e = min(start, end), max(start, end)
                
                line = f"{subject}\t{s}\t{e}\n"
                fh.write(line)
                if not (subject in vis.BLAST_FEATURE_MAPPING and vis.BLAST_FEATURE_MAPPING[subject] == "tdna"):
                    fh_no.write(line)
                    
    ref_file = f"{references_base}.fasta"
    ref_file_no_tdna = f"{references_base}_no_tdna.fasta"
    
    run_cmd(f"bedtools getfasta -fi {allfasta} -bed {bedfile} -fo {ref_file}")
    run_cmd(f"bedtools getfasta -fi {allfasta} -bed {bedfile_no_tdna} -fo {ref_file_no_tdna}")
    
    return contigs, ref_file, ref_file_no_tdna

def create_contigfile_from_assembly(assemblyfile, blastfile, outdir, lineid, extension=50000):
    tdna_ranges = []
    with open(blastfile, 'r') as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < 12: continue
            if float(parts[2]) > 95.0 and int(parts[3]) > 50:
                qstart, qend = sorted([int(parts[6]), int(parts[7])])
                tdna_ranges.append([parts[0], qstart, qend])

    tdna_ranges.sort(key=lambda x: (x[0], x[1]))
    unique_tdnas = []
    
    for current in tdna_ranges:
        if not unique_tdnas:
            unique_tdnas.append(current)
            continue
        prev = unique_tdnas[-1]
        if prev[0] == current[0] and current[1] - prev[2] <= 10000:
            prev[2] = max(prev[2], current[2])
        else:
            unique_tdnas.append(current)

    run_cmd(f"samtools faidx {assemblyfile}")
    contigs_length = {}
    if os.path.exists(f"{assemblyfile}.fai"):
        with open(f"{assemblyfile}.fai", 'r') as fh:
            for line in fh:
                parts = line.split()
                contigs_length[parts[0]] = int(parts[1])
                
    bedfile = f"{blastfile}.bed"
    with open(bedfile, "w") as fh:
        for subj, start, end in unique_tdnas:
            s = max(0, start - extension)
            e = min(contigs_length.get(subj, end + extension), end + extension)
            fh.write(f"{subj}\t{s}\t{e}\n")
            
    contigfile = os.path.join(outdir, f"{lineid}_contigs.fasta")
    run_cmd(f"bedtools getfasta -fi {assemblyfile} -bed {bedfile} -fo {contigfile}")
    run_cmd(f"sed -i 's/:/_/g' {contigfile}")
    
    return contigfile

# --- DRY HELPER FOR ANALYSIS ---

def process_and_annotate(target_file, lineid, allfasta, outdir, webdir, statistics, stat_key, is_fastq=False, fastq_for_mapping=None):
    if not os.path.isfile(target_file):
        return None, None
        
    prefix = stat_key.split("_")[0]
    blast_vs_allfasta = os.path.join(outdir, f"{lineid}_{prefix}_vs_allfasta.bls")
    run_blast(target_file, blast_vs_allfasta, allfasta)
    
    ref_base = os.path.join(outdir, f"{lineid}_{prefix}_with_flanking")
    contigs, ref_file, ref_no_tdna = get_annotation_from_blast_result(target_file, blast_vs_allfasta, allfasta, ref_base)
    statistics[stat_key] = contigs
    
    mappings_by_contig = None
    if not is_fastq and fastq_for_mapping:
        bed_file = os.path.join(outdir, f"{lineid}_reads_vs_{prefix}.bed")
        run_minimap2(fastq_for_mapping, target_file, bed_file)
        mappings_by_contig = get_mappings_from_bedfile(fastq_for_mapping, target_file, bed_file)

    os.makedirs(webdir, exist_ok=True)
    for c, data in contigs.items(): 
        c_name = c.replace(":", "_")
        img_name = os.path.join(webdir, f"{c_name}.png")
        
        if BE_VERBOSE: print(f"Creating image {img_name} for {lineid}.")
        mapping_data = mappings_by_contig.get(c) if mappings_by_contig else None
        vis.draw_insertion(img_name, data, mapping_data)
        
        statistics[stat_key][c]["image"] = os.path.basename(img_name)
        
    return ref_file, ref_no_tdna

# --- MAIN PIPELINE FUNCTIONS ---

def analyse_insertion_assembly(lineid, fastqfile, tdnafile, allfasta, outdir, webdir, assembly):
    os.makedirs(outdir, exist_ok=True)
    statistics = {"fastq_stats": get_fastq_stats(fastqfile)}
    
    blastfile = os.path.join(outdir, f"{lineid}.bls")
    run_blast(assembly, blastfile, tdnafile)
    
    contigfile = create_contigfile_from_assembly(assembly, blastfile, outdir, lineid)
    
    ref_file, ref_no_tdna = process_and_annotate(
        contigfile, lineid, allfasta, outdir, webdir, statistics, "contigs_and_blast_results", fastq_for_mapping=fastqfile
    )
    
    html_summary(lineid, statistics, webdir)
    return contigfile, ref_file, ref_no_tdna      

def analyse_insertion(lineid, fastqfile, tdnafile, allfasta, outdir, webdir, filter_reads=None):
    os.makedirs(outdir, exist_ok=True)
    statistics = {"fastq_stats": get_fastq_stats(fastqfile)}
    
    blastfile = os.path.join(outdir, f"{lineid}.bls")
    run_blast(fastqfile, blastfile, tdnafile)
    idfile = get_id_file(blastfile, filter_reads)
    
    filtered_fastq = os.path.join(outdir, f"{lineid}.fastq")
    filter_sequences(idfile, fastqfile, filtered_fastq, WRITE_FASTQ_WITHOUT_TDNA)
    statistics["filtered_fastq_stats"] = get_fastq_stats(filtered_fastq)
    
    contigfile, unassembledfile = run_canu(filtered_fastq, statistics, outdir)
    
    assemblytigfile = f"{os.path.splitext(contigfile)[0]}.layout.tigInfo"
    if os.path.isfile(assemblytigfile):
        with open(assemblytigfile, "r") as fh:
            statistics["assembly_statistics"] = [line.split() for line in fh]

    # DRY: Alle drei Annotations-Blöcke (Contigs, Unassembled, Reads) wurden in diese 3 Zeilen komprimiert
    process_and_annotate(contigfile, lineid, allfasta, outdir, webdir, statistics, "contigs_and_blast_results", fastq_for_mapping=filtered_fastq)
    process_and_annotate(unassembledfile, lineid, allfasta, outdir, webdir, statistics, "unassembled_and_blast_results", fastq_for_mapping=filtered_fastq)
    
    fastafile = create_fasta_from_fastq(filtered_fastq)
    ref_file, ref_no_tdna = process_and_annotate(fastafile, lineid, allfasta, outdir, webdir, statistics, "reads_and_blast_results", is_fastq=True)

    html_summary(lineid, statistics, webdir)
    return contigfile, ref_file, ref_no_tdna, idfile
