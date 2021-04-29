#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 12:30:47 2019

@author: nkleinbo
"""

#whether to be verbose or not:
BE_VERBOSE = True;
#no of CPUs to use at max:
NO_CPUS = 8;
#BLAST parameters for annotation:
BLAST_PARAMS = " -word_size 7 -perc_identity 70 -evalue 1e-50"
#BLAST parameters for finding reads:
#BLAST_PARAMS = " -word_size 7 -perc_identity 70 -evalue 1e-50"
#expected genome size:
GENOME_SIZE=145000000
#write fastq files WITHOUT T-DNA:
WRITE_FASTQ_WITHOUT_TDNA = False;
#allowed overlap of BLAST results for contig annotation:
ALLOWED_OVERLAP = 20
#allowed identity of BLAST results for contig annotation:
ALLOWED_IDENTITY = 80
#use corrected reads as input:
CORRECTED_READS = False
#bed quality value cutoff for mapping
MAPPING_QUAL_CUTOFF = 60
#mapping length cutoff
MAPPING_LENGTH_CUTOFF = 200
#path to stylesheet:
STYLESHEET = "style/style.css"
#Rerun BLASTS:
RERUN_BLAST = False
#Rerun Assemblies:
RERUN_ASSEMBLY = False
#Rerun Assemblies:
RERUN_MAPPING = False
#Recreate filtered fastqs
RECREATE_FASTQ = False
#RERUN ALL
RERUN_ALL = True

#Image options:

WIDTH = 1600
HEIGHT_WITHOUT_MAPPINGS = 100
OFFSET = 50
OFFSET_TOP = 50
OFFSET_MAPPINGS = HEIGHT_WITHOUT_MAPPINGS
MAPPING_WIDTH = 2
MAPPPING_COLOUR = (0,0,0)
MAPPPING_COLOUR_READ = (200,200,200)
LINEWIDTH = 5
FEATUREWIDTH = 16
FONTSIZE = 12
ARROW_OFFSET = 4
FONT = "fonts/TruenoRg.otf"

BLAST_FEATURE_MAPPING = {
    "Agrobacterium_circular": "agrobac",
    "Agrobacterium_plasmid": "agrobac",
    "Agrobacterium_linear": "agrobac",
    "At_Chr1": "plant",
    "At_Chr2": "plant",
    "At_Chr3": "plant",
    "At_Chr4": "plant",
    "At_Chr5": "plant",
    "Chr1_AthCol0": "plant",
    "Chr2_AthCol0": "plant",
    "Chr3_AthCol0": "plant",
    "Chr4_AthCol0": "plant",
    "Chr5_AthCol0": "plant",
    "Mitochondria": "mitochondria",
    "Chloroplast": "chloroplast",
    "T732": "adapter",
    "LR32": "adapter",
    "pAC106_tdna": "tdna",
    "pAC106_vec": "vector",
    "pAC161_tdna": "tdna",
    "pAC161_vec": "vector",
    "pGABI1_tdna": "tdna",
    "pGABI1_vec": "vector",
    "pADIS1_tdna": "tdna",
    "pADIS1_vec": "vector",
    "pROK2_vec": "vector",
    "pRok2_tdna": "tdna",
    "LB": "LB",
    "RB": "RB"
}


COLOURS = {
    "tdna": (200,0,0,1),
    "vector": (150,0,0,1),
    "agrobac": (150,0,150,1),
    "adapter": (100,0,0,1),
    "plant": (0,200,0,1),
    "chloroplast": (0,155,0,1),
    "mitochondria": (0,0,255,1),
    "LB": (255,0,0,1),
    "RB": (255,0,0,1),
    "other": (155,155,155,1)
}