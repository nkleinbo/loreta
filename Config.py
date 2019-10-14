#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 12:30:47 2019

@author: nkleinbo
"""

#wether to be verbose or not:
BE_VERBOSE = True;
#no of CPUs to use at max:
NO_CPUS = 64;
#BLAST parameters:
BLAST_PARAMS = " -perc_identity 85"
#expected genome size:
GENOME_SIZE=145000000
#write fastq files WITHOUT T-DNA:
WRITE_FASTQ_WITHOUT_TDNA = True;

#Image options:

WIDTH = 800
HEIGHT = 200
LINEWIDTH = 5
FEATUREWIDTH = 10

BLAST_FEATURE_MAPPING = {
    "Agrobacterium_circular": "agrobac",
    "Agrobacterium_plasmid": "agrobac",
    "Agrobacterium_linear": "agrobac",
    "At_Chr1": "plant",
    "At_Chr2": "plant",
    "At_Chr3": "plant",
    "At_Chr4": "plant",
    "At_Chr5": "plant",
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
    "tdna": (200,0,0),
    "vector": (150,0,0),
    "adapter": (100,0,0),
    "plant": (0,255,0),
    "chlorplast": (0,155,0),
    "mitochondria": (0,0,255),
    "LB": (255,0,0),
    "RB": (255,0,0)
}