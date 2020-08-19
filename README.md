# tdna_nanopore

Nanopore T-DNA Anlaysis Tool
============================

This tool extracts reads containing T-DNA out of a full genome sequencing run of Oxford Nanopore data. 
It performs basically the following steps:
- BLAST all reads vs. T-DNA sequences of interest
- run an assembly using canu with these reads
- annotate the resulting contigs using the reference sequence of the host species
- annotate all reads as well
- if available, you can search for insertions in a finished complete genome assembly (see Example call)
- create a summary html page with png images showing the annotation

These results can then be inspected. For simple T-DNA insertions, the insertions are usually fully explained by the assembly and/or individual reads.
Difficult cases might be harder to interpret, for example inverted T-DNA repeats. One repeat is usually of good quality (the downstream one), the other of bad quality. In order to resolve these cases, reads running in both directions need to be identified. 
Other difficult cases are larger inversions, duplications or translocations, in most of these cases a look on the complete sequencing run will be necessary.

Requirements
------------
In order to run the tool, you need to clone the repository::

  git clone https://github.com/nkleinbo/tdna_nanopore
  
In addition, you need some additional tools like BLAST etc. If you are a researcher in Germany, you can apply for a SimpleVM project in the de.NBI cloud, there is a "Nanopore Workbench" image, which contains all necessary tools to run the script.
Otherwise you will need to install the following::

  python3:
  sudo apt install python3
  
  python3 package for image drawing:
  pip3 install Pillow
  
  BLAST:
  sudo apt install ....
  
  Minimap2:
  ...
  
  seqtk:
  ...
  
  canu:
  ...
  
  EMBOSS:
  ...
  
  bedtools:
  ...
  
  samtools:
  ...
  
And, finally, you will need reference files. If you are analysing GABI-Kat or SALK lines, you can use the reference files we provide here: www.www.www
You need two files:
- one file containing all T-DNA sequences of interest, this file is used for identification of reads containing T-DNA
- one file containing all sequences you expect, in the files mentioned above: T-DNA, the vector backbone, A. thaliana. Mitochondria, Chloroplast and A. tumefaciens. This file is used to annotate the generated contigs and reads.
  
  
Example call
------------

The call is as follows (fill in paths to your files/directories)

  python3 run_all.py 
  -f <path to your data path, each line needs a seperate folder in here>
  -o <output directory for result files like fasta files and assemblies>
  -t <location of the fasta file with your T-DNA sequences>
  -a <fasta file with all your references for annotation>
  -w <path to a directory, where the html files should be created>
  -g <path to an assembly, if you have a precomputed assembly of your sequencing run>


Configuration
-------------
  
There are several places where you can configure the behaviour of the script. 

### Reference sequences ###
If you are not using our T-DNA and reference file above, you will need to add the names of your references in the dictionary in visulisation.py::

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

This is to ensure, that each annotated sequence gets the correct colour. Unknown features will get a gray colour. If You need to add all sequence names from your fasta annotation file. For example, if there is a sequence called "genome" which is your plant genome, you need to add "genome": "plant" in this dictionary. If you want to use different colours then the predefined, you can add it in this dictionary::


    COLOURS = {
        "tdna": (200,0,0,1),
        "vector": (150,0,0,1),
        "agrobac": (150,0,150,1),
        "adapter": (100,0,0,1),
        "plant": (0,200,0,1),
        "chloroplast": (0,155,0,1),
        "mitochondria": (0,0,255,1),
        "LB": (255,0,0,1),
        "RB": (255,0,0,1)
    }

... and use its key in the BLAST_FEATURE_MAPPING array.

### Tool behaviour ###

The behaviour of the tool can be modified in several ways. The parameters are constants in the file analyse_insertion.py

You can set the verbose mode and the number of CPUs/threads, that the tool can use::

    BE_VERBOSE = True;
    NO_CPUS = 64;
    
The BLAST parameters, if you don't find what you expect. However, the default values should work fine for most cases::
    
    BLAST_PARAMS = " -perc_identity 85 -evalue 1e-50"

If you are working with something else then A. thaliana, you should change the genome size::

    GENOME_SIZE=145000000

If you want to analyse the reads that do not contain T-DNA you can choose, that the script writes such a file (ending with .notdna.fastq)::

    WRITE_FASTQ_WITHOUT_TDNA = False;

The allowed overlap between BLAST results for the annotation of contigs and the minimum identity for BLAST results::

    ALLOWED_OVERLAP = 10
    ALLOWED_IDENTITY = 90

You can use corrected reads as input (mostly important for canu parameters)::

    CORRECTED_READS = False

Mapping Cutoffs::

    MAPPING_QUAL_CUTOFF = 60
    MAPPING_LENGTH_CUTOFF = 200

Then you can choose, whether you want to rerun intermediate steps of the script in every call, RERUN_ALL overwrites all other states, if True::

    RERUN_BLAST = True
    RERUN_ASSEMBLY = True
    RERUN_MAPPING = True
    RECREATE_FASTQ = True
    RERUN_ALL = True


### Visualisation ###

You can change the generated images by changing the following constants in visualisation.py::

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

However, better leave it as it is. :)