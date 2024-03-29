# loreta

Long Read(-based) T-DNA (insertion) Analysis 
============================

This tool extracts reads containing T-DNA out of a full genome sequencing run of Oxford Nanopore Technologies (ONT) data. 
It performs basically the following steps:
- BLAST all reads vs. T-DNA sequences of interest
- run an assembly using canu with these reads
- annotate the resulting contigs using the reference sequence of the host species
- annotate features in all reads
- if available, you can search for insertions in a complete genome assembly (see 'Example call')
- create a summary HTML page with PNG images showing the annotation

These results can then be inspected. For simple T-DNA insertions, the insertions are usually fully explained by the assembly and/or individual reads.
Difficult cases might be harder to interpret, for example inverted T-DNA repeats. One repeat is usually of good quality (the downstream one), the other of bad quality. In order to resolve these cases, reads running in both directions need to be identified. 
Other difficult cases are larger inversions, duplications or translocations, in most of these cases a look on the complete sequencing run will be necessary.

Requirements
------------
In order to run the tool, you need to clone the repository:
```bash
git clone https://github.com/nkleinbo/loreta
```

  
Some additional tools are required. If you are a researcher in Germany, you can apply for a SimpleVM project in the de.NBI cloud, there is a "Nanopore Workbench" image, which contains all necessary tools to run the script.


Otherwise you will need to install the following::

* [Python](https://www.python.org/) 3.7 or later
* [Pillow](https://pypi.org/project/Pillow/)
* [BLAST+](https://www.ncbi.nlm.nih.gov/books/NBK279690/)
* [Samtools](http://www.htslib.org/) 
* [Minimap2](https://github.com/lh3/minimap2)
* [seqtk](https://github.com/lh3/seqtk)
* [canu](https://github.com/marbl/canu)
* [EMBOSS](http://emboss.sourceforge.net/)
* [bedtools](https://github.com/arq5x/bedtools2)
* [gnu parallel](https://www.gnu.org/software/parallel/)  
  
Finally, you will need files containing reference sequences. If you are analysing [GABI-Kat](https://www.gabi-kat.de/) or [SALK](http://signal.salk.edu/tdna_FAQs.html) lines, you can use the reference files provided in the folder 'references' (gunzip them).
You need two files:
- one file containing all T-DNA sequences of interest. This file is used for identification of reads containing T-DNA
- one file containing all sequences you expect. In the files mentioned above: T-DNA sequences, the vector backbone, *A. thaliana* genome sequence, chondrome (mitochondria), plastome (chloroplast), and *A. tumefaciens* genome sequence. This file is used to annotate the generated contigs and reads.

 
  
Example call
------------

The call is as follows (fill in paths to your files/directories)
```bash
python3 run_all.py 
  -f <path to your data directory, each plant line needs a seperate FASTQ file in this directory>
  -o <output directory for result files like FASTA files and assemblies>
  -t <location of the FASTA file with your T-DNA sequences>
  -a <FASTA file with all your references for annotation>
  -w <path to a directory, where the html files should be created>
  -g <path to an assembly, if you have a precomputed assembly of your sequencing run>

```
The script will search for FASTQ files in the data dir (make sure, there is nothing else in it) and run an analysis for each FASTQ file. The prefix of the FASTQ file will serve as folder name for the directories in the result and web directory.
If you provide assemblies with the -g option, make sure, that your assemblies are all stored as FASTA files in the given directory and that the prefix of the FASTA file matches your FASTQ files. For example: If you have a FASTA assembly file "040A11.fasta" in the assembly directory, you need a corresponding "040A11.fastq" file in the FASTQ directory.

Note that, if you run with -g option, no assemblies of filtered reads are computed.

Configuration
-------------
  
The bevaviour of the script can be modified in the file Config.py. First of all, make sure your references are defined (see next paragraph).

### Reference sequences ###
If you are not using our T-DNA and reference sequence file above, you will need to add the names of your references in the dictionary in Config.py::

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

This is to ensure, that each annotated sequence gets the correct colour. Unknown features will get a gray colour. If You need to add all sequence names from your fasta annotation file. For example, if there is a sequence called "genome" which is your plant genome sequence, you need to add "genome": "plant" in this dictionary. If you want to use different colours then the predefined, you can add it in this dictionary::


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
        "other: (155,155,155,1)
    }

... and use its key in the BLAST_FEATURE_MAPPING array.

### Tool behaviour ###

The behaviour of the tool can be modified in several ways.

You can set the verbose mode and the number of CPUs/threads, that the tool can use::

    BE_VERBOSE = True;
    NO_CPUS = 64;
    
The BLAST parameters, if you do not find what you expect. However, the default values should work fine for most cases::
    
    BLAST_PARAMS = " -perc_identity 85 -evalue 1e-50"

If you are working with something else then *A. thaliana*, you should change the genome size::

    GENOME_SIZE=145000000

If you want to analyse the reads that do not contain T-DNA sequences you can choose, that the script writes such a file (ending with .notdna.fastq)::

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


### Run in a docker container

There is a Dockerfile available, so you can run the tool inside a docker container as well. To do so, first install docker:

```bash
sudo apt install docker.io
```
You will need the docker file, so just clone the repository:

```bash
git clone https://github.com/nkleinbo/loreta
```
Then build the docker image:
```bash
cd loreta
sudo docker build -t loreta .
```

Then start the tool, you will need to mount all local directories that contain data into the container. In this example, all our data is in ~/data/ and is mounted to /data in the docker image:
```bash
sudo docker run --name loreta_container -v ~/data/:/data/ --rm -i -t loreta bash
```
This will start an interactive shell inside the container with your ~/data/ (containing a folder "input_data" with your input data) directory mounted to /data. You can then go on running the tool:
```bash
cd loreta
python3 run_all.py -f /data/input_data/ -o /data/docker_results/ -t references/tdna.fas -a references/all_fasta.fas -w /data/docker_html/
```

You could do the sam non interactively with:
```bash
sudo docker run --name loreta_container -v ~/data/:/data/ --rm -i -t loreta python3 loreta/run_all.py -f /data/input_data/ -o /data/docker_results/ -t loreta/references/tdna.fas -a loreta/references/all_fasta.fas -w /data/docker_html/
```

