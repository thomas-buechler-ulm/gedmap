# Preamble

GED-MAP is a prototype for pangemonic read mapping.
This software is part of the article *Efficient short read mapping to a pangenome that is represented by a graph of ED strings* by Thomas Büchler, Jannik Olbrich and Enno Ohlebusch. (Submitted to Bioinformatics)

# Requirements

A modern, C++11 ready compiler such as g++ version 4.9 or higher or clang version 3.2 or higher .

Installation of the SDSL library. (https://github.com/simongog/sdsl-lite)

# Installation

	git clone https://github.com/thomas-buechler-ulm/gedmap/
	cd gedmap
	make gedmap

First, if sdsl is not installed in your home directory, edit the PATHS file and set the correct path to sdsl-lite directory, before executing 'make'.
This produces the executable 'gedmap'.

## Short description

The binary 'gedmap' contains 4 programs:
- gedmap parse: Parses a FASTA and a VCF file into and GEDS (+ adjacency file)
- gedmap index: Calculated a minimizer index from of the GEDS (+ adjacency file) 
- gedmap align: Aligns reads of a FASTQ file to the GEDS with help of the index.
- gedmap sample: Sampled reads from a GEDS

Run 'gedmap parse -h', etc. to get the information about this program.
The parameter for all programs are also described on the bottom of this page.

## A little example
A little example is given in the directory with this name.  (It also contains a README file.)
To test your installation you can run:

	make a_little_example



# Requierements for the experiments

Installations of:
- hisat2 (version 2.2.1) from http://daehwankimlab.github.io/hisat2/manual/)
- samtools
- vg ( version v1.37.0-11-g2f6837d33) from https://github.com/vgteam/vg


# Rerun experiments

Set correct paths in the PATHS file. (Not nescessary if all software tools are installed in your home directory.)

	make data
	make all_eperiments

The first command downloads all human sequence data into the directory 'data'.
The second command runs the three mapping tools. 
The indexes and alignments will be calculated and stored in 'exp/gedmap',  'exp/hista2' and  'exp/giraffe'.
Samples will be generated and stored in 'exp/sample'.
The program 'eval_sam' and 'com_sams' (also generated by 'make all_eperiments') evaluate the sam files to obtain the results given in the paper.

At the end the directories 'exp/gedmap',  'exp/hista2' and  'exp/giraffe'  will contain following files: 'index.time.out', 'map.time.out' and 'eval.out'
The first and second are generated by the time program and contain the information about time and space consumption.
The files 'eval.out' and 'compare_sam_ged_hisat2' contain the information about the accuracy

# Parameter

### gedmap parse:
This  parses a FA file and a VCF file to an EDS graph.

3 Arguments expected:

	[1] filename of FA
	[2] filename of VCF
	[3] output filename

Optional parameters: 

	-nocnv, do not include CNV entries. (Then only and EDS is generated)
	-tmp tmp_dir , to set tmp direcoty (DEFAULT: tmp_dir=/tmp)

### gedmap index:           

This calculates a minimizer index of a given EDS graph.
1 Arguments expected:

	[1] filename of GEDS

Optional parameters: 

	-a fname    file name of the adijacency file
	-2fa fname  file of the 2fa file
	-k k        kmer/minimizer size k (DEFAULT: k=15)
	-w w        window size w (DEFAULT: w=7)
	-o fname    file name of minimizer index (DEFAULT: fname= geds_fname.min)
	-n n        maximum number of N in minimizer (DEFAULT: n=2)
	-tc x       maximum number of threads used (DEFAULT: use as many as avaiable)
	-trim x     removes sets greater than x from index (DEFAULT: x=1000) (use x=0 for no trim)
                                                       
### gedmap align:

This algings reads to the given GEDS and minimizer index

3 Arguments expected:

	[1] filename of GEDS
	[2] filename of FASTQ
	[3] filename of MINI
	
Optional parameters: 

	-o        fname, output will be stored in file fname (DEFAULT = [2].sam)
	-2fa     .2fa-file , if given this is used to transform GEDS-positions to FA positions
	-a fname  file name of the adijacency file
	-rc       reversed complement of pattern will be searched, too
	-oa       only aligned reads will be reported
	-io       output reads in the same order as in the input (may be a bit slower and with higher memory)
	-mc x     minimizer count, x minimizers will be looked up per read(DEFAULT: x=100)
	-ws x     window size of hotspot (DEFAULT: x=500)
	-wh x     minimum minimizer score (DEFAULT: x=1)
	-dd x     doubt distance, when best alignment has a distance x or greater, go to next round (DEFAULT: x=3)
	-mac x    max number of alignments completely calculated (DEFAULT: x=5)
	-mat x    max number of alignments tried to calculate (DEFAULT: x=10)
	-mao x    max number of alignments in output (DEFAULT: x=1)
	-d x      max distance in alignment (DEFAULT: x=50)
	-tmp tmp_dir, to change DEFAULT tmp direcoty from /tmp to tmp_dir
	-tc x     maximum number of threads used per index copy (DEFAULT: uses as many as avaiable)
                                   
### gedmap sample:
This samples reads from the given GEDS.
1 Arguments expected:

	[1] filename of GEDS

Optional parameters: 

	-c x   sample x reads (count, DEFAULT x=100)
	-l x   reads length = x (length, DEFAULT x=100)
	-e x   probability of a base to be false = 1/x, (error rate, DEFAULT: x=100)
	-e_d x probability of a base indel = 1/x, (insertion rate, DEFAULT: x=1000)
	-e_i x probability of a base indel = 1/x, (deletion rate, DEFAULT: x=1000)
	-rc x  probability of a read to be a reverse complement = 1/x, (rev complement rate, DEFAULT: x=2)
	-s x   seed for rng=x, (DEFAULT: x=0)
	-o fname       output will be written to file fname (DEFAULT:  fname=[1].sample.fastq)
	-2fa 2fa-file  if given this is used to transform GEDS-positions to FA positions
	-a fname       file name of the adijacency file
	-ae fname      file name of the adijacency file, only sample reads that go over edges in the graph

for all probabilities above: if x=0 then the probability is 0, other wise the probability is 1/x: 100% for x = 1, 10% for x=10, etc.


