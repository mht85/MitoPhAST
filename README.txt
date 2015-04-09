TITLE: MitoPhAST - Mitogenome Phylogenetic Analysis of Sequences Tool


VERSION
=======
1.0


CONTACT
=======
Tan Mun Hua <tan.mun.hua@monash.edu>


INTRODUCTION
============
The increased rate at which complete mitogenomes are being sequenced and their increasing use for phylogenetic studies have resulted in a bioinformatics bottleneck in preparing and utilising such data for phylogenetic data analysis. We present, MitoPhAST, an automated tool that:
1) Identifies annotated protein-coding gene features and generates a standardized, concatenated and partitioned amino acid alignment directly from complete/partial GenBank/EMBL-format mitogenome flat files
2) Generates a maximum likelihood phylogenetic tree using optimized protein models
3) Reports various mitochondrial genes and sequence information in a simple table format.


PREREQUISITES
=============
Both install and automated run scripts were tested on Ubuntu 12.04 64-bit and MacOS 10.7.5 systems.
Please ensure systems are equipped with working versions of:

- Perl (http://www.perl.org/get.html)
- BioPerl (http://www.bioperl.org/wiki/Getting_BioPerl)
- Python (https://www.python.org)
- BioPython (http://biopython.org/wiki/Main_Page)


INSTALLATION
============
1) Decompress the zipped folder, e.g.:
unzip MitoPhAST-master.zip
2) Change into the created directory, e.g.:
cd MitoPhAST-master
3) Run installation script, e.g.:
./install.sh
(for MAC users, run ./install_mac.sh instead)

MitoPhAST compiles raxmlHPC-PTHREADS-SSE3 by default. To install a different flavor of RAxML, manually compile the target flavor and soft-link the binary file.
Example:
make -f Makefile.SSE3.gcc && rm *.o
ln -s raxmlHPC-SSE3 raxml


USAGE
=====

Example: ./run_MitoPhASTv1.sh -i example/input/ -n HG942366,HG799086 -m -I -B 100 -R 50 -T 15

-i <path> input directory containing Genbank and/or EMBL files [at least one of -i and -n is required]
-n <string>       comma-separated string of mitogenome accession IDs to download from NCBI [at least one of -i or -n is required]
-m                optional argument to exclude genes which are absent in some files (default: include all genes, missing genes are substituted with hyphens)
-I                turns off interactive prompt for user to check for missing genes (default: interactive) If prompt is disabled, program will only provide warning to standard output
-S                stops program after supermatrix construction, turns off model estimation and ML analysis
-B <int>          number of bootstrap replicates (default: autoMRE)
-R <int>          number of ML trees generated if standard bootstrapping is used - one tree out of <int> with best likelihood score is produced (default: 1)
-T <int>          number of threads (default: 2, must be greater than 1)
-r                run ML analysis with rapid bootstrapping (default: standard bootstrapping)

OUTPUT
======
An output folder is automatically generated with the input folder as prefix (e.g. $PATH/example/input-out).

+ input-out
  |- download
  |- input_tmp
  |- phylogeny
     |- align_trim
     |- ml_analysis
     |- model_select
  |- profiles
     |- hmmer_out
     |- fpscan_out
  |- sequences
  |- summaries


Briefly:
 
Extracted sequences can be found in the sequences folder with .cds.fa suffixes.
FASTA and PHYLIP files after alignment can be found in the align_trim folder as all_13pcg.fasta and all_13pcg.phy.
Final tree from RAxML analysis with standard/rapid bootstrapping can be found in the ml_analysis folder in the MLtree.renamed.tre.


In detail:

$PATH/example/input-out/download/
---------------------------------
Mitogenome files are downloaded from NCBI to this folder based on accession numbers provided to the -n option.

$PATH/example/input-out/input_tmp/
----------------------------------
Genbank/EMBL files are copied to this folder and worked on to ensure no modifications were made to the original files.

$PATH/example/input-out/summaries/
----------------------------------
A summary report is generated for each GenBank/EMBL file. The number of PCGs, rRNAs, tRNAs and control regions are counted, followed by a table containing the names, strand, start/stop positions, lengths and intergenic nucleotides for each gene.

$PATH/example/input-out/sequences/
----------------------------------
Amino acid sequences for all 13 mitochondrial PCGs are extracted from each GenBank/EMBL file.

$PATH/example/input-out/profiles/
---------------------------------
To compensate for non-standard nomenclature in gene names in Genbank/EMBL files, extracted amino acid sequences are assigned as specific PCGs based on PFAM/PRINT profiles specific to each PCG. It is recommended that users manually check these files to ensure sequences are correctly assigned. Profiling has worked for majority of the sequences we have tested. If you encounter a gene which is not assigned properly, the workaround is to ensure the gene is named properly according to standards (e.g. ATP6, ATP8, COX1, COX2, COX3, CYTB, ND1, ND2, ND3, ND4, ND4L, ND5, ND6).

$PATH/example/input-out/phylogeny/align_trim/
---------------------------------------------
Sequences from each organism are grouped according to the 13 mitochondrial PCGs. Sequences within each PCG group are aligned with Clustal-Omega and ambiguously-aligned regions are trimmed with trimAl (-automated1 option). If the -no-missing option is specified, any PCG group with sequences absent in one or more organisms are completely excluded from the phylogenetic analysis. PCG groups which are completely empty after alignment and trimming are also excluded. Remaining PCG sequences are concatenated for each organism. Final FASTA and PHYLIP files are generated. A resulting partition file is generated for use by RAxML.

$PATH/example/input-out/phylogeny/model_select
--------------------------------------------
ProtTest is used to find the best-fitting model for the alignment.

$PATH/example/input-out/phylogeny/ml_analysis/
----------------------------------------------
RAxML is used to run a maximum-likelihood analysis on the resulting PHYLIP file with PCGs. Analysis is run with standard/rapid bootstrapping. The resulting tree can be found in the MLtree*.renamed.tre files.


Gene ID list
============
In the event where a protein sequence fails to be identified by its PFAM/PRINT domain, MitoPhAST relies on gene tags in the GenBank/EMBL file to group protein sequences.
The following are gene IDs recognized by MitoPhAST (case insensitive):
ATP6 / ATPASE6 / ATPASE_6
ATP8 / ATPASE8 / ATPASE_8
COX1 / COI / CO1
COX2 / COII / CO2
COX3 / COIII / CO3
COB / CYTB / CYT_B
NAD1 / ND1 / NADH1
NAD2 / ND2 / NADH2
NAD3 / ND3 / NADH3
NAD4 / ND4 / NADH4
NAD4L / ND4L / NADH4L
NAD5 / ND5 / NADH5
NAD6 / ND6 / NADH6


Versions of software/programs used
==================================
Argtable2 v2.13
Clustal Omega v1.2.0
FingerPRINTScan v3.596
HMMER v3.0
ProtTest3 v3.4
RAxML v8.0.25
trimAl v1.2rev59


GitHub link
===========
https://github.com/munhua/MitoPhAST
