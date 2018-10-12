TITLE: MitoPhAST - Mitogenome Phylogenetic Analysis of Sequences Tool

=======
VERSION
=======
3.0

=======
CONTACT
=======
MunHua <mun.tan@deakin.edu.au>

============
INTRODUCTION
============
The increased rate at which complete mitogenomes are being sequenced and their increasing use for phylogenetic studies have resulted in a bioinformatics bottleneck in preparing and utilising such data for phylogenetic data analysis. We present, MitoPhAST, an automated tool that:
1) Identifies annotated protein-coding gene (PCG) and ribosomal RNA (rRNA) features directly from complete/partial GenBank/EMBL-format mitogenome files and generates a standardized, concatenated and partitioned nucleotide/amino acid alignments.
2) Incorporates nuclear gene sequences (if supplied by user) into phylogenetic analysis.
2) Generates a maximum likelihood phylogenetic tree using optimized substitution models performed by IQ-TREE.
3) Reports various mitochondrial genes and sequence information in a simple table format.
4) Extracts mitochondrial gene order (MGO) from GenBank files and clusters identical MGOs into groups

=============
PREREQUISITES
=============
Both install and automated run scripts were tested on Ubuntu 16.04 64-bit and MacOS 10.7.5 systems.
Please ensure systems are equipped with working versions of:

- Perl (http://www.perl.org/get.html)
- BioPerl (http://www.bioperl.org/wiki/Getting_BioPerl)
- Bio::DB::EUtilities
- LWP::Protocol::https
- R

============
INSTALLATION
============
1) Decompress the zipped folder, e.g.:
unzip MitoPhAST-master.zip
2) Change into the created directory, e.g.:
cd MitoPhAST-master
3) Run installation script, e.g.:
./install.sh
(for MAC users, run ./install_mac.sh instead)

=====
USAGE
=====
./run_pipeline.pl -h

Usage example: ./run_pipeline.pl -in mitofiles -acc NC_001645.1,NC_008795.1 -out analysis_out --summarize --build_trees -mito_gcode 1 -mode custom -other_pcg H3,H4 -other_pcg_gcode 1 -chars nt --gene_order -ori COX1

Options:

	INPUT (at least one REQUIRED):
	-in <folder>		input directory containing all GenBank and/or EMBL files
	-acc <string>		string of accession numbers to download from NCBI (comma-separated)
	-list <file>		file of accession numbers to download from NCBI (one accession/line)
	
	OUTPUT:
	-out <folder>		output directory to store all results and intermediate files (REQUIRED)
				
	ANALYSIS (at least one REQUIRED):
	--summarize		to generate summaries for each mitogenome
	--extract_only		to only extract and group sequences from GenBank and/or EMBL files (no need to invoke if running --build_trees)
	--build_trees		to run phylogenetic analysis
	--gene_order		to run gene order comparison analysis

	IF --build_trees:
		-mito_gcode <int>		genetic code for mitogenome sequences (e.g. 1=Standard, 2=Vertebrate Mitochondrial, 5=Invertebrate Mitochondrial, etc) (REQUIRED)
		-mode <string>			indicate genes used for phylogenetic analyis, default=pcg (modes = pcg/pcg12s16s/12s16s/all/custom)
		-pcg <string>			if -mode 'custom', list PCG genes (comma-separated) for phylogenetic analysis (or '-pcg all' for all mitogenome PCGs)
		-rrna <string>			if -mode 'custom', list rRNA genes (comma-separated) for phylogenetic analysis (or '-rrna all' for all mitogenome rRNAs)
		-other_pcg <string>		if -mode 'custom', list other non-mitogenome PCGs (comma-separated) for phylogenetic analysis (see README)
		-other_pcg_gcode <int>		genetic code for genes listed at -other_pcg (e.g. 1=Standard, 2=Vertebrate Mitochondrial, 5=Invertebrate Mitochondrial, etc)
		-other_nonpcg <string>		if -mode 'custom', list other non-mitogenome non-PCGs (comma-separated) for phylogenetic analysis (see README)
		-chars <string>			indicate type of character used for phylogenetic analysis, default=both (chars = aa/nt/both)
		-bb <on/off>			to run ultrafast bootstrap - UFBoot (Minh et al, 2013), default=on
		-bbrep <int>			number of ultrafast bootstrap replicates, default=1000
		-alrt <on/off>			to run SH-alrt test - SH-aLRT (Guindon et al, 2010), default=on
		-alrtrep <int>			number of SH-alrt replicates, default=1000

	IF --gene_order:
		-ori <string>			re-orient linear gene orders to this gene (REQUIRED)
																																											
	OTHER:
	-cpu <int>		number of threads, default=1
	-interactive <on/off>	control interactive messages (default=on, recommended)
	-h/help			for this help message

================
MITOGENOME INPUT
================
MitoPhAST automatically detects GenBank/EMBL files located in the input directory specified at -in.
Users can also provide comma-separated string (-acc) or a list (-list) of accession numbers to be downloaded from NCBI

=============
NUCLEAR INPUT
=============
Users can also perform phylogenetic analysis on combined mitogenome + nuclear datasets.

1) For nuclear PCGs (-other_pcg), have two fasta files per nuclear gene in the input directory (-in).
   The fasta file needs to be named all_<gene>.aa.fa or all_<gene>.nt.fa for amino acid and nucleotide sequences, respectively.
   For nuclear non-coding genes (-other_nonpcg), you will have only one fasta file per nuclear gene, named all_<gene>.nt.fa.

	   E.g. -other_pcg H3,H4 will need files all_H3.aa.fa, all_H3.nt.fa, all_H4.aa.fa, all_H4.nt.fa
	        -other_nonpcg 18S,28S will need files all_18S.nt.fa, all_28S.nt.fa

2) The identifier for each fasta entry will contain accession number and species name in the correct format (>Accession_species).
   This helps the pipeline match nuclear sequences to the correct mitogenome sequences extracted from GenBank/EMBL files.
	
	   E.g. in all_H3.aa.fa
		>NMVJ53370.1_Ibacus_alticrenatus
		MARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVKKPHRYRPGTVALREIRRYQKSTELLIRKLPFQRLVREIAQDFKTDLRFQSSAVMALQEASEAYLVGLFEDTNLCAIHAKRVTIMPKDIQLARRIRGERA
		>NMVJ55596.1_Puerulus_angulatus
		MARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVKKPHRYRPGTVALREIRRYQKSTELLIRKLPFQRLVREIAQDFKTDLRFQSSAVMALQEASEAYLVGLFEDTNLCAIHAKRVTIMPKDIQLARRIRGERA

3) Specify genetic code for nuclear PCGs (typically -other_pcg_gcode is set to 1 for the Standard code)

===============
MGO ORIENTATION
===============
To compare mitochondrial gene orders, linear orders need to be re-oriented to a common gene.
Users can specify the gene of choice at -ori with the following options:
ATP6, ATP8, COX1, COX2, COX3, CYTB, ND1, ND2, ND3, ND4, ND4L, ND5, ND6, rrnL, rrnS, trnA, trnC, trnD, trnE, trnF, trnG, trnH, trnI, trnK, trnL, trnL, trnM, trnN, trnP, trnQ, trnR, trnS, trnS, trnT, trnV, trnW, trnY

=================
BRIEF DESCRIPTION
=================

- Nucleotide and amino acid sequences for 13 PCG and 2 rRNA genes are extracted from each mitogenome GenBank/EMBL file.
- To verify annotation and to compensate for non-standard nomenclature in gene names in existing GenBank/EMBL files, extracted mitogenome PCG sequences are first searched against PFAM/PRINT profiles specific to each of the 13 PCGs. It is recommended that users manually check these files to ensure sequences are correctly assigned. Profiling has worked for the majority of the sequences we have tested. If you find a gene that is unassigned to any group, the workaround is to ensure the gene is named according to standards (see section "Gene ID list" for details). Users will be flagged by the pipeline if it encounters any anomalies.
- Nuclear PCG genes do not undergo any verification process.
- Verified nucleotide or amino acid sequences are then aligned using MAFFT and TranslatorX, respectively, both trimmed with Gblocks to remove ambiguously-aligned regions.
- 'Clean' alignments are concatenated with FASconCAT-G into supermatrices.
- Supermatrices, along with partition information, are provided to IQ-TREE for Maximum-likelihood analysis. PCG (amino acid) and non-coding (nucleotide) alignments are partitioned according to individual genes, whereas PCG nucleotide alignments are partitioned according to individual genes and codon positions (first, second, third).
- The resulting trees with SH-aLRT/ultrafast bootstrap support values constructed by IQ-TREE can be found in the output directory.
- The pipeline also retains mitochondrial gene order (MGO) information from GenBank/EMBL files and performs comparisons within the dataset, providing users fasta format files that can be used as input for CREx and TREx as well as MGOs visualized in PDF files to facilitate further comparative analysis.

======
OUTPUT
======
An output folder is automatically generated based on name specified at -out <folder>

+ <output directory>
  |- results
     |- 1.mitofiles
     |- 2.summaries
     |- 3.sequences
     |- 4.trees
     |- 5.gene_orders
  |- work
     |- extract
     |- gene_orders
     |- phylogenetics
     |- summaries


$PATH/<OUTPUT_DIR>/results/1.mitofiles/
---------------------------------
Contains all analyzed mitogenome files, including those downloaded from NCBI based on accession numbers provided to the -acc or -list options.

$PATH/<OUTPUT_DIR>/results/2.summaries/
---------------------------------
Contains summary reports generated for each GenBank/EMBL file. The number of PCGs, rRNAs, tRNAs and control regions are counted, followed by a table containing the names, strand, start/stop positions, lengths and intergenic nucleotides for each gene.

$PATH/<OUTPUT_DIR>/results/3.sequences/
----------------------------------
Contains Nucleotide and Amino acid sequences for 13 PCGs and 2 rRNAs compiled from analyzed mitogenome files.

$PATH/<OUTPUT_DIR>/results/4.trees/
----------------------------------
Contains Fasta and Phylip files used for phylogenetic tree construction. Also contains maximum likelihood trees (Newick) with SH-aLRT/UFboot support values at nodes for every alignment.

$PATH/<OUTPUT_DIR>/results/5.gene_orders
----------------------------------
Contains grouping of mitochondrial gene orders (MGOs) in FASTA and PDF output files

$PATH/<OUTPUT_DIR>/work/
----------------------------------
Contains intermediate files including extracted gene sequences from each mitogenome file, raw and trimmed multiple sequence alignments of each gene (nt and aa), output files from IQ-TREE, etc.

============
Gene ID list
============
In the event where a protein-coding gene sequence fails to be identified by its PFAM/PRINT domain, MitoPhAST relies on gene tags in the GenBank/EMBL file to group PCG sequences.
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

==================================
Versions of software/programs used
==================================
FASconCAT-G v1.02
FingerPRINTScan v3.596
Gblocks v0.91b
HMMER v3.1b2
IQ-TREE v1.5.5
MAFFT v7.394
BLAST+ v2.6.0

===========
GitHub link
===========
https://github.com/mht85/MitoPhAST
