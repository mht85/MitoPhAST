#!/usr/bin/perl

##########################
#Email mun.tan@deakin.edu.au
#Github mht85
#Date 08/08/17
#This script automates the pipeline and all scripts according to user configuration
#Example ./run_pipeline.pl -in mitofiles -acc NC_001645.1,NC_008795.1 -out analysis_out --summarize --build_trees -mito_gcode 1 -mode pcg -chars nt --gene_order -ori COX1
##########################

use strict;
use warnings;
use Getopt::Long;
use Cwd 'abs_path';
use Bio::SeqIO;


### USAGE STATEMENT ###
my ($indir, $inacc, $inlist, $outdir, $run_summarize, $run_extract, $run_trees, $run_orders, $mito_gcode, $mode, $chars, $custom_pcg, $custom_rrna, $custom_otherpcg, $custom_otherpcg_gcode, $custom_othernonpcg, $ori, $cpu, $usage, $help, $bb, $bbrep, $alrt, $alrtrep, $interactive);

$cpu = 1;
$mode = "pcg";
$chars = "aa";
$bb = "on";
$bbrep = "1000";
$alrt = "on";
$alrtrep = "1000";
$interactive = "on";

GetOptions ("in=s" => \$indir,
	    "acc=s" => \$inacc,
	    "list=s" => \$inlist,
	    "out=s" => \$outdir,
	    "summarize" => \$run_summarize,
	    "extract_only" => \$run_extract,
	    "build_trees" => \$run_trees,
	    "gene_order" => \$run_orders,
	    "mito_gcode=s" => \$mito_gcode,
	    "mode=s" => \$mode,
	    "chars=s" => \$chars,
	    "bb=s" => \$bb,
	    "bbrep=s" => \$bbrep,
	    "alrt=s" => \$alrt,
	    "alrtrep=s" => \$alrtrep,
	    "pcg=s" => \$custom_pcg,
	    "rrna=s" => \$custom_rrna,
	    "other_pcg=s" => \$custom_otherpcg,
	    "other_pcg_gcode=s" => \$custom_otherpcg_gcode,
	    "other_nonpcg=s" => \$custom_othernonpcg,
	    "ori=s" => \$ori,
	    "cpu=s" => \$cpu,
	    "interactive=s" => \$interactive,
	    "h" => \$help,
	    "help" => \$help);

$usage =
"Usage example: ./run_pipeline.pl -in data -out analysis_out --summarize --build_trees -mito_gcode 5 -mode custom -pcg all -other_pcg H3,H4 -other_pcg_gcode 1 -chars nt --gene_order -ori COX1

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
			-pcg <string>		if -mode 'custom', list PCG genes (comma-separated) for phylogenetic analysis (or '-pcg all' for all mitogenome PCGs)
			-rrna <string>		if -mode 'custom', list rRNA genes (comma-separated) for phylogenetic analysis (or '-rrna all' for all mitogenome rRNAs)
			-other_pcg <string>	if -mode 'custom', list other non-mitogenome PCGs (comma-separated) for phylogenetic analysis (see README)
			-other_pcg_gcode <int>	genetic code for genes listed at -other_pcg (e.g. 1=Standard, 2=Vertebrate Mitochondrial, 5=Invertebrate Mitochondrial, etc)
			-other_nonpcg <string>	if -mode 'custom', list other non-mitogenome non-PCGs (comma-separated) for phylogenetic analysis (see README)
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
	-h/help			for this help message";

if ($help) {
        die ("\n$usage\n\n");
}

########## CHECK ARGUMENTS ##############################################################################################

unless (($indir || $inacc || $inlist) && $outdir && ($run_summarize || $run_extract || $run_trees || $run_orders)) {
        die ("\nError: missing arguments, ./run_pipeline.pl -h for usage\n\n");
}

if ($run_trees) {
	unless ($mito_gcode) {
		die ("\nError: set -mito_gcode to run --build_trees\n\n");
	}
	if ($mode eq "custom") {
		unless ($custom_pcg || $custom_rrna || $custom_otherpcg || $custom_othernonpcg) {
			die ("\nError: must have at least ONE of -pcg, -rrna, -other_pcg or -other_nonpcg to run -mode custom\n\n");
		}
		if ($custom_otherpcg && !defined ($custom_otherpcg_gcode)) {
			if ($chars eq "aa") {
				$custom_otherpcg_gcode = 1;
			} else {
				die ("\nERROR: genetic code needs to be set for genes listed in -other_pcg\n\n");
			}
		}
	}
}

if ($run_orders) {
	unless ($ori) {
		die ("\nError: set -ori to run --gene_order\n\n");
	}
}


########## SET PATHS ####################################################################################################
my ($bindir, $maindir, $dbdir, $exedir, $rdir);

($bindir = abs_path($0)) =~ s/run_pipeline.pl/bin/g;
($maindir = $bindir) =~ s/bin//g;
($dbdir = $bindir) =~ s/bin/db/g;
($exedir = $bindir) =~ s/bin/exe/g;
$rdir = "$exedir/Rpackage";

########## SET OUTPUT DIRECTORY #########################################################################################
my ($ans, $resultdir, $workdir);

if (-e $outdir && -d $outdir) {
	print "Output folder already exists. Would you like to overwrite the whole directory? (Y/N)";
	$ans = <STDIN>;
	chomp $ans;
	if ($ans !~ /^Y$/i) {
		exit;
	}
}
system("rm -rf $outdir");
$resultdir = "$outdir/results";
$workdir = "$outdir/work";
mkdir $outdir;
mkdir $resultdir;
mkdir $workdir;

########## PREPARE INPUT ################################################################################################
my ($mitodir, $fh1, @filecount, $init_count, $dldir, @custom_otherpcglist, @custom_othernonlist);

$mitodir = "$resultdir/1.mitofiles";
mkdir $mitodir;

if ($inacc || $inlist) {
	print "\n### MITOGENOME COMPILATION ###\n";
	$dldir = "$workdir/downloads";
	mkdir "$dldir";
	if ($inacc) {
		if ($inlist) {
			print "Command: $bindir/get_ncbi.pl -acc $inacc -list $inlist -out $dldir\n\n";
			system("$bindir/get_ncbi.pl -acc $inacc -list $inlist -out $dldir");
		} else {
			print "Command: $bindir/get_ncbi.pl -acc $inacc -out $dldir\n\n";
			system("$bindir/get_ncbi.pl -acc $inacc -out $dldir");
		}
	} elsif ($inlist) {
		print "Command: $bindir/get_ncbi.pl -list $inlist -out $dldir\n\n";
		system("$bindir/get_ncbi.pl -list $inlist -out $dldir");
	}
	system("cp $dldir/1.downloads/* $mitodir/");
	if (! -e "$dldir/download.success") {
		die ("\nError: Pipeline terminated - download step failed\n");
	}
}

if ($indir) {
	system("cp $indir/* $mitodir/");
}

opendir(INDIR, $mitodir) or die $!;
while ($fh1 = readdir(INDIR)) {
	if ($fh1 =~ /\.(gb|gbk|gbff|gbf|embl|eml)$/) {
		push (@filecount, $fh1);
	}
}
closedir(INDIR);
$init_count = scalar@filecount;

print "\n[ $init_count ] mitogenome files can be found in $mitodir\n";

@filecount = ();

########## SUMMARIZE ####################################################################################################
my ($sumworkdir, $sumresdir, $new_count);

if ($run_summarize) {
	$sumworkdir = "$workdir/summaries";
	$sumresdir = "$resultdir/2.summaries";
	mkdir $sumworkdir;
	mkdir $sumresdir;
	print "\n### SUMMARIZE ###\n";
	print "Command: $bindir/summarize.pl -in $mitodir -out $sumworkdir\n\n";
	system("$bindir/summarize.pl -in $mitodir -out $sumworkdir");
	system("cp $sumworkdir/1.summaries/*.summary.txt $sumresdir/");
	system("cp $sumworkdir/1.summaries/summary_report.txt $sumresdir/");
	opendir(INDIR, $sumresdir) or die $!;
	while ($fh1 = readdir(INDIR)) {
		if ($fh1 =~ /\.summary\.txt$/) {
			push (@filecount, $fh1);
		}
	}
	closedir(INDIR);
	$new_count = scalar@filecount;
	print "\n[ $new_count ] summary files can be found in $sumresdir\n";
	if (! -e "$sumworkdir/summary.success") {
		die ("\nError: Pipeline terminated - summary step failed\n");
	}
#	if ($new_count != $init_count) {
#		die ("\nError: Pipeline terminated - failed to generate summaries for every mitogenome file\n\n");
#	}
	@filecount = ();
	$new_count = "";
}
if (! defined ($run_extract || $run_trees || $run_orders)) {
	printEnd();
	exit;
}

########## EXTRACT SEQUENCES ############################################################################################
my ($seqworkdir, $seqresdir, @custom_othernonpcg, $custom_othernonpcgfile, @custom_otherpcg, $custom_otherpcgfile, $gene, @seqcount, $seqcount);

if ($run_extract || $run_trees || $run_orders) {
	$seqworkdir = "$workdir/extract";
	$seqresdir = "$resultdir/3.sequences";
	mkdir $seqworkdir;
	mkdir $seqresdir;
	print "\n### SEQUENCE EXTRACT ###\n";
	print "Command: $bindir/extract_seqs.pl -in $mitodir -out $seqworkdir -cpu $cpu -interactive $interactive\n\n";
	system("$bindir/extract_seqs.pl -in $mitodir -out $seqworkdir -cpu $cpu -interactive $interactive");
	system("cp $seqworkdir/1.sequences/all_*.aa.fa $seqresdir");
	system("cp $seqworkdir/1.sequences/all_*.nt.fa $seqresdir");
	if ($custom_otherpcg) {
		@custom_otherpcg = split(/\,/, $custom_otherpcg);
		foreach $custom_otherpcgfile (@custom_otherpcg) {
			system("cp $indir/all_$custom_otherpcgfile.nt.fa $seqresdir");
			system("cp $indir/all_$custom_otherpcgfile.aa.fa $seqresdir");
		}
	}
	if ($custom_othernonpcg) {
		@custom_othernonpcg = split(/\,/, $custom_othernonpcg);
		foreach $custom_othernonpcgfile (@custom_othernonpcg) {
			system("cp $indir/all_$custom_othernonpcgfile.nt.fa $seqresdir");
		}
	}
	opendir(INDIR, $seqresdir) or die $!;
	while ($fh1 = readdir(INDIR)) {
		if ($fh1 =~ /\.aa\.fa/ || $fh1 =~ /\.nt\.fa/) {
			($gene = $fh1) =~ s/(all_|.aa.fa|.nt.fa)//g;
			open INFILE, "$seqresdir/$fh1" or die $!;
			while (<INFILE>) {
				if ($_ =~ /^>/) {
					push (@seqcount, $_);
				}
			}
			$seqcount = scalar@seqcount;
			print "[ $seqcount ] sequences extracted for gene $gene\n";
#			if ($seqcount != $init_count) {
#				print "WARNING: The pipeline failed to extract sequences for gene $gene for every mitogenome file.\n";
#				print "Would you like to proceed? (Y/N)";
#				$ans2 = <STDIN>;
#				if ($ans2 !~ /^Y$/i) {
#					exit;
#				}
#			}
		}
		@seqcount = ();
		$seqcount = "";
	}
	print "\nExtracted sequences can be found in $seqresdir\n";
	if (! -e "$seqworkdir/extract.success") {
		die ("\nError: Pipeline terminated - sequence extraction step failed\n");
	}
}
if (! defined ($run_trees || $run_orders)) {
	printEnd();
	exit;
}


########## PHYLOGENETICS ################################################################################################
my ($treeworkdir, $treeresdir, $custom_string);

if ($run_trees) {
	$treeworkdir = "$workdir/phylogenetics";
	$treeresdir = "$resultdir/4.trees";
	mkdir $treeworkdir;
	mkdir $treeresdir;
	print "\n### PHYLOGENETIC ANALYSIS ###\n";
	print "Mode: $mode\n";
	print "Characters: $chars\n";
	if ($mode =~ /^(pcg|pcg12s16s|12s16s|all)$/) {
		print "Command: $bindir/build_trees.pl -in $seqresdir -out $treeworkdir -mito_gcode $mito_gcode -mode $mode -chars $chars -bb $bb -bbrep $bbrep -alrt $alrt -alrtrep $alrtrep -cpu $cpu -interactive $interactive\n\n";
		system("$bindir/build_trees.pl -in $seqresdir -out $treeworkdir -mito_gcode $mito_gcode -mode $mode -chars $chars -bb $bb -bbrep $bbrep -alrt $alrt -alrtrep $alrtrep -cpu $cpu -interactive $interactive");
	} elsif ($mode =~ /custom/) {
		if ($custom_pcg) {
			$custom_string .= " -pcg $custom_pcg";
		}
		if ($custom_rrna) {
			$custom_string .= " -rrna $custom_rrna";
		}
		if ($custom_otherpcg) {
			$custom_string .= " -other_pcg $custom_otherpcg -other_pcg_gcode $custom_otherpcg_gcode";
		}
		if ($custom_othernonpcg) {
			$custom_string .= " -other_nonpcg $custom_othernonpcg";
		}
		print "Command: $bindir/build_trees.pl -in $seqresdir -out $treeworkdir -mito_gcode $mito_gcode -mode $mode $custom_string -chars $chars -bb $bb -bbrep $bbrep -alrt $alrt -alrtrep $alrtrep -cpu $cpu -interactive $interactive\n\n";
		system("$bindir/build_trees.pl -in $seqresdir -out $treeworkdir -mito_gcode $mito_gcode -mode $mode $custom_string -chars $chars -bb $bb -bbrep $bbrep -alrt $alrt -alrtrep $alrtrep -cpu $cpu -interactive $interactive");
	}
	system("cp $treeworkdir/2.buildtree/*.tre $treeresdir");
	system("cp $treeworkdir/2.buildtree/*.fasta $treeresdir");
	system("cp $treeworkdir/2.buildtree/*.phy $treeresdir");
	print "\nAlignment and tree files can be found in $treeresdir\n";
#	opendir(INDIR, $treeresdir) or die $!;
#	while ($fh1 = readdir(INDIR)) {
#		if ($fh1 =~ /\.tre$/) {
#			push (@filecount, $fh1);
#		}
#	}
#	closedir(INDIR);
#	$new_count = scalar@filecount;
#	if ($new_count == 0) {
#		die ("\nError: Pipeline terminated - phylogenetic analysis failed\n\n");
#	} else {
#		print "\nAlignment and tree files can be found in $treeresdir\n";
#	}
#	@filecount = ();
	if (! -e "$treeworkdir/phylogenetics.success") {
		die ("\nError: Pipeline terminated - tree building step failed\n");
	}
}
if (! defined ($run_orders)) {
	printEnd();
	exit;
}


########## GENE ORDER ###################################################################################################
my ($goresdir, $goworkdir);

if ($run_orders) {
	$goresdir = "$resultdir/5.gene_orders";
	$goworkdir = "$workdir/gene_orders";
	mkdir $goresdir;
	mkdir $goworkdir;
	system("cp $seqworkdir/1.sequences/sample_gene_order.txt $goresdir");
	print "\n### GENE ORDER ANALYSIS ###\n";
	print "Linearized gene orders re-oriented to begin with $ori gene\n";
	print "Command: $bindir/check_order.pl -in $goresdir/sample_gene_order.txt -out $goworkdir -ori $ori\n\n";
	system("$bindir/check_order.pl -in $goresdir/sample_gene_order.txt -out $goworkdir -ori $ori");
#	opendir(INDIR, $goworkdir) or die $!;
#	while ($fh1 = readdir(INDIR)) {
#		if ($fh1 =~ /order\./) {
#			push (@filecount, $fh1);
#		}
#	}
#	closedir(INDIR);
#	$new_count = scalar@filecount;
##	if ($new_count != 5) {
#	if ($new_count != 4) {
#		die ("\nError: Pipeline terminated - gene order comparison step failed\n\n");
#	} else {
	system("cp $goworkdir/1.orders/order.*.txt $goresdir");
	system("cp $goworkdir/1.orders/order.*.pdf $goresdir");
	print "Results from gene order analysis can be found in $goresdir\n";
#	}
#	@filecount = ();
#	$new_count = "";
	if (! -e "$goworkdir/order.success") {
		die ("\nError: Pipeline terminated - gene order step failed\n");
	}
}
printEnd();


### SUBROUTINE TO PRINT END MESSAGE

sub printEnd {
	print "\n\n### PIPELINE END ###\n";
	print "\nFinal result files are available in at $resultdir\n";
	print "\n### END ###\n";
}
