#!/usr/bin/perl

##########################
#Email mun.tan@deakin.edu.au
#GitHub mht85
#Date 23/06/17
#This script runs nucleotide- and amino acid-based phylogenetic analyses using various genes
#Example ./build_trees.pl -in <folder> -out <folder> -mito_gcode <int> -pcg ATP6,ATP8,COX1,COX2,COX3,CYTB -rrna rrnS -other_pcg AAA,Bb -other_pcg_gcode 1 -other_nonpcg CcC,dD -chars nt -cpu 10
##########################

use strict;
use warnings;
use Getopt::Long;
use Cwd 'abs_path';
use Cwd;

### USAGE STATEMENT ###
my ($indir, $outdir, $usage, $help, $gcode, $mode, $pcg, $rrna, $otherpcg, $otherpcg_gcode, $othernonpcg, $char, $cpu, $bb, $bbrep, $alrt, $alrtrep, $interactive);
$mode = "pcg";
$cpu = "1";
$char = "aa";
$bb = "on";
$bbrep = "1000";
$alrt = "on";
$alrtrep = "1000";
$interactive = "on";

GetOptions ("in=s" => \$indir,
	    "out=s" => \$outdir,
	    "h" => \$help,
	    "help" => \$help,
	    "mito_gcode=s" => \$gcode,
	    "mode=s" => \$mode,
	    "pcg=s" => \$pcg,
	    "rrna=s" => \$rrna,
	    "other_pcg=s" => \$otherpcg,
	    "other_pcg_gcode=s" => \$otherpcg_gcode,
	    "other_nonpcg=s" => \$othernonpcg,
	    "chars=s" => \$char,
	    "bb=s" => \$bb,
	    "bbrep=s" => \$bbrep,
	    "alrt=s" => \$alrt,
	    "alrtrep=s" => \$alrtrep,
	    "cpu=s" => \$cpu,
	    "interactive=s" => \$interactive);

$usage =
"usage: ./build_trees.pl -in <folder> -out <folder> -mito_gcode <int> -mode custom -pcg ATP6,ATP8,COX1,COX2,COX3,CYTB -rrna rrnS -other_pcg AAA,Bb -other_pcg_gcode 1 -other_nonpcg CcC,dD -chars nt -cpu 10
options:
        -h/help			for this help message
        -in <folder>		input directory containing input data (one FASTA file per gene) (REQUIRED, see README for non-mitogenome genes)
        -out <folder>		output directory (REQUIRED)
	-mito_gcode <int>	genetic code (e.g. 1=Standard, 2=Vertebrate Mitochondrial, 5=Invertebrate Mitochondrial) (REQUIRED, but is only applied on mitogenome PCGs)
	-mode <mode>		indicate genes used for phylogenetic analyis, default=pcg (modes = pcg/pcg12s16s/12s16s/all/custom)
		-pcg <list>		if -mode 'custom', list PCG genes (comma-separated) for phylogenetic analysis (or '-pcg all' for all mitogenome PCGs)
		-rrna <list>		if -mode 'custom', list rRNA genes (comma-separated) for phylogenetic analysis (or '-rrna all' for all mitogenome rRNAs)
		-other_pcg <list>	if -mode 'custom', list other non-mitogenome PCGs (comma-separated) for phylogenetic analysis (see README)
		-other_pcg_gcode <int>	genetic code for genes listed at -other_pcg
		-other_nonpcg <list>	if -mode 'custom', list other non-mitogenome non-PCGs (comma-separated) for phylogenetic analysis (see README)
	-chars <chars>		indicate type of character used for phylogenetic analysis, default=aa (chars = aa/nt/both)
	-bb <on/off>		to run ultrafast bootstrap - UFBoot (Minh et al, 2013), default=on
	-bbrep <int>		number of ultrafast bootstrap replicates, default=1000
	-alrt <on/off>		to run SH-alrt test - SH-aLRT (Guindon et al, 2010), default=on
	-alrtrep <int>		number of SH-alrt replicates, default=1000
        -cpu <int>		number of threads, default=1
	-interactive <on/off>	control interactive messages (default=on, recommended)";

if ($help) {
	die ("\n$usage\n\n");
}

unless ($indir && $outdir && $gcode) {
	die ("\nERROR: missing arguments, ./build_trees.pl -h for usage\n\n");
}

if ($mode ne "custom" && ($pcg || $rrna || $otherpcg || $othernonpcg)) {
	die ("\nERROR: -pcg, -rrna, -other_pcg and -other_nonpcg options only available when -mode custom\n\n");
} elsif ($mode eq "custom" && (!defined ($pcg || $rrna || $otherpcg || $othernonpcg))) {
	die ("\nERROR: -mode custom but no PCG/rRNA/genes listed for analysis\n\n");
}

if ($otherpcg && (!defined $otherpcg_gcode)) {
	die ("\nERROR: genetic code needs to be set for genes listed in -other_pcg\n\n");
}

### SET PATHS ###
my ($bindir, $maindir, $dbdir, $exedir);

($bindir = abs_path($0)) =~ s/build_trees.pl//g;
($maindir = $bindir) =~ s/bin//g;
($dbdir = $bindir) =~ s/bin/db/g;
($exedir = $bindir) =~ s/bin/exe/g;


### SET OUTPUT DIRECTORY ###
my ($ans, $aligndir, $treedir);
#if (-e $outdir && -d $outdir) {
#        print "Output folder already exists. Would you like to overwrite the whole directory? (Y/N)";
#        $ans = <STDIN>;
#        chomp $ans;
#        if ($ans !~ /^Y$/i) {
#                exit 0;
#        }
#}
#system("rm -rf $outdir");
#mkdir $outdir;
mkdir "$outdir";
mkdir "$outdir/1.align";
mkdir "$outdir/2.buildtree";
$outdir = abs_path("$outdir");
$aligndir = abs_path("$outdir/1.align");
$treedir = abs_path("$outdir/2.buildtree");
mkdir "$aligndir/pcgaa";
mkdir "$aligndir/pcgnt";
mkdir "$aligndir/nonpcg";


### EXECUTE SUBROUTINES ###

my @pcglist = ("ATP6", "ATP8", "COX1", "COX2", "COX3", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6", "CYTB");
my @rrnlist = ("rrnL", "rrnS");
my (%pcglist, %rrnlist);
$pcglist{$_} = 1 foreach @pcglist;
$rrnlist{$_} = 1 foreach @rrnlist;

my (@configs, $genes_pcg, $genes_nucpcg, $genes_nonpcg, @genes_pcg, @genes_nucpcg, @genes_nonpcg, @chars, @nexus);
@configs = save_config();
print "Save run configurations: OK\n";
($genes_pcg, $genes_nucpcg, $genes_nonpcg) = check_genes();
@genes_pcg = @$genes_pcg;
@genes_nucpcg = @$genes_nucpcg;
@genes_nonpcg = @$genes_nonpcg;
print "Check genes: OK\n";
@chars = check_chars();
print "Check characters: OK\n";
check_files_pcg(@genes_pcg);
check_files_nucpcg(@genes_nucpcg);
check_files_nonpcg(@genes_nonpcg);
print "Check files: OK\n";
get_all_taxa();
check_taxa();
print "\nCheck taxa: OK\n";
print "Aligning sequences...";
align_pcg();
align_nonpcg();
print "done.\n";
print "Concatenating sequences...";
foreach (@configs) {
	@nexus = make_phylip($_);
}
print "done.\n";
print "Running IQTREE...\n\n";
foreach (@nexus) {
	run_iqtree($_);
}
#print "done.\n";
#print "Final constructed trees can be found in $outdir\n";
system("touch $outdir/phylogenetics.success");
####################################
## SUBROUTINES
####################################

# SUBROUTINE: Save configurations
my (@save_config);
sub save_config {
	if ($char eq "both") {
		push (@save_config, "$mode-aa");
		push (@save_config, "$mode-nt");
	} else {
		push (@save_config, "$mode-$char");
	}
	return @save_config;
}

# SUBROUTINE: Check gene names are OK and save list
my (@save_pcg, @save_nucpcg, @save_nonpcg, @tmp);
sub check_genes {
	if ($mode eq "all" || $mode eq "pcg12s16s") {
		@save_pcg = ("COX1", "COX2", "COX3", "ATP6", "ATP8", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6", "CYTB");
		@save_nonpcg = ("rrnS", "rrnL");
	} elsif ($mode eq "pcg") {
		@save_pcg = ("COX1", "COX2", "COX3", "ATP6", "ATP8", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6", "CYTB");
	} elsif ($mode eq "12s16s") {
		@save_nonpcg = ("rrnS", "rrnL");
	} elsif ($mode eq "custom") {
		if ($pcg) {
			if ($pcg eq "all") {
				push (@save_pcg, "COX1", "COX2", "COX3", "ATP6", "ATP8", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6", "CYTB");
			} else {
				@tmp = split(/\,/, $pcg);
				foreach (@tmp) {
					if (exists $pcglist{$_}) {
						push (@save_pcg, $_);
					} else {
						die ("\nERROR: -pcg $_ is not recognized\n\n");
					}
				}
				@tmp = ();
			}
		}
		if ($rrna) {
			if ($rrna eq "all") {
				push (@save_nonpcg, "rrnS", "rrnL");
			} else {
	                        @tmp = split(/\,/, $rrna);
				foreach (@tmp) {
					if (exists $rrnlist{$_}) {
						push (@save_nonpcg, $_);
					} else {
						die ("\nERROR: -rrna $_ is not recognized\n\n");
					}
				}
	                        @tmp = ();
			}
		}
		if ($otherpcg) {
			@tmp = split(/\,/, $otherpcg);
			foreach (@tmp) {
#				$_ =~ s/^.*\/all_//g;
#				$_ =~ s/.(nt|aa).fa//g;
				push (@save_nucpcg, $_);
			}
		}
		if ($othernonpcg) {
			@tmp = split(/\,/, $othernonpcg);
			foreach (@tmp) {
#				$_ =~ s/^.*\/all_//g;
#				$_ =~ s/.nt.fa//g;
				push (@save_nonpcg, $_);
			}
		}
	} else {
		die ("\nERROR: -mode set to unknown mode. Check usage below\n\n$usage\n\n");
	}
	return (\@save_pcg, \@save_nucpcg, \@save_nonpcg);
}

# SUBROUTINE: Check run characters are OK and save list
my (@save_chars);
sub check_chars {
	if ($mode eq "12s16s") {
		if ($char eq "both" || $char eq "aa") {
			die ("\nERROR: no amino acid information since -mode set to '$mode'\n\n");
		}
	} elsif ($mode eq "custom" && !defined($pcg) && !defined($otherpcg)) {
		if ($char eq "both" || $char eq "aa") {
			die ("\nERROR: no amino acid information since no PCGs listed at -pcg or -other_pcg\n\n");
		}
	} else {
		if ($char eq "both") {
			@save_chars = ("aa", "nt");
		} elsif ($char eq "aa") {
			@save_chars = ("aa");
		} elsif ($char eq "nt") {
			@save_chars = ("nt");
		} else {
			die ("\nERROR: -chars set to unknown mode. Check usage below\n\n$usage\n\n");
		}
	}
	return @save_chars;
}

# SUBROUTINE: Ensure necessary files are available for every gene
my (@tocheck, $gene, $eachchar, $filename, %filehash, $type, $filename_aa);
sub check_files_pcg {
	(@tocheck) = @_;
	foreach $gene (@tocheck) {
		foreach $eachchar (@chars) {
			$filename = "$indir/all_$gene.$eachchar.fa";
			report_files($filename, "pcg");
		}
	}
}

sub check_files_nucpcg {
	(@tocheck) = @_;
	foreach $gene (@tocheck) {
		foreach $eachchar (@chars) {
			$filename = "$indir/all_$gene.$eachchar.fa";
			report_files($filename, "nucpcg");
		}
	}
}

sub check_files_nonpcg {
	(@tocheck) = @_;
	foreach $gene (@tocheck) {
		$filename = "$indir/all_$gene.nt.fa";
		report_files($filename, "nonpcg");
	}
}

sub report_files {
	($filename, $type) = @_;
	if (! -e $filename) {
		die ("\nERROR: $filename does not exist!\n\n");
	} else {
		$filehash{$filename} = "$type";
		if ($char eq "nt" && ($type eq "pcg" || $type eq "nucpcg")) {
			($filename_aa = $filename) =~ s/.nt.fa/.aa.fa/g;
			if (! -e $filename_aa) {
				die ("\nERROR: $filename_aa does not exist! Pipeline needs an accompanying AA file even for NT analysis\n\n");
			}
		}
	}
}

# SUBROUTINE: Ensure each taxa is represented in every gene file
my ($file, %alltaxa);
sub get_all_taxa {
	foreach $file (keys %filehash) {
		open IN, $file;
		while (<IN>) {
			chomp;
			if ($_ =~ /^>/) {
				$_ =~ s/-(ATP|COX|ND|CYTB|rrn).*$//g;
				$_ =~ s/^>//g;
				if (! exists $alltaxa{$_}) {
					$alltaxa{$_} = 1;
				}
			}
		}
	}
}

my (%checktaxa, $taxa, @notaxa, %notaxa);
sub check_taxa {
	foreach $file (keys %filehash) {
		open IN, $file;
		while (<IN>) {
			chomp;
			if ($_ =~ /^>/) {
				$_ =~ s/-(ATP|COX|ND|CYTB|rrn).*$//g;
				$_ =~ s/^>//g;
				if (! exists $checktaxa{$_}) {
					$checktaxa{$_} = 1;
				}
			}
		}
		foreach $taxa (keys %alltaxa) {
			if (! exists $checktaxa{$taxa}) {
				$gene = $file;
				$gene =~ s/(.*\/all_|\.\w+.fa)//g;
				if (! exists $notaxa{"$taxa\t$gene"}) {
					push (@notaxa, "$taxa\t$gene");
					$notaxa{"$taxa\t$gene"} = 1;
				}
			}
		}
		%checktaxa = ();
	}
	if (scalar @notaxa > 0) {
		print "\n>>>WARNING: The following sequences are missing:\n";
		print "$_\n" foreach @notaxa;
		print "Sequences listed above are missing.";
		if ($interactive eq "on") {
			print "Proceed? (Y/N)";
			my $answer = <STDIN>;
			chomp $answer;
			if ($answer !~ /(Y|y)/) {
				exit 0;
			}
		}
	}
}

# SUBROUTINE: Align & trim PCG (aa, nt)
my (@pcg_files, @nucpcg_files, $file_aln, $seqcount, $minflankpos, $out_aln, $out_trim, $out_aln_aa, $file_aln_aa, $file_aln_aa_tmp, $file_aln_aa_stop, $seq);
sub align_pcg {
	@pcg_files = grep { $filehash{$_} eq 'pcg' } keys %filehash;
	@nucpcg_files = grep { $filehash{$_} eq 'nucpcg' } keys %filehash;
	foreach $file_aln (@pcg_files) {
		open ALN_FILE, $file_aln or die $!;
		while (<ALN_FILE>) {
			if ($_ =~ /^>/) {
				$seqcount += 1;
			}
		}
		$minflankpos = int(0.55*$seqcount + 0.5);
		($out_aln = (split(/\//, $file_aln))[-1]) =~ s/.fa$/.align/g;
		if ($file_aln =~ /.aa.fa$/) {
			system("$exedir/mafft-7.394-with-extensions/bin/mafft --localpair --thread $cpu --adjustdirection $file_aln > $aligndir/pcgaa/$out_aln 2>$aligndir/pcgaa/mafft.log");
			system("sed -i.bak 's/\_R\_//g' $aligndir/pcgaa/$out_aln");
			system("rm $aligndir/pcgaa/$out_aln.bak");
			system("$exedir/Gblocks_0.91b/Gblocks $aligndir/pcgaa/$out_aln -t=p -b2=$minflankpos -b3=8 -b4=5 -b5=h -e=.fas >> $aligndir/pcgaa/gblocks.log 2>&1");
			system("sed -i.bak 's/ //g' $aligndir/pcgaa/$out_aln.fas");
			system("rm $aligndir/pcgaa/$out_aln.fas.bak");
		} elsif ($file_aln =~ /.nt.fa$/) {
			($file_aln_aa = $file_aln) =~ s/.nt.fa/.aa.fa/g;
			($file_aln_aa_tmp = (split(/\//, $file_aln_aa))[-1]) =~ s/.fa$/.stop.fa/g;
			$file_aln_aa_stop = "$aligndir/pcgnt/$file_aln_aa_tmp";
			open FILE_ALN_AA, "$file_aln_aa" or die $!;
			open FILE_ALN_AA2, ">$file_aln_aa_stop" or die $!;
			while (<FILE_ALN_AA>) {
				chomp;
				if ($_ !~ /^>/) {
					$seq .= $_;
				} else {
					if ($seq) {
						print FILE_ALN_AA2 $seq,"X\n";
						undef $seq;
					}
					print FILE_ALN_AA2 "$_\n";
				}
			}
			print FILE_ALN_AA2 $seq,"X\n";
			undef $seq;
			($out_aln_aa = (split(/\//, $file_aln_aa_stop))[-1]) =~ s/.fa$/.align/g;
			system("$exedir/mafft-7.394-with-extensions/bin/mafft --localpair --thread $cpu --adjustdirection $file_aln_aa_stop > $aligndir/pcgnt/$out_aln_aa 2>$aligndir/pcgnt/mafft.log");
			system("sed -i.bak 's/\_R\_//g' $aligndir/pcgnt/$out_aln_aa");
			system("rm $aligndir/pcgnt/$out_aln_aa.bak");
			system("perl $exedir/TranslatorX_v1.1/translatorx_vLocal.edit.pl -i $file_aln -o $aligndir/pcgnt/$out_aln -a $aligndir/pcgnt/$out_aln_aa -c $gcode -p F -g \"-b2=$minflankpos -b3=8 -b4=5 -b5=h\" >> $aligndir/pcgnt/translatorx.log 2>&1");
#			system("perl $exedir/pal2nal.v14/pal2nal.pl $aligndir/pcgnt/$out_aln.aa_cleanali.fasta $aligndir/pcgnt/$out_aln.nt_cleanali.fasta -output fasta -nogap -codontable $gcode > $aligndir/pcgnt/$out_aln.nt_cleanali.fas 2> $aligndir/pcgnt/pal2nal.log");
			system("ln -s $out_aln.nt_cleanali.fasta $aligndir/pcgnt/$out_aln.nt_cleanali.fas");
		}
		$seqcount = "0";
		$minflankpos = "";
	}
	foreach $file_aln (@nucpcg_files) {
		open ALN_FILE, "$file_aln" or die $!;
		while (<ALN_FILE>) {
			if ($_ =~ /^>/) {
				$seqcount += 1;
			}
		}
		$minflankpos = int(0.55*$seqcount + 0.5);
		($out_aln = (split(/\//, $file_aln))[-1]) =~ s/.fa$/.align/g;
		if ($file_aln =~ /.aa.fa$/) {
			system("$exedir/mafft-7.394-with-extensions/bin/mafft --localpair --thread $cpu --adjustdirection $file_aln > $aligndir/pcgaa/$out_aln 2>$aligndir/pcgaa/mafft.log");
			system("sed -i.bak 's/\_R\_//g' $aligndir/pcgaa/$out_aln");
			system("rm $aligndir/pcgaa/$out_aln.bak");
			system("$exedir/Gblocks_0.91b/Gblocks $aligndir/pcgaa/$out_aln -t=p -b2=$minflankpos -b3=8 -b4=5 -b5=h -e=.fas >> $aligndir/pcgaa/gblocks.log 2>&1");
			system("sed -i.bak 's/ //g' $aligndir/pcgaa/$out_aln.fas");
			system("rm $aligndir/pcgaa/$out_aln.fas.bak");
		} elsif ($file_aln =~ /.nt.fa$/) {
			($file_aln_aa = $file_aln) =~ s/.nt.fa/.aa.fa/g;
			($file_aln_aa_tmp = (split(/\//, $file_aln_aa))[-1]) =~ s/.fa$/.stop.fa/g;
			$file_aln_aa_stop = "$aligndir/pcgnt/$file_aln_aa_tmp";
			open FILE_ALN_AA, "$file_aln_aa" or die $!;
			open FILE_ALN_AA2, ">$file_aln_aa_stop" or die $!;
			while (<FILE_ALN_AA>) {
				chomp;
				if ($_ !~ /^>/) {
					$seq .= $_;
				} else {
					if ($seq) {
						print FILE_ALN_AA2 $seq,"X\n";
						undef $seq;
					}
					print FILE_ALN_AA2 "$_\n";
				}
			}
			print FILE_ALN_AA2 $seq,"X\n";
			undef $seq;
			($out_aln_aa = (split(/\//, $file_aln_aa_stop))[-1]) =~ s/.fa$/.align/g;
			system("$exedir/mafft-7.394-with-extensions/bin/mafft --localpair --thread $cpu --adjustdirection $file_aln_aa > $aligndir/pcgnt/$out_aln_aa 2>$aligndir/pcgnt/mafft.log");
			system("sed -i.bak 's/\_R\_//g' $aligndir/pcgnt/$out_aln_aa");
			system("rm $aligndir/pcgnt/$out_aln_aa.bak");
			system("perl $exedir/TranslatorX_v1.1/translatorx_vLocal.edit.pl -i $file_aln -o $aligndir/pcgnt/$out_aln -a $aligndir/pcgnt/$out_aln_aa -c $gcode -p F -g \"-b2=$minflankpos -b3=8 -b4=5 -b5=h\" >> $aligndir/pcgnt/translatorx.log 2>&1");
#			system("perl $exedir/pal2nal.v14/pal2nal.pl $aligndir/pcgnt/$out_aln.aa_cleanali.fasta $aligndir/pcgnt/$out_aln.nt_cleanali.fasta -output fasta -nogap -codontable $gcode > $aligndir/pcgnt/$out_aln.nt_cleanali.fas 2> $aligndir/pcgnt/pal2nal.log");
			system("ln -s $out_aln.nt_cleanali.fasta $aligndir/pcgnt/$out_aln.nt_cleanali.fas");
		}
		$seqcount = "0";
		$minflankpos = "";
	}
}

# SUBROUTINE: Align & trim non-PCG (nt)
my (@nonpcg_files);
sub align_nonpcg {
	@nonpcg_files = grep { $filehash{$_} eq 'nonpcg' } keys %filehash;
	foreach $file_aln (@nonpcg_files) {
		open ALN_FILE, "$file_aln" or die $!;
		while (<ALN_FILE>) {
			if ($_ =~ /^>/) {
				$seqcount += 1;
			}
		}
		$minflankpos = int(0.55*$seqcount + 0.5);
		($out_aln = (split(/\//, $file_aln))[-1]) =~ s/.fa$/.align/g;
		$out_trim = "$out_aln.trim";
		system("$exedir/mafft-7.394-with-extensions/bin/mafft --localpair --thread $cpu --adjustdirection $file_aln > $aligndir/nonpcg/$out_aln 2> $aligndir/nonpcg/mafft.log");
		system("sed -i.bak 's/\_R\_//g' $aligndir/nonpcg/$out_aln");
		system("rm $aligndir/nonpcg/$out_aln.bak");
		system("$exedir/Gblocks_0.91b/Gblocks $aligndir/nonpcg/$out_aln -t=p -b2=$minflankpos -b3=8 -b4=5 -b5=h -e=.fas >> $aligndir/nonpcg/gblocks.log 2>&1");
		system("sed -i.bak 's/ //g' $aligndir/nonpcg/$out_aln.fas");
		system("rm $aligndir/nonpcg/$out_aln.fas.bak");
		$seqcount = "0";
		$minflankpos = "";
	}
}

# SUBROUTINE: Substitute _R_ attached to name when mafft adjusts direction
my ($subfile);
sub sub_direction {
	($subfile) = @_;
	open SUB, "$subfile" or die $!;
	open SUBOUT, ">$subfile.sub";
	while (<SUB>) {
		chomp;
		$_ =~ s/\_R\_//g;
		print SUBOUT "$_\n";
	}
	system("mv $subfile.sub $subfile");
	close(SUB);
	close(SUBOUT);
}

my ($subfile2);
sub sub_header {
	($subfile2) = @_;
	open SUB2, "$subfile2" or die $!;
	open SUB2OUT, ">$subfile2.sub";
	while (<SUB2>) {
		chomp;
		$_ =~ s/_1$//g;
		print SUB2OUT "$_\n";
	}
	system("mv $subfile2.sub $subfile2");
	close(SUB2);
	close(SUB2OUT);
}

# SUBROUTINE: Concatenate and convert to phylip format
my ($config, $mode2, $char2, $char3, $pcgphy, $pcgnex, $nonpcgphy, $nonpcgnex, $mixnex);
sub make_phylip {
	($config) = @_;
	($mode2, $char2) = split(/\-/, $config);
	if ($char2 eq "aa") {
		$char3 = "AA";
	} elsif ($char2 eq "nt") {
		$char3 = "DNA";
	}
	$pcgphy = "PCG$char2.phy";
	$pcgnex = "PCG$char2.nex";
	$nonpcgphy = "nonPCG.phy";
	$nonpcgnex = "nonPCG.nex";
	$mixnex = "PCG$char2"."_nonPCG.nex";
	if ($mode2 eq "all") {
		concat("$aligndir/pcg$char2/", "$pcgphy");
		concat("$aligndir/nonpcg/", "$nonpcgphy");
		print_header("$treedir/$pcgnex");
		get_pcg_parts("$aligndir/pcg$char2/FcC_info.xls", "$pcgphy", "$treedir/$pcgnex", "$char3");					# PCG (aa/nt)
		print_end("$treedir/$pcgnex");
		print_header("$treedir/$mixnex");
		get_pcg_parts("$aligndir/pcg$char2/FcC_info.xls", "$pcgphy", "$treedir/$mixnex", "$char3");					# PCG and nonPCG (aa/nt)
		get_nonpcg_parts("$aligndir/nonpcg/FcC_info.xls", "$nonpcgphy", "$treedir/$mixnex", "DNA");
		print_end("$treedir/$mixnex");
		print_header("$treedir/$nonpcgnex");
		get_nonpcg_parts("$aligndir/nonpcg/FcC_info.xls", "$nonpcgphy", "$treedir/$nonpcgnex", "DNA");					# nonPCG (nt)
		print_end("$treedir/$nonpcgnex");
		push (@nexus, "$pcgnex");
		push (@nexus, "$mixnex");
		push (@nexus, "$nonpcgnex");
	} elsif (($config =~ "^pcg-") || (($mode2 eq "custom") && (defined($pcg || $otherpcg)) && (! defined($rrna || $othernonpcg)))) {
		concat("$aligndir/pcg$char2/", "$pcgphy");
		print_header("$treedir/$pcgnex");
		get_pcg_parts("$aligndir/pcg$char2/FcC_info.xls", "$pcgphy", "$treedir/$pcgnex", "$char3");					# PCG (aa/nt)
		print_end("$treedir/$pcgnex");
		push (@nexus, "$pcgnex");
	} elsif (($config =~ "^pcg12s16s-") || (($mode2 eq "custom") && (defined($pcg || $otherpcg)) && (defined($rrna || $othernonpcg)))) {
		concat("$aligndir/pcg$char2/", "$pcgphy");
		concat("$aligndir/nonpcg/", "$nonpcgphy");
		print_header("$treedir/$mixnex");
		get_pcg_parts("$aligndir/pcg$char2/FcC_info.xls", "$pcgphy", "$treedir/$mixnex", "$char3");					# PCG and nonPCG (aa/nt)
		get_nonpcg_parts("$aligndir/nonpcg/FcC_info.xls", "$nonpcgphy", "$treedir/$mixnex", "DNA");	
		print_end("$treedir/$mixnex");
		push (@nexus, "$mixnex");
	} elsif (($config =~ "^12s16s") || (($mode2 eq "custom") && (!defined($pcg && $otherpcg)) && (defined($rrna || $othernonpcg)))) {
		concat("$aligndir/nonpcg/", "$nonpcgphy");
		print_header("$treedir/$nonpcgnex");
		get_nonpcg_parts("$aligndir/nonpcg/FcC_info.xls", "$nonpcgphy", "$treedir/$nonpcgnex", "DNA");					# nonPCG (nt)
		print_end("$treedir/$nonpcgnex");
		push (@nexus, "$nonpcgnex");
	}
	return @nexus;
}

my ($target_dir, $phylipfile, $fastafile);
sub concat {
	($target_dir, $phylipfile) = @_;
	($fastafile = $phylipfile) =~ s/phy/fasta/g;
	chdir("$target_dir") or die $!;
	if (! -e "FcC_supermatrix.phy" || ! -e "FcC_supermatrix.fas") {
		system("perl $exedir/FASconCAT-G_v1.02/FASconCAT-G_v1.02.pl -p -p -s > fasconcat.log 2>&1");
		if (! -e "FcC_supermatrix.phy" || ! -e "FcC_supermatrix.fas") {
			die ("\nAn error has occurred during concatenation, run terminated\n")
		}
	}
	system("cp FcC_supermatrix.phy $treedir/$phylipfile");
	system("cp FcC_supermatrix.fas $treedir/$fastafile");
	chdir("../../../../");
}


my ($outpart);
sub print_header {
	($outpart) = @_;
	open OUT, ">$outpart";
	print OUT "#NEXUS\nbegin sets;\n";
}

my ($inpart, $phyfile, $chartype, $name, $c1_start, $c2_start, $c3_start, @c1, @c2, @c3, $c1, $c2, $c3, $start, $stop);
sub get_pcg_parts {
	($inpart, $phyfile, $outpart, $chartype) = @_;
	open PART1, "$inpart" or die $!;
	open OUT, ">>$outpart";
	while (<PART1>) {
		chomp;
		if ($_ =~ /^all_/) {
			$name = (split(/\t/, $_))[0];
			$name = (split(/\./, $name))[0];
			$name =~ s/all_//g;
			if ($chartype eq "DNA") {
				$c1_start = (split(/\t/, $_))[1];
				$c2_start = $c1_start + 1;
				$c3_start = $c2_start + 1;
				$stop = (split(/\t/, $_))[2];
				print OUT "\t\tcharset $name.c1=$phyfile:$chartype, $c1_start-$stop\\3;\n";
				print OUT "\t\tcharset $name.c2=$phyfile:$chartype, $c2_start-$stop\\3;\n";
				print OUT "\t\tcharset $name.c3=$phyfile:$chartype, $c3_start-$stop\\3;\n";
			} elsif ($chartype eq "AA") {
				$start = (split(/\t/, $_))[1];
				$stop = (split(/\t/, $_))[2];
				print OUT "\t\tcharset $name=$phyfile:$chartype, $start-$stop;\n";
			}
		}
	}
	close(PART1);
}

sub get_nonpcg_parts {
	($inpart, $phyfile, $outpart, $chartype) = @_;
	open PART1, "$inpart" or die $!;
	open OUT, ">>$outpart";
	while (<PART1>) {
		chomp;
		if ($_ =~ /^all_/) {
			$name = (split(/\t/, $_))[0];
			$name = (split(/\./, $name))[0];
			$name =~ s/all_//g;
			$start = (split(/\t/, $_))[1];
			$stop = (split(/\t/, $_))[2];
			print OUT "\t\tcharset $name=$phyfile:$chartype, $start-$stop;\n";
		}
	}
	close(PART1);
}

sub print_end {
	open OUT, ">>$outpart";
	print OUT "\nend;\n";
}

# SUBROUTINE: Partition-finding and ML analysis with IQ tree
my ($nexus, $treefile, $cpurun, $lines);
sub run_iqtree {
	($nexus) = @_;
	($treefile = $nexus) =~ s/.nex/.tre/g;
	chdir("$treedir");
#	system("cp ../*.phy .");
	if ($nexus =~ /^PCGnt/) {
		$cpurun = scalar@save_pcg + scalar@save_nucpcg + scalar@save_nonpcg;
	} else {
		open NEX, $nexus or die $!;
		$cpurun = "0";
		while (<NEX>) {
			if ($_ =~ /charset/) {
				$cpurun++;
			}
		}
	}
	if ($cpu >= $cpurun) {
		$cpu = $cpurun;
	}
	if ($bb eq "on") {
		if ($alrt eq "on") {
			system("$exedir/iqtree-omp-1.5.5/bin/iqtree-omp -spp $nexus -nt $cpu -bb $bbrep -alrt $alrtrep -m TESTMERGE > iqtree.log 2>&1");	
		} else {
			system("$exedir/iqtree-omp-1.5.5/bin/iqtree-omp -spp $nexus -nt $cpu -bb $bbrep -m TESTMERGE > iqtree.log 2>&1");
		}
	} elsif ($bb eq "off") {
		if ($alrt eq "on") {
			system("$exedir/iqtree-omp-1.5.5/bin/iqtree-omp -spp $nexus -nt $cpu -alrt $alrtrep -m TESTMERGE > iqtree.log 2>&1");
		} else {
			system("$exedir/iqtree-omp-1.5.5/bin/iqtree-omp -spp $nexus -nt $cpu -m TESTMERGE > iqtree.log 2>&1");
		}
	}
	if (-e "$nexus.treefile") {
		system("ln -s $nexus.treefile $treefile");
		print "$treefile constructed\n";
	} else {
		exit;
	}
	chdir("$outdir");
}





