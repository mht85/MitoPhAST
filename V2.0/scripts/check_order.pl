#!/usr/bin/perl

##########################
#Email mun.tan@deakin.edu.au
#GitHub mht85
#Date 23/06/17
#This script runs compared gene orders and reports shared vs unique gene orders
#Example ./check_order.pl -in <ffile> -out <folder> -ori <gene to orient to>
##########################

use strict;
use warnings;
use Getopt::Long;
use Cwd 'abs_path';

# Input FASTA format:
# >SPECIES|ACCESSION
# GENE_ORDER (tab-separated)

### USAGE STATEMENT ###
my ($infile, $outdir, $help, $ori, $usage, $fail);
$fail = "";

GetOptions ("in=s" => \$infile,
            "out=s" => \$outdir,
            "h" => \$help,
            "ori=s" => \$ori);

$usage =
"usage: ./check_order.pl -in <file> -out <folder> -ori <gene to orient to>
options:
	-h			for this help message
	-in <file>		linear gene orders in FASTA format (REQUIRED)
	-out <folder>		output directory (REQUIRED)
	-ori <gene>		re-orient linear gene orders to this gene (REQUIRED)";

if ($help) {
        die ("\n$usage\n\n");
}

unless ($infile && $outdir && $ori) {
        die ("\nERROR: missing arguments, ./check_order.pl -h for usage\n\n");
}


### SET PATHS ###
my ($bindir, $maindir, $exedir, $rdir);

($bindir = abs_path($0)) =~ s/check_order.pl//g;
($maindir = $bindir) =~ s/scripts//g;
($exedir = $bindir) =~ s/scripts/programs/g;
$rdir = "$exedir/Rpackage";


### SET OUTPUT DIRECTORY ###
#my ($ans);
#if (-e $outdir && -d $outdir) {
#	print "Output folder already exists. Would you like to overwrite the whole directory? (Y/N)";
#	$ans = <STDIN>;
#	chomp $ans;
#	if ($ans !~ /^Y$/i) {
#		exit 0;
#	}
#}
#system("rm -rf $outdir");
#mkdir $outdir;


### EXECUTE SUBROUTINES ###

open IN, "$infile" or die $!;

while (<IN>) {
	chomp;
	make_hash($_);
}
close(IN);


if ($fail ne "fail") {
	cluster_go();
	report();
	plot("$outdir/order.patterntable.txt", "$outdir/order.patterns.pdf");
	plot("$outdir/order.speciestable.txt", "$outdir/order.species.pdf");
}

unlink("plot_order.Rout");
unlink("$outdir/order.patterntable.txt");
unlink("$outdir/order.speciestable.txt");


####################################
## SUBROUTINES
####################################

# SUBROUTINE: Orient linear gene orders to specified gene, invert gene order such that oriented gene is on forward strand and save in hash

my ($signal1, $signal2, $species, $acc, $oldline, @old_order, $gene, @new_order, @temporder1, @temporder2, $newline, %hash1, %check, $read);

sub make_hash {
	($read) = @_;
	$signal1 = "no";
	$signal2 = "fwd";
	$read =~ s/\t+$//g;
	if ($read =~ /^>/) {
		($species, $acc) = split(/\|/, $read);
		$species =~ s/>//g;
	} else {
		$oldline = $read;
		@old_order = split(/\t/, $oldline);
		foreach (@old_order) {
			$gene = "$_";
			if ($gene =~ /(-$ori|$ori)$/) {
				$signal1 = "yes";
				if ($gene =~ /^-/) {
					$signal2 = "rev";
					$gene = invert_strand($gene);
				}
				if (scalar @new_order == 0) {
					push (@new_order, $gene);
				} else {
					$fail = "fail";
					die ("\nERROR: Gene specified in -ori exists in duplicates.\nProcess terminated\n\n");
				}
			} else {
				if ($signal1 eq "no") {
					push (@temporder1, $gene);
				} elsif ($signal1 eq "yes") {
					push (@temporder2, $gene);
				}
			}
		}
		if ($signal2 eq "fwd") {
			push (@temporder2, @temporder1);
			push (@new_order, @temporder2);
		} elsif ($signal2 eq "rev") {
			@temporder1 = reverse @temporder1;
			@temporder2 = reverse @temporder2;
			push (@temporder1, @temporder2);
			foreach $gene (@temporder1) {
				$gene = invert_strand($gene);
				push (@new_order, $gene);
			}
		}
		if (scalar @new_order > 0) {
			$newline = join(" ", @new_order);
			$hash1{$species."|".$acc} = $newline;
			if (not exists $check{$newline}) {
				$check{$newline} = 1;
			}
		}
		@new_order = ();
		@temporder1 = ();
		@temporder2 = ();
		$species = "";
		$acc = "";
	}
}

# SUBROUTINE: Invert strand designation (input gene)

my $to_inv;

sub invert_strand {
	($to_inv) = @_;
	if ($to_inv !~ /^-/) {
		$to_inv = "-".$to_inv;
	} else {
		$to_inv =~ s/^-//g;
	}
	return $to_inv;
}


# SUBROUTINE: Cluster gene orders

my (%hash2, $go, $total, @members, $members, %hash3);

sub cluster_go {
	foreach $go (sort keys %check) {
		@members = grep { $hash1{$_} eq $go } keys %hash1;
		@members = sort { $a cmp $b } @members;
		$total = scalar(@members);
		$members = join(";", @members);
		$hash3{$members."%%".$go} = $total;
	}
}


# SUBROUTINE: Report gene orders in various output formats

sub report {
	my $i = 0;
	open OUT1, ">$outdir/order.species.txt";
	open OUT2, ">$outdir/order.patterns.txt";
	open OUT3, ">$outdir/order.detailed.txt";
	open OUT4, ">$outdir/order.patterntable.txt";
	open OUT5, ">$outdir/order.speciestable.txt";
	print OUT4 "Pattern_Num\tNum_Species\tGene_Order\n";
	print OUT5 "Species\tNum_Species\tGene_Order\n";
	foreach (sort { $hash3{$b} <=> $hash3{$a} } keys %hash3) {
		$i += 1;
		($members, $go) = split("%%", $_);
		@members = split(";", $members);
		$total = $hash3{$_};
		print OUT1 ">Pattern_$i ## $total species\n";
		foreach (@members) {
			print OUT1 "---$_\n";
		}
		print OUT1 "\n";
		print OUT2 ">Pattern_$i ## $total species ## $members\n";
		print OUT2 "$go\n";
		foreach (@members) {
			print OUT3 ">$_ ## Pattern_$i ## $total species\n";
			print OUT3 "$go\n";
			print OUT5 "$_ ... Pattern_$i\t$total\t$go\n";
		}
		print OUT4 "Pattern_$i\t$total\t$go\n";
	}
}

my ($plotin, $plotout);
sub plot {
	($plotin, $plotout) = @_;
	system("R CMD BATCH --vanilla '--args $plotin $plotout $rdir' $bindir/plot_order.R");
}
