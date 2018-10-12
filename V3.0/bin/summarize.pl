#!/usr/bin/perl

##########################
#Email mun.tan@deakin.edu.au
#Github mht85
#Date 22/06/17
#This script provides a summary for each mitogenome GenBank/EMBL file
#Example ./summarize.pl -in ../example_data -out summary
##########################

use strict;
use warnings;
use Getopt::Long;
use Cwd 'abs_path';
use Bio::SeqIO;


### USAGE STATEMENT ###
my ($indir, $outdir, $usage, $help);

GetOptions ("in=s" => \$indir,
	    "out=s" => \$outdir,
            "h" => \$help,
	    "help" => \$help);

$usage =
"usage: ./summarize.pl -in <folder> -out <folder>
options:
        -h/help		for this help message
        -in <folder>	input directory containing all GenBank or EMBL files (REQUIRED)
        -out <folder>	output directory to print all summaries (REQUIRED)";

if ($help) {
        die ("\n$usage\n\n");
}

unless ($indir && $outdir) {
        die ("\nError: missing arguments, ./summarize.pl -h for usage\n\n");
}


### SET PATHS ###
my ($bindir, $maindir, $sumdir);

($bindir = abs_path($0)) =~ s/summarize.pl//g;
($maindir = $bindir) =~ s/bin//g;
$sumdir = "$outdir/1.summaries";
mkdir $sumdir;

### SET OUTPUT DIRECTORY ###
#my ($ans);
#if (-e $outdir && -d $outdir) {
#	print "Output folder already exists. Would you like to overwrite the whole directory? (Y/N)";
#	$ans = <STDIN>;
#	chomp $ans;
#	if ($ans !~ /^Y$/i) {
#		exit;
#	}
#}
#system("rm -rf $outdir");
#mkdir $outdir;


### EXECUTE SUBROUTINES ###
my ($fh1, $in_file, $out_file);
opendir(IN_DIR, $indir) or die $!;
while ($fh1 = readdir(IN_DIR)) {
	$in_file = "$indir/$fh1";
	($out_file = $in_file) =~ s/^.*\///;
	next unless (-f "$in_file" && $in_file =~ m/\.(gb|gbk|gbff|gbf|embl|eml)$/);
	summarize();
}
closedir(IN_DIR);

my (@reportlist, @flaglist, $ans);
open OUT_REPORT,  ">$sumdir/summary_report.txt";
print OUT_REPORT "AccessionID\tSpecies\tLength\t%A\t%C\t%G\t%T\t#CDS\t#tRNA\t#rRNA\t#CR\n";
foreach (@reportlist) {
	print OUT_REPORT "$_\n";
}

if (scalar@flaglist > 0) {
        print "\nWARNING: Mitogenomes for the following samples have unexpected number of genes for CDS (13), tRNA (22) and rRNA (2)\n";
	foreach (@flaglist) {
	        print "$_\n";
	}
	print "\n";
#	print "Proceed? (Y/N)";
#	$ans = <STDIN>;
#	chomp $ans;
#	if ($ans !~ /^Y$/i) {
#		exit;
#	}
}

system("touch $outdir/summary.success");



####################################
## SUBROUTINES
####################################

# SUBROUTINE: Prints summary for each file

my ($printname, $fmt, $seqio, $seqobj, $acc, $species, $length, $dna, $a, $c, $g, $t, @features, $feat, $start, $stop, $strand, $len, @CDS, @tRNA, @rRNA, @CR, $name, @array, $prior, $name2, $strand2, $start2, $stop2, $len2, $inter, $name3, $strand3, $start3, $stop3, $len3, @toprint, $note, $misc_obj, $numCDS, $numtRNA, $numrRNA, $numCR, $loc, @pos, $pos1, $pos2);

sub summarize {
	$prior="";
	($printname = $in_file) =~ s/^.*\///g;
	print "Printing summary for $printname\n";
	open OUT_SUM, ">$sumdir/$out_file.summary.txt";
	if ($in_file =~ /\.(gb|gbk|gbff|gbf)$/) {
		$fmt = "genbank";
	} elsif ($in_file =~ /\.(embl|eml)$/) {
		$fmt = "embl";
	}
	$seqio = Bio::SeqIO->new(-file => $in_file, -format => $fmt);
	$seqobj = $seqio->next_seq;
	if ($seqobj->accession ne "unknown") {
		if ($seqobj->seq_version) {
			$acc = $seqobj->accession.".".$seqobj->seq_version;      # Accession.version
		} else {
			$acc = $seqobj->accession.".V0";
		}
	} else {
		$loc = $seqobj->display_id;
		$acc = $loc.".V0";
	}
	$species = $seqobj->species->node_name;
	$length = $seqobj->length;
	$dna = $seqobj->seq;
	$dna =~ tr/acgt/ACGT/;
	$a = sprintf("%.2f", (($dna=~tr/A//)/$length)*100);
	$c = sprintf("%.2f", (($dna=~tr/C//)/$length)*100);
	$g = sprintf("%.2f", (($dna=~tr/G//)/$length)*100);
	$t = sprintf("%.2f", (($dna=~tr/T//)/$length)*100);
	print OUT_SUM "AccessionID:\t",$acc,"\nSpecies:\t",$species,"\nLength:\t",$length,"\n\n%A: ",$a,"\n%C: ",$c,"\n%G: ",$g,"\n%T: ",$t,"\n%AT: ",$a+$t,"\n%GC: ",$g+$c,"\n";
	@features = $seqobj->get_SeqFeatures();
	foreach $feat (@features) {
		if ( $feat->location->isa('Bio::Location::SplitLocationI')) {
			@pos = $feat->location->sub_Location;
			$pos1 = shift @pos;
			$start = $pos1->start;
			$pos2 = pop @pos;
			$stop = $pos2->end;
			$len = $length - $start + 1 + $stop - 1 + 1;
		} else {
			$start = $feat->start;
			$stop = $feat->end;
			$len = $stop - $start + 1;
		}
		$strand = $feat->strand;
		$strand =~ s/\-1/\-/g;
		$strand =~ s/1/\+/g;
		if ($feat->primary_tag eq 'CDS' || $feat->primary_tag eq 'tRNA' || $feat->primary_tag eq 'rRNA' || $feat->primary_tag eq 'misc_feature' || $feat->primary_tag eq 'D-loop') {
			if ($feat->primary_tag eq 'CDS') {
				push (@CDS, "yes");
			} elsif ($feat->primary_tag eq 'tRNA') {
				push (@tRNA, "yes");
			} elsif ($feat->primary_tag eq 'rRNA') {
				push (@rRNA, "yes");
			} elsif ($feat->primary_tag eq 'misc_feature' || $feat->primary_tag eq 'D-loop') {
				foreach ($feat) {
					($misc_obj) = $feat;
					if ($misc_obj->has_tag('note')) {
						($note) = $misc_obj->get_tag_values('note');
					} elsif ($misc_obj->has_tag('product')) {
						($note) = $misc_obj->get_tag_values('product');
					}
					if ($note =~ /(control|A\+T)/i) {
						push (@CR, "yes");
					}
				}
			}
			if ($feat->has_tag('gene')) {
				for $name ($feat->get_tag_values('gene')) {
					push (@array, $name."\t".$strand."\t".$start."\t".$stop."\t".$len);
				}
			} elsif ($feat->has_tag('product')) {
				for $name ($feat->get_tag_values('product')) {
					push (@array, $name."\t".$strand."\t".$start."\t".$stop."\t".$len);
				}
			} elsif ($feat->has_tag('note')) {
				for $name ($feat->get_tag_values('note')) {
					push (@array, $name."\t".$strand."\t".$start."\t".$stop."\t".$len);
				}
			} else {
				push (@array, "misc_feature\t".$strand."\t".$start."\t".$stop."\t".$len);
			}
		}
	}
	$numCDS = scalar@CDS;
	chomp $numCDS;
	$numtRNA = scalar@tRNA;
	chomp $numtRNA;
	$numrRNA = scalar@rRNA;
	chomp $numrRNA;
	$numCR = scalar@CR;
	chomp $numCR;
	print OUT_SUM "#PCG: ",$numCDS,"\n#tRNA: ",$numtRNA,"\n#rRNA: ",$numrRNA,"\n#possible CR: ",$numCR,"\n\n";
	print OUT_SUM "Gene\tStrand\tPosition\tLength (bp)\tIntergenic_nucleotides (bp)\n";
	foreach (@array) {
		if ($prior =~ /\w+/) {
			($name2, $strand2, $start2, $stop2, $len2) = split(/\t/, $_);
			$inter = $start2 - $prior - 1;
			push (@toprint, $name2."\t".$strand2."\t".$start2."-".$stop2."\t".$len2."\t".$inter);
			$prior = $stop2;
		} else {
			($name3, $strand3, $start3, $stop3, $len3) = split(/\t/, $_);
			$prior = $stop3;
		}
	}
	$inter = $length - $prior;
	unshift (@toprint, $name3."\t".$strand3."\t".$start3."-".$stop3."\t".$len3."\t".$inter);
	foreach (@toprint) {
		print OUT_SUM $_,"\n";
	}
	$prior="";
	if (scalar@CDS ne "13" || scalar@tRNA ne "22" || scalar@rRNA ne "2") {
		push (@flaglist, "$acc\t$species\t$numCDS CDS, $numtRNA tRNA, $numrRNA, rRNA");
	}
	push (@reportlist, "$acc\t$species\t$length\t$a\t$c\t$g\t$t\t$numCDS\t$numtRNA\t$numrRNA\t$numCR");
	@array=();
	@toprint=();
	@CDS=();
	@tRNA=();
	@rRNA=();
	@CR=();
}
