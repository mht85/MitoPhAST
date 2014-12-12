#!usr/bin/perl
#use strict;
#use warnings;
use Getopt::Long;
use Bio::SeqIO;
use Bio::Species;
use Bio::SeqFeatureI;
use autodie;

####################################
# Command line options and Usage

GetOptions ("input_dir=s" => \$indir,
	    "output_dir=s" => \$outdir)
or die("Error in command line arguments\n");

$usage = "Usage: summarize.pl <input directory path>
options:
        -input_dir <path> input directory containing all GenBank or EMBL files
	-output_dir <path> output directory to print all summaries";

unless ($indir && $outdir) {
	die ("$usage\n");
}

####################################

# Run subroutines

opendir(IN_DIR, $indir) or die $!;

while ($fh1 = readdir(IN_DIR)) {
	$in_file = "$indir/$fh1";
	($out_file = $in_file) =~ s/^.*\///;
	next unless (-f "$in_file" && $in_file =~ m/\.(gb|gbk|gbff|gbf|embl|eml)$/);
	if (! -d $outdir ) {
		mkdir $outdir;
	}
	summarize();
}
closedir(IN_DIR);

sub summarize {
#	print "Printing summary for $in_file\n";
	open OUT_SUM, ">$outdir/$out_file.summary.txt";

	if ($in_file =~ /\.(gb|gbk|gbff|gbf)$/) {
		$fmt = "genbank";
	} elsif ($in_file =~ /\.(embl|eml)$/) {
		$fmt = "embl";
	}

	$seqio = Bio::SeqIO->new(-file => $in_file, -format => $fmt);
	$seqobj = $seqio->next_seq;
	$acc = $seqobj->accession;
	$species = $seqobj->species->node_name;
	$length = $seqobj->length;
	$dna = $seqobj->seq;
	$dna =~ tr/acgt/ACGT/;
	$a = sprintf("%.2f", (($dna=~tr/A//)/$length)*100);
	$c = sprintf("%.2f", (($dna=~tr/C//)/$length)*100);
	$g = sprintf("%.2f", (($dna=~tr/G//)/$length)*100);
	$t = sprintf("%.2f", (($dna=~tr/T//)/$length)*100);
	print OUT_SUM "AccessionID:\t",$acc,"\nSpecies:\t",$species,"\nLength:\t",$length,"\n\n%A:\t",$a,"\n%C:\t",$c,"\n%G:\t",$g,"\n%T:\t",$t,"\n%AT:\t",$a+$t,"\n%GC:\t",$g+$c,"\n";
	@features = $seqobj->get_SeqFeatures();
	foreach $feat (@features) {
		$start = $feat->start;
		$stop = $feat->end;
		$strand = $feat->strand;
		$strand =~ s/\-1/\-/g;
		$strand =~ s/1/\+/g;
		$len = $stop-$start+1;
		if ($feat->primary_tag eq 'CDS' || $feat->primary_tag eq 'tRNA' || $feat->primary_tag eq 'rRNA' || $feat->primary_tag eq 'misc_feature') {
			if ($feat->primary_tag eq 'CDS') {
				push (@CDS, "yes");
			} elsif ($feat->primary_tag eq 'tRNA') {
				push (@tRNA, "yes");
			} elsif ($feat->primary_tag eq 'rRNA') {
				push (@rRNA, "yes");
			} elsif ($feat->primary_tag eq 'misc_feature') {
				push (@CR, "yes");
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
	print OUT_SUM "#PCG: ",scalar@CDS,"\n#tRNA: ",scalar@tRNA,"\n#rRNA: ",scalar@rRNA,"\n#possible CR: ",scalar@CR,"\n\n";
	print OUT_SUM "Gene\tStrand\tPosition\tLength (bp)\tIntergenic_nucleotides (bp)\n";
	@CDS=();
	@tRNA=();
	@rRNA=();
	@CR=();
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
	@array=();
	@toprint=();
}
