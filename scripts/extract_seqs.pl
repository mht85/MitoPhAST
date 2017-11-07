#!usr/bin/perl
#use strict;
#use warnings;
use Getopt::Long;
use Bio::Seq;
use Bio::SeqIO;
use Bio::Species;
use autodie;

####################################
# Command line options and Usage

GetOptions ("input_dir=s" => \$indir,
	    "output_dir=s" => \$outdir),
or die("Error in command line arguments\n");

$usage = "Usage: extract_seqs.pl -input_dir <input directory path> -output_dir <output directory path>
options:
        -input_dir <path> input directory containing all GenBank or EMBL files
	-output_dir <path> output directory containing all extracted PCG and rRNA sequences";

unless ($indir && $outdir) {
	die ("$usage\n");
}

####################################


opendir(IN_DIR, $indir) or die $!;

while ($fh1 = readdir(IN_DIR)) {
	$in_file = "$indir/$fh1";
	($out_file = $in_file) =~ s/^.*\///;
	next unless (-f "$in_file" && $in_file =~ m/\.(gb|gbk|gbff|gbf|embl|eml)$/);
	if (! -d $outdir ) {
                mkdir $outdir;
        }
	extract_seqs();
}
closedir(IN_DIR);

sub extract_seqs {
#	print "Extracting sequences for $in_file\n";
	open OUT_CDS, ">$outdir/$out_file.cds.fa";

	if ($in_file =~ /\.(gb|gbk|gbff|gbf)$/) {
		$fmt = "genbank";
	} elsif ($in_file =~ /\.(embl|eml)$/) {
		$fmt = "embl";
	}

	$seqio = Bio::SeqIO->new(-file => $in_file, -format => $fmt);
	$seqobj = $seqio->next_seq;
	$acc = $seqobj->accession;
	$species = $seqobj->species->node_name;
	for $feat ($seqobj->get_SeqFeatures) {
		if ($feat->primary_tag eq "CDS") {
			for $seq ($feat->get_tag_values('translation')) {
				if ($feat->has_tag('gene')) {
					for $val ($feat->get_tag_values('gene')) {
						$val =~ s/\s+/_/g;
						$val =~ s/\-/_/g;
						print OUT_CDS ">$acc-$val\n$seq\n";
						$val = "";
					}
				} elsif ($feat->has_tag('product')) {
					for $val ($feat->get_tag_values('product')) {
						$val =~ s/\s+/_/g;
						$val =~ s/\-/_/g;
						print OUT_CDS ">$acc-$val\n$seq\n";
						$val = "";
					}
				}
			}
		}
	}
}
