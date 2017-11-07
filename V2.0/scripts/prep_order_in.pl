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

### Make tRNA table
my (%trnhash);
%trnhash = (
	Ala => 'A',
	Arg => 'R',
	Asn => 'N',
	Asp => 'D',
	Cys => 'C',
	Gln => 'Q',
	Glu => 'E',
	Gly => 'G',
	His => 'H',
	Ile => 'I',
	Leu => 'L',
	Lys => 'K',
	Met => 'M',
	Phe => 'F',
	Pro => 'P',
	Ser => 'S',
	Thr => 'T',
	Trp => 'W',
	Tyr => 'Y',
	Val => 'V',
);

opendir(IN_DIR, $indir) or die $!;

while ($fh1 = readdir(IN_DIR)) {
	$in_file = "$indir/$fh1";
	($out_file = $in_file) =~ s/^.*\///;
	next unless (-f "$in_file" && $in_file =~ m/\.(gb|gbk|gbff|gbf|embl|eml)$/);
	if (! -d $outdir ) {
                mkdir $outdir;
        }
	extract_genes();
}
closedir(IN_DIR);
gene_order();

sub extract_genes {
#	print "Extracting sequences for $in_file\n";
#	open OUT_CDS, ">$outdir/$out_file.cds.fa";
	if ($in_file =~ /\.(gb|gbk|gbff|gbf)$/) {
		$fmt = "genbank";
	} elsif ($in_file =~ /\.(embl|eml)$/) {
		$fmt = "embl";
	}
	$seqio = Bio::SeqIO->new(-file => $in_file, -format => $fmt);
	$seqobj = $seqio->next_seq;
	$acc = $seqobj->accession;
	$species = $seqobj->species->node_name;
	$species =~ s/ /_/g;
	$species =~ s/'//g;
	$species =~ s/\(/_/g;
	$species =~ s/\)/_/g;
	for $feat ($seqobj->get_SeqFeatures) {
		# Extract PCG and save as key in array %allcds
		if ($feat->primary_tag eq 'CDS') {
			my ($cds_obj, $start, $end, $gene, $strand);
			foreach ($feat) {
				($cds_obj) = $feat;
				$start = $cds_obj->start;
				$end = $cds_obj->end;
				$strand = $cds_obj->strand;
				if ($cds_obj->has_tag('gene')) {
					($gene) = $cds_obj->get_tag_values('gene');
				} elsif ($cds_obj->has_tag('product')) {
					($gene) = $cds_obj->get_tag_values('product');
				}
				$gene =~ s/ /_/g;
				$gene =~ tr/a-z/A-Z/;
				$gene =~ s/COIII/COX3/g;
				$gene =~ s/COII/COX2/g;
				$gene =~ s/COI/COX1/g;
				$gene =~ s/CYTB/COB/g;
				$gene =~ s/NADH_DEHYDROGENASE_SUBUNIT_/NAD/g;
				$gene =~ s/ND/NAD/g;
				push (@annotarray, "$acc\#$species\#$start\#$end\#$strand\#$gene");
				$start=$end=$gene=$strand="";
			}
		# Extract rRNA and save as key in array %allrrn
		} elsif ($feat->primary_tag eq 'rRNA') {
			my ($rrn_obj, $start, $end, $nt, $gene, $note, $rrnline, $strand);
			foreach ($feat) {
				($rrn_obj) = $feat;
				$start = $rrn_obj->start;
				$end = $rrn_obj->end;
				$strand = $rrn_obj->strand;
				if ($rrn_obj->has_tag('gene')) {
					($gene) = $rrn_obj->get_tag_values('gene');
				} elsif ($rrn_obj->has_tag('note')) {
					($gene) = $rrn_obj->get_tag_values('note');
				} elsif ($rrn_obj->has_tag('product')) {
					($gene) = $rrn_obj->get_tag_values('product');
				}
				$gene =~ s/ /_/g;
				$gene =~ s/(12S.*|small.*|s-.*)/rrnS/g;
				$gene =~ s/(16S.*|large.*|l-.*)/rrnL/g;
				$gene =~ s/rnr1/rrnS/ig;
				$gene =~ s/rnr2/rrnL/ig;
				push (@annotarray, "$acc\#$species\#$start\#$end\#$strand\#$gene");
				$start=$end=$gene=$strand="";
			}
		} elsif ($feat->primary_tag eq 'tRNA') {
			my ($trn_obj, $start, $end, $nt, $gene, $note, $trnline, $strand);
			foreach ($feat) {
				($trn_obj) = $feat;
				$start = $trn_obj->start;
				$end = $trn_obj->end;
				$strand = $trn_obj->strand;
				if ($trn_obj->has_tag('gene')) {
					($gene) = $trn_obj->get_tag_values('gene');
				} elsif ($trn_obj->has_tag('product')) {
					($gene) = $trn_obj->get_tag_values('product');
				} elsif ($trn_obj->has_tag('note')) {
					($gene) = $trn_obj->get_tag_values('note');
				}
				$gene =~ s/\((\w+|\-+)\)$//g;
				$gene =~ s/\d+$//g;
				foreach (keys %trnhash) {
					$gene =~ s/$_/$trnhash{$_}/ig;	
				}
				$gene =~ s/tRNA-/trn/g;
				$gene =~ s/TRN/trn/g;
				push (@annotarray, "$acc\#$species\#$start\#$end\#$strand\#$gene");
				$start=$end=$gene=$strand="";
			}
		}
	}
}

# SUBROUTINE: Save gene order in tab files
my ($trn, $header, $header2, %check, $order, @all);
sub gene_order {
	open GENE, ">$outdir/sample_gene_order.txt";
	foreach (@annotarray) {
		$id = (split(/\#/, $_))[0];
		($species = (split(/\#/, $_))[1]) =~ s/_/ /g;
		$header = "$species|$id";
		$strand = (split(/\#/, $_))[4];
		$gene = (split(/\#/, $_))[5];
		if ($strand eq "1") {
			$strand = "";
		} elsif ($strand eq "-1") {
			$strand = "-";
		}
		if (! exists $check{$header}) {
			if (@all > 0) {
				$order = join("\t", @all);
				print GENE "$order\n";
			}
			print GENE ">$header\n";
			@all = ();
			$check{$header} = 1;
		}
		if ($_ =~ /\#trn/) {
			$trn = (split(/\#/, $_))[5];
			push (@all, $strand.$trn);
		} else {
			push (@all, $strand.$gene);
		}
	}

$order = join("\t", @all);
print GENE "$order\n";
@all = ();
$check{$header} = 1;
}
