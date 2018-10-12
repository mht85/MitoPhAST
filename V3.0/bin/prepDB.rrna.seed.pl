#!/usr/bin/perl

#########################y
#Email mun.tan@deakin.edu.au
#Github mht85
#Date 21/06/17
#This script extracts rRNA nucleotide sequences and runs BLASTN to identify rRNAs that contain homology to seed rRNA sequences
#Example ./prepDB.rrna.seed.pl -in 090617_ncbi_metazoa_mito_refseq.gb -cpu 10
##########################

use strict;
use warnings;
use Getopt::Long;
use Cwd;
use Bio::SeqIO;


### USAGE STATEMENT ###
my ($in, $usage, $help, $cpu);
$cpu = 1;

GetOptions ("in=s" => \$in,
	    "h" => \$help,
	    "cpu" => \$cpu);

$usage =
"usage: ./prepDB.rrna.seed.pl -in <file>
options:
	-h		for this help message
	-in <file>	input file (one file) containing all metazoa mitogenomes on NCBI's RefSeq database (see README.txt)
	-cpu <int>	number of threads";

if ($help) {
	print "$usage\n";
	exit;
}

unless ($in) {
	die ("Error: missing input file, ./prepDB.rrna.seed.pl -h for usage\n");
}


### SET PATHS ###
my ($bindir, $maindir, $dbdir, $exedir);

$bindir = getcwd;
($maindir = getcwd) =~ s/bin//g;
($dbdir = getcwd) =~ s/bin/db/g;
($exedir = getcwd) =~ s/bin/exe/g;


### EXTRACT SEQUENCES ###
my ($seqio, $seqobj, $feat, $locus, $test, $id, $species, $rrn_obj, $start, $end, $gene, $nt, %allrrn, $note, $rrnline, $strand);
my ($domquery, $domhit, @pfamhit, @genegroup, %genehit, $pfam, $gene2);

system("$exedir/ncbi-blast-2.6.0+/bin/makeblastdb -in db.rrna.seed.fa -dbtype nucl");
print "Extracting sequences from $in\n";

$seqio = Bio::SeqIO->new('-file'=>"$in", '-format'=>"genbank");
while ($seqobj = $seqio->next_seq) {
	if ($seqobj->accession ne "unknown") {
		if ($seqobj->seq_version) {
			$id = $seqobj->accession.".".$seqobj->seq_version;	# Accession.version
		} else {
			$id = $seqobj->accession.".V";				# No version number
		}
	} else {
		$locus = $seqobj->display_id;
		$id = $locus.".V";						# No accession/version number
	}
	$species = $seqobj->species->node_name;					# Species name
	$species =~ s/ /_/g;
	$species =~ s/'//g;
	$species =~ s/(\(|\))/_/g;
	for $feat ($seqobj->get_SeqFeatures) {
		# Extract rRNA and save as key in array %allrrn
		if ($feat->primary_tag eq 'rRNA') {
			foreach ($feat) {
				($rrn_obj) = $feat;
				$start = $rrn_obj->start;
				$end = $rrn_obj->end;
				$strand = $rrn_obj->strand;
				$nt = $rrn_obj->spliced_seq->seq;
				if ($rrn_obj->has_tag('gene')) {
					($gene) = $rrn_obj->get_tag_values('gene');
				} elsif ($rrn_obj->has_tag('product')) {
					($gene) = $rrn_obj->get_tag_values('product');
				} elsif ($rrn_obj->has_tag('note')) {
					($gene) = $rrn_obj->get_tag_values('note');
				}
				$gene =~ s/ /_/g;
				$rrnline = "$id\#$species\#$start\#$end\#$strand\#$gene";
				$allrrn{$rrnline} = $nt;
				open OUT_TMP2, ">>seq.rrna.tmp";
				print OUT_TMP2 ">$rrnline\n";
				print OUT_TMP2 "$nt\n";
				close OUT_TMP2;
				system("$exedir/ncbi-blast-2.6.0+/bin/blastn -query seq.rrna.tmp -db db.rrna.seed.fa -evalue 1e-20 -outfmt 6 -out seq.rrna.tmp.out -num_threads $cpu -max_target_seqs 1");
				if (-z "seq.rrna.tmp.out") {
					open OUT_UNCLASS, ">>db.rrna.unclassified.fa";
					print OUT_UNCLASS ">$rrnline\n";
					print OUT_UNCLASS "$allrrn{$rrnline}\n";
					close OUT_UNCLASS;
				} else {
					open IN, "seq.rrna.tmp.out";
					my ($query, $hit, %done, @done);
					while (<IN>) {
						$query = (split(/\t/, $_))[0];
						$hit = (split(/\t/, $_))[1];
						if (! exists $done{$query}) {
							$len = length($allrrn{$query});
							open OUT_CLASS, ">>db.$hit.fa";
							print OUT_CLASS ">$query\n";
							print OUT_CLASS "$allrrn{$query}\n";
							close OUT_CLASS;
							open OUT_ALLCLASS, ">>db.rrna.classified.fa";
							print OUT_ALLCLASS ">$hit-$len\n";
							print OUT_ALLCLASS "$allrrn{$query}\n";
							close OUT_ALLCLASS;
							push (@done, $query);
						}
						$done{$_} = 1 foreach @done;
					}
				}
				system("rm seq.rrna.tmp seq.rrna.tmp.out");
				$start=$end=$nt=$gene="";
			}
		}
	}
}

system("rm db.rrna.seed.fa.nhr db.rrna.seed.fa.nin db.rrna.seed.fa.nsq");
