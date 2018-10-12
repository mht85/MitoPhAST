#!/usr/bin/perl

#########################y
#Email mun.tan@deakin.edu.au
#GitHub mht85
#Date 21/06/17
#This script runs after prepDB.rrna.seed.pl to identify rRNAs that share homology to previously classified rRNA.
#Example ./prepDB.rrna.blast.pl -in db.rrna.unclassified.fa -db db.rrna.classified.fa -cpu 10
##########################

use strict;
use warnings;
use Getopt::Long;
use Cwd;
use Bio::SeqIO;


### USAGE STATEMENT ###
my ($in, $db, $usage, $help, $cpu);
$cpu = 1;

GetOptions ("in=s" => \$in,
	    "db=s" => \$db,
	    "h" => \$help,
	    "cpu=s" => \$cpu);

$usage =
"usage: ./prepDB.blast.pl -in <file>
options:
	-h		for this help message
	-in <file>	FASTA file with unclassified PCGs
	-db <file>	FASTA file with curated and classified PCGs
	-cpu <int>	number of threads";

if ($help) {
	print "$usage\n";
	exit;
}

unless ($in && $db) {
	die ("Error: missing input files, ./prepDB.blast.pl -h for usage\n");
}


### SET PATHS ###
my ($bindir, $maindir, $dbdir, $exedir);

$bindir = getcwd;
($maindir = getcwd) =~ s/bin//g;
($dbdir = getcwd) =~ s/bin/db/g;
($exedir = getcwd) =~ s/bin/exe/g;


### BLAST SEQUENCES ###
print "BLASTing sequences from $in to $db\n";

system("$exedir/ncbi-blast-2.6.0+/bin/makeblastdb -in $db -dbtype nucl");
system("$exedir/ncbi-blast-2.6.0+/bin/blastn -query $in -db $db -evalue 1e-5 -outfmt 6 -out $in.out -num_threads $cpu -max_target_seqs 1");

open IN, "$in";
my ($header, $sequence, %seqhash);
while (<IN>) {
	chomp;
	if ($_ =~ /^>/) {
		$header = $_;
		$header =~ s/>//g;
	} else {
		$sequence = $_;
		$seqhash{$header} = $sequence;
	}
}
$seqhash{$header} = $sequence;

open IN2, "$in.out";
my ($query, $hit, %rescued, @done, %done);
while (<IN2>) {
	$query = (split(/\t/, $_))[0];
	$hit = (split(/\t/, $_))[1];
	if (! exists $done{$query}) {
		$rescued{$query} = 1;
		$len = length($seqhash{$query});
		open OUT_CLASS, ">>db.$hit.fa";
		print OUT_CLASS ">$query\n";
		print OUT_CLASS "$seqhash{$query}\n";
		close OUT_CLASS;
		open OUT_ALLCLASS, ">>$db";
		print OUT_ALLCLASS ">$hit-$len\n";
		print OUT_ALLCLASS "$seqhash{$query}\n";
		close OUT_ALLCLASS;
		push (@done, $query);
	}
	$done{$_} = 1 foreach @done;
}

system("rm $in");
system("rm $db.nhr $db.nin $db.nsq");
system("rm $in.out");

foreach (keys %seqhash) {
	if (! exists $rescued{$_}) {
		open OUT_UNCLASS, ">>$in";
		print OUT_UNCLASS ">$_\n";
		print OUT_UNCLASS "$seqhash{$_}\n";
		close OUT_UNCLASS;
	}
}
