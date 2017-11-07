#!usr/bin/perl
#use strict;
#use warnings;
use Getopt::Long;
use Bio::SeqIO;
use Bio::Species;
use autodie;

####################################
# Command line options and Usage

GetOptions ("input_dir=s" => \$indir,
	    "no-missing" => \$nomissing,
	    "species_list=s" => \$species,
	    "output_gap_dir=s" => \$outdir)
or die("Error in command line arguments\n");

$usage = "Usage: add_gaps.pl -input_dir <path> -output_gap_dir <path>
options:
	-input_dir <path> input directory containing all trimAl-trimmed alignments
	-species_list <file> file containing list of accession IDs of individuals
	-no-missing include only genes found in all individuals
	-output_gap_dir <path> output directory containing sequences with added gaps for each gene";

unless ($indir && $species && $outdir) {
	die ("$usage\n");
}

####################################

species_list();

opendir(IN_DIR, $indir) or die $!;

while ($fh1 = readdir(IN_DIR)) {
	$in_file = "$indir/$fh1";
	next unless (-f "$in_file" && $in_file =~ m/\-trimal.fasta$/);
	if ( ! -d $outdir ) {
                mkdir $outdir;
        }
	# test if file is empty
	open GBfile, "$in_file" or die $!;
	while (<GBfile>) {
		chomp;
		if ($_ !~ /^>.*$/ && $_ !~ /^$/) {
			push (@GBfile, $_);
		}
	}
	if (scalar@GBfile > 0) {
		dna_hash();
		id_missing(@acc);
		if (! defined $nomissing) {
			add_gaps();
		}
		print_seqs();
	}
	@GBfile=();
	@missing=();
	@present=();
	%seqhash=();
	%seqhash2=();
}
closedir(IN_DIR);
#create_phylip();
#create_partition_config();

sub species_list {
	open SPECIES, "$species" or die $!;
	while (<SPECIES>) {
		chomp;
		($acc, $name) = split(/\t/, $_);
		push (@acc, $acc);
	}
}


sub dna_hash {
	open SEQ, "$in_file" or die $!;
	while (<SEQ>) {
		chomp;
		if ($_ =~ /^>/) {	
			if (scalar @dna > 0) {
				$dna = join("", @dna);
				$seqhash{$header} = $dna;
				@dna = ();
			}
			($header = $_) =~ s/>//g;
		} else {
			$_ =~ s/\s//g;
			push (@dna, $_);
		}	
	}
	$dna = join("", @dna);
	$seqhash{$header} = $dna;
	@dna = ();
}

sub id_missing {
	foreach $key (keys %seqhash) {
		$value = $seqhash{$key};
		($key2, $gene) = split(/\-/, $key);
		$seqhash2{$key2} = $value;
	}
	foreach $acc (@acc) {
		if (! exists $seqhash2{$acc}) {
			push (@missing, $acc."-".$gene);
		} else {
			push (@present, $acc."-".$gene);
		}
	}
}

sub add_gaps {
	foreach (keys %seqhash2) {
		$seq = $seqhash2{$_};
		$length = length($seq);
	}
	for $miss (@missing) {
		$seqhash{$miss} = "-"x$length;
	}
}

sub print_seqs {
	if (! defined $nomissing) {
		open OUT, ">$outdir/$fh1-2";
		foreach (sort keys %seqhash) {
			print OUT ">$_\n$seqhash{$_}\n";
		}
	} elsif (defined $nomissing) {
		if (scalar@missing == 0) {
			open OUT, ">$outdir/$fh1-2";
			foreach (sort keys %seqhash) {
				print OUT ">$_\n$seqhash{$_}\n";
			}
		}
	}
}

