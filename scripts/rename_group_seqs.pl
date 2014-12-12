#!usr/bin/perl
#use strict;
#use warnings;
use Getopt::Long;
use Bio::SeqIO;
use Bio::Species;
use autodie;

####################################
# Command line options and Usage

GetOptions ("input_profiles_dir=s" => \$indir,
	    "input_sequences_dir=s" => \$indir2,
	    "output_renamed_dir=s" => \$outdir,
	    "output_grouped_dir=s" => \$outdir2)
or die("Error in command line arguments\n");

$usage = "Usage: rename_group_seqs.pl -input_profiles_dir <path> -input_sequences_dir <path> -output_renamed_dir <path> -output_grouped_dir <path>
options:
        -input_profiles_dir <path> input directory containing all PFAM/PRINT profile assignments
	-input_sequences_dir <path> input directory containing sequences from all individuals
	-output_renamed_dir <path> output directory containing renamed sequence files
	-output_grouped_dir <path> output directory containing grouped sequence files";

unless ($indir && $indir2 && $outdir && $outdir2) {
	die ("$usage\n");
}

####################################

opendir(IN_DIR, $indir) or die $!;

while ($fh1 = readdir(IN_DIR)) {
	$in_file = "$indir/$fh1";
	$fh1 =~ s/\.final/\.fa/g;
	$in_file2 = "$indir2/$fh1";
	($out_file = $in_file2) =~ s/^.*\///;
	$out_file =~ s/\.fa/\.renamed\.fa/g;
	next unless (-f "$in_file" && $in_file =~ m/\.final$/);
	if ( ! -d $outdir ) {
                mkdir $outdir;
        }
	if ( ! -d $outdir2 ) {
		mkdir $outdir2;
	}
	rename_seqs();
}
closedir(IN_DIR);
print "\n";
group_seqs("ATP8","ATP6","COX1","COX2","COX3","COB","NAD1","NAD2","NAD3","NAD4","NAD4L","NAD5","NAD6");




sub rename_seqs {
#	print "Renaming sequences for $fh1\n";

	open IN, "$in_file" or die $!;
	open IN2, "$in_file2" or die $!;
	open OUT_RENAME, ">$outdir/$out_file";

	while ($profiles=<IN>) {
		chomp $profiles;
		($old, $new) = split(/\t/, $profiles);
		($keep, $null) = split(/\-/, $old);
		$new = $keep."-".$new;
		$hashprofile{$old} = $new;
	}
	while ($tosub=<IN2>) {
		chomp $tosub;
		if ($tosub =~ /^>/) {
			$tosub2 = $tosub;
			$tosub2 =~ s/>//g;
			if (exists $hashprofile{$tosub2}) {
				$tosub =~ s/$tosub2/$hashprofile{$tosub2}/g;
			} elsif ($tosub2 =~ /(ATP6|ATP8|COX1|COX2|COX3|COB|NAD1|NAD2|NAD3|NAD4|NAD4L|NAD5|NAD6)/i) {
				$tosub =~ tr/[a-z]/[A-Z]/;
			} elsif ($tosub2 =~ /ATPASE(_)/i) {
				$tosub =~ s/ATPASE(_)/ATP/ig;
			} elsif ($tosub2 =~ /CO(X)III/i) {
				$tosub =~ s/CO(X)III/COX3/ig;
			} elsif ($tosub2 =~ /CO(X)II/i) {
				$tosub =~ s/CO(X)II/COX2/ig;
			} elsif ($tosub2 =~ /CO(X)I/i) {
				$tosub =~ s/CO(X)I/COX1/ig;
			} elsif ($tosub2 =~ /CO3/i) {
				$tosub=~ s/CO3/COX3/ig;
			} elsif ($tosub2 =~ /CO2/i) {
				$tosub=~ s/CO2/COX2/ig;
			} elsif ($tosub2 =~ /CO1/i) {
				$tosub=~ s/CO1/COX1/ig;
			} elsif ($tosub2 =~ /ND/i) {
				$tosub =~ s/ND/NAD/ig;
			} elsif ($tosub2 =~ /NADH/i) {
				$tosub =~ s/NADH/NAD/ig;
			} elsif ($tosub2 =~ /L/i) {
				$tosub =~ s/L/L/ig;
			} elsif ($tosub2 =~ /CYT(_)B/i) {
				$tosub =~ s/CYT(_)B/COB/ig;
			} else {
				$tosub = "";
			}
			if ($tosub =~ /\w+/) {
				$header = $tosub;
				print OUT_RENAME "$tosub\n";
			}
		} else {
			if ($header =~ /\w+/) {
				$sequence = $tosub;
				$nexthash{$header} = $sequence;
				$header = "";
				$sequence = "";
				print OUT_RENAME "$tosub\n";
			}
		}
	}
	$nexthash{$header} = $sequence;
}

sub group_seqs {
	foreach $gene (@_) {
		open OUT_GENE, ">$outdir2/all_$gene.cds.fa";
		foreach $key (keys %nexthash) {
			if ($key =~ /$gene\b/) {
				print OUT_GENE "$key\n$nexthash{$key}\n";
			}
		}
	}
}
