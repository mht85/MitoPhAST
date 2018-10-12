#!/usr/bin/perl


##########################
#Email mun.tan@deakin.edu.au
#Github mht85
#Date 08/08/17
#This script downloads GenBank files from NCBI when given accession numbers
#Example ./get_ncbi.pl -acc AAAAAAAA,BBBBBBBB -list <file> -out <folder>
##########################

use strict;
use warnings;
use Getopt::Long;
use Cwd 'abs_path';
use Bio::DB::EUtilities;
use Bio::SeqIO;

### USAGE STATEMENT ###
my ($inacc, $inlist, $outdir, $usage, $help);

GetOptions ("acc=s" => \$inacc,
	    "list=s" => \$inlist,
	    "out=s" => \$outdir,
	    "help" => \$help,
            "h" => \$help);

$usage =
"usage: ./get_ncbi.pl -acc KT984196,KU500619 -list <file> -out <folder>
options:
        -h/help		for this help message
        -acc <string>	string of accession numbers to download, comma-separated string (at least one of -acc or -list REQUIRED)
	-list <file>	file of accession numbers to download, one accession/line (at least one of -acc or -list REQUIRED)
        -out <folder>	output directory to store downloaded files (REQUIRED)";

if ($help) {
        die ("\n$usage\n\n");
}

unless (($inacc || $inlist) && $outdir) {
        die ("\nError: missing arguments, ./get_ncbi.pl -h for usage\n\n");
}


### SET PATHS ###
my ($bindir, $maindir, $dldir);

($bindir = abs_path($0)) =~ s/get_ncbi.pl//g;
($maindir = $bindir) =~ s/bin//g;
$dldir = "$outdir/1.downloads";
mkdir $dldir;

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


### EXECUTE TASK ###

my (@inacc, @allacc, $id, $gbfetch, $outfile);

if ($inacc) {
	chomp $inacc;
	@inacc = split(/\,/, $inacc);
	push (@allacc, @inacc);
}

if ($inlist) {
	open IN, "$inlist" or die $!;
	while (<IN>) {
		chomp;
		push (@allacc, $_);
	}
	close(IN);
}

foreach $id (@allacc) {
	print "Downloading entry for accession $id\n";
        $gbfetch = Bio::DB::EUtilities->new(-eutil => "efetch",
                                            -db      => "nucleotide",
                                            -rettype => "gbwithparts",
					    -email => "null",
                                            -id      => "$id");
	$id =~ s/\.\d+//g;
	$outfile = "$dldir/$id.gb";
	$gbfetch->get_Response(-file => $outfile);
	if (! -e $outfile) {
		die ("\nError: failed to download file for accession $id\n\n");
	}
}

system("touch $outdir/download.success");

