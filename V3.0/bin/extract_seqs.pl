#!/usr/bin/perl

##########################
#Email mun.tan@deakin.edu.au
#GitHub mht85
#Date 22/06/17
#This script extracts PCG and rRNA sequences from GenBank/EMBL files and is clustered according to homology using BLASTP/N to sequence database in <db> folder
#Example ./extract_seqs.pl -in <folder> -out <folder> -cpu 10
##########################

use strict;
use warnings;
use Getopt::Long;
use Cwd 'abs_path';
use Bio::SeqIO;


### USAGE STATEMENT ###
my ($indir, $outdir, $usage, $help, $cpu, $interactive);
$cpu = 1;
$interactive = "on";

GetOptions ("in=s" => \$indir,
	    "out=s" => \$outdir,
	    "help" => \$help,
	    "h" => \$help,
	    "cpu=s" => \$cpu,
	    "interactive=s" => \$interactive);

$usage =
"usage: ./extract_seqs.pl -in <file>
options:
	-h/help			for this help message
	-in <folder>		input directory containing all GenBank/EMBL files (REQUIRED)
	-out <folder>		output directory (REQUIRED)
	-cpu <int>		number of threads, default=1
	-interactive <on/off>	control interactive messages (default=on, recommended)";

if ($help) {
	die ("\n$usage\n\n");
}

unless ($indir && $outdir) {
	die ("\nError: missing input files, ./extract_seqs.pl -h for usage\n\n");
}


### SET PATHS ###
my ($bindir, $maindir, $dbdir, $exedir, $seqdir, $profdir);

($bindir = abs_path($0)) =~ s/extract_seqs.pl//g;
($maindir = $bindir) =~ s/bin//g;
($dbdir = $bindir) =~ s/bin/db/g;
($exedir = $bindir) =~ s/bin/exe/g;
$seqdir = "$outdir/1.sequences";
$profdir = "$outdir/2.profiles";
mkdir $seqdir;
mkdir $profdir;

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


### Collection of PFAM profiles and e-values
my (%pfamhash);
%pfamhash = (
	'PF00032.hmm' => '1e-5',
	'PF00115.hmm' => '1e-5',
	'PF00116.hmm' => '1e-5',
	'PF00119.hmm' => '1e-5',
	'PF00146.hmm' => '1e-5',
	'PF00499.hmm' => '1e-5',
	'PF00507.hmm' => '1e-5',
	'PF00510.hmm' => '1e-5',
	'PF00895.hmm' => '1e-3',
	'PF00420.hmm' => '1e-2',
);


### Collection of PRINTS profiles and e-values
my (%printhash);
%printhash = (
	'PR01434.pval' => '1e-5',
	'PR01436.pval' => '1e-5',
	'PR01437.pval' => '1e-5',
);


### Gene name substitution
my (%genenamehash);
%genenamehash = (
	'ATP-synt_A' => 'ATP6',
	'ATP-synt_8' => 'ATP8',
	'COX1' => 'COX1',
	'COX2' => 'COX2',
	'COX3' => 'COX3',
	'Cytochrom_B_C' => 'CYTB',
	'NADHdh' => 'ND1',
	'NADHDHGNASE2' => 'ND2',
	'Oxidored_q4' => 'ND3',
	'PRINT_ND4' => 'ND4',
	'Oxidored_q2' => 'ND4L',
	'PRINT_ND5' => 'ND5',
	'Oxidored_q3' => 'ND6',
);


### Check gene names
my (%genes);
%genes = (
	'ATP6' => 'ATP6/ATPASE6/ATPASE_6',
	'ATP8' => 'ATP8/ATPASE8/ATPASE_8',
	'COX1' => 'COX1/COI/CO1',
	'COX2' => 'COX2/COII/CO2',
	'COX3' => 'COX3/COIII/CO3',
	'CYTB' => 'CYTB/COB/CYT_B',
	'ND1' => 'ND1/NAD1/NADH1',
	'ND2' => 'ND2/NAD2/NADH2',
	'ND3' => 'ND3/NAD3/NADH3',
	'ND4' => 'ND4/NAD4/NADH4',
	'ND4L' => 'ND4L/NAD4L/NADH4L',
	'ND5' => 'ND5/NAD5/NADH5',
	'ND6' => 'ND6/NAD6/NADH6',
	'rrnL' => 'rrnL/16S_RIBOSOMAL_RNA',
	'rrnS' => 'rrnS/12S_RIBOSOMAL_RNA',
);


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
my ($file, $format, $in_file, $id);

opendir(IN_DIR, $indir) or die $!;
while ($file = readdir(IN_DIR)) {
	if ($file !~ /\.(gb|gbk|gbf|gbff|embl|eml)$/) {
		next;
	} elsif ($file =~ /\.(gb|gbk|gbf|gbff)$/) {
		$format = "genbank";
	} elsif ($file =~ /\.(embl|eml)$/) {
		$format = "embl";
	}
	$in_file = "$indir/$file";
	extract_seqs($in_file, $format);
	validatePCG($in_file);
	validaterRNA($in_file);
	findMissing($in_file);
}
validateName();
printSeqs();
gene_order();
closedir(IN_DIR);
#cluster_seqs();
#gene_order();
system("touch $outdir/extract.success");

####################################
## SUBROUTINES
####################################

# SUBROUTINE: Extract PCG and rRNA sequences from each Genbank/EMBL file

my ($in_sub, $printname, %allaa, %allnt, %allrrn, @annotarray, $acc, @pos, $pos1, $pos2);
sub extract_seqs {
	($in_sub, $format) = @_;
	($printname = $in_sub) =~ s/^.*\///g;
	print "Extracting sequences for $printname\n";
	open OUT_AA, ">$seqdir/$printname.pcg.aa.fa";
	open OUT_NT, ">$seqdir/$printname.pcg.nt.fa";
	open OUT_RRN, ">$seqdir/$printname.rrna.fa";
	my ($seqio, $seqobj, $feat, $loc, $test, $species);
	$seqio = Bio::SeqIO->new(-file => "$in_sub", '-format' => "$format");
	while ($seqobj = $seqio->next_seq) {
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
		$species = $seqobj->species->node_name;                 # Species name
		$species =~ s/ /_/g;
		$species =~ s/'//g;
		$species =~ s/\(/_/g;
		$species =~ s/\)/_/g;
		for $feat ($seqobj->get_SeqFeatures) {
			# Extract PCG and save as key in array %allcds
			if ($feat->primary_tag eq 'CDS') {
				my ($cds_obj, $start, $end, $codon_start, $nt, $cut, $transl_except, $transl_except_pos, @transl_except_pos, $add, $gene, $aa, $note, $aaline, $ntline, $strand);
				foreach ($feat) {
					($cds_obj) = $feat;
					if ( $cds_obj->location->isa('Bio::Location::SplitLocationI')) {
						@pos = $cds_obj->location->sub_Location;
						$pos1 = shift @pos;
						$start = $pos1->start;
						$pos2 = pop @pos;
						$end = $pos2->end;
					} else {
						$start = $cds_obj->start;
						$end = $cds_obj->end;
					}
					($codon_start) = $cds_obj->get_tag_values('codon_start');
					$strand = $cds_obj->strand;
					$nt = $cds_obj->spliced_seq->seq;
					$nt = substr($nt, $codon_start-1);
					if ($cds_obj->has_tag('transl_except')) {
						($transl_except) = $cds_obj->get_tag_values('transl_except');
						$transl_except_pos = (split(/\,/, $transl_except))[0];
						$transl_except_pos =~ s/(pos|complement|\(|\)|\:)//g;
						@transl_except_pos = split(/\.\./, $transl_except_pos);
						if (scalar @transl_except_pos == 1) {
							$add = "2";
						} else {
							$add = "1";
						}
						$nt .= "A"x$add;
					}
					$gene=$aa=$note="";
					if ($cds_obj->has_tag('gene')) {
						($gene) = $cds_obj->get_tag_values('gene');
					} elsif ($cds_obj->has_tag('product')) {
						($gene) = $cds_obj->get_tag_values('product');
					}
					$gene =~ s/ /_/g;
					($aa) = $cds_obj->get_tag_values('translation') if $cds_obj->has_tag('translation');
					$cut = int(length($nt)/3) - 1;
					$aa = substr($aa, "-$cut");
					($note) = $cds_obj->get_tag_values('note') if ($cds_obj->has_tag('note'));
					$aaline = "$acc\#$species\#$start\#$end\#$strand\#$gene";
					$ntline = "$acc\#$species\#$start\#$end\#$strand\#$gene";
					push (@annotarray, "$acc\#$species\#$start\#$end\#$strand");
					$allaa{$aaline} = "$aa";
					$allnt{$ntline} = "$nt";
					print OUT_AA ">$aaline\n";
					print OUT_AA "$aa\n";
					print OUT_NT ">$ntline\n";
					print OUT_NT "$nt\n";
					$start=$end=$nt=$gene=$aa=$note=$aaline=$ntline=$strand="";
				}
			# Extract rRNA and save as key in array %allrrn
			} elsif ($feat->primary_tag eq 'rRNA') {
				my ($rrn_obj, $start, $end, $nt, $gene, $note, $rrnline, $strand);
				foreach ($feat) {
					($rrn_obj) = $feat;
					if ( $rrn_obj->location->isa('Bio::Location::SplitLocationI')) {
						@pos = $rrn_obj->location->sub_Location;
						$pos1 = shift @pos;
						$start = $pos1->start;
						$pos2 = pop @pos;
						$end = $pos2->end;
					} else {
						$start = $rrn_obj->start;
						$end = $rrn_obj->end;
					}
					$strand = $rrn_obj->strand;
					$nt = $rrn_obj->spliced_seq->seq;
					if ($rrn_obj->has_tag('gene')) {
						($gene) = $rrn_obj->get_tag_values('gene');
					} elsif ($rrn_obj->has_tag('note')) {
						($gene) = $rrn_obj->get_tag_values('note');
					} elsif ($rrn_obj->has_tag('product')) {
						($gene) = $rrn_obj->get_tag_values('product');
					}
					$gene =~ s/ /_/g;
					push (@annotarray, "$acc\#$species\#$start\#$end\#$strand");
					$rrnline = "$acc\#$species\#$start\#$end\#$strand\#$gene";
					$allrrn{$rrnline} = "$nt";
					print OUT_RRN ">$rrnline\n";
					print OUT_RRN "$nt\n";
					$start=$end=$nt=$gene="";
				}
			} elsif ($feat->primary_tag eq 'tRNA') {
				my ($trn_obj, $start, $end, $nt, $gene, $note, $trnline, $strand);
				foreach ($feat) {
					($trn_obj) = $feat;
					if ( $trn_obj->location->isa('Bio::Location::SplitLocationI')) {
						@pos = $trn_obj->location->sub_Location;
						$pos1 = shift @pos;
						$start = $pos1->start;
						$pos2 = pop @pos;
						$end = $pos2->end;
					} else {
						$start = $trn_obj->start;
						$end = $trn_obj->end;
					}
					$strand = $trn_obj->strand;
					$nt = $trn_obj->spliced_seq->seq;
					if ($trn_obj->has_tag('gene')) {
						($gene) = $trn_obj->get_tag_values('gene');
					} elsif ($trn_obj->has_tag('product')) {
						($gene) = $trn_obj->get_tag_values('product');
					} elsif ($trn_obj->has_tag('note')) {
						($gene) = $trn_obj->get_tag_values('note');
					}
					$gene =~ s/\((\w+|\-+)\)$//g;
					$gene =~ s/\-(a|c|g|u|t)+$//g;
					$gene =~ s/\d+$//g;
					foreach (keys %trnhash) {
						$gene =~ s/$_/$trnhash{$_}/ig;
					}
					$gene =~ s/tRNA-/trn/g;
					push (@annotarray, "$acc\#$species\#$start\#$end\#$strand\t$gene");
				}
#			} elsif ($feat->primary_tag eq 'misc_feature' || $feat->primary_tag eq 'D-loop') {
#				my ($misc_obj, $start, $end, $nt, $gene, $note, $miscline, $strand);
#				foreach ($feat) {
#					($misc_obj) = $feat;
##					$start = $misc_obj->start;
##					$end = $misc_obj->end;
#					$strand = $misc_obj->strand;
#					$nt = $misc_obj->seq->seq;
#					if ($misc_obj->has_tag('note')) {
#						($gene) = $misc_obj->get_tag_values('note');
#					} elsif ($misc_obj->has_tag('product')) {
#						($gene) = $misc_obj->get_tag_values('product');
#					}
#					if ($gene =~ /(control|A\+T)/i) {
#						push (@annotarray, "$acc\#$species\#$start\#$end\#$strand\#CR");
#					}
#				}
			}
		}
	}
}

# SUBROUTINE: Validation of PCGs based on Pfam/PRINTS
my ($hmmfile, $pfam, $eval, $next, $pfamquery, $pfamhit, @selected, $printfile, $print, $printhit, $printeval, %printcollect, $printquery, @printcollect, $printtop, @printtop);
sub validatePCG {
#	print "Identifying protein sequences based on PFAM/PRINTS profiles\n";
	($in_sub) = @_;
	($printname = $in_sub) =~ s/^.*\///g;
	for $hmmfile (keys %pfamhash) {
		($pfam = $hmmfile) =~ s/.hmm//g;
		$eval = $pfamhash{$hmmfile};
		system("$exedir/hmmer-3.1b2-linux-intel-x86_64/bin/hmmsearch --cpu $cpu --tblout $profdir/$pfam.$printname.out -E $eval $dbdir/$hmmfile $seqdir/$printname.pcg.aa.fa >> $profdir/hmmer.log 2>> $profdir/hmmer.err");
		open PFAM_OUT, "$profdir/$pfam.$printname.out";
		while (<PFAM_OUT>) {
			if ($_ =~ /E-value/) {
				$next = <PFAM_OUT>;
				$next = <PFAM_OUT>;
				if ($next !~ /^#/) {
					$pfamquery = (split (/\s+/, $next))[0];
					$pfamhit = (split (/\s+/, $next))[2];
					$pfamhit = $genenamehash{$pfamhit};
					push (@selected, "$pfamquery\t$pfamhit");
				}
			}
		}
		close (PFAM_OUT);
	}
	for $printfile (keys %printhash) {
		($print = $printfile) =~ s/.pval//g;
		$eval = $printhash{$printfile};
		system("$exedir/FingerPRINTScan_3596/fingerPRINTScan $dbdir/$printfile $seqdir/$printname.pcg.aa.fa -e $eval -o 8 -R > $profdir/$print.$printname.out 2>> $profdir/fpscan.log");
		open PRINTS_OUT, "$profdir/$print.$printname.out";
		while (<PRINTS_OUT>) {
			if ($_ =~ /^Sn\;/) {
				if ($printhit && $printeval) {
					$printcollect{$printeval} = $printquery."HIT".$printhit;
					$printhit = "";
					$printeval = "";
				}
				$printquery = (split(/\s+/, $_))[1];
			} elsif ($_ =~ /^1TBH/) {
				$printhit = (split(/\s+/, $_))[1];
				$printhit = $genenamehash{$printhit};
				$printeval = (split(/\s+/, $_))[2];
			}
		}
		if ($printeval) {
			$printcollect{$printeval} = $printquery."HIT".$printhit;
			$printhit = "";
			$printeval = "";
		}
		$printquery = "";
		foreach (sort {$a <=> $b} keys %printcollect) {
			push (@printcollect, $printcollect{$_});
		}
		if (scalar @printcollect > 0) {
			$printtop = shift @printcollect;
			@printtop = split(/HIT/, $printtop);
			push (@selected, "$printtop[0]\t$printtop[1]");
		}
		%printcollect=();
		@printcollect=();
	}
}

# SUBROUTINE: Validation of rRNA based on homology
my (@blastout, $pid, $aln, $query, $sub, $hit, $len, %annothash, %exists, $species, $start, $end, $strand, $gene, $seq);
sub validaterRNA {
#	print "Identifying rRNA sequences based on BLAST homology\n";
	($in_sub) = @_;
	($printname = $in_sub) =~ s/^.*\///g;
	@blastout = `$exedir/ncbi-blast-2.6.0+/bin/blastn -query $seqdir/$printname.rrna.fa -db $dbdir/db.rrna.classified.db -evalue 1e-5 -outfmt 6 -max_target_seqs 1 -num_threads $cpu 2>$profdir/blast.log`;
	foreach (@blastout) {
		$pid = (split(/\t/, $_))[2];
		$aln = (split(/\t/, $_))[3];
		$query = (split(/\t/, $_))[0];
		$sub = (split(/\t/, $_))[1];
		$hit = (split(/\-/, $sub))[0];
		$len = (split(/\-/, $sub))[1];
		($id, $species, $start, $end, $strand, $gene) = split(/\#/, $query);
		$annothash{"$id\#$species\#$start\#$end\#$strand"} = $hit;
#		$cov = $aln/$len;
		if (! exists $exists{$query}) {
			$exists{$query} = 1;
#			open OUT_INDRRN, ">>$seqdir/$printname.$hit.nt.fa";
			($id, $species, $start, $end, $strand, $gene, $seq) = split(/\#/, $allrrn{$query});
			push (@selected, "$query\t$hit");
#			print OUT_INDRRN ">",$acc,"_",$species,"\n",$seq,"\n";
			$id=$species=$start=$end=$seq=$strand="";
#			close(OUT_INDRRN);
		}
	}
}

# SUBROUTINE: Check for missing genes
my (%withannot, @missinggenes, $missingstr, %selected, %collectmissing, @query);
sub findMissing {
	($in_sub) = @_;
	($printname = $in_sub) =~ s/^.*\///g;
	foreach (@selected) {
		($query, $hit) = split(/\t/, $_);
		@query = split(/\#/, $query);
		pop @query;
		$acc = $query[0];
		$species = $query[1];
		$query = join("#", @query);
		$withannot{$query} = $hit;
		$selected{$hit} = $query;
	}
	foreach (keys %genes) {
		if (! exists $selected{$_}) {
			push (@missinggenes, $_);
		}
	}
	if (scalar @missinggenes > 0) {
		$missingstr = join("\,", sort @missinggenes);
		$collectmissing{$acc."\#".$species} = $missingstr;
	}
	@missinggenes=();
	%selected=();
	@selected=();
}

# SUBROUTINE: Validation of PCGs based on name
my ($ans, $missgene, @altgenes, %altgenes, $genename, @notfound, %notfound, $notfoundstr, %notfoundmissing);
sub validateName {
	if (%collectmissing) {
		print "\n>>>WARNING: PFAM/PRINT/rRNA searches failed to identify genes for the following samples:\n";
		foreach (sort keys %collectmissing) {
			($acc, $species) = split(/\#/, $_);
			print "Missing genes for $species ($acc): $collectmissing{$_}\n";
		}
		print ">>>END\n";
		print "\nThese gene sequences may or may not be present in the GenBank/EMBL files, you may want to manually check to confirm this.\n";
		print "In the next step, MitoPhAST will attempt to identify missing sequences based on gene names provided in the input files.\n";
		if ($interactive eq "on") {
			print "Proceed? (Y/N)";
			$ans = <STDIN>;
			chomp $ans;
			if ($ans !~ /^Y$/i) {
				exit;
			}
		}
		foreach (sort keys %collectmissing) {
			($acc, $species) = split(/\#/, $_);
			print "\nIdentifying gene sequences for $species ($acc) based on gene IDs in GenBank/EMBL files\n";
			@missinggenes = split(/\,/, $collectmissing{$_});
			foreach $missgene (@missinggenes) {
				@altgenes = split(/\//, $genes{$missgene});
				$altgenes{$_} = 1 foreach @altgenes;
				foreach (keys %allaa) {
					if ($_ =~ /^$acc\#/) {
						$genename = (split(/\#/, $_))[-1];
						$genename =~ tr/a-z/A-Z/;
						if ($genename eq $missgene) {
							if (exists $altgenes{$genename}) {
#								push (@found, $missgene);
								@query = split(/\#/, $_);
								pop @query;
								$query = join("#", @query);
								$withannot{$query} = $missgene;
							} else {
								if (! exists $notfound{$missgene}) {
									push (@notfound, $missgene);
									$notfound{$_} = 1 foreach @notfound;
								}
							}
						}
					}
				}
				foreach (keys %allrrn) {
					if ($_ =~ /^$acc\#/) {
						$genename = (split(/\#/, $_))[-1];
						$genename =~ tr/a-z/A-Z/;
						if ($genename eq $missgene) {
							if (exists $altgenes{$genename}) {
#								push (@found, $missgene);
								@query = split(/\#/, $_);
								pop @query;
								$query = join("#", @query);
								$withannot{$query} = $missgene;
							} else {
								if (! exists $notfound{$missgene}) {
									push (@notfound, $missgene);
									$notfound{$_} = 1 foreach @notfound;
								}
							}
						}
					}
				}
			}
#			if (scalar @found > 0) {
#				$foundstr = join(",", @found);
#				$foundmissing{$accession} = $foundstr;
#				@found = ();
#			}
			if (scalar @notfound > 0) {
				$notfoundstr = join(",", @notfound);
				$notfoundmissing{$acc} = $notfoundstr;
				@notfound = ();
			}
		}
		if (%notfoundmissing) {
			print "\n>>>WARNING: Identification by gene ID still failed to identify some genes for the following samples:\n";
			foreach (sort keys %notfoundmissing) {
				print "Missing genes for $_: $notfoundmissing{$_}\n";
			}
			print ">>>END\n";
			print "\nThese gene sequences may or may not be present in the GenBank/EMBL files, you may want to manually check to confirm this.\n";
			print "If these sequences are present, it is likely that the gene tag (for CDS) is missing or it does not conform to IDs recognized by MitoPhAST (see README.txt under Gene ID list section)\n";
			print "If you choose to continue, absent genes will be treated as missing/unknown.\n";
			if ($interactive eq "on") {
				print "Proceed? (Y/N)";
				$ans = <STDIN>;
				chomp $ans;
				if ($ans !~ /^Y$/i) {
					exit;
				}
			}	
		} else {
			print "Gene ID successfully identified all missing genes for all samples.\n";
		}
	} else {
		print "\nPFAM/PRINTS/rRNA searches successfully identified all genes for all samples.\n";
	}
}


# SUBROUTINE: Print sequences to file according to genes
my ($printgene, @string, $string, $newstring);
sub printSeqs {
	foreach $printgene (sort keys %genes) {
		if ($printgene !~ /^rrn/) {
			open PRINTAA, ">$seqdir/all_$printgene.aa.fa";
			open PRINTNT, ">$seqdir/all_$printgene.nt.fa";
			foreach (sort keys %allaa) {
				@string = split(/\#/, $_);
				pop @string;
				$string = join("#", @string);
				$newstring = $string[0]."_".$string[1];
				if ($withannot{$string} eq $printgene) {
					print PRINTAA ">$newstring\n$allaa{$_}\n";
					print PRINTNT ">$newstring\n$allnt{$_}\n";
#				} elsif ($withannot{$string} =~ /----------/) {
#					print PRINTAA ">$string\n----------\n";
#					print PRINTNT ">$string\n----------\n";
				}
			}
		} else {
			open PRINTRRN, ">$seqdir/all_$printgene.nt.fa";
			foreach (sort keys %allrrn) {
				@string = split(/\#/, $_);
				pop @string;
				$string = join("#", @string);
				$newstring = $string[0]."_".$string[1];
				if ($withannot{$string} eq $printgene) {
					print PRINTRRN ">$newstring\n$allrrn{$_}\n";
#				} elsif ($withannot{$string} =~ /----------/) {
#					print PRINTRRN ">$string\n----------\n";
				}
			}
		}
	}
}


# SUBROUTINE: Save gene order in tab files
my (@split, $header, %check, @all, $order);
sub gene_order {
	open GENE, ">$seqdir/sample_gene_order.txt";
	foreach (@annotarray) {
		@split = split(/\t/, $_);
		$acc = (split(/\#/, $split[0]))[0];
		$species = (split(/\#/, $split[0]))[1];
		$header = "$species|$acc";
		$strand = (split(/\#/, $split[0]))[4];
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
		if (scalar @split == 1) {
			if (exists $withannot{$split[0]}) {
				push (@all, $strand.$withannot{$split[0]});
			}
		} elsif (scalar @split == 2) {
			push (@all, $strand.$split[1]);
		}
	}
	$order = join("\t", @all);
	print GENE "$order\n";
	@all = ();
}
