#!/bin/bash

threads="1"
bb="1000"
alrt="1000"

function checkargs {
if [[ $OPTARG =~ ^-[i/n/m/B/T/I/S]$ ]]
then
	echo "$OPTARG is an invalid argument to -$opt"
	exit
fi
}

while getopts "i:n:T:O:b:a:mIS" opt
do
case $opt in
i)
	checkargs
	in_dir_tmp=`echo $OPTARG | sed 's/\/$//g'`;;
n) 
	checkargs
	download=`echo $OPTARG | sed 's/\,/ /g'`;;
m)
	nomiss="nomiss";;
T)
	checkargs
	threads=$OPTARG;;
I)
	noninteractive="int";;
S)
	supermatstop="stop";;
b)
	checkargs
	bb=`echo $OPTARG`;;
a)
	checkargs
	alrt=`echo $OPTARG`;;
O)
	checkargs
	ORI=`echo $OPTARG`;;
:)
	echo "option -$OPTARG requires an argument"
	exit
esac
done


if [ -z "$in_dir_tmp" ]
then
	if [ -z "$download" ]
	then
		echo "MitoPhAST usage example:"
		echo
		echo "run_MitoPhASTv1.sh -i example/input/ -n HG942366,HG799086 -T 15 -O COX1"
		echo
		echo "-i <path>	input directory containing Genbank and/or EMBL files [at least one of -i and -n is required]"
		echo "-n <string>	comma-separated string of mitogenome accession IDs to download from NCBI [at least one of -i or -n is required]"
		echo "-m		optional argument to exclude genes which are absent in some files (default: include all genes, missing genes are substituted with hyphens)"
		echo "-I		turns off interactive prompt for user to check for missing genes (default: interactive) If prompt is disabled, program will only provide warning to standard output."
		echo "-S		stops program after supermatrix construction, turns off model estimation and ML analysis"
		echo "-b <int>	number of ultrafast bootstrap replicates run on IQ-TREE (default: 1000, must be >= 1000)"
		echo "-a <int>	number of SH-aLRT replicates run on IQ-TREE (default: 1000)"
		echo "-T <int>	number of threads (default: 1)"
		echo "-O <string>	gene to orient to for linear gene order comparison analysis (if not provided, analysis is skipped), see README for allowed genes"
		echo
		exit
	fi
fi

if [ ! -z "$in_dir_tmp" ]
then
	echo "Using mitogenome files from $in_dir_tmp/"
	echo
fi
if [ ! -z "$download" ]
then
	echo "Mitogenome files with accession IDs $download will be downloaded."
	echo
fi
if [ ! -z "$nomiss" ]
then
	echo "Genes absent in some organisms will be excluded from analysis completely! (-m option)"
	echo
fi
if [ -z "$supermatstop" ]
then
	echo "Number of ultrafast bootstrap replicates: $bb"
	echo "Number of SH-aLRT replicates: $alrt"
	echo "Number of threads used: $threads"
	echo
fi

path=$(cd "$(dirname "$0")"; pwd)
scripts=$path/scripts
programs=$path/programs
cd $in_dir_tmp
in_dir=`pwd`
cd $path
out_dir=$in_dir-out

if [ ! -d $out_dir ]; then
	mkdir $out_dir
	mkdir $out_dir/downloads/
	mkdir $out_dir/input_tmp/
	mkdir $out_dir/sequences/
	mkdir $out_dir/summaries/
	mkdir $out_dir/profiles/
	mkdir $out_dir/gene_orders/
	mkdir $out_dir/phylogeny/
	mkdir $out_dir/phylogeny/align_trim/
	mkdir $out_dir/phylogeny/ml_analysis/
	echo "Output directory $out_dir created"
	echo
else
	echo "Directory $out_dir already exists!"
	echo
	exit
fi

########## PULLING GENBANK FILES BASED ON ACCESSION ##########

if [ ! -z "$download" ]
then
	cd $out_dir/downloads/
	echo "Downloading mitogenomes for $download to $out_dir/downloads/"
	echo
	python $scripts/pull_genomes_as_genbank.py -a $download > download.log 2>&1
fi
cd $path

########## SUMMARIZING INFORMATION FROM GENBANK/EMBL FILES AND EXTRACTING SEQUENCES ##########

if [ ! -z "$in_dir_tmp" ]
then
	for file in $in_dir/*
	do
	name=`echo $file | sed 's/^.*\///g' | sed 's/ /_/g'`
	noacc=`grep "^ACCESSION" $file | awk '{if (NF == "1") print "none"}'`
	noacc2=`grep "^AC" $file | perl -pe 's/\;//g' | awk '{if (NF == "1") print "none"}'`
	if [ "$noacc" = "none" ]
	then
		name2=`echo "PSEUDO"$RANDOM`
		sed 's/D-loop/misc_feature/g' $file | sed "s/LOCUS       /LOCUS       $name2/g" | sed "s/ACCESSION/ACCESSION   $name2/g" > $out_dir/input_tmp/$name
	elif [ "$noacc2" = "none" ]
	then
		name2=`echo "PSEUDO"$RANDOM`
		sed 's/D-loop/misc_feature/g' $file | sed "s/^ID   /ID   $name2;/g" | sed "s/^AC   /AC   $name2/g" > $out_dir/input_tmp/$name
	else
		sed 's/D-loop/misc_feature/g' $file > $out_dir/input_tmp/$name
	fi
	done
fi

if [ ! -z "$download" ]
then
	for file in $out_dir/downloads/*gb*
	do
	name=`echo $file | sed 's/^.*\///g' | sed 's/ /_/g'`
	sed 's/D-loop/misc_feature/g' $file > $out_dir/input_tmp/$name
	done
fi

echo "Printing summary files to $out_dir/summaries/"
echo
perl $scripts/summarize.pl -input_dir $out_dir/input_tmp/ -output_dir $out_dir/summaries

if [ ! -z "$ORI" ]
then
	echo "Clustering gene orders, results in $out_dir/gene_orders/"
	echo
	perl $scripts/prep_order_in.pl -input_dir $out_dir/input_tmp/ -output_dir $out_dir/gene_orders
	perl $scripts/check_order.pl -in $out_dir/gene_orders/sample_gene_order.txt -out $out_dir/gene_orders/ -ori $ORI
	
fi

echo "Extracting sequences to $out_dir/sequences/"
echo

perl $scripts/extract_seqs.pl -input_dir $out_dir/input_tmp/ -output_dir $out_dir/sequences
cat $out_dir/summaries/*txt | grep -P "(AccessionID|Species)" | perl -pe 's/\n/\t/g' | perl -pe 's/\tAccessionID/\nAccessionID/g' | awk -F"\t" '{print $2"\t"$4" ["$2"]"}' > $out_dir/summaries/specieslist

########## GROUPING SEQUENCES ACCORDING TO PFAM AND PRINT ##########

echo "Identifying protein sequences based on PFAM/PRINTS profiles"
echo
cd $out_dir/profiles/
mkdir hmmer_out
mkdir fpscan_out

for seqfilepath in $out_dir/sequences/*.cds.fa
do
	seqfile=`echo $seqfilepath | sed 's/^.*\///g'`

# PFAM domain identification

	for hmmfile in PF00032.hmm PF00115.hmm PF00116.hmm PF00119.hmm PF00146.hmm PF00499.hmm PF00510.hmm PF0507.hmm
	do
		$programs/hmmer-3.0/src/hmmsearch --cpu $threads --tblout $out_dir/profiles/hmmer_out/$hmmfile.$seqfile.out -E 1e-5 $path/hmmer_print_profiles/$hmmfile $seqfilepath >> $out_dir/profiles/hmmer_out/hmmer.log 2>> $out_dir/profiles/hmmer_out/hmmer.err
	done

	for hmmfile in PF00895.hmm
	do
		$programs/hmmer-3.0/src/hmmsearch --cpu $threads --tblout $out_dir/profiles/hmmer_out/$hmmfile.$seqfile.out -E 1e-3 $path/hmmer_print_profiles/$hmmfile $seqfilepath >> $out_dir/profiles/hmmer_out/hmmer.log 2>> $out_dir/profiles/hmmer_out/hmmer.err
	done

	for hmmfile in PF00420.hmm
	do
		$programs/hmmer-3.0/src/hmmsearch --cpu $threads --tblout $out_dir/profiles/hmmer_out/$hmmfile.$seqfile.out -E 1e-2 $path/hmmer_print_profiles/$hmmfile $seqfilepath >> $out_dir/profiles/hmmer_out/hmmer.log 2>> $out_dir/profiles/hmmer_out/hmmer.err
	done

	cat $out_dir/profiles/hmmer_out/*.$seqfile.out | grep -A2 "target" | grep -v "^#" | grep -vP "^--$" | awk '{print $1"\t"$3}' > $seqfile.select

# PRINTS domain identification

	for printsfile in PR01434.pval PR01436.pval PR01437.pval
	do
		$programs/FingerPRINTScan_3596/fingerPRINTScan $path/hmmer_print_profiles/$printsfile $seqfilepath -e 1e-5 -o 8 -R > $out_dir/profiles/fpscan_out/$printsfile.$seqfile.out 2>> $out_dir/profiles/fpscan_out/fpscan.log
	done

	grep "PRINT_ND5" $out_dir/profiles/fpscan_out/PR01434.pval.$seqfile.out | sort -gk3 | head -1 | fgrep -B3 -wf - $out_dir/profiles/fpscan_out/PR01434.pval.$seqfile.out | head -1 | awk '{print $2"\tNAD5"}' > $seqfile.select2
	grep "NADHDHGNASE2" $out_dir/profiles/fpscan_out/PR01436.pval.$seqfile.out | sort -gk3 | head -1 | fgrep -B3 -wf - $out_dir/profiles/fpscan_out/PR01436.pval.$seqfile.out | head -1 | awk '{print $2"\tNAD2"}' >> $seqfile.select2
	grep "PRINT_ND4" $out_dir/profiles/fpscan_out/PR01437.pval.$seqfile.out | sort -gk3 | head -1 | fgrep -B3 -wf - $out_dir/profiles/fpscan_out/PR01437.pval.$seqfile.out | head -1 | awk '{print $2"\tNAD4"}' >> $seqfile.select2

# Combine PFAM & PRINTS results

	cat $seqfile.select $seqfile.select2 | sed 's/Cytochrom_B_C/COB/g' | sed 's/ATP-synt_A/ATP6/g' | sed 's/NADHdh/NAD1/g' | sed 's/Oxidored_q2/NAD4L/g' | sed 's/Oxidored_q3/NAD6/g' | sed 's/ATP-synt_8/ATP8/g' | sed 's/Oxidored_q4/NAD3/g' | sort -k2 | uniq > ${seqfile%%fa}final

	rm $seqfile.select*
done

# Check for protein sequences missed by PFAM/PRINTS identification

cat /dev/null > $out_dir/profiles/missed1

for file in $out_dir/profiles/*.cds.final
do
	file2=`echo $file | awk -F"/" '{print $NF}' | perl -pe 's/\.\w+\.cds\.final//g'`
	miss1=`awk -F"\t" '{print $2}' $file | grep -P "(ATP6|ATP8|COX1|COX2|COX3|COB|NAD1|NAD2|NAD3|NAD4|NAD4L|NAD5|NAD6)" | sort | uniq | cat - $path/hmmer_print_profiles/pcg.list | sort | uniq -c | awk '{if ($1 != 2) print $2}' | perl -pe 's/\n/\, /g' | perl -pe 's/\, $//g'`
	if [ ! -z "$miss1" ]
	then
		echo -e "$file2\t(missing: $miss1)" >> $out_dir/profiles/missed1
	fi
done

if [ -s $out_dir/profiles/missed1 ]
then
	echo "####################"
	echo "WARNING: PFAM/PRINTS profiles failed to identify some protein sequences for the following samples:"
	echo "$(cat $out_dir/profiles/missed1)"
	echo "####################"
	echo
	echo "These protein sequences may or may not be present in the GenBank/EMBL files, you may want to manually check to confirm this."
	echo "In the next step, MitoPhAST will attempt to identify missing sequences based on gene names provided in the input files."
	if [ -z "$noninteractive" ]
	then
		echo -n "Would you like to proceed (Y/N)?"
		read answer
		if test "$answer" != "Y" -a "$answer" != "y";
		then
			exit 0;
		fi
	fi
else
	echo "PFAM/PRINTS profiles successfully identified all protein sequences for all samples."
fi
echo

########## GROUPING SEQUENCES ACCORDING TO ANNOTATED GENE ID ##########

# Substitute names in sequences and group according to genes

echo "Identifying protein sequences based on gene IDs in GenBank/EMBL files"
perl $scripts/rename_group_seqs.pl -input_profiles_dir $out_dir/profiles/ -input_sequences_dir $out_dir/sequences/ -output_renamed_dir $out_dir/phylogeny/align_trim/ -output_grouped_dir $out_dir/phylogeny/align_trim/


if [ -s $out_dir/profiles/missed1 ]
then
# Output sequences 'rescued' by this second round identification
	cat $out_dir/phylogeny/align_trim/*.cds.renamed.fa | perl -pe 's/\n/\t/g' | perl -pe 's/\t>/\n>/g' > $out_dir/profiles/temp.fasta
	cat /dev/null > $out_dir/profiles/rescued_proteins.fa

	while read line
	do
		missedacc=`echo "$line" | awk -F"\t" '{print $1}' | perl -pe 's/\.\d+$//g'`
		missedgene=`echo "$line" | awk -F"\t" '{print $2}' | sed 's/(missing: //g' | sed 's/)//g' | perl -pe 's/\, / /g'`
		for var in $missedgene
		do
			fgrep -w ">$missedacc-$var" $out_dir/profiles/temp.fasta | perl -pe 's/\t/\n/g' >> $out_dir/profiles/rescued_proteins.fa
		done
	done<"$out_dir/profiles/missed1"
	rm $out_dir/profiles/temp.fasta

# Check for protein sequences missed by this second round identification
	cat /dev/null > $out_dir/profiles/missed2

	for file in $out_dir/phylogeny/align_trim/*.cds.renamed.fa
	do
	file2=`echo $file | awk -F"/" '{print $NF}' | perl -pe 's/\.\w+\.cds\.renamed.fa//g'`
	miss2=`grep ">" $file | awk -F"-" '{print $NF}' | grep -P "(ATP6|ATP8|COX1|COX2|COX3|COB|NAD1|NAD2|NAD3|NAD4|NAD4L|NAD5|NAD6)" | sort | uniq | cat - $path/hmmer_print_profiles/pcg.list | sort | uniq -c | awk '{if ($1 != 2) print $2}' | perl -pe 's/\n/\, /g' | perl -pe 's/\, $//g'`
	if [ ! -z "$miss2" ]
	then
		echo -e "$file2\t(missing: $miss2)" >> $out_dir/profiles/missed2
	fi
	done	
	if [ -s $out_dir/profiles/missed2 ]
	then
		echo "####################"
		echo "WARNING: Identification by gene ID still failed to identify some protein sequences for the following samples:"
		echo "$(cat $out_dir/profiles/missed2)"
		echo "####################"	
		echo
		echo "These protein sequences may or may not be present in the GenBank/EMBL files, you may want to manually check to confirm this."
		echo "If these sequences are present, it is likely that the gene tag (for CDS) is missing or it does not conform to IDs recognized by MitoPhAST (see README.txt under Gene ID list section)"
		echo "Protein sequences which were successfully identified by gene IDs can be found in $out_dir/profiles/rescued_proteins.fa"
		if [ -z "$noninteractive" ]
		then
			echo -n "Would you like to proceed to phylogenetic analysis (Y/N)?"
			read answer
			if test "$answer" != "Y" -a "$answer" != "y";
			then
				exit 0;
			fi
		fi
	else
	        echo "Identification by gene ID successfully identified all protein sequences for all samples."
		echo "Protein sequences which were successfully identified by gene IDs can be found in $out_dir/profiles/rescued_proteins.fa"
	fi
	echo
fi


########## ALIGNMENT AND TRIMMING ##########

echo "Multiple sequence alignment with Clustal Omega and trimming with trimAl"
echo

cd $out_dir/phylogeny/align_trim/
cat all_*.fa | grep ">" | awk -F"-" '{print $1}' | sed 's/>//g' | sort | uniq > acclist
for file in all_*.cds.fa
do
	$programs/clustal-omega/clustalo -i $file -o $file.clustalo.fasta >> clustalo.log
	$programs/trimAl/source/trimal -in $file.clustalo.fasta -out $file.clustalo-trimal.fasta -automated1 >> trimal.log
done

if [ ! -z $nomiss ]
then
	perl $scripts/add_gaps.pl -input_dir $out_dir/phylogeny/align_trim/ -species_list $out_dir/summaries/specieslist -no-missing -output_gap_dir $out_dir/phylogeny/align_trim/
else
	perl $scripts/add_gaps.pl -input_dir $out_dir/phylogeny/align_trim/ -species_list $out_dir/summaries/specieslist -output_gap_dir $out_dir/phylogeny/align_trim/
fi

while read line
do
	echo ">$line" > $line.cat.fasta
	cat all_*trimal.fasta-2 | perl -pe 's/\n/\t/g' | perl -pe 's/\t>/\n>/g' | fgrep ">$line-" | perl -pe 's/\t/\n/g' | grep -v ">" | perl -pe 's/\n//g' | awk '{print $0}' >> $line.cat.fasta
done<"acclist"

cat *.cat.fasta > all_13pcg.fasta
awk '$1=substr($1"       ",1)' FS="\n" OFS= RS=\> all_13pcg.fasta |awk 'NF>0{s=length($2);t=t"\n"$0}END{print NR-1,s"\n"t}' > all_13pcg.phy

if [ ! -z "$supermatstop" ]
then
	echo "----------PROGRAM RUN COMPLETED----------"
	echo
	echo "Supermatrix files all_13pcg.fasta and all_13pcg.phy can be found in $out_dir/phylogeny/align_trim/"
	echo
	if [ -s $out_dir/profiles/missed1 ]
	then
		echo "NOTE: PFAM/PRINTS profiles failed to identify some protein sequences for the following samples:"
		echo "$(cat $out_dir/profiles/missed1)"
		echo
		if [ -s $out_dir/profiles/missed2 ]
		then
			echo "NOTE: Identification by gene ID still failed to identify some protein sequences for the following samples:"
			echo "$(cat $out_dir/profiles/missed2)"
			echo
                fi
		echo "See above messages or README.txt for more information."
        fi
	echo "-----------------------------------------"
	echo
	rm $out_dir/profiles/missed1
	rm $out_dir/profiles/missed2
	exit
fi


########## PARTITIONING ##########

cd $out_dir/phylogeny/ml_analysis
cp ../align_trim/all_13pcg.fasta .
cp ../align_trim/all_13pcg.phy .

if [ ! -s all_13pcg.fasta ]; then
        echo "No fasta file found for protein data, script terminated!"
        exit
elif [ ! -s all_13pcg.phy ]; then
        echo "No phylip alignment found for protein data, script terminated!"
        exit
fi

echo "Generating Nexus file"
echo

cat /dev/null > partitions.tmp.txt

for f in ../align_trim/all_*trimal.fasta-2
do
	name=`echo $f | sed 's/^.*\///g' | sed 's/all_//g' | perl -pe 's/.cds.fa.clustalo-trimal.fasta-2//g'`
	len=`grep -v ">" $f | head -1 | awk '{print length}'`
	echo "$name $len" >> partitions.tmp.txt
done

echo "#NEXUS" > all_13pcg.nex
echo "begin sets;" >> all_13pcg.nex
echo "" >> all_13pcg.nex

while read line
do
	gene=`echo $line | awk '{print $1}'`
	length=`echo $line | awk '{print $2}'`
	if [ -z "$start" ]
	then
		stop=$((1 + length - 1))
		echo "	charset $gene=all_13pcg.phy:AA, 1-$stop;" >> all_13pcg.nex
		start=$((stop + 1))
	else
		stop=$((start + length - 1))
		echo "	charset $gene=all_13pcg.phy:AA, $start-$stop;" >> all_13pcg.nex
		start=$((stop + 1))
	fi
done<"partitions.tmp.txt"

echo "" >> all_13pcg.nex
echo "end;" >> all_13pcg.nex
unset start

########## IQ-TREE ##########

echo "Running IQ-TREE"
echo

cd ../ml_analysis

cp $out_dir/phylogeny/align_trim/all_13pcg.phy $out_dir/phylogeny/ml_analysis/

# Running IQ-TREE

if [ "$threads" -ge "13" ]
then
	$programs/iqtree-omp-1.5.5/bin/iqtree-omp -spp all_13pcg.nex -nt 13 -bb $bb -alrt $alrt -m TESTMERGE > iqtree.log 2>&1
else
	$programs/iqtree-omp-1.5.5/bin/iqtree-omp -spp all_13pcg.nex -nt $threads -bb $bb -alrt $alrt -m TESTMERGE > iqtree.log 2>&1
fi

if [ ! -s all_13pcg.nex.treefile ]
then
	echo "IQ-TREE failed, script terminated!"
	exit
fi

cp all_13pcg.nex.treefile MLtree.renamed.tre

while read line
do
	tosub=`echo "$line" | awk -F"\t" '{print $1}'`
	subwith=`echo "$line" | awk -F"\t" '{print $2}'`
	sed "s/$tosub\:/$subwith\:/g" MLtree.renamed.tre > tempMLfile
	mv tempMLfile MLtree.renamed.tre
done<"$out_dir/summaries/specieslist"


########## Reiterate missing genes by PFAM/PRINTS and Name ##########

echo "----------PROGRAM RUN COMPLETED----------"
echo
echo "ML tree file with bootstrap values can be found in file $out_dir/phylogeny/ml_analysis/MLtree.renamed.tre"
echo
if [ -s $out_dir/profiles/missed1 ]
then
	echo "NOTE: PFAM/PRINTS profiles failed to identify some protein sequences for the following samples:"
	echo "$(cat $out_dir/profiles/missed1)"
	echo
	if [ -s $out_dir/profiles/missed2 ]
	then
		echo "NOTE: Identification by gene ID still failed to identify some protein sequences for the following samples:"
		echo "$(cat $out_dir/profiles/missed2)"
		echo
		rm $out_dir/profiles/missed2
	fi
	echo "See above messages or README.txt for more information."
	rm $out_dir/profiles/missed1
fi
echo "-----------------------------------------"
echo
