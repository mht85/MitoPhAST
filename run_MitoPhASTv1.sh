#!/bin/bash

bootstrap="autoMRE"
mlrep="1"
threads="2"

function checkargs {
if [[ $OPTARG =~ ^-[i/n/m/B/R/T/I/S/r]$ ]]
then
	echo "$OPTARG is an invalid argument to -$opt"
	exit
fi
}

while getopts "i:n:B:R:T:mISr" opt
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
B)
	checkargs
	bootstrap=$OPTARG;;
R)
	checkargs
	mlrep=$OPTARG;;
T)
	checkargs
	threads=$OPTARG;;
I)
	noninteractive="int";;
S)
	supermatstop="stop";;
r)
	rapidBS="rapid";;
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
		echo "run_MitoPhASTv1.sh -i example/input/ -n HG942366,HG799086 -m -I -B 1000 -R 5 -T 15"
		echo
		echo "-i <path>	input directory containing Genbank and/or EMBL files [at least one of -i and -n is required]"
		echo "-n <string>	comma-separated string of mitogenome accession IDs to download from NCBI [at least one of -i or -n is required]"
		echo "-m		optional argument to exclude genes which are absent in some files (default: include all genes, missing genes are substituted with hyphens)"
		echo "-I		turns off interactive prompt for user to check for missing genes (default: interactive) If prompt is disabled, program will only provide warning to standard output."
		echo "-S		stops program after supermatrix construction, turns off model estimation and ML analysis" 
		echo "-B <int>	number of bootstrap replicates (default: autoMRE)"
		echo "-R <int>	number of ML trees generated if standard bootstrapping is used - one tree out of <int> with best likelihood score is produced (default: 1"
		echo "-T <int>	number of threads (default: 2, must be greater than 1)"
		echo "-r		run ML analysis with rapid bootstrapping (default: standard bootstrapping)"
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
	echo "Number of bootstrap replicates: $bootstrap"
	echo "Number of ML trees generated (if standard bootstrapping is used): $mlrep"
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
	mkdir $out_dir/phylogeny/
	mkdir $out_dir/phylogeny/align_trim/
	mkdir $out_dir/phylogeny/model_select/
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


########## PARTITIONING AND MODEL SELECTION ##########

cd $out_dir/phylogeny/model_select/
cp ../align_trim/all_13pcg.fasta .
cp ../align_trim/all_13pcg.phy .

if [ ! -s all_13pcg.fasta ]; then
        echo "No fasta file found for protein data, script terminated!"
        exit
elif [ ! -s all_13pcg.phy ]; then
        echo "No phylip alignment found for protein data, script terminated!"
        exit
fi

echo "Running ProtTest"
echo

cd $programs/prottest-3.4-20140123/
java -Xmx2G -jar prottest-3.4.jar -i $out_dir/phylogeny/align_trim/all_13pcg.fasta -o $out_dir/phylogeny/model_select/all_13pcg.prottest.out -all-distributions -F -AIC -BIC -AICC -all -S 2 -threads $threads > $out_dir/phylogeny/model_select/prottest.log 2> $out_dir/phylogeny/model_select/prottest.err
prottest=`fgrep "Best model according to AIC:" $out_dir/phylogeny/model_select/all_13pcg.prottest.out | awk -F": " '{print $2}'`
gamma=`echo $prottest | awk '{if ($0 ~ /\+G/) print "GAMMA"; else print "CAT"}'`
freq=`echo $prottest | awk '{if ($0 ~ /\+F/) print "F"}'`
model=`echo $prottest | awk -F"+" '{print $1}' | perl -pe 'tr/[a-z]/[A-Z]/'`

cd $out_dir/phylogeny/model_select/

echo "Partitioning data"
echo

cat /dev/null > partitions.tmp.txt

for f in ../align_trim/all_*trimal.fasta-2
do
name=`echo $f | sed 's/^.*\///g' | sed 's/all_//g' | perl -pe 's/.cds.fa.clustalo-trimal.fasta-2//g'`
len=`grep -v ">" $f | head -1 | awk '{print length}'`
echo "$name $len" >> partitions.tmp.txt
done

while read line
do
gene=`echo $line | awk '{print $1}'`
length=`echo $line | awk '{print $2}'`
if [ -z "$start" ]; then
	stop=$((1 + length - 1))
	echo "$model$freq, $gene = 1-$stop" > partitions.txt
	start=$((stop + 1))
else
	stop=$((start + length - 1))
	echo "$model$freq, $gene = $start-$stop" >> partitions.txt
	start=$((stop + 1))
fi
done<"partitions.tmp.txt"

unset start
cp partitions.txt ../ml_analysis/

########## RAxML ##########

cd ../ml_analysis

cp $out_dir/phylogeny/align_trim/all_13pcg.phy $out_dir/phylogeny/ml_analysis/
cp $out_dir/phylogeny/model_select/partitions.txt $out_dir/phylogeny/ml_analysis/

# Running RAxML with rapid bootstrapping

if [ ! -z "$rapidBS" ]
then
	p_rand=$RANDOM
	bx_rand=$RANDOM

	mkdir $out_dir/phylogeny/ml_analysis/tree
	cd $out_dir/phylogeny/ml_analysis/tree

	echo "Running RAxML with rapid bootstrapping with command: "
	echo "$programs/standard-RAxML-master/raxml -f a -m PROT$gamma$model$freq -q $out_dir/phylogeny/ml_analysis/partitions.txt -p $p_rand -x $bx_rand -# $bootstrap -s $out_dir/phylogeny/ml_analysis/all_13pcg.phy -n R1 -T $threads > $out_dir/phylogeny/ml_analysis/tree/raxml.log 2>&1"
	echo

	$programs/standard-RAxML-master/raxml -f a -m PROT$gamma$model$freq -q $out_dir/phylogeny/ml_analysis/partitions.txt -p $p_rand -x $bx_rand -# $bootstrap -s $out_dir/phylogeny/ml_analysis/all_13pcg.phy -n R1 -T $threads > $out_dir/phylogeny/ml_analysis/tree/raxml.log 2>&1

	if [ ! -s RAxML_bipartitions.R1 ]; then
		echo "RAxML failed, script terminated!"
		exit
	fi

	cp RAxML_bipartitions.R1 MLtree.renamed.tre

	while read line
	do
		tosub=`echo "$line" | awk -F"\t" '{print $1}'`
		subwith=`echo "$line" | awk -F"\t" '{print $2}'`
		sed "s/$tosub\:/$subwith\:/g" MLtree.renamed.tre > tempMLfile
		mv tempMLfile MLtree.renamed.tre
	done<"$out_dir/summaries/specieslist"
fi

# Running RAxML with standard bootstrapping

if [ -z "$rapidBS" ]
then
	mkdir $out_dir/phylogeny/ml_analysis/tree/
	cd $out_dir/phylogeny/ml_analysis/tree/
	touch RAxML_bestTree.S1-all

	# Generate ML trees
	for ((rep=1; rep<=$mlrep; rep++))
	do
		p_rand=$RANDOM
		echo "Running RAxML tree$rep with commands: "
		echo "$programs/standard-RAxML-master/raxml -m PROT$gamma$model$freq -q $out_dir/phylogeny/ml_analysis/partitions.txt -p $p_rand -s $out_dir/phylogeny/ml_analysis/all_13pcg.phy -n S1-$rep -T $threads > $out_dir/phylogeny/ml_analysis/tree/raxml.S1-$rep.log 2>&1"
		$programs/standard-RAxML-master/raxml -m PROT$gamma$model$freq -q $out_dir/phylogeny/ml_analysis/partitions.txt -p $p_rand -s $out_dir/phylogeny/ml_analysis/all_13pcg.phy -n S1-$rep -T $threads > $out_dir/phylogeny/ml_analysis/tree/raxml.S1-$rep.log 2>&1
		cat RAxML_bestTree.S1-$rep >> RAxML_bestTree.S1-all
	done

	# Find tree with best likelihood
	$programs/standard-RAxML-master/raxml -f N -m PROT$gamma$model$freq -q $out_dir/phylogeny/ml_analysis/partitions.txt -s $out_dir/phylogeny/ml_analysis/all_13pcg.phy -n fN -z RAxML_bestTree.S1-all -T $threads > raxml.S1-fN.log 2>&1
	best=`grep -P "^\d+" RAxML_info.fN | head -1 | awk '{print $1+1}'`
	echo
	echo "Tree number $best out of $mlrep has best likelihood score. See $out_dir/phylogeny/ml_analysis/tree/RAxML_info.fN for more details"
	echo

	# Standard bootstrapping
	p_rand=$RANDOM
	bx_rand=$RANDOM	
	echo "$programs/standard-RAxML-master/raxml -m PROT$gamma$model$freq -q $out_dir/phylogeny/ml_analysis/partitions.txt -p $p_rand -s $out_dir/phylogeny/ml_analysis/all_13pcg.phy -n S2-bs$bootstrap -b $bx_rand -# $bootstrap -T $threads > raxml.S2-bs$bootstrap.log 2>&1"
	echo
	$programs/standard-RAxML-master/raxml -m PROT$gamma$model$freq -q $out_dir/phylogeny/ml_analysis/partitions.txt -p $p_rand -s $out_dir/phylogeny/ml_analysis/all_13pcg.phy -n S2-bs$bootstrap -b $bx_rand -# $bootstrap -T $threads > raxml.S2-bs$bootstrap.log 2>&1

	# Map bootstrap to best tree
	p_rand=$RANDOM
	echo "$programs/standard-RAxML-master/raxml -f b -m PROT$gamma$model$freq -t RAxML_bestTree.S1-$best -p $p_rand -z RAxML_bootstrap.S2-bs$bootstrap -n S3-t$best-bs$bootstrap -T $threads > raxml.S3-t$best-bs$bootstrap.log 2>&1"
	echo
	$programs/standard-RAxML-master/raxml -f b -m PROT$gamma$model$freq -t RAxML_bestTree.S1-$best -p $p_rand -z RAxML_bootstrap.S2-bs$bootstrap -n S3-t$best-bs$bootstrap -T $threads > raxml.S3-t$best-bs$bootstrap.log 2>&1

	if [ ! -s RAxML_bipartitions.S3-t$best-bs$bootstrap ]; then
		echo "RAxML with standard bootstrapping failed, script terminated!"
		exit
	fi
	cp RAxML_bipartitions.S3-t$best-bs$bootstrap MLtree.renamed.tre

	while read line
	do
		tosub=`echo "$line" | awk -F"\t" '{print $1}'`
		subwith=`echo "$line" | awk -F"\t" '{print $2}'`
		sed "s/$tosub\:/$subwith\:/g" MLtree.renamed.tre > tempMLfile
		mv tempMLfile MLtree.renamed.tre
	done<"$out_dir/summaries/specieslist"
fi

########## Reiterate missing genes by PFAM/PRINTS and Name ##########

echo "----------PROGRAM RUN COMPLETED----------"
echo
echo "ML tree file with bootstrap values can be found in file $out_dir/phylogeny/ml_analysis/tree/MLtree.renamed.tre"
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
