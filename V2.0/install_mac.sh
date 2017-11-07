path=`pwd`
programs=$path/programs/

cd $programs
rm -rf hmmer-3.0
rm -rf argtable2-13
rm -rf clustal-omega/clustalo
rm -rf clustal-omega/clustal-omega-1.2.0
rm -rf iqtree-omp-1.5.5
rm -rf trimAl
rm -rf FingerPRINTScan_3596

# Unpack hmmer
tar -zxvf hmmer-3.0.tar.gz
cd hmmer-3.0
./configure --prefix=$programs/hmmer-3.0/ && make && make install
cd ..

if [ ! -s $programs/hmmer-3.0/bin/hmmsearch ]
then
echo "ERROR: hmmer installation failed"
exit
fi

# Unpack FingerPRINTScan
mkdir $programs/FingerPRINTScan_3596
ln -s $programs/fingerPRINTScan $programs/FingerPRINTScan_3596/fingerPRINTScan

if [ ! -s $programs/FingerPRINTScan_3596/fingerPRINTScan ]
then
echo "ERROR: FingerPRINTScan installation failed"
exit
fi

# Unpack Argtable2
tar -zxvf argtable2-13.tar.gz
cd argtable2-13
./configure --prefix=$programs/argtable2-13 && make && make install
cd ..

if [ ! -s $programs/argtable2-13/include/argtable2.h ]
then
echo "ERROR: argtable installation failed"
exit
elif [ ! -s $programs/argtable2-13/lib/libargtable2.a ]
then
echo "ERROR: argtable installation failed"
exit
fi

# Unpack ClustalOmega
cd clustal-omega
ln -s $programs/clustal-omega/clustal-omega-1.2.0-macosx clustalo
cd ..

# Unpack IQ-TREE
unzip iqtree-omp-1.5.5-MacOSX.zip
mv iqtree-omp-1.5.5-MacOSX iqtree-omp-1.5.5

if [ ! -s $programs/iqtree-omp-1.5.5/bin/iqtree-omp ]
then
echo "ERROR: IQ-TREE installation failed"
exit
fi

# Unpack trimAl and compile
tar -zxvf trimal.v1.2rev59.tar.gz
cd trimAl/source/
make

if [ ! -s $programs/trimAl/source/trimal ]
then
echo "ERROR: trimal installation failed"
exit
fi
