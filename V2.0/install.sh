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

# Install hmmer
tar -zxvf hmmer-3.0.tar.gz
cd hmmer-3.0
./configure --prefix=$programs/hmmer-3.0/ && make && make install
cd ..

if [ ! -s $programs/hmmer-3.0/bin/hmmsearch ]
then
echo "ERROR: hmmer installation failed"
exit
fi

# Install FingerPRINTScan
tar -zxvf FingerPRINTScan_3596.tgz
cd FingerPRINTScan_3596
./configure --prefix=$programs/FingerPRINTScan_3596 && make && make install
cd ..

if [ ! -s $programs/FingerPRINTScan_3596/fingerPRINTScan ]
then
echo "ERROR: FingerPRINTScan installation failed"
exit
fi

# Install Argtable2
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

# Install ClustalOmega
cd clustal-omega
tar -zxvf clustal-omega-1.2.0.tar.gz
cd clustal-omega-1.2.0/
./configure --prefix=$programs/clustal-omega/clustal-omega-1.2.0/ CFLAGS="-I$programs/argtable2-13/include/" LDFLAGS="-L$programs/argtable2-13/lib" && make && make install
cd ..

if [ ! -s $programs/clustal-omega/clustal-omega-1.2.0/bin/clustalo ]
then
echo "ERROR: clustal omega installation failed"
exit
fi

ln -s $programs/clustal-omega/clustal-omega-1.2.0/bin/clustalo clustalo
cd ..

# Install IQ-TREE
tar -zxvf iqtree-omp-1.5.5-Linux.tar.gz
mv iqtree-omp-1.5.5-Linux iqtree-omp-1.5.5

if [ ! -s $programs/iqtree-omp-1.5.5/bin/iqtree-omp ]
then
echo "ERROR: IQ-TREE installation failed"
exit
fi

# Install trimAl
tar -zxvf trimal.v1.2rev59.tar.gz
cd trimAl/source/
make

if [ ! -s $programs/trimAl/source/trimal ]
then
echo "ERROR: trimal installation failed"
exit
fi