path=`pwd`
programs=$path/programs/

cd $programs
rm -rf hmmer-3.0
rm -rf argtable2-13
rm -rf clustal-omega/clustalo
rm -rf clustal-omega/clustal-omega-1.2.0
rm -rf prottest-3.4-20140123
rm -rf standard-RAxML-master
rm -rf trimal-trimAl_1.4

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
tar -zxvf FingerPRINTScan_3596.tgz
cd FingerPRINTScan_3596
./configure --prefix=$programs/FingerPRINTScan_3596 && make && make install
cd ..

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

# Unpack ProtTest
tar -zxvf prottest-3.4-20140123.tar.gz

if [ ! -s $programs/prottest-3.4-20140123/prottest-3.4.jar ]
then
echo "ERROR: prottest installation failed"
exit
fi

# Unpack RAxML and compile
unzip standard-RAxML-master.zip
cd standard-RAxML-master
make -f Makefile.SSE3.PTHREADS.gcc && rm *.o
ln -s raxmlHPC-PTHREADS-SSE3 raxml
cd ..

if [ ! -s $programs/standard-RAxML-master/raxml ]
then
echo "ERROR: raxml installation failed"
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
