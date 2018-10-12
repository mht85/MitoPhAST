path=`pwd`
exe=`echo "$path/exe"`

# clean up previous installations
cd $exe
rm -rf FASconCAT-G_v1.02
rm -rf FingerPRINTScan_3596
rm -rf Gblocks_0.91b
rm -rf hmmer-3.1b2-linux-intel-x86_64
rm -rf iqtree-omp-1.5.5
rm -rf mafft-7.222-with-extensions
rm -rf ncbi-blast-2.6.0+

# install FASconCAT
mkdir FASconCAT-G_v1.02
cd FASconCAT-G_v1.02
ln -s ../FASconCAT-G.zip
unzip FASconCAT-G.zip
cd ..

if [ ! -s $exe/FASconCAT-G_v1.02/FASconCAT-G_v1.02.pl ]
then
echo "ERROR: FASconCAT-G installation failed"
exit
fi

# install FingerPRINTScan
tar -zxvf FingerPRINTScan_3596.tgz
cd FingerPRINTScan_3596
./configure --prefix=$exe/FingerPRINTScan_3596/ && make && make install
cd ..

if [ ! -s $exe/FingerPRINTScan_3596/bin/fingerPRINTScan ]
then
echo "ERROR: FingerPRINTScan installation failed"
exit
fi

# install Gblocks
tar -xvf Gblocks_Linux64_0.91b.tar

if [ ! -s $exe/Gblocks_0.91b/Gblocks ]
then
echo "ERROR: Gblocks installation failed"
exit
fi

# install hmmer
tar -zxvf hmmer-3.1b2-linux-intel-x86_64.tar.gz
cd hmmer-3.1b2-linux-intel-x86_64
./configure --prefix=$exe/hmmer-3.1b2-linux-intel-x86_64/ && make && make install
cd ..

if [ ! -s $exe/hmmer-3.1b2-linux-intel-x86_64/bin/hmmsearch ]
then
echo "ERROR: HMMer installation failed"
exit
fi

# install IQ-TREE
tar -zxvf iqtree-omp-1.5.5-Linux.tar.gz
mv iqtree-omp-1.5.5-Linux iqtree-omp-1.5.5

if [ ! -s $exe/iqtree-omp-1.5.5/bin/iqtree-omp ]
then
echo "ERROR: IQ-TREE installation failed"
exit
fi

# install MAFFT
tar -zxvf mafft-7.394-with-extensions-src.tgz
cd mafft-7.394-with-extensions/core
cp Makefile Makefile.backup
sed -i "s|/usr/local|$exe/mafft-7.394-with-extensions|g" Makefile
make && make install
cd ../../

if [ ! -s $exe/mafft-7.394-with-extensions/bin/mafft ]
then
echo "ERROR: MAFFT installation failed"
exit
fi

# install NCBI BLAST
tar -zxvf ncbi-blast-2.6.0+.tar.gz

if [ ! -s $exe/ncbi-blast-2.6.0+/bin/blastn ]
then
echo "ERROR: BLAST installation failed"
exit
fi

echo "Installation completed successfully."
