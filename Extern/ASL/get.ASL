#!/bin/sh

echo " "
echo "Running script for downloading the source code for the ASL"
echo " "

rm -f solvers.tar
echo "Downloading the source code from www.netlib.org..."
wget ftp://www.netlib.org/ampl/solvers.tar

echo "Unpacking the source code..."
tar xf solvers.tar

echo "Uncompressing the source files..."
gunzip -fr solvers

echo "Deleting the tar file..."
rm solvers.tar

echo " "
echo "Done downloading the source code for ASL."
echo " "

