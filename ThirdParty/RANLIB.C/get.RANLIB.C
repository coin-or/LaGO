#!/bin/sh

wgetcount=`which wget 2>/dev/null | wc -w`
if test ! $wgetcount = 1; then
  echo "Utility wget not found in your PATH."
  exit -1
fi

echo " "
echo "Running script for downloading the source code for RANLIB.C"
echo " "

rm -f ranlib.c.tar.gz
echo "Downloading the source code from www.netlib.org..."
wget http://www.netlib.org/random/ranlib.c.tar.gz

rm -rf ranlib.c
echo "Unpacking the source code..."
tar xzf ranlib.c.tar.gz

rm -f ranlib.c.tar.gz

echo " "
echo "Done downloading the RANLIB.C source."
echo " "
