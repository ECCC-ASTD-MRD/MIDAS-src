#!/bin/bash
#

cd aai
./compile_aai.sh
cd ../

echo ''
echo ' ----------------------------------------------------------- '
read -rsp $' Press enter to compile next program...\n ----------------------------------------------------------- \n \n'

cd ominusf
./compile_ominusf.sh
cd ../

echo ''
echo ' ----------------------------------------------------------- '
read -rsp $' Press enter to compile next program...\n ----------------------------------------------------------- \n \n'

cd ensmanip
./compile_ensmanip.sh
cd ../

echo ''
echo ' ----------------------------------------------------------- '
read -rsp $' Press enter to compile next program...\n ----------------------------------------------------------- \n \n'

cd diagbmatrix
./compile_diagbmatrix.sh
cd ../

echo ''
echo ' ----------------------------------------------------------- '
read -rsp $' Press enter to compile next program...\n ----------------------------------------------------------- \n \n'

cd randompert
./compile_randompert.sh
cd ../

echo ''
echo ' ----------------------------------------------------------- '
read -rsp $' Press enter to compile next program...\n ----------------------------------------------------------- \n \n'

cd calcstats
./compile_calcstats.sh
cd ../

echo ''
echo ' ----------------------------------------------------------- '
read -rsp $' Press enter to compile next program...\n ----------------------------------------------------------- \n \n'

./compile_oavar.sh

