#!/bin/bash -
#===============================================================================
#
#          FILE: run_POG.sh
#
#         USAGE: ./run_POG.sh [options]
#
#   DESCRIPTION: bash srcipt for running POG code
#===============================================================================
function print_usage {
  echo "Bash srcipt for running POG code"
  echo "Usage"
  echo "   1) ./run_POG.sh      => print this help message;"
  echo "   2) ./run_POG.sh -tec => produce Tecplot output;"
  echo "   3) ./run_POG.sh -vtk => produce VTK output."
  exit 1
}
if [ $# -eq 0 ] ; then
  print_usage
fi
rm -rf output input                                                                                         # cleaning working directory
mkdir output                                                                                                # creating output directiry
ln -fs ../off/output input                                                                                  # link to OFF output
# last simulation iteration
for file in $( ls ./input/*b001*N_*.sol ); do
  nlast=`basename $file | awk -F . '{print $3}' | awk -F _ '{print $2}' | sed 's/^0*//'`
done
n=`printf "%10.10d" $nlast`
tout=$1" yes"                                                                                        # type of output
./POG -m ./input/sod.g01.b001.geo -s ./input/sod.g01.b001-N_$n.sol -o output/sod.g01.b001-N_$n $tout #-cell # post-processing block 1
./POG -m ./input/sod.g01.b002.geo -s ./input/sod.g01.b002-N_$n.sol -o output/sod.g01.b002-N_$n $tout #-cell # post-processing block 2
exit 0
