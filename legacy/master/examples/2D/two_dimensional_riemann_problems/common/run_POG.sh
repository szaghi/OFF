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
  echo "   3) ./run_POG.sh -vtk   => produce VTK output;"
  echo "   4) ./run_POG.sh -gnu   => produce gnuplot output;"
  echo "   5) ./run_POG.sh -clean => clean current directory."
  exit 1
}
if [ $# -eq 0 ] ; then
  print_usage
fi
if [ "$1" == "-clean" ] ; then
  rm -rf output input        # cleaning working directory
  ln -fs ../off/output input # link to OFF output
else
  if [ ! -d ../off/output ]; then
    echo "OFF output dir not found!"
    echo "Before run POG you must run OFF:"
    echo "cd ../off/ ; run_OFF.sh -mpi[or -no_mpi] ; cd -"
    exit 1
  fi
  if [ ! -f POG ]; then
    echo "POG code has not been generated!"
    echo "Go into the root directory and compile POG"
    exit 1
  fi
  rm -rf output input        # cleaning working directory
  ln -fs ../off/output input # link to OFF output
  mkdir -p output            # creating outout directory
  ulimit -s unlimited        # increasing the stack size
  for file in $( ls ./input/*b001*N_*.sol ); do
    nlast=`basename $file | awk -F . '{print $3}' | awk -F _ '{print $2}' | sed 's/^0*//'`
  done
  n=`printf "%10.10d" $nlast`
  if [ "$1" = "-tec" ]; then
    tout="-tec yes -vtk no -gnu no"
  elif [ "$1" = "-vtk" ]; then
    tout="-tec no -vtk yes -gnu no"
  elif [ "$1" = "-gnu" ]; then
    tout="-tec no -vtk no -gnu yes"
  else
    echo "Unknown switch $1"
    print_usage
  fi
  ./POG -m ./input/2DRP.g01.b001.geo -s ./input/2DRP.g01.b001-N_$n.sol -o output/2DRP.g01.b001 $tout # post-processing block 1
  ./POG -m ./input/2DRP.g01.b002.geo -s ./input/2DRP.g01.b002-N_$n.sol -o output/2DRP.g01.b002 $tout # post-processing block 2
  ./POG -m ./input/2DRP.g01.b003.geo -s ./input/2DRP.g01.b003-N_$n.sol -o output/2DRP.g01.b003 $tout # post-processing block 2
  ./POG -m ./input/2DRP.g01.b004.geo -s ./input/2DRP.g01.b004-N_$n.sol -o output/2DRP.g01.b004 $tout # post-processing block 4
fi
exit 0
