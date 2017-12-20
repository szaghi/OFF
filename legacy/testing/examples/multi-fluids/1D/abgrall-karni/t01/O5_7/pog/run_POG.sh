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
  echo "   1) ./run_POG.sh        => print this help message;"
  echo "   2) ./run_POG.sh -tec   => produce Tecplot output;"
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
  for file in $( ls ./input/*N_* ); do
    flast=`basename $file`
  done
  if [ "$1" = "-tec" ]; then
    tout="-tec yes -vtk no -gnu no"
  elif [ "$1" = "-vtk" ]; then
    tout="-tec no -vtk yes -gnu no"
  elif [ "$1" = "-gnu" ]; then
    tout="-tec no -vtk no -gnu yes -ascii"
  else
    echo "Unknown switch $1"
    print_usage
  fi
  froot=`ls ./input/*.geo | xargs basename -s .geo`
  ./POG -m ./input/$froot".geo" -s ./input/$flast -o output/$froot $tout
  if [ "$1" = "-gnu" ]; then
    gnuplot r.gnu
  fi
fi
exit 0
