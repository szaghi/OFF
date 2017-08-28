#!/bin/bash -
#===============================================================================
#
#          FILE: run_IBM.sh
#
#         USAGE: ./run_IBM.sh [options]
#
#   DESCRIPTION: bash srcipt for running IBM code
#===============================================================================
function print_usage {
  echo "Bash srcipt for running IBM code"
  echo "Usage"
  echo "   1) ./run_IBM.sh        => print this help message;"
  echo "   2) ./run_IBM.sh -out   => produce output for simulation;"
  echo "   3) ./run_IBM.sh -clean => clean current directory."
  exit 1
}
if [ $# -eq 0 ] ; then
  print_usage
fi
if [ "$1" == "-clean" ] ; then
  rm -rf output # cleaning working directory
elif [ "$1" == "-out" ] ; then
  if [ ! -f IBM ]; then
    echo "IBM code has not been generated!"
    echo "Go into the root directory and compile IBM"
    exit 1
  fi
  rm -rf output                            # cleaning working directory
  mkdir -p output                          # creating outout directory
  ./IBM ibm_options.dat                    # running IBM with options defined in "ibm_options.dat"
  cp ./input/solver.dat ./output/          # copying solver options
  cp ./input/procmap-no_mpi.dat ./output/  # copyng procmap without MPI
  cp ./input/procmap-mpi.dat ./output/     # copyng procmap with    MPI
else
  echo "Unknown switch $1"
  print_usage
fi
exit 0
