#!/bin/bash -
#===============================================================================
#
#          FILE: run_OFF.sh
#
#         USAGE: ./run_OFF.sh [options]
#
#   DESCRIPTION: bash srcipt for running OFF code
#===============================================================================
function print_usage {
  echo "Bash srcipt for running OFF code"
  echo "Usage"
  echo "   1) ./run_OFF.sh         => print this help message;"
  echo "   2) ./run_OFF.sh -no_mpi => run simulation without MPI (link to procmap-no_mpi.dat);"
  echo "   3) ./run_OFF.sh -mpi    => run simulation with    MPI (link to procmap-mpi.dat);"
  echo "   4) ./run_OFF.sh -clean  => clean current directory."
  exit 1
}
if [ $# -eq 0 ] ; then
  print_usage
fi
if [ "$1" == "-clean" ] ; then
  rm -rf output lockfile input # cleaning working directory
else
  if [ ! -d ../ibm/output ]; then
    echo "IBM output dir not found!"
    echo "Before run OFF you must run IBM:"
    echo "cd ../ibm/ ; run_IBM.sh -out ; cd -"
    exit 1
  fi
  if [ ! -f OFF ]; then
    echo "OFF code has not been generated!"
    echo "Go into the root directory and compile OFF"
    exit 1
  fi
  rm -rf output lockfile input # cleaning working directory
  ln -fs ../ibm/output input   # link to IBM output
  mkdir -p output              # creating outout directory
  if [ "$1" = "-no_mpi"  ] ; then
    cd input ; ln -fs procmap-no_mpi.opt procmap.opt ; cd ../ # linking procmap.dat to map without MPI
    ./OFF off_options.dat                                     # running OFF without MPI
  elif [ "$1" = "-mpi"  ] ; then
    cd input ; ln -fs procmap-mpi.opt procmap.opt ; cd ../ # linking procmap.dat to map with MPI
    mpirun -n 2 ./OFF off_options.dat                      # running OFF with MPI
  else
    echo "Unknown switch $1"
    print_usage
  fi
fi
exit 0
