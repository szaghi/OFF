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
  echo "   3) ./run_OFF.sh -mpi    => run simulation with    MPI (link to procmap-mpi.dat)."
  exit 1
}
if [ $# -eq 0 ] ; then
  print_usage
fi
rm -rf output lockfile input # cleaning working directory
ln -fs ../ibm/output input # link to IBM output
if [ "$1" = "-no_mpi"  ] ; then
  cd input
  ln -fs procmap-no_mpi.dat procmap.dat
  cd ../
  ./OFF off_options.dat # running OFF without MPI
elif [ "$1" = "-mpi"  ] ; then
  cd input
  ln -fs procmap-mpi.dat procmap.dat
  cd ../
  mpirun -n 2 ./OFF off_options.dat # running OFF with MPI
fi
exit 0
