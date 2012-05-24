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
  echo "   1) ./run_IBM.sh         => print this help message;"
  echo "   2) ./run_IBM.sh -no_mpi => produce output for simulation without MPI (link to procmap-no_mpi.dat);"
  echo "   3) ./run_IBM.sh -mpi    => produce output for simulation with    MPI (link to procmap-mpi.dat)."
  exit 1
}
if [ $# -eq 0 ] ; then
  print_usage
fi
rm -rf output         # cleaning working directory
./IBM ibm_options.dat # running IBM with options defined in "ibm_options.dat"
# making simbolic links
cd output
ln -fs ../solver.dat .
if [ "$1" = "-no_mpi"  ] ; then
  ln -fs ../procmap-no_mpi.dat procmap.dat
elif [ "$1" = "-mpi"  ] ; then
  ln -fs ../procmap-mpi.dat procmap.dat
fi
cd ../
exit 0
