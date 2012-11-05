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
  echo "   2) ./run_IBM.sh -out => produce output for simulation."
  exit 1
}
if [ $# -eq 0 ] ; then
  print_usage
fi
rm -rf output         # cleaning working directory
mkdir output          # creating outout directory
./IBM ibm_options.dat # running IBM with options defined in "ibm_options.dat"
# making simbolic links
cd output
ln -fs ../solver.dat .
ln -fs ../procmap-no_mpi.dat .
ln -fs ../procmap-mpi.dat .
cd ../
exit 0
