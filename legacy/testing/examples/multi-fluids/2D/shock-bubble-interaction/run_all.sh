#!/bin/bash -
#===============================================================================
#
#          FILE: run_all.sh
#
#         USAGE: ./run_all.sh [options]
#
#   DESCRIPTION: bash srcipt for running the present test
#===============================================================================
function print_usage {
  echo "Bash srcipt for running IBM code"
  echo "Usage"
  echo "   1) ./run_IBM.sh        => print this help message;"
  echo "   2) ./run_IBM.sh -run   => run IBM, OFF and POG;"
  echo "   3) ./run_IBM.sh -clean => clean the test."
  exit 1
}
if [ $# -eq 0 ] ; then
  print_usage
fi
if [ "$1" == "-clean" ] ; then
  cd ibm
  ./run_IBM.sh -clean > /dev/null 2>&1
  cd ..
  cd off
  ./run_OFF.sh -clean > /dev/null 2>&1
  cd ..
  cd pog
  ./run_POG.sh -clean > /dev/null 2>&1
  cd ..
elif [ "$1" == "-run" ] ; then
  cd ibm
  ./run_IBM.sh -out > /dev/null 2>&1
  cd ..
  cd off
  ./run_OFF.sh -no_mpi > /dev/null 2>&1
  cd ..
  cd pog
  ./run_POG.sh -tec > /dev/null 2>&1
  lay_export rupg.lay
  cd ..
else
  echo "Unknown switch $1"
  print_usage
fi
exit 0
