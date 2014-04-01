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
  echo "   1) ./run_POG.sh             => print this help message;"
  echo "   2) ./run_POG.sh -tec [opts] => produce Tecplot output;"
  echo "   3) ./run_POG.sh -vtk [opts] => produce VTK output;"
  echo "   4) ./run_POG.sh -gnu [opts] => produce gnuplot output;"
  echo "   5) ./run_POG.sh -clean      => clean current directory."
  echo " Options are:"
  echo "   -all               => post process all files;"
  echo "   -n #n1 #n2 [#nout] => post process files from #n1 to #n2 (optionally with frequency #nout)."
  echo " Examples:"
  echo "   1) ./run_POG.sh -tec           => post process the last solution file;"
  echo "   2) ./run_POG.sh -tec -all      => post process all solution files;"
  echo "   3) ./run_POG.sh -tec -n 10 320 => post process solution files from 10 to 320."
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
  rm -rf input               # cleaning working directory
  ln -fs ../off/output input # link to OFF output
  mkdir -p output            # creating outout directory
  # setting POG options
  if [ "$1" = "-tec" ]; then
    tout="-tec yes -vtk no -gnu no -mirrorY"
  elif [ "$1" = "-vtk" ]; then
    tout="-tec no -vtk yes -gnu no"
  elif [ "$1" = "-gnu" ]; then
    tout="-tec no -vtk no -gnu yes -ascii"
  else
    echo "Unknown switch $1"
    print_usage
  fi
  froot=`ls ./input/*.geo | xargs basename -s .geo`
  if [ $# -ge 2 ] ; then
    if [ "$2" == "-all" ] ; then
      for file in $( ls ./input/*N_* ); do
        flast=`basename $file`
        nf=`echo  $flast | awk -F _ '{print $2}'`
        ./POG -m ./input/$froot".geo" -s ./input/$flast -o output/$froot-$nf $tout
        cd output
        ln -fs $froot-$nf* $froot".plt"
        cd ../
      done
    elif [ "$2" == "-n" ] ; then
      if [ $# -eq 4 ] ; then
        n1=$3
        n2=$4
        nout=1
      elif [ $# -eq 5 ] ; then
        n1=$3
        n2=$4
        nout=$5
      else
        echo "With -n switch the arguments must be #n1 #n2 [#nout]"
        print_usage
      fi
      for file in $( ls ./input/*N_* ); do
        flast=`basename $file`
        nf=`echo  $flast | awk -F _ '{print $2}'`
        n=`echo  $flast | awk -F _ '{print $2}' | bc -l`
        mod=$(($n % $nout))
        if [ $n -ge $n1 ] && [ $n -le $n2 ] && [ $mod -eq 0 ] ; then
          ./POG -m ./input/$froot".geo" -s ./input/$flast -o output/$froot-$nf $tout
          cd output
          ln -fs $froot-$nf* $froot".plt"
          cd ../
        fi
      done
    else
      echo "Unknown switch $2"
      print_usage
    fi
  else
    for file in $( ls ./input/*N_* ); do
      flast=`basename $file`
      nf=`echo  $flast | awk -F _ '{print $2}'`
    done
    ./POG -m ./input/$froot".geo" -s ./input/$flast -o output/$froot-$nf $tout
    cd output
    ln -fs $froot-$nf* $froot".plt"
    cd ../
  fi
fi
exit 0
