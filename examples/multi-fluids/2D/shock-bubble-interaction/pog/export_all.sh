#!/bin/bash -
mkdir -p pngs
for file in $( ls ./output/*-*.plt ); do
  flast=`basename $file`
  froot=`echo  $flast | awk -F - '{print $1}'`
  nf=`basename $flast .plt | awk -F - '{print $2}'`
  cd output
  ln -fs $flast $froot".plt"
  cd ../
  ./lay_export sch.lay
  mv sch.png pngs/$froot-$nf-sch.png
done
