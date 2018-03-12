#!/bin/bash
for d in $(find ./src/third_party/*/ -maxdepth 0 -type d)
do
   cd $d
   echo $PWD
   git checkout master
   git pull
   cd -
done
