#! /bin/bash

ipara=$(awk '/iparallel/{print $1}' AFSSH.inp)

rm hops.out
rm fort.17
rm fort.18
rm fort.40
for((i=1;i<=$ipara;i++));do
   cat $i/hops.out >> hops.out
   cat $i/fort.17 >> fort.17
   cat $i/fort.18 >> fort.18
   cat $i/fort.40 >> fort.40
done

