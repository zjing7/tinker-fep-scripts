#!/bin/bash

barname=solvwaterboxproddyn.bar
barout=bar2.out
flist=`cat fep-list.solv.txt | sed "s|$|/$barname|g"`
python ../../../scripts/t_bar.py $flist -s -o solv.t_bar.txt > solv.t_bar.log

flist=`cat fep-list.gas.txt | sed "s|$|/$barname|g"`
python ../../../scripts/t_bar.py $flist -s -o gas.t_bar.txt > gas.t_bar.log

flist=`cat fep-list.solv.txt | sed "s|$|/$barout|g"`
grep "Free Energy via BAR Iteration" $flist | awk '{print $1,$7,$9; s+=$7;v+=$9*$9;} END {print "Total",s,sqrt(v)}' > solv.bariter.txt
flist=`cat fep-list.gas.txt | sed "s|$|/$barout|g"`
grep "Free Energy via BAR Iteration" $flist | awk '{print $1,$7,$9; s+=$7;v+=$9*$9;} END {print "Total",s,sqrt(v)}' > gas.bariter.txt
