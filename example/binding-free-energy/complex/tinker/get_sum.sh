#!/bin/bash

grep "Free Energy via BAR Iteration" fep-solv-*.out | awk '{print $1,$7,$9; s+=$7;v+=$9*$9;} END {print "Total",s,sqrt(v)}' > table_bariter.txt
rsync -a table_bariter.txt ../
