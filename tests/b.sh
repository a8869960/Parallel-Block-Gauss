#!/bin/bash
# 4 runs
for ((m=1; m<=4; m++)) ; do for((p=1; p<=6; p++)) ; do echo $'\n' "n=4 m=$m p=$p r=4 s=0 b.txt ------------------" ; ./a.out 4 $m $p 4 0 matr/b.txt ; done ; done
