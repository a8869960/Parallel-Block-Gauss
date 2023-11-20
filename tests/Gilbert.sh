#!/bin/bash
# 145 runs
for ((n=12; n<=12;n++)) ; do for ((m=3;m<=$n;m+=1)) ; do for((p=1; p<=$m+1; p++)) ; do echo $'\n' "n=$n m=$m p=$p r=3 s=4 --------------" ; ./a.out $n $m $p 3 4 ; done ; done ; done
