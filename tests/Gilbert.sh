#!/bin/bash
# 145 runs
for ((n=2; n<=12;n++)) ; do for ((m=3;m<=$n;m+=1)) ; do echo $'\n' "n=$n m=$m r=3 s=4 --------------" ; ./a.out $n $m 3 4 ; done ; done
