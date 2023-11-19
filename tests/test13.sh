#!/bin/bash
# 145 runs
for ((n=2; n<=12;n++)) ; do for ((m=3;m<=$n;m+=3)) ; do for ((s=1;s<=4;s++)) ; do echo "n=$n m=$m r=3 s=$s" ; ./a.out $n $m 3 $s ; done ; done ; done
