#!/bin/bash

for ((i=8; i<=$1; i++))
do
   args="${i} 9 123 $(( i > 16 ? 9 : i-7 )) 2 99"
   echo "args: t flags seed N cutoff target == " $args >> subsetsum-mem${i}.txt
   cmd="../SCIPR-build/scipr_pcp/PCP --gtest --gtest_filter=cs2Brex.MEMcoNPsubsetsum --extra-args "
   echo $cmd "${args}"
   $cmd "${args}" >> subsetsum-mem${i}.txt
done
