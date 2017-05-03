#!/bin/bash

function calcn {
	n=1
	local m=1
	while [ $((2**m * (m * 9 + 7))) -lt $((2**i)) ]
	do
	n=$m
	m=$[$m+1]
	done
}

for ((i=$1; i<=$2; i++))
do
   calcn
   args="${i} 9 123 ${n} 99"
   echo "args: t flags seed N target == " $args >> subsetsum${i}.txt
   cmd="../SCIPR-build/scipr_pcp/PCP --gtest --gtest_filter=cs2Brex.coNPsubsetsum --extra-args "
   echo $cmd "${args}"
   $cmd "${args}" >> subsetsum${i}.txt
done
