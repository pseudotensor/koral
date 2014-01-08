#!/bin/bash
for i in dumps/0/res*head; do cp $i dumps; done
for i in dumps/0/res*dat; do cp $i dumps; done
typeset -i i END
let END=$1 i=1
while ((i<END)); do
    cd dumps/$i
    pwd
    for j in res*dat; do cat $j >> ../$j; done
    cd ../..
    let i++
done
#for i in $(seq 1 $1);
#do
#echo $i;


