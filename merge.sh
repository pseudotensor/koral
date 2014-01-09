#!/bin/bash
#merges all the res... and avg... files from per core folders and puts them in dumps/
for i in dumps/0/res*head; do cp $i dumps; done
for i in dumps/0/res*dat; do cp $i dumps; done
typeset -i i END
let END=$1 i=1
while ((i<END)); do
    cd dumps/$i
    for j in res*dat; do cat $j >> ../$j; done
    cd ../..
    let i++
done
for k in dumps/0/avg*head; do cp $k dumps; done
for k in dumps/0/avg*dat; do cp $k dumps; done
typeset -i i END
let END=$1 i=1
while ((i<END)); do
    cd "dumps/$i"
    pwd
    for j in avg*dat; do cat $j >> ../$j; done
    cd ../..
    let i++
done





