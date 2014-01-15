#!/bin/bash
#input: ./merge_res.sh #cores #start #end #iter
#merges the res... files from per core folders and puts them in dumps/
typeset -i i j END

let END=$3 i=$2
while ((i<=END)); do
    cp dumps/0/res`printf %04d ${i}`.head dumps
    cp dumps/0/res`printf %04d ${i}`.dat dumps
    let i=i+$4
done
let ENDF=$1 i=1
while ((i<ENDF)); do
    echo $i
    cd dumps/$i
    let END=$3 j=$2
    while ((j<=END)); do
    cat res`printf %04d ${j}`.dat >> ../res`printf %04d ${j}`.dat
    let j=j+$4
    done
    cd ../..
    let i++
done




