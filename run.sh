#!/bin/bash

I=1
# 1: semi-Dirac; 2: massive; 3: massless
./a.out $I 1 -1 0 0.01 1 0.001 10 10

if [ -e "D$I.mat" ]
then
    rm D$I.mat
fi
txt2mat D$I.txt D$I D$I.mat
rm D$I.txt

F=$(echo "10/10000*435" | bc -l)
echo $F
