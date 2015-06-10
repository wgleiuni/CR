#!/bin/bash

I=1
# 1: semi-Dirac; 2: massive; 3: massless
./a.out $I 2 -1 0 1.7 1 0.001 100 100000

if [ -e "D$I.mat" ]
then
    rm D$I.mat
fi

if [ -e "L$I.mat" ]
then
    rm L$I.mat
fi

txt2mat D$I.txt D$I D$I.mat
txt2mat L$I.txt L$I L$I.mat
rm D$I.txt L$I.txt
