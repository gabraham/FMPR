#!/bin/bash

N=$(cat X.txt | wc -l)
p=$(head -1 X.txt | wc -w)
K=$(head -1 Y.txt | wc -w)
~/Downloads/gflasso/gflasso Y.txt X.txt \
   $N $p $K config.txt gc Ynetwork.txt \
   gflasso
