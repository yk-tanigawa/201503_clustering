#!/bin/sh

make

seed=0
out_file="../results_Hamerly.txt"

for n in 1000 10000 100000
do
    for d in 2 10 100
    do
	for k in 10 100 1000
	do
	    file="../data/n${n}_d${d}.dat"
	    ./Hamerly ${n} ${d} ${k} ${file} >> ${out_file}
	done
    done
done

make clean
