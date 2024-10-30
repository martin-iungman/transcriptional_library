#!/bin/bash
for Gate in 2 3 4 5 6 7
do
	for Lib in 1 3
	do
		echo ${Gate}-${Lib} >> tmp.txt
	done
done
