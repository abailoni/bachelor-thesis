#!/bin/bash


i=1
input[1]=$2
input[2]=$3

if [ input[0] ]
then
	echo "Ciao!"
	let input[0]=4
fi
#while [ "${input[i-1]}" != "." ] 
# do
# 	echo "Specificare file input dati $i:"
# 	read input[i]
# 	let i=i+1
# done
# echo "Ciao ;)"
