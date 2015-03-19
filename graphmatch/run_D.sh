#!/bin/bash -l
#$ -S /bin/bash

# if [ "$#" -lt 2 ]; 
# then
#   echo "Usage: $0 " 
#   exit 1
# fi


datasets=("kdd" "sigir" "sigmod" "icdm" "cikm" "icde")
dataList1=("kdd"  "sigir" "sigmod")
dataList2=("icdm" "cikm"  "icde" )

T=10
Epsilon=0.01
topK=100

for i in 0 1 2 
do
	data1=${dataList1[i]}
	data2=${dataList2[i]}
	python initRun.py $data1"_"$data2 
done


for D in 2 5 20 50 100 200
do
	echo "=======$D======="

	for data in ${datasets[@]}
	do
		./main_release $data $T $D $Epsilon
		# python randomproject.py $data $T $D $Epsilon
	done
	for i in 0 1 2 
	do
		
		data1=${dataList1[i]}
		data2=${dataList2[i]}
		python insertText.py $data1"_"$data2 $D
		python calrecallCoauthorD.py  $data1 $data2 "pathvec" $topK	
	done
done

