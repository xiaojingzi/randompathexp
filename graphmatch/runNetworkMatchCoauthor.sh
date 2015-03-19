#!/bin/bash -l
#$ -S /bin/bash

# if [ "$#" -lt 2 ]; 
# then
#   echo "Usage: $0 " 
#   exit 1
# fi
######## for network matching on coauthor data set ##########
######## only test ourmethod, refex, basic features, and random method ======

datasets=("kdd" "sigir" "sigmod" "icdm" "cikm" "icde")
dataList1=("kdd" "sigir" "sigmod")
dataList2=("icdm" "cikm" "icde")

topK=100
T=5
D=50
Epsilon=0.01

echo "=======our method======="

for data in ${datasets[@]}
do
	# ./main_release $data $T $D $Epsilon
	python randomproject.py $data $T $D $Epsilon
done

for i in 0 1 2
do
	data1=${dataList1[i]}
	data2=${dataList2[i]}
	python calrecallCoauthor.py  $data1 $data2 "pathvec" $topK
done


echo "=======refex======"
for data in ${datasets[@]}
do
	python refex.py $data
done

for i in 0 1 2
do
	data1=${dataList1[i]}
	data2=${dataList2[i]}
	python calrecallCoauthor.py  $data1 $data2 "ReFex" $topK
done



echo "=======degree, clustering coefficient, closeness, betweenness, pagerank, triad number======="
for data in ${datasets[@]}
do
	python basicfeatures2.py $data
done
for i in 0 1 2
do
	data1=${dataList1[$i]}
	data2=${dataList2[$i]}
	echo "~~~~~$data1~~~~$data2~~~~"
	echo "calcualte recall of "$data1"_"$data2".Degree"
	python calrecallCoauthor.py  $data1 $data2 "Degree" $topK
	echo "calcualte recall of "$data1"_"$data2".ClusCoeff"
	python calrecallCoauthor.py  $data1 $data2 "ClusCoeff" $topK
	echo "calcualte recall of "$data1"_"$data2".Closeness"
	python calrecallCoauthor.py  $data1 $data2 "Closeness" $topK
	echo "calcualte recall of "$data1"_"$data2".Betweenness"
	python calrecallCoauthor.py  $data1 $data2 "Betweenness" $topK
	echo "calcualte recall of "$data1"_"$data2".Pagerank"
	python calrecallCoauthor.py  $data1 $data2 "Pagerank" $topK
	echo "calcualte recall of "$data1"_"$data2".Trianglenumber"
	python calrecallCoauthor.py  $data1 $data2 "Trianglenumber" $topK
done

echo "=======random======="
for i in 0 1 2
do
	data1=${dataList1[i]}
	data2=${dataList2[i]}
	python calrecallRandom.py $data1 $data2 $topK
done

