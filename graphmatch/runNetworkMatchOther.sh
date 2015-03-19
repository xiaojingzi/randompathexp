#!/bin/bash -l
#$ -S /bin/bash

# if [ "$#" -lt 2 ]; 
# then
#   echo "Usage: $0 " 
#   exit 1
# fi

######## only test ourmethod, refex, basic features, and random method


datasets=("twitter")
maxnodes=(112411)


topK=100
T=5
D=50
Epsilon=0.005
topDegree=1000

# echo "=======our method======="
# for data in ${datasets[@]}
# do
# 	./main_release $data"1" $T $D $Epsilon
# 	./main_release $data"2" $T $D $Epsilon
# 	max=${maxnodes[i]}
# 	./ann/bin/ann_sample -d $D -max $max -nn $topK -qf "result/"$data"1.pathvec" -df "result/"$data"2.pathvec"  -rf "result/"$data"1_"$data"2.pathvec"
# 	python calrecallOther.py $data "pathvec" $topDegree

# done

echo "=======refex======"
D=18

for data in ${datasets[@]}
do
	python refex.py $data"1"
	python refex.py $data"2"

	max=${maxnodes[i]}
	./ann/bin/ann_sample -d $D -max $max -nn $topK -qf "result/"$data"1.ReFex" -df "result/"$data"2.ReFex"  -rf "result/"$data"1_"$data"2.ReFex"
	python calrecallOther.py  $data "ReFex" $topDegree
done



echo "=======degree, clustering coefficient, closeness, betweenness, pagerank, triad number======="
for data in ${datasets[@]}
do
	python basicfeatures2.py $data"1"
	python basicfeatures2.py $data"2"
done

D=1
for i in 0
do
	data1=${dataList1[$i]}
	data2=${dataList2[$i]}
	max=${maxnodes[i]}
	echo "~~~~~$data1~~~~$data2~~~~"
	echo "degree"
	./ann/bin/ann_sample -d $D -max $max -nn $topK -qf $data1"1.Degree" -df $data"2.Degree"  -rf $data"1_"$data"2.Degree"
	python calrecallOther.py  $data "Degree" $topDegree	
	echo "clustering"
	./ann/bin/ann_sample -d $D -max $max -nn $topK -qf $data"1.ClusCoeff" -df $data"2.ClusCoeff"  -rf $data"1_"$data"2.ClusCoeff"
	python calrecallOther.py  $data "ClusCoeff" $topDegree	
	echo "Closeness"
	./ann/bin/ann_sample -d $D -max $max -nn $topK -qf $data"1.Closeness" -df $data"2.Closeness"  -rf $data"1_"$data"2.Closeness"
	python calrecallOther.py  $data "Closeness" $topDegree	
	echo "betweenness"
	./ann/bin/ann_sample -d $D -max $max -nn $topK -qf $data"1.Betweenness" -df $data"2.Betweenness"  -rf $data"1_"$data"2.Betweenness"
	python calrecallOther.py  $data "Betweenness" $topDegree	
	echo "pagerank"
	./ann/bin/ann_sample -d $D -max $max -nn $topK -qf $data"1.Pagerank" -df $data"2.Pagerank"  -rf $data"1_"$data"2.Pagerank"
	python calrecallOther.py  $data "Pagerank" $topDegree	
	echo "Trianglenumber"
	./ann/bin/ann_sample -d $D -max $max -nn $topK -qf $data"1.Trianglenumber" -df $data"2.Trianglenumber"  -rf $data"1_"$data"2.Trianglenumber"
	python calrecallOther.py  $data "Trianglenumber" $topDegree	
done

echo "=======random======="
for i in 0 
do
	data1=${dataList1[i]}
	data2=${dataList2[i]}
	python calrecallRandom.py $data1 $data2 $topK
done
