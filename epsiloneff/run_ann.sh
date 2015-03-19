#!/bin/bash -l
#$ -S /bin/bash

# if [ "$#" -lt 2 ]; 
# then
#   echo "Usage: $0 " 
#   exit 1
# fi


# datasets=("mobile2" "twitter2")
# maxs=(194526 112411)
dataset="tencent"
maxs=(6523 25844 48837 169209 230103 443070 702049 2767344 5355507 26033969 51640620 103675609 355591065)
E=(10000 50000 100000 500000 1000000 5000000 10000000 50000000 100000000 500000000 1000000000 2000000000 5958853072) 
#Epsilons=(0.06 0.03 0.01 0.006 0.003 0.001 0.0006 0.0003 0.0001 0.00003)
#dpsilons=(600 300 100 60 30 10 6 3 1)

topK=5
T=5
D=50



# for i in 0 1 2 3 4 5 6
# do
# 	data=$dataset$i
# 	echo $data
# 	python findCommonNeighbors.py $data
# 	python caldegreestructurehole.py $data
# done
# python initRun.py "neighborsim" "Epsilon"
# python initRun.py "structuresim" "Epsilon"
	
	for i in 9 10 11 12 
	do
		data=$dataset$i
		max=1000000
		echo "~~~~~", $data, $max, "~~~~~"
		Epsilon=0`echo "scale = 10; sqrt( 1.0/${E[i]})" |bc`
		echo $Epsilon
		displayEpsilon=`echo "$Epsilon*1000000" |bc`
		echo $displayEpsilon
		displayEpsilon=`/usr/bin/printf "%.0f\n" $displayEpsilon`
		echo $displayEpsilon
		rdsextr/res/tencent/ann/bin/ann_sample -d $D -max $max -nn $topK -df "result/"$data"_"$T"_"$D"_"$displayEpsilon".pathvec" -qf "result/"$data"_"$T"_"$D"_"$displayEpsilon".pathvec" -rf "result/"$data"_"$T"_"$D"_"$displayEpsilon".topK.pathvec" 
	done

