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
maxs=(6523 25844 48837 169209 230103 443070 702049 2767344)
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
for j in 0 1 2 3 4 5 6 7 8
do
	Epsilon=${Epsilons[j]}
	displayEpsilon=${dpsilons[j]}
	echo "=======$Epsilon=======$displayEpsilon===="
	# python insertText.py $Epsilon "neighborsim" "Epsilon"
	# python insertText.py $Epsilon "structuresim" "Epsilon"
	
	for i in 0 1 2 3 4 5 6 7 9 10 11
	do
		data=$dataset$i
		max=${maxs[i]}
		echo "~~~~~", $data, $max, "~~~~~"
		./main_release $data $T $D $Epsilon
		cp "result/"$data"_"$T"_"$D"_"$displayEpsilon".pathsim" "result/"$data"_"$T"_"$D"_"$displayEpsilon".topK.pathsim"
		# python calNeighborSimScore.py $data $T $D $Epsilon "pathsim" $topK "Epsilon" $max

		#./ann/bin/ann_sample -d $D -max $max -nn $topK -df "result/"$data"_"$T"_"$D"_"$displayEpsilon".pathvec" -qf "result/"$data"_"$T"_"$D"_"$displayEpsilon".pathvec" -rf "result/"$data"_"$T"_"$D"_"$displayEpsilon".topK.pathvec" 
		# python calStructureSimScore.py $data $T $D $Epsilon "pathvec" $topK "Epsilon" $max
	done
done

