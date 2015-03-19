#!/bin/bash -l
#$ -S /bin/bash

# if [ "$#" -lt 2 ]; 
# then
#   echo "Usage: $0 " 
#   exit 1
# fi


datasets=("mobile2" "twitter2" "kdd" "sigir" "sigmod" "icdm" "cikm" "icde")
maxs=(194526 112411 2867 2607 2851 3548 2616 2559)
Epsilons=(0.001 0.001 0.01 0.01 0.01 0.01 0.01 0.01)
dpsilons=(10 10 100 100 100 100 100 100)
D=50
topK=50
for i in  1 #1 2 3 4 5 6 7
do
	data=${datasets[i]}
	python initRun.py $data "neighborsim" "T"
	python initRun.py $data "structuresim" "T"

done

for i in 1 #1 2 3 4 5 6 7
do

	data=${datasets[i]}
	max=${maxs[i]}
	Epsilon=${Epsilons[i]}
	displayEpsilon=${dpsilons[i]}
	echo $data, $max, $Epsilon, $displayEpsilon
	for T in  2  5 10 20 50 100
	do
		echo "=======$T======="
		./main_release $data $T $D $Epsilon
		cp "result/"$data"_"$T"_"$D"_"$displayEpsilon".pathsim" "result/"$data"_"$T"_"$D"_"$displayEpsilon".topK.pathsim"
		python insertText.py $data $T "neighborsim" "T"
		python calNeighborSimScore.py $data $T $D $Epsilon "pathsim" $topK "T"

		./ann/bin/ann_sample -d $D -max $max -nn $topK -df "result/"$data"_"$T"_"$D"_"$displayEpsilon".pathvec" -qf "result/"$data"_"$T"_"$D"_"$displayEpsilon".pathvec" -rf "result/"$data"_"$T"_"$D"_"$displayEpsilon".topK.pathvec" 
		python insertText.py $data $T "structuresim" "T"
		python calStructureSimScore.py $data $T $D $Epsilon "pathvec" $topK "T"
	done
done

