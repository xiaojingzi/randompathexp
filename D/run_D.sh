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
displayEpsilons=(10 10 100 100 100 100 100 100)
topK=50
T=5

for i in 0 1 2 3 4 5 6 7
do
	data=${datasets[i]}
	python initRun.py $data "structuresim" "D"

done

for i in 0 1 2 3 4 5 6 7
do
	data=${datasets[i]}
	max=${maxs[i]}
	Epsilon=${Epsilons[i]}
	displayEpsilon=${displayEpsilons[i]}
	echo $data, $max, $Epsilon, $displayEpsilon
		
	for D in 2 5 10 20 50 100
	do
		echo "=======$D======="
		./main_release $data $T $D $Epsilon
		./ann/bin/ann_sample -d $D -max $max -nn $topK -df "result/"$data"_"$T"_"$D"_"$displayEpsilon".pathvec" -qf "result/"$data"_"$T"_"$D"_"$displayEpsilon".pathvec" -rf "result/"$data"_"$T"_"$D"_"$displayEpsilon".topK.pathvec" 
		python insertText.py $data $D "structuresim" "D"
		python calStructureSimScore.py $data $T $D $Epsilon "pathvec" $topK "D"
	done
done

