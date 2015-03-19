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
Klist=(5 10 20 30 40 50 60 70 80 90 100)
#dpsilons=(600 300 100 60 30 10 6 3 1)

topK=5
T=5
D=50



for j in 0 1 2 3 4 5 6 7 8 9 10
do
	K=${Klist[j]} 
        echo "==========$K========="	
	for i in 0 
	do
		data=$dataset$i
		max=${maxs[i]}
		echo "~~~~~", $data, $max, "~~~~~"
		Epsilon=0`echo "scale = 10; sqrt( 1.0/${E[i]})" |bc`
		echo $Epsilon
		displayEpsilon=`echo "$Epsilon*1000000" |bc`
		displayEpsilon=`/usr/bin/printf "%.0f\n" $displayEpsilon`
		echo $displayEpsilon
		./rdsextr2/build/RdSExtr $data $T $K $Epsilon 1
		#cp "result/"$data"_"$T"_"$D"_"$displayEpsilon".pathsim" "result/"$data"_"$T"_"$D"_"$displayEpsilon".topK.pathsim"
		# python calNeighborSimScore.py $data $T $D $Epsilon "pathsim" $topK "Epsilon" $max

		#./ann/bin/ann_sample -d $D -max $max -nn $topK -df "result/"$data"_"$T"_"$D"_"$displayEpsilon".pathvec" -qf "result/"$data"_"$T"_"$D"_"$displayEpsilon".pathvec" -rf "result/"$data"_"$T"_"$D"_"$displayEpsilon".topK.pathvec" 
		# python calStructureSimScore.py $data $T $D $Epsilon "pathvec" $topK "Epsilon" $max
	done
done

