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
maxs=(6523 25844 48837 169209 230103 443070 702049 2767344 5355507)
topK=5
T=5
D=18



for i in 0 1 2 3 4 5 
do
	data=$dataset$i
	max=${maxs[i]}
	echo "~~~~~", $data, $max, "~~~~~"
	#python refex.py $data
	#./ann/bin/ann_sample -d $D -max $max -nn $topK -df "result/"$data".ReFex" -qf "result/"$data".ReFex" -rf "result/"$data".topK.ReFex"
	python calReFexStructureSimScore.py $data "ReFex" $topK 
done


