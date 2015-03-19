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


for i in 0 1 2 3 4 5 6 7 8
do
	data=$dataset$i
	
	echo "~~~~~", $data, "~~~~~"
	python convertTopSimRank.py $data
	java -jar -Xmx300G  topsimrank.jar "data/"$data".SimRankGraph" "result/"$data".topK.SimRank"
done

