
dataList=("kdd" "icdm" "sigir" "cikm" "sigmod" "icde" "twitter2" "mobile2")
topK=50
for i in  6 7
do
    data=${dataList[i]}
    echo "data=$data"
    python convertTopSimRank.py $data

    # #==========simrank ===========
    java -jar -Xmx200G  topsimrank.jar "data/"$data".SimRankGraph" "result/"$data".topK.SimRank"
    python calNeighborSimScore.py $data "SimRank" $topK

done

for i in 6 7 
do
    data=${dataList[i]}
    python calStructureSimScore.py $data "SimRank" $topK
done
