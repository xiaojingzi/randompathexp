
dataList=("kdd" "icdm" "sigir" "cikm" "sigmod" "icde")
# dataList=("twitter2")
maxs=(2867 2607 2851 3548 2616 2559)
# maxs=(112411)
T=10
D=50
D2=18
Epsilon=0.01
topK=50
for i in  0 1 2 3 4 5 
do
    data=${dataList[i]}
    max=${maxs[i]}
    echo "data=$data, max=$max"
    # python findCommonNeighbors.py $data

    #=========pathsim ==========
    # ./main_release $data $T $D $Epsilon
    # python randomproject.py $data $T $D $Epsilon
    # cp "result/"$data".pathsim" "result/"$data".topK.pathsim"
    python calNeighborSimScore.py $data "pathsim"


    # #=========pathVec ===========
    # ./ann/bin/ann_sample -d $D -max $max -nn $topK -df "result/"$data".pathvec" -qf "result/"$data".pathvec" -rf "result/"$data".topK.pathvec"
    python calNeighborSimScore.py $data "pathvec"
    

    # #=========refex =============
    # python refex.py $data
    # ./ann/bin/ann_sample -d $D2 -max $max -nn $topK -df "result/"$data".ReFex" -qf "result/"$data".ReFex" -rf "result/"$data".topK.ReFex"
    python calNeighborSimScore.py $data "ReFex"
    

    #=========rolesim==============
    #   # convet the graph format for the rolesim code
    # python convertGraphFormat.py $data
    # ./rolesim/main -g  "data/"$data".RoleSimGraph" -o "result/"$data".RoleSim"
    # python getTopK.py $data $topK "RoleSim"
    python calNeighborSimScore.py $data "RoleSim"
        

    # #==========simrank ===========
    # python simrank.py $data
    # python getTopK.py $data $topK "SimRank"
    python calNeighborSimScore.py $data "SimRank"

done

for i in 0 1 2 3 4 5 
do
    data=${dataList[i]}
    echo $data
    python calStructureSimScore.py $data "pathsim"
    python calStructureSimScore.py $data "pathvec"
    python calStructureSimScore.py $data "ReFex"
    python calStructureSimScore.py $data "RoleSim"
    python calStructureSimScore.py $data "SimRank"
done
