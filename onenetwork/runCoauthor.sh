
dataList=("kdd" "icdm" "sigir" "cikm" "sigmod" "icde")
# dataList=("twitter2")
maxs=(2867 2607 2851 3548 2616 2559)
# maxs=(112411)
T=5
D=50
D2=18
Epsilon=0.01
topK=50
for i in  0 1 2 3 4 5 
do
    data=${dataList[i]}
    max=${maxs[i]}
    echo "data=$data, max=$max"
    # # python findCommonNeighbors.py $data

    #=========pathsim ==========
    # ./main_release $data $T $D $Epsilon
    # python randomproject.py $data $T $D $Epsilon
    # cp "result/"$data".pathsim" "result/"$data".topK.pathsim"
    #python calNeighborSimScore.py $data "pathsim" $topK


    # #=========pathVec ===========
    # ./ann/bin/ann_sample -d $D -max $max -nn $topK -df "result/"$data".pathvec" -qf "result/"$data".pathvec" -rf "result/"$data".topK.pathvec"
    #python calNeighborSimScore.py $data "pathvec" $topK
    

    # #=========refex =============
    # python refex.py $data
    # ./ann/bin/ann_sample -d $D2 -max $max -nn $topK -df "result/"$data".ReFex" -qf "result/"$data".ReFex" -rf "result/"$data".topK.ReFex"
    #python calNeighborSimScore.py $data "ReFex" $topK
    

    #=========rolesim==============
      # convet the graph format for the rolesim code
    # python convertGraphFormat.py $data
    # ./rolesim/main -g  "data/"$data".RoleSimGraph" -m "3-Iceberg" -o "result/"$data".RoleSim"
    # python getTopK.py $data $topK "RoleSim"
    #python calNeighborSimScore.py $data "RoleSim" $topK
        

    # #==========simrank ===========
    # python simrank.py $data
    # python getTopK.py $data $topK "SimRank"
    #python calNeighborSimScore.py $data "SimRank" $topK

   ## ===========random walk with restart ====
   python generate_rwr2.py $data  0.9
   python getTopK.py $data $topK "rwr"
   python calNeighborSimScore.py $data "rwr" $topK
done

for i in 0 1 2 3 4 5 
do
    data=${dataList[i]}
    #python calStructureSimScore.py $data "pathsim" $topK
    #python calStructureSimScore.py $data "pathvec" $topK
    #python calStructureSimScore.py $data "ReFex" $topK
    #python calStructureSimScore.py $data "RoleSim" $topK
    #python calStructureSimScore.py $data "SimRank" $topK
    python calStructureSimScore.py $data "rwr" $topK
done
