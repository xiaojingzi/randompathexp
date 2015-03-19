datasets=("kdd" "icdm" "sigir" "cikm" "sigmod" "icde" "twitter2" "mobile2")

for dataset in ${datasets[@]}
do
	python transfer2rwrformat.py $dataset /home/jing/randompath/onenetwork/data/$dataset
done
