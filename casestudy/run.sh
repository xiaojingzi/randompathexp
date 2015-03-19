T=5
D=50
Epsilon=0.01
data="netscience"
./main_release $data $T $D $Epsilon
python calculatesimilarity.py $data "newman,_m"
python calculatesimilarity.py $data "rinzel,_j"
python calculatesimilarity.py $data "robert,_f"
python transfer2gml.py $data "newman,_m"
python transfer2gml.py $data "rinzel,_j"
python transfer2gml.py $data "robert,_f"