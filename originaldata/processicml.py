import sys
from collections import defaultdict


names = []
name2id = {}
edgeMap = defaultdict(int)


for line in open("icml.txt", "r"):
	authorList = line.strip().split("\\\t")
	print authorList

	for i in range(len(authorList)):
		name1 = authorList[i].strip().lower()
		if name1 not in name2id:
			name2id[name1] = len(names) 
			names.append(name1) 
		id1 = name2id[name1]
		for j in range(i+1, len(authorList)):
			name2 = authorList[j].strip().lower()
			if name2 not in name2id:
				name2id[name2] = len(names)
				names.append(name2)
			id2 = name2id[name2]

			edgeName1 = str(id1)+"_"+str(id2)
			edgeName2 = str(id2)+"_"+str(id1)
			print edgeName1, edgeName2
			if edgeName1 in edgeMap:
				edgeMap[edgeName1] +=1
			elif edgeName2 in edgeMap:
				edgeMap[edgeName2] +=1
			else:
				edgeMap[edgeName1] +=1

print len(names), len(edgeMap)
print "output" 
# save names to file
outputNameFile = "icml.dict" 
f = open(outputNameFile, "w")
for i in range(len(names)):
	f.write(names[i]+"\t"+str(i)+"\n")
f.close()

outputNetworkFile = "icml.graph"
f = open(outputNetworkFile, "w")
f.write(str(len(names))+"\t"+str(len(edgeMap))+"\n");
for edge, weight in edgeMap.items():
	nodes = edge.split("_")
	f.write(str(nodes[0])+"\t"+str(nodes[1])+"\t"+str(weight)+"\n")
f.close()