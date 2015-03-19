from collections import defaultdict
f = open("mobile.txt", "r")
gf1 = open("mobile1.graph", "w")
gf2 = open("mobile2.graph", "w")

df = open("mobile.dict", "w")
nodedict = {}
nodedict1 = {}
nodedict2 = {}
edges1 = defaultdict(float)
edges2 = defaultdict(float)
for line in f:
	elements = line.strip().split(",")
	id1 = elements[0]
	id2 = elements[1]
	key = id1+"_"+id2
	week = int(elements[2])
	weight = float(elements[3])

	if week == 1 or week == 2:
		edges1[key] += weight
		nodedict1[id1] = id1
		nodedict1[id2] = id2

	elif week == 3 or week == 4:
		edges2[key] += weight
		nodedict2[id1] = id1
		nodedict2[id2] = id2

	nodedict[id1] = id1
	nodedict[id2] = id2

print "dict len = ", len(nodedict)
print "dict1 len = ", len(nodedict1)
print "dict2 len = ", len(nodedict2)
print "edges in week 1 & 2 = ", len(edges1)
print "edges in week 3 & 4 = ", len(edges2)
gf1.write(str(len(nodedict))+"\t"+str(len(edges1))+"\n")
for key,value in edges1.items():
	keyList = key.split("_")
	gf1.write(keyList[0]+"\t"+keyList[1]+"\t"+str(value)+"\n")
gf1.close()

gf2.write(str(len(nodedict))+"\t"+str(len(edges2))+"\n")
for key,value in edges2.items():
	keyList = key.split("_")
	gf2.write(keyList[0]+"\t"+keyList[1]+"\t"+str(value)+"\n")
gf2.close()
	

for k, v in nodedict.items():
	df.write(k+"\t"+v+"\n")
df.close()
