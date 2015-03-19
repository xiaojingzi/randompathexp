import re

fnode = open("netscience_node.gml", "r")
fdict = open("netscience.dict","w")
lines = fnode.readlines()
for i in range(len(lines)/5):
	idStr = lines[i*5+2]
	id = idStr.strip().split(" ")[1]
	labelStr = lines[i*5+3]
	label = re.split('"', labelStr)[1].lower()
	fdict.write(label+"\t"+id+"\n")

fdict.close()
fdict.close()

fedge = open("netscience_edge.gml", "r")
fgraph = open("netscience.graph","w")
lines = fedge.readlines()
for i in range(len(lines)/6):
	id1Str = lines[i*6+2]
	id1 = id1Str.strip().split(" ")[1]

	id2Str = lines[i*6+3]
	id2 = id2Str.strip().split(" ")[1]

	valueStr = lines[i*6+4]
	value = valueStr.strip().split(" ")[1]

	fgraph.write(id1+"\t"+id2+"\t"+value+"\n")
fgraph.close()
