
import sys
import networkx as nx
import re

if len(sys.argv) < 3:
	print "Usage: transfer2gml dataset username"
	exit()

dataset = sys.argv[1]
username = sys.argv[2]
value_file = "result/"+username
network_file = "data/"+dataset+".graph"
map_file = "data/"+dataset+".dict"


G = nx.Graph()

# read node map
name2id = {}
for line in open(map_file, "r"):
    elements = re.compile('\t').split(line.strip())
    name = elements[0]
    id = int(elements[1])
    name2id[name] = id
    G.add_node(id, name=name)

node_num = G.number_of_nodes()

# read network
i = 0
for line in open(network_file, "r"):
	if i == 0:
		i+=1
		continue
	i +=1
	elements = re.compile('\t').split(line.strip())
	node1 = int(elements[0])
	node2 = int(elements[1])
	w = float(elements[2])
	G.add_edge(node1, node2, weight=w)
	G.add_edge(node2, node1, weight=w) 
    
# read value
nodeValues = []
for line in open(value_file, "r"):
	nodeValue = float(line.strip())
	nodeValues.append(nodeValue)

graphFile = "result/"+username+".gml"
f = open(graphFile, "w")

f.write("graph ["+"\n")
for node in G.nodes(data=True):
	f.write("\tnode\n\t[\n")
	f.write("\t\tid "+str(node[0])+"\n")
	f.write("\t\tvalue "+str(nodeValues[node[0]])+"\n") 
	f.write('\t\tlabel "'+str(node[1]['name'])+'"'+"\n")
	f.write("\t]\n")

for edge in G.edges(data=True):
	f.write("\tedge\n\t[\n")
	f.write("\t\tsource "+str(edge[0])+"\n")
	f.write("\t\ttarget "+str(edge[1])+"\n")
	f.write("\t\tvalue "+str(edge[2]['weight'])+"\n")
	f.write("\t]\n")
f.write("]")

f.close()
