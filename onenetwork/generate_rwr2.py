
import sys
import re
import numpy as np
import numpy.linalg as linalg
import networkx as nx
import time

if len(sys.argv) < 3:
    print "usage: data c"
    exit()

data = sys.argv[1]
c = float(sys.argv[2])
network_file = "data/"+data+".graph"
map_file ="data/"+ data+".dict"


G = nx.Graph()

for line in open(map_file, "r"):
    elements = re.compile('\t').split(line.strip())
    if len(elements) == 2:
        name = elements[0]
        id = int(elements[1])
    elif len(elements) == 1:
        name = "no name"
        id = int(elements[0])
    G.add_node(id, name=name, weight=0.0)
N = G.number_of_nodes()


i =0
for line in open(network_file, "r"):
    if i ==0:
        es = line.split('\t')
    else:
        elements = re.compile('\t').split(line.strip())
        node1 = int(elements[0])
        node2 = int(elements[1])
        w = float(elements[2])
        G.add_edge(node1, node2, weight=w)
        G.add_edge(node2, node1, weight=w)
    i +=1
M = G.number_of_edges()
print "node number = ", N , ", edge number = ", M
    
# conduct random walk with restart until convergence
v = np.diag([1.0]*N)
t1 = time.time()
start = time.time()
for i in range(N):
    
    count = 0
    # start = time.time()
    while True:
        
        v_old = v[i].copy()
        v[i] = 0
        
        for node in G.nodes():
            # reach probability in last step is zero
            if v_old[node] == 0:
                continue
            # walk to neighbors with probability 1-c
            for neighbor in G.neighbors(node):
                v[i][neighbor] += (1-c) * v_old[node]*G[node][neighbor]['weight'] / G.degree(node, weight='weight')
            
            # walk to start node i with probability c
            v[i][i] += c * v_old[node]

        diff = np.linalg.norm(v[i] - v_old, ord=1)
        count += 1
        if (diff < np.expm1(1e-5)):
            break
    # end = time.time()
    if i % 1 == 0:
        	end = time.time()
		print "Time of ",i," is ", (end -start)
		start = end

end = time.time()
print "time=",(end - t1)
output_file = "result/"+data+".rwr"
f = open(output_file, "w")
d = v.shape[1]
for i in range(N):
    for j in range(N):
        #if  not np.isnan(v[i][j]):
        f.write(str(v[i][j])+" ")
    f.write("\n")
f.close()
