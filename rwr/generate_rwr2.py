
import sys
import re
import numpy as np
import numpy.linalg as linalg
import networkx as nx
import time

if len(sys.argv) < 3:
    print "usage: datafolder c"
    exit()

datafolder = sys.argv[1]
c = float(sys.argv[2])
network_file = datafolder+"/nwdata"
map_file = datafolder+"/node.dict"


G = nx.Graph()

# read node map
for line in open(map_file, "r"):
    elements = re.compile('\t').split(line.strip())
    name = elements[0]
    id = int(elements[1])
    G.add_node(id, name=name)

# read network
for line in open(network_file, "r"):
    elements = re.compile('\t').split(line.strip())
    node1 = int(elements[0])
    node2 = int(elements[1])
    weight = float(elements[2])
    G.add_edge(node1, node2, weight=weight)
    G.add_edge(node2, node1, weight=weight)

N = G.number_of_nodes()
M = G.number_of_edges()
    
# conduct random walk with restart until convergence
v = np.diag([1.0]*N)
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
    if i % 100 == 0:
        print i

end = time.time()
print "time=",(end - start)
# output rwr score
output_file = "result/"+datafolder+"/rwr"
np.savez(output_file, x=v)