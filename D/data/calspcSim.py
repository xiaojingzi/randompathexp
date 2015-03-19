import sys
import networkx as nx
import time
import re
import numpy as np

def init( data_param=None):
    global dataset, network_file, map_file
    dataset = data_param

    network_file= dataset+".graph"
    map_file= dataset+".dict"
    

def loadNode():
    global G, map_file, N
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
    # print "load nodes done."
    N = G.number_of_nodes()


def loadLink():
    global G, newtork_file, N, M
    e = 0

    for line in open(network_file, "r"):
        elements = re.compile('\t').split(line.strip())
        node1 = int(elements[0])
        node2 = int(elements[1])
        w = float(elements[2])
        G.add_edge(node1, node2, weight=w)
        G.node[node1]['weight'] += w
        G.add_edge(node2, node1, weight=w)
        G.node[node2]['weight'] += w
        e +=1
        # if e % 10000 == 0:
            # print e, " edge"

    # print "load edges dones."
    M = G.number_of_edges()
    # print "#nodes = ", N, ", #edges = ", M


def getspcSim(data=None):
    t1 = time.time()
    init(data)
    loadNode()

    # print "Time of loading nodes = ", (t2 - t1)
    loadLink()
    # t3 = time.time()
    # print "Time of loading edges = ", (t3 - t2)
    input_file = dataset +'.spcid'
    spcIDdict = {}
    f = open(input_file, 'r')
    for line in f:
        line = line.strip()
        id = int(line)
        spcIDdict[id] = 1
    f.close()
    print spcIDdict.keys()
    spcSim = np.zeros(shape=(N,N))
    for u in G.nodes():
        for v in G.nodes():
            if u in spcIDdict and v in spcIDdict:
                spcSim[u][v] = 1
            elif (not u in spcIDdict) and (not v in spcIDdict):
                spcSim[u][v] = 1
    np.save(data + 'spcSim',spcSim)
    for u in range(5):
        for v in range(5):
            print spcSim[u][v],
        print '\n'

if len(sys.argv) < 2:
    print "usage: dataset"
    print "output: dataset + 'spcSim'"
    exit()
data=sys.argv[1]

# print "dataset=", data
getspcSim(data=data)
