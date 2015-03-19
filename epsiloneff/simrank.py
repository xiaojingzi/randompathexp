import sys
import pickle
import networkx as nx
import re
import time
from collections import defaultdict
import copy

def simrank(G, r=0.9, max_iter=100):
  # init. vars
  sim_old = defaultdict(list)
  sim = defaultdict(list)
  for n in G.nodes():
    sim[n] = defaultdict(int)
    sim[n][n] = 1
    sim_old[n] = defaultdict(int)
    sim_old[n][n] = 0


  # recursively calculate simrank
  for iter_ctr in range(max_iter):
    if _is_converge(sim, sim_old):
      break
    sim_old = copy.deepcopy(sim)
    for u in G.nodes():
      for v in G.nodes():
        if u == v:
          continue
        if len(G.neighbors(u)) == 0 or len(G.neighbors(v)) == 0:
          continue
        s_uv = 0.0
        for n_u in G.neighbors(u):
          for n_v in G.neighbors(v):
            s_uv += sim_old[n_u][n_v]
        sim[u][v] = (r * s_uv / (len(G.neighbors(u)) * len(G.neighbors(v))))
      # if u % 100 == 0:
        # print u

  return sim

def _is_converge(s1, s2, eps=1e-4):
  for i in s1.keys():
    for j in s1[i].keys():
      if abs(s1[i][j] - s2[i][j]) >= eps:
        return False
  return True


def main(map_file, network_file, sim_file):
  
  t1 = time.time()
  G = nx.Graph()

  # read node map
  for line in open(map_file, "r"):
    elements = re.compile('\t').split(line.strip())
    if len(elements) == 2:
      name = elements[0]
      id = int(elements[1])
    elif len(elements) == 1:
      name = "no name"
      id = int(elements[0])
    G.add_node(id, name=name)
  # print "load nodes done."

  # read network
  i = 0

  for line in open(network_file, "r"):
    if i == 0:
      i +=1 
      continue
    i+=1
    elements = re.compile('\t').split(line.strip())
    node1 = int(elements[0])
    node2 = int(elements[1])
    G.add_edge(node1, node2)
    G.add_edge(node2, node1)
  # print "load edges dones."

  N = G.number_of_nodes()
  M = G.number_of_edges()
  t2 = time.time()
  print "Time of loading = ", (t2-t1)

  # print "#nodes = ", N, ", #edges = ", M

  # calcualte simrank on the graph
  start = time.time()
  sim = simrank(G)
  end = time.time()
  print "simrank time cost = ", (end - start)

  t3 = time.time()
  # print sim
  f = open(sim_file, "w")
  for i in range(N):
    for j in sim[i]:
      f.write(str(j) + "\t")
    f.write('\n')
  f.close()
  #pickle.dump(sim, open(sim_file, "wb"))

  t4 = time.time()
  print "Time of saving = ", (t4 - t3)

if __name__ == '__main__':
  if len(sys.argv) < 2:
      print "usage: simrank.py dataset"
      exit()
  dataset = sys.argv[1]

  map_file = "data/"+dataset + ".dict"
  network_file = "data/"+dataset + ".graph"
  sim_file =  "result/"+dataset +".SimRank"
  main(map_file, network_file, sim_file)
