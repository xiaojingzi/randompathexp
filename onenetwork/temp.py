import sys

if len(sys.argv) < 2:
    print "usage: dataset"
    exit()
dataset = sys.argv[1]

filen = "result/"+dataset+"top50.simrank"
f = open(filen, "r")
wf = open("result/"+dataset+".topK.SimRank", "w")
for line in f:
	elements = line.split(" ")
	for element in elements:
		a = element.split(":")
		wf.write(a[0]+" ")
wf.close()
