import sys
if len(sys.argv) < 3:
    print "usage: simtype testparameter"
    exit()
simtype = sys.argv[1]
testparameter = sys.argv[2]

f = open("result/"+testparameter+'.'+simtype+'.score', 'w')
f.close()
