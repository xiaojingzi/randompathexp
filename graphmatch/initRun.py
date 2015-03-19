import sys
if len(sys.argv) < 2:
    print "usage: data"
    exit()
data = sys.argv[1]


print 'start to deal with ' + data
f = open("result/"+data + '_D.recall', 'w')
f.close()
