import sys
if len(sys.argv) < 3:
    print "usage: data D"
    exit()
data = sys.argv[1]
D = sys.argv[2]
f = open("result/"+data + '_D.recall', 'a')
f.write('#'+D + '\n')
f.close()

