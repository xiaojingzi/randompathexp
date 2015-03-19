import sys
if len(sys.argv) < 4:
    print "usage: parametervalue simtype testparameter"
    exit()
parametervalue = sys.argv[1]
simtype = sys.argv[2]
testparameter = sys.argv[3]
f = open("result/"+testparameter+'.'+simtype+".score", 'a')
f.write('#'+parametervalue + '\n')
f.close()

