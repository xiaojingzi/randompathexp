CC	= g++
CPPFLAGS= -Wno-deprecated -O2 -c -pg
LDFLAGS	= -O2 -pg
SOURCES	= IcebergSimMap.cpp Util.cpp SimMatrix.cpp GraphRoleUtil.cpp Graph.cpp PM/misc.cpp PM/PMinterface.cpp PM/PMinit.cpp PM/PMduals.cpp PM/PMexpand.cpp PM/PMshrink.cpp PM/PMmain.cpp PM/MinCost/MinCost.cpp  main.cpp
OBJECTS	= $(SOURCES:.cpp=.o)
EXECUTABLE=main

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE) : $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o : 
	$(CC) $(CPPFLAGS) $< -o $@

clean:
	-rm -f *.o PM/*.o PM/*/*.o
