OBJECTS=TSMA.o ccdeRepart.o base.o wdn.o

wdn: $(OBJECTS)
	g++ -std=c++11 -o wdn $(OBJECTS) libtoolkit.a -lm

wdn.o: wdn.cpp base.h
	g++ -std=c++11 -c wdn.cpp

base.o: base.h toolkit.h base.cpp
	g++ -std=c++11 -c base.cpp

TSMA.o: TSMA.h base.h TSMA.cpp
	g++ -std=c++11 -c TSMA.cpp

ccdeRepart.o: ccdeRepart.h base.h ccdeRepart.cpp
	g++ -std=c++11 -c ccdeRepart.cpp

.PHONY: clean
clean:
	rm -f wdn $(OBJECTS) *.rep *.par

