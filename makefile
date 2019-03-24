CC = g++

CFLAGS  = -Wall -g -O3 -march=native -I/usr/local/include -std=c++11 -I../../include
LDFLAGS = -lm -lgsl -lgflags -lglog -lgslcblas -L/usr/local/lib
BIN = embedder
OBJ = hyperboliclinear.o hyperbolic.o graph.o random.o springembedder.o geometricds.o powerlawCommon.o NLEHelper.o

main: main.o
	$(CC) $(CFLAGS) $(DEFINE) -o $(BIN) main.o $(OBJ) $(LDFLAGS)

main.o: hyperboliclinear.o hyperbolic.o graph.o springembedder.o NLEHelper.o
NLEHelper.o: plfit/powerlawCommon.o
springembedder.o: geometricds.o
hyperboliclinear.o: hyperbolic.o graph.o random.o geometricds.o
graph.o: random.o
hyperbolic.o: graph.o random.o


%.o: %.cpp
	$(CC) $(CFLAGS) -c $<

.PHONY: clean
clean:
	rm -rf $(BIN) $(OBJ) main.o