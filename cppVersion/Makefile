CC  = g++

LDFLAGS = `root-config --libs` -Wall -Wextra -g -O0 -fsanitize=address
CCFLAGS = `root-config --cflags` -g -c -std=c++11 -Wall -Wextra -O0 -ffast-math -fsanitize=address -I./include/

all: algoBenchmark

algoBenchmark: tmp/algoBenchmark.o tmp/Event.o tmp/GenParticles.o tmp/RecHits.o tmp/SimClusters.o tmp/ImagineAlgo.o tmp/Hexel.o tmp/BasicCluster.o
	$(CC) $^ -o $@ $(LDFLAGS)

tmp/algoBenchmark.o: algoBenchmark.cpp
	$(CC) $^ -o $@ $(CCFLAGS)

tmp/%.o: src/%.cpp
	$(CC) $^ -o $@ $(CCFLAGS)

clean:
	rm -f tmp/*.o algoBenchmark
