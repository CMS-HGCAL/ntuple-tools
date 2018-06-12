CC  = g++

LDFLAGS = `root-config --libs` -Wall -Wextra -g -O0 -fsanitize=address
CCFLAGS = `root-config --cflags` -g -c -std=c++11 -Wall -Wextra -O0 -ffast-math -fsanitize=address

all: algoBenchmark

algoBenchmark: tmp/algoBenchmark.o tmp/HGEvent.o tmp/HGGenParticles.o tmp/HGRecHits.o tmp/HGSimClusters.o tmp/HGCalImagineAlgo.o tmp/Hexel.o tmp/BasicCluster.o
	$(CC) $^ -o $@ $(LDFLAGS)

tmp/%.o: %.cpp
	$(CC) $^ -o $@ $(CCFLAGS)

clean:
	rm -f tmp/*.o algoBenchmark
