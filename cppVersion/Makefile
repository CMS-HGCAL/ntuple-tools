CC  = g++

LDFLAGS = `root-config --libs` -Wall -Wextra -g -O0 -fsanitize=address
CCFLAGS = `root-config --cflags` -g -c -std=c++11 -Wall -Wextra -O0 -ffast-math -fsanitize=address

all: algoBenchmark

algoBenchmark: algoBenchmark.o HGEvent.o HGGenParticles.o HGRecHits.o HGSimClusters.o
	$(CC) $^ -o $@ $(LDFLAGS)

%.o: %.cpp
	$(CC) $^ -o $@ $(CCFLAGS)

clean:
	rm -f *.o algoBenchmark
