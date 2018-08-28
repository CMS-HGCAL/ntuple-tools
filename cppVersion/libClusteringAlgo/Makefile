CC  = g++

current_dir = $(shell pwd)

LIB_DIR = ${current_dir}/lib
OBJ_DIR = tmp
SRC_DIR = src
INC_DIR = include

LDFLAGS = `root-config --libs` -Wall -Wextra -g -O0
CCFLAGS = `root-config --cflags` -g -c -std=c++1y -Wall -Wextra -O0 -I${INC_DIR}

all: ${LIB_DIR}/libAlgo.so

MODULES = ${OBJ_DIR}/Event.o ${OBJ_DIR}/GenParticles.o ${OBJ_DIR}/RecHits.o ${OBJ_DIR}/SimClusters.o ${OBJ_DIR}/ImagingAlgo.o ${OBJ_DIR}/Hexel.o ${OBJ_DIR}/BasicCluster.o ${OBJ_DIR}/ConfigurationManager.o ${OBJ_DIR}/Clusters2D.o ${OBJ_DIR}/ClusterMatcher.o ${OBJ_DIR}/MatchedClusters.o

${LIB_DIR}/libAlgo.so: $(MODULES)
	$(CC) -shared -o $@ $^ $(LDFLAGS)

${OBJ_DIR}/%.o: ${SRC_DIR}/%.cpp
	$(CC) $^ -o $@ $(CCFLAGS)

clean:
	rm -f ${OBJ_DIR}/*.o ${LIB_DIR}/libAlgo.so
