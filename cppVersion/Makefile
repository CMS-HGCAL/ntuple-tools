CC  = g++

LDFLAGS = `root-config --libs` -Wall -Wextra -g -O0
CCFLAGS = `root-config --cflags` -g -c -std=c++1y -Wall -Wextra -O0 -I./include/

TMP_DIR = tmp

all: algoBenchmark createQualityPlots

algoBenchmark: ${TMP_DIR}/algoBenchmark.o ${TMP_DIR}/Event.o ${TMP_DIR}/GenParticles.o ${TMP_DIR}/RecHits.o ${TMP_DIR}/SimClusters.o ${TMP_DIR}/ImagingAlgo.o ${TMP_DIR}/Hexel.o ${TMP_DIR}/BasicCluster.o ${TMP_DIR}/ConfigurationManager.o
	$(CC) $^ -o $@ $(LDFLAGS)

createQualityPlots: ${TMP_DIR}/createQualityPlots.o ${TMP_DIR}/Event.o ${TMP_DIR}/GenParticles.o ${TMP_DIR}/RecHits.o ${TMP_DIR}/SimClusters.o ${TMP_DIR}/ImagingAlgo.o ${TMP_DIR}/Hexel.o ${TMP_DIR}/BasicCluster.o ${TMP_DIR}/ConfigurationManager.o
	$(CC) $^ -o $@ $(LDFLAGS)


${TMP_DIR}/algoBenchmark.o: algoBenchmark.cpp
	@mkdir -p $(@D)
	$(CC) $^ -o $@ $(CCFLAGS)

${TMP_DIR}/createQualityPlots.o: createQualityPlots.cpp
	@mkdir -p $(@D)
	$(CC) $^ -o $@ $(CCFLAGS)

${TMP_DIR}/%.o: src/%.cpp
	$(CC) $^ -o $@ $(CCFLAGS)

clean:
	rm -f ${TMP_DIR}/*.o algoBenchmark createQualityPlots
