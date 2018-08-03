CC  = g++

current_dir = $(shell pwd)

OBJ_DIR = tmp
SRC_DIR = src

ALGO_DIR = ${current_dir}/libClusteringAlgo
ALGO_INC = ${ALGO_DIR}/include
ALGO_LIB = ${ALGO_DIR}/lib

GEN_DIR = ${current_dir}/libChromosome
GEN_LIB = ${GEN_DIR}/lib
GEN_INC = ${GEN_DIR}/include

EXECS = algoBenchmark createQualityPlots matchingBenchmark geneticOptimizer rootOptimizer

LDFLAGS = `root-config --libs` -Wall -Wextra -g -O0 -L${ALGO_LIB} -L${GEN_LIB}
CCFLAGS = `root-config --cflags` -g -c -std=c++1y -Wall -Wextra -O0 -I${ALGO_INC} -I${GEN_INC}

all: ${EXECS}

algoBenchmark: ${OBJ_DIR}/algoBenchmark.o ${ALGO_LIB}/libAlgo.so
	$(CC) $^ -o $@ $(LDFLAGS) -lAlgo -Wl,-rpath,${ALGO_LIB}/libAlgo.so

createQualityPlots: ${OBJ_DIR}/createQualityPlots.o ${ALGO_LIB}/libAlgo.so
	$(CC) $^ -o $@ $(LDFLAGS) -lAlgo -Wl,-rpath,${ALGO_LIB}/libAlgo.so

matchingBenchmark: ${OBJ_DIR}/matchingBenchmark.o ${ALGO_LIB}/libAlgo.so
	$(CC) $^ -o $@ $(LDFLAGS) -lAlgo -Wl,-rpath,${ALGO_LIB}/libAlgo.so

geneticOptimizer: ${OBJ_DIR}/geneticOptimizer.o ${GEN_LIB}/libChromo.so
	$(CC) $^ -o $@ $(LDFLAGS) -lChromo -Wl,-rpath,${GEN_LIB}/libChromo.so

rootOptimizer: ${OBJ_DIR}/rootOptimizer.o ${ALGO_LIB}/libAlgo.so
	$(CC) $^ -o $@ $(LDFLAGS) -lAlgo -Wl,-rpath,${ALGO_LIB}/libAlgo.so -lMinuit

${ALGO_LIB}/libAlgo.so:
	+$(MAKE) -C ${ALGO_DIR}

${GEN_LIB}/libChromo.so:
	+$(MAKE) -C ${GEN_DIR}

${OBJ_DIR}/%.o: ${SRC_DIR}/%.cpp
	$(CC) $^ -o $@ $(CCFLAGS)

clean:
	rm -f ${OBJ_DIR}/*.o ${EXECS}
	cd ${ALGO_DIR} && $(MAKE) clean
	cd ${GEN_DIR} && $(MAKE) clean
