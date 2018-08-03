CC  = g++

current_dir = $(shell pwd)
main_dir = $(shell cd ..;pwd)

OBJ_DIR = tmp
SRC_DIR = src
INC_DIR = include
LIB_DIR = ${current_dir}/lib

ALGO_LIB = ${main_dir}/libClusteringAlgo/lib
ALGO_INC = ${main_dir}/libClusteringAlgo/include

LDFLAGS = `root-config --libs` -Wall -Wextra -g -O0 -L${ALGO_LIB} -lAlgo -Wl,-rpath,${ALGO_LIB}/libAlgo.so
CCFLAGS = `root-config --cflags` -g -c -std=c++1y -Wall -Wextra -O0 -I${INC_DIR} -I${ALGO_INC}

all: ${LIB_DIR}/libChromo.so

MODULES = ${OBJ_DIR}/Chromosome.o

${LIB_DIR}/libChromo.so: $(MODULES) ${ALGO_LIB}/libAlgo.so
	$(CC) -shared -o $@ $^ $(LDFLAGS)

${OBJ_DIR}/%.o: ${SRC_DIR}/%.cpp
	$(CC) $^ -o $@ $(CCFLAGS)

clean:
	rm -f ${OBJ_DIR}/*.o ${LIB_DIR}/libChromo.so
