CC  = g++

LDFLAGS = `root-config --libs` -Wall -Wextra -g -O0
CCFLAGS = `root-config --cflags` -g -c -std=c++1y -Wall -Wextra -O0 -I./include/

TMP_DIR = tmp

all: geneticOptimizer

geneticOptimizer: ${TMP_DIR}/geneticOptimizer.o ${TMP_DIR}/Chromosome.o
	$(CC) $^ -o $@ $(LDFLAGS)

${TMP_DIR}/geneticOptimizer.o: geneticOptimizer.cpp
	@mkdir -p $(@D)
	$(CC) $^ -o $@ $(CCFLAGS)

${TMP_DIR}/%.o: src/%.cpp
	$(CC) $^ -o $@ $(CCFLAGS)

clean:
	rm -f ${TMP_DIR}/*.o geneticOptimizer
