CC  = g++

LDFLAGS = `root-config --libs` -Wall -Wextra -g -O0 -lMinuit
CCFLAGS = `root-config --cflags` -g -c -std=c++1y -Wall -Wextra -O0 -I./include/

TMP_DIR = tmp

all: rootOptimizer

rootOptimizer: ${TMP_DIR}/rootOptimizer.o
	$(CC) $^ -o $@ $(LDFLAGS)

${TMP_DIR}/rootOptimizer.o: rootOptimizer.cpp
	@mkdir -p $(@D)
	$(CC) $^ -o $@ $(CCFLAGS)

${TMP_DIR}/%.o: src/%.cpp
	$(CC) $^ -o $@ $(CCFLAGS)

clean:
	rm -f ${TMP_DIR}/*.o rootOptimizer
