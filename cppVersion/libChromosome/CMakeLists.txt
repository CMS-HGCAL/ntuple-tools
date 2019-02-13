add_library(libChromosome src/Chromosome.cpp
                          src/ChromosomeProcessor.cpp)

set_property(TARGET libChromosome PROPERTY CXX_STANDARD 17)
