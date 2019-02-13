cmake_minimum_required (VERSION 2.6)
project (HGCalClustering)

find_package(ROOT REQUIRED COMPONENTS Minuit)
include(${ROOT_USE_FILE})

add_executable(createQualityPlots src/createQualityPlots.cpp)
add_executable(algoBenchmark src/algoBenchmark.cpp)
add_executable(rootOptimizer src/rootOptimizer.cpp)
add_executable(geneticOptimizer src/geneticOptimizer.cpp)
add_executable(analyzeTestbeam src/analyzeTestbeam.cpp)

set_property(TARGET createQualityPlots PROPERTY CXX_STANDARD 17)
set_property(TARGET algoBenchmark PROPERTY CXX_STANDARD 17)
set_property(TARGET rootOptimizer PROPERTY CXX_STANDARD 17)
set_property(TARGET geneticOptimizer PROPERTY CXX_STANDARD 17)
set_property(TARGET analyzeTestbeam PROPERTY CXX_STANDARD 17)


include_directories("${PROJECT_SOURCE_DIR}/libChromosome/include")
include_directories("${PROJECT_SOURCE_DIR}/libClusteringAlgo/include")

add_subdirectory(libChromosome)
add_subdirectory(libClusteringAlgo)

target_link_libraries(createQualityPlots libChromosome libClusteringAlgo ${ROOT_LIBRARIES})
target_link_libraries(algoBenchmark libChromosome libClusteringAlgo ${ROOT_LIBRARIES})
target_link_libraries(rootOptimizer libChromosome libClusteringAlgo ${ROOT_LIBRARIES})
target_link_libraries(geneticOptimizer libChromosome libClusteringAlgo ${ROOT_LIBRARIES})
target_link_libraries(analyzeTestbeam libChromosome libClusteringAlgo ${ROOT_LIBRARIES})
