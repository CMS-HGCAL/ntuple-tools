cmake_minimum_required (VERSION 2.6)
project (HGCalClustering)

find_package(ROOT REQUIRED COMPONENTS Minuit)
include(${ROOT_USE_FILE})

set_property(TARGET ROOT::Core PROPERTY
    INTERFACE_INCLUDE_DIRECTORIES "${ROOT_INCLUDE_DIRS}")

add_library(ROOT::Flags_CXX IMPORTED INTERFACE)

# Fix for ROOT_CXX_FLAGS not actually being a CMake list
separate_arguments(ROOT_CXX_FLAGS)
set_property(TARGET ROOT::Flags_CXX APPEND PROPERTY
    INTERFACE_COMPILE_OPTIONS ${ROOT_CXX_FLAGS})

# Add definitions
separate_arguments(ROOT_DEFINITIONS)
foreach(_flag ${ROOT_EXE_LINKER_FLAG_LIST})
    # Remove -D or /D if present
    string(REGEX REPLACE [=[^[-//]D]=] "" _flag ${_flag})
    set_property(TARGET ROOT::Flags APPEND PROPERTY INTERFACE_LINK_LIBRARIES ${_flag})
endforeach()

# This also fixes a bug in the linker flags
string(REPLACE "-L " "-L" ROOT_EXE_LINKER_FLAGS "${ROOT_EXE_LINKER_FLAGS}")
separate_arguments(ROOT_EXE_LINKER_FLAGS)
set_property(TARGET ROOT::Flags_CXX APPEND PROPERTY
    INTERFACE_LINK_LIBRARIES ${ROOT_EXE_LINKER_FLAGS})

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

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
