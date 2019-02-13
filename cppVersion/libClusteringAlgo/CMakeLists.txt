add_library(libClusteringAlgo src/BasicCluster.cpp
                              src/ConfigurationManager.cpp
                              src/Hexel.cpp
                              src/RecHits.cpp
                              src/ClusterMatcher.cpp
                              src/Event.cpp
                              src/ImagingAlgo.cpp
                              src/SimClusters.cpp
                              src/Clusters2D.cpp
                              src/GenParticles.cpp
                              src/MatchedClusters.cpp
                              src/TestbeamTrack.cpp
                              )

set_property(TARGET libClusteringAlgo PROPERTY CXX_STANDARD 17)
