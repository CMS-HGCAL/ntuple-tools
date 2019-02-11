//
//  ClusterMatcher.hpp
//
//  Created by Jeremi Niedziela on 05/07/2018.
//

#ifndef ClusterMatcher_hpp
#define ClusterMatcher_hpp

#include "BasicCluster.hpp"
#include "RecHits.hpp"
#include "Helpers.hpp"
#include "MatchedClusters.hpp"

class ClusterMatcher {
public:
  ClusterMatcher();
  ~ClusterMatcher();
  
  /// Fills a vector of matched rec and sim clusters finding those that share the most hits (by detID)
  /// \param matched Vector that will be filled with rec and sim clusters
  /// \param recHitsPerCluster Vector of rec clusters in given layer
  /// \param simHitsPerCluster Vector or sim clusters in given layer
  /// \param draw Should matched clusters in this layer be plotted
  void MatchClustersByDetID(std::vector<std::shared_ptr<MatchedClusters>> &matched,
                            std::vector<std::unique_ptr<RecHits>> &recHitsPerCluster,
                            std::vector<std::unique_ptr<RecHits>> &simHitsPerCluster,
                            bool draw=false);
  
  /// Fills a vector of matched rec and sim clusters finding the nearest rec cluster for given sim cluster
  /// \param matched Vector that will be filled with rec and sim clusters
  /// \param recHitsPerCluster Vector of rec clusters in given layer
  /// \param simHitsPerCluster Vector or sim clusters in given layer
  void MatchClustersClosest(std::vector<std::shared_ptr<MatchedClusters>> &matched,
                            std::vector<std::unique_ptr<RecHits>> &recHitsPerCluster,
                            std::vector<std::unique_ptr<RecHits>> &simHitsPerCluster);
  
  /// Fills a vector of unmatched rec and sim clusters (simply assigns all sim clusters to all rec each rec cluster)
  /// \param matched Vector that will be filled with rec and sim clusters
  /// \param recHitsPerCluster Vector of rec clusters in given layer
  /// \param simHitsPerCluster Vector or sim clusters in given layer
  void MatchClustersAllToAll(std::vector<std::shared_ptr<MatchedClusters>> &matched,
                             std::vector<std::unique_ptr<RecHits>> &recHitsPerCluster,
                             std::vector<std::unique_ptr<RecHits>> &simHitsPerCluster);
  
  /// Plots sim hits, rec clusters and rec hits that belong to them from a vactor of matched clusters.
  /// If there are two sim hits in the same position, the second one will be drawn as opened circle on top
  /// of a closed one.
  /// \param matched Vector of matched clusters to be drawn
  void DrawMatched(std::vector<std::shared_ptr<MatchedClusters>> &matched);
};


#endif /* ClusterMatcher_hpp */
