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
  void MatchClustersByDetID(std::vector<MatchedClusters*> &matched,
                            std::vector<std::unique_ptr<RecHits>> &recHitsPerCluster,
                            std::vector<std::unique_ptr<RecHits>> &simHitsPerCluster);
  
  /// Fills a vector of matched rec and sim clusters finding the nearest rec cluster for given sim cluster
  /// \param matched Vector that will be filled with rec and sim clusters
  /// \param recHitsPerCluster Vector of rec clusters
  /// \param simHitsPerCluster Vector or sim clusters
  /// \param layer Layer index
  void MatchClustersClosest(std::vector<MatchedClusters*> &matched,
                            std::vector<RecHits*> &recHitsPerCluster,
                            std::vector<RecHits*> &simHitsPerCluster,
                            int layer);
  
  /// Fills a vector of unmatched rec and sim clusters (simply assigns all sim clusters to all rec each rec cluster)
  /// \param matched Vector that will be filled with rec and sim clusters
  /// \param recHitsPerCluster Vector of rec clusters
  /// \param simHitsPerCluster Vector or sim clusters
  /// \param layer Layer index
  void MatchClustersAllToAll(std::vector<MatchedClusters*> &matched,
                             std::vector<RecHits*> &recHitsPerCluster,
                             std::vector<RecHits*> &simHitsPerCluster,
                             int layer);
};


#endif /* ClusterMatcher_hpp */
