//
//  ClusterMatcher.hpp
//
//  Created by Jeremi Niedziela on 05/07/2018.
//

#ifndef ClusterMatcher_hpp
#define ClusterMatcher_hpp

#include "BasicCluster.hpp"
#include "RecHits.hpp"

struct MatchedClusters;

class ClusterMatcher {
public:
  ClusterMatcher();
  ~ClusterMatcher();
  
  void MatchClustersByDetID(std::vector<MatchedClusters*> &matched, std::vector<RecHits*> &recHitsPerCluster, std::vector<RecHits*> &simHitsPerCluster, int layer);
  
  /// Fills a vector of matched rec and sim clusters finding the nearest rec cluster for given sim cluster
  /// \param matched Vector that will be filled with rec and sim clusters
  /// \param recHitsPerCluster Vector of rec clusters
  /// \param simHitsPerCluster Vector or sim clusters
  /// \param layer Layer index
  void MatchClustersClosest(std::vector<MatchedClusters*> &matched, std::vector<RecHits*> &recHitsPerCluster, std::vector<RecHits*> &simHitsPerCluster, int layer);
  
  /// Fills a vector of unmatched rec and sim clusters (simply assigns all sim clusters to all rec each rec cluster)
  /// \param matched Vector that will be filled with rec and sim clusters
  /// \param recHitsPerCluster Vector of rec clusters
  /// \param simHitsPerCluster Vector or sim clusters
  /// \param layer Layer index
  void MatchClustersAllToAll(std::vector<MatchedClusters*> &matched, std::vector<RecHits*> &recHitsPerCluster, std::vector<RecHits*> &simHitsPerCluster, int layer);
  
  
private:
  BasicCluster* GetBasicClusterFromRecHits(std::unique_ptr<RecHits> &hits);
};


/// Struct containing one rec cluster and a vector of sim clusters matched to it
struct MatchedClusters {
  
  MatchedClusters(){
    recCluster = nullptr;
    simClusters = new std::vector<BasicCluster*>;
  }
  
  BasicCluster *recCluster;
  std::vector<BasicCluster*> *simClusters;
  
  inline double GetTotalSimEnergy(){
    double energy = 0;
    for(BasicCluster *cluster : *simClusters){
      energy += cluster->GetEnergy();
    }
    return energy;
  }
};

#endif /* ClusterMatcher_hpp */
