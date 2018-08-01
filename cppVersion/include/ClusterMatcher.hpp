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
    recClusters = new std::vector<BasicCluster*>;
    simClusters = new std::vector<BasicCluster*>;
    recHits = std::unique_ptr<RecHits>(new RecHits());
    simHits = std::unique_ptr<RecHits>(new RecHits());
  }
  
  std::vector<BasicCluster*> *recClusters;
  std::vector<BasicCluster*> *simClusters;
  
  std::vector<int> recIndices;
  std::vector<int> simIndices;
  
  std::vector<unsigned int> recDetIDs;
  std::vector<unsigned int> simDetIDs;
  
  std::vector<double> recEnergies;
  std::vector<double> simEnergies;
  
  std::unique_ptr<RecHits> recHits;
  std::unique_ptr<RecHits> simHits;
  
  inline BasicCluster* GetBasicClusterFromRecHits(std::unique_ptr<RecHits> &hits)
  {
    double recEnergy = hits->GetTotalEnergy();
    double xMax   = hits->GetXmax();
    double xMin   = hits->GetXmin();
    double yMax   = hits->GetYmax();
    double yMin   = hits->GetYmin();
    
    double clusterX = (xMax+xMin)/2.;
    double clusterY = (yMax+yMin)/2.;
    double clusterEta = hits->GetCenterEta();
    double clusterR = std::max((xMax-xMin)/2.,(yMax-yMin)/2.);
    
    BasicCluster *basicCluster = new BasicCluster(recEnergy,clusterX,clusterY,0,clusterEta,clusterR);
    return basicCluster;
  }
  
  inline BasicCluster* GetMergedRecCluster(){
    return GetBasicClusterFromRecHits(recHits);
  }
  
  inline BasicCluster* GetMergedSimCluster(){
    return GetBasicClusterFromRecHits(simHits);
  }
  
  inline void AddRecDetIDs(std::vector<unsigned int> &IDs){
    recDetIDs.insert(recDetIDs.end(), IDs.begin(), IDs.end());
  }
  
  inline void AddSimDetIDs(std::vector<unsigned int> &IDs){
    simDetIDs.insert(simDetIDs.end(), IDs.begin(), IDs.end());
  }
  
  inline void AddRecEnergies(std::vector<double> &e){
    recEnergies.insert(recEnergies.end(), e.begin(), e.end());
  }
  
  inline void AddSimEnergies(std::vector<double> &e){
    simEnergies.insert(simEnergies.end(), e.begin(), e.end());
  }
  
  inline double GetSharedFraction(){
    std::vector<unsigned int> common;
    std::set_intersection(simDetIDs.begin(), simDetIDs.end(),
                          recDetIDs.begin(), recDetIDs.end(),
                          std::back_inserter(common));
    
    double energySumShared=0;
    double energySumTotal=0;
    
    for(int i=0;i<common.size();i++){
      auto indx = std::find(simDetIDs.begin(),simDetIDs.end(), common[i]);
      energySumShared += simEnergies[indx-simDetIDs.begin()];
    }
    for(int i=0;i<simEnergies.size();i++){
      energySumTotal += simEnergies[i];
    }
    return energySumShared/energySumTotal;
  }
  
  inline double GetTotalRecEnergy(){
    double energy = 0;
    for(BasicCluster *cluster : *recClusters){
      energy += cluster->GetEnergy();
    }
    return energy;
  }
  
  inline double GetTotalSimEnergy(){
    double energy = 0;
    for(BasicCluster *cluster : *simClusters){
      energy += cluster->GetEnergy();
    }
    return energy;
  }
};

#endif /* ClusterMatcher_hpp */
