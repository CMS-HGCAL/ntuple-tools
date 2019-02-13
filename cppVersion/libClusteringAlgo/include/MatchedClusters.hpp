//
//  MatchedClusters.hpp
//
//  Created by Jeremi Niedziela on 27/08/2018.
//

#ifndef MatchedClusters_hpp
#define MatchedClusters_hpp

#include "BasicCluster.hpp"
#include "RecHits.hpp"
#include "Helpers.hpp"

#include <vector>

/// Class containing vector of rec cluster and a vector of sim clusters matched together
class MatchedClusters {
  
public:
  /// Default constructor
  MatchedClusters();
  
  /// Default destructor
  ~MatchedClusters();
  
  void Print();
  
  void AddRecCluster(int index,std::unique_ptr<RecHits> &hits);
  void AddSimCluster(int index,std::unique_ptr<RecHits> &hits);
  
  void Merge(std::shared_ptr<MatchedClusters> &clusters);
  
  std::shared_ptr<BasicCluster> GetRecClusterByIndex(int index);
  std::shared_ptr<BasicCluster> GetSimClusterByIndex(int index);
  
  double GetSharedFraction();
  double GetSharedFractionWithRecHits(std::vector<unsigned int> &detIDs);
  inline std::vector<unsigned int> GetRecDetIDs(){return recDetIDs;}
  inline std::vector<unsigned int> GetSimDetIDs(){return simDetIDs;}
  
  bool    HasRecClusters();
  double  GetRecX();
  double  GetRecY();
  double  GetRecRadius();
  double  GetRecEta();
  double  GetRecEnergy();
  
  inline int GetFirstRecIndex(){return recIndices[0];}
  inline void GetRecHits(std::unique_ptr<RecHits> &hits){hits->AddHits(recHits);}
  inline void GetSimHits(std::unique_ptr<RecHits> &hits){hits->AddHits(simHits);}
  
  bool    HasSimClusters();
  double  GetSimX();
  double  GetSimY();
  double  GetSimRadius();
  double  GetSimEta();
  double  GetSimEnergy();
  
  inline int GetNsimClusters(){return (int)simClusters->size();}
  inline int GetNrecClusters(){return (int)recClusters->size();}
  
  bool ContainsSimCluster(int simClusterIndex);
  
private:
  std::unique_ptr<std::vector<std::shared_ptr<BasicCluster>>> simClusters;  ///< Vector of simulated clusters
  std::vector<int> simIndices;
  std::vector<unsigned int> simDetIDs;
  std::unique_ptr<RecHits> simHits;
  
  std::unique_ptr<std::vector<std::shared_ptr<BasicCluster>>> recClusters;  ///< Vector of reconstructed clusters
  std::vector<int> recIndices;
  std::vector<unsigned int> recDetIDs;
  std::unique_ptr<RecHits> recHits;
  
  std::shared_ptr<BasicCluster> GetBasicClusterFromRecHits(std::unique_ptr<RecHits> &hits);
};

#endif /* MatchedClusters_hpp */
