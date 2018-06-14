//
//  HGSimClusters.hpp
//
//  Created by Jeremi Niedziela on 07/06/2018.
//

#ifndef HGSimClusters_hpp
#define HGSimClusters_hpp

#include <TTree.h>

///  SimCluster keeps basic information about simulated clusters
class SimClusters {
public:
  /// Default constructor
  /// \param _tree Input tree from which sim clusters will be read
  SimClusters(TTree *_tree);
  ~SimClusters();
  
  /// Returns number of clusters
  inline int N(){return (int)eta->size();}
  
  /// Returns detIDs of hits grouped by cluster
  inline std::vector<std::vector<unsigned int>>* GetHits(){return hits;}
  
private:
  std::vector<float> *eta;    ///< Pseudorapidity values for each of the clusters
  std::vector<float> *phi;    ///< Polar anlges values for each of the clusters
  std::vector<float> *energy; ///< Energy values for each of the clusters (GeV)
  std::vector<std::vector<unsigned int>> *hits; /// Vector of detIDs of hits belonging to each of the clusters
  
};

#endif /* HGSimClusters_hpp */
