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
  
  /// Prints basic information about gen particle at given index
  void Print(int index);
  
  /// Deletes all objects and leaves collection with size zero
  void Clean();
  
  /// Returns number of clusters
  inline int N(){return (int)eta->size();}
  
  /// Returns energy of i-th cluster
  inline float GetEnergy(int i){return energy->at(i);}
  
  /// Returns pseudorapidity of i-th cluster
  inline float GetEta(int i){return eta->at(i);}
  
  /// Returns polar angle of i-th cluster
  inline float GetPhi(int i){return phi->at(i);}
  
  /// Returns transverse momentum of i-th cluster
  inline float GetPt(int i){return pt->at(i);}
  
  /// Returns detIDs of hits grouped by cluster
  inline std::vector<std::vector<unsigned int>>* GetHits(){return hits;}
  
  /// Returns energy fractions for hits grouped by cluster
  inline std::vector<std::vector<float>>* GetFractions(){return fractions;}
  
  /// Returns layer indices of hits in given cluster
  inline std::vector<unsigned int> GetLayersInCluster(int i){return layers->at(i);}
  
  int GetNsimClustersInLayer(unsigned int iLayer);
private:
  std::vector<float> *eta;    ///< Pseudorapidity values for each of the clusters
  std::vector<float> *phi;    ///< Polar anlges values for each of the clusters
  std::vector<float> *energy; ///< Energy values for each of the clusters (GeV)
  std::vector<float> *pt;     ///< Transverse momentum values for each of the clusters (GeV)
  std::vector<std::vector<unsigned int>> *hits;   ///< Vector of detIDs of hits belonging to each of the clusters
  std::vector<std::vector<unsigned int>> *layers; ///< Vector of detIDs of hits belonging to each of the clusters
  std::vector<std::vector<float>> *fractions;     ///< Vector of energy fractions to be applied to each hit
};

#endif /* HGSimClusters_hpp */
