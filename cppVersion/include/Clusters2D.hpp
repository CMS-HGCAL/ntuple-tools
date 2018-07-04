//
//  Clusters2D.hpp
//
//  Created by Jeremi Niedziela on 04/07/2018.
//

#ifndef Clusters2D_hpp
#define Clusters2D_hpp

#include <TTree.h>

#include <vector>

/// Collection of generated particles
class Clusters2D{
public:
  /// Default constructor
  Clusters2D();
  
  /// Default constructor
  /// \param _tree Input tree from which generated partiles will be read
  Clusters2D(TTree *_tree);
  ~Clusters2D();
  
  /// Returns number of stored RecHit objects
  inline int N(){return (int)energy->size();}
  
  /// Returns pseudorapidity of i-th cluster
  inline float GetEta(int i){return eta->at(i);}
  
  /// Returns polar angle of i-th cluster
  inline float GetPhi(int i){return phi->at(i);}
  
  /// Returns energy of i-th cluster
  inline float GetEnergy(int i){return energy->at(i);}
  
  /// Returns layer index of i-th cluster
  inline int GetLayer(int i){return layer->at(i);}
  
  /// Returns a set of clusters laying in given layer
  void GetClustersInLayer(std::unique_ptr<Clusters2D> &clustersInLayer, int _layer);
  
private:
  std::vector<float> *eta;        ///< Pseudorapidity values of hit
  std::vector<float> *phi;        ///< Polar angle values of hit
  std::vector<float> *energy;     ///< Energy values of hit (GeV)
  std::vector<int> *layer;      ///< Layer index
  
  void AddCluster(float _eta, float _phi, float _energy, int _layer);
};

#endif /* Clusters2D_hpp */
