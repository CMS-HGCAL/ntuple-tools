//
//  HGGenParticle.hpp
//
//  Created by Jeremi Niedziela on 07/06/2018.
//

#ifndef HGGenParticle_hpp
#define HGGenParticle_hpp

#include <TTree.h>

#include <vector>

/// Collection of generated particles
class GenParticles{
public:
  /// Default constructor
  /// \param _tree Input tree from which generated partiles will be read
  GenParticles(TTree *_tree);
  ~GenParticles();
  
  void Clean();
  
  /// Returns pointer to vector telling for each gen particle if it reached the EE
  inline std::vector<int>* GetReachedEE(){return reachedEE;}
  
private:
  std::vector<int> *reachedEE; ///< Info whether or not particle reached EE
};

#endif /* HGGenParticle_hpp */
