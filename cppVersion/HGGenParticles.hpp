//
//  HGGenParticle.hpp
//
//  Created by Jeremi Niedziela on 07/06/2018.
//

#ifndef HGGenParticle_hpp
#define HGGenParticle_hpp

#include <TTree.h>

#include <vector>

class HGGenParticles {
public:
  HGGenParticles(TTree *_tree);
  ~HGGenParticles();
  
  std::vector<int> *reachedEE;
  
private:
//  TTree *tree;
  
};

#endif /* HGGenParticle_hpp */
