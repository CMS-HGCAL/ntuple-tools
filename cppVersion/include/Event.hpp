//
//  HGCalEvent.hpp
//
//  Created by Jeremi Niedziela on 07/06/2018.
//

#ifndef HGCalEvent_hpp
#define HGCalEvent_hpp

#include "GenParticles.hpp"
#include "RecHits.hpp"
#include "SimClusters.hpp"

#include <TTree.h>

class Event {
public:
  Event(TTree* _tree);
  ~Event();
  
  inline void GoToEvent(int event){tree->GetEvent(event);}
  
  GenParticles *genParticles;
  RecHits *recHits;
  SimClusters *simClusters;
private:
  TTree *tree;
  
};

#endif /* HGCalEvent_hpp */
