//
//  HGCalEvent.hpp
//
//  Created by Jeremi Niedziela on 07/06/2018.
//

#ifndef HGCalEvent_hpp
#define HGCalEvent_hpp

#include "HGGenParticles.hpp"
#include "HGRecHits.hpp"
#include "HGSimClusters.hpp"

#include <TTree.h>

class HGEvent {
public:
  HGEvent(TTree* _tree);
  ~HGEvent();
  
  inline void GoToEvent(int event){tree->GetEvent(event);}
  
  HGGenParticles *genParticles;
  HGRecHits *recHits;
  HGSimClusters *simClusters;
private:
  TTree *tree; // do I need to keep this tree here?
  
};

#endif /* HGCalEvent_hpp */
