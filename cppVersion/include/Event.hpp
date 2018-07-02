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

/// Event class containing pointers to generated particles, sim clusters and rec hits
class Event {
public:
  /// Default constructor, requires a tree from which particles, clusters and hits will be read
  Event(TTree* _tree);
  ~Event();
  
  /// Jump to event
  inline void GoToEvent(int event){
    tree->GetEvent(event);
  }
  
  /// Returns collection of generated particles
  inline std::shared_ptr<GenParticles> GetGenParticles(){return genParticles;}
  
  /// Returns collection of rec hits
  inline std::shared_ptr<RecHits> GetRecHits(){return recHits;}
  
  /// Returns collection of sim clusters
  inline std::shared_ptr<SimClusters> GetSimClusters(){return simClusters;}
  
private:
  TTree *tree;  ///< Pointer to tree containing HGCal events
  
  std::shared_ptr<GenParticles> genParticles; ///< Collection of generated particles
  std::shared_ptr<RecHits> recHits;           ///< Collection of rec hits
  std::shared_ptr<SimClusters> simClusters;   ///< Collection of sim clusters
};

#endif /* HGCalEvent_hpp */
