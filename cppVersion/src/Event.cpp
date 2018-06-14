//
//  HGCalEvent.cpps
//
//  Created by Jeremi Niedziela on 07/06/2018.
//

#include "Event.hpp"

using namespace std;

Event::Event(TTree* _tree) : tree(_tree)
{
  genParticles = new GenParticles(_tree);
  recHits = shared_ptr<RecHits>(new RecHits(_tree));
  simClusters = new SimClusters(_tree);
}

Event::~Event()
{
  delete genParticles;
  delete simClusters;
}
