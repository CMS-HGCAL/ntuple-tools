//
//  HGCalEvent.cpps
//
//  Created by Jeremi Niedziela on 07/06/2018.
//

#include "Event.hpp"

using namespace std;

Event::Event(TTree* _tree) : tree(_tree)
{
  genParticles = shared_ptr<GenParticles>(new GenParticles(_tree));
  recHits = shared_ptr<RecHits>(new RecHits(_tree));
  simClusters = shared_ptr<SimClusters>(new SimClusters(_tree));
  clusters2D = shared_ptr<Clusters2D>(new Clusters2D(_tree));
}

Event::~Event()
{
}
