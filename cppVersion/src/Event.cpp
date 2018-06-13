//
//  HGCalEvent.cpp
//  xHGCalClustering
//
//  Created by Jeremi Niedziela on 07/06/2018.
//  Copyright © 2018 Jeremi Niedziela. All rights reserved.
//

#include "Event.hpp"

Event::Event(TTree* _tree) : tree(_tree)
{
  genParticles = new GenParticles(_tree);
  recHits = new RecHits(_tree);
  simClusters = new SimClusters(_tree);
}

Event::~Event()
{
  delete genParticles;
  delete recHits;
  delete simClusters;
}
