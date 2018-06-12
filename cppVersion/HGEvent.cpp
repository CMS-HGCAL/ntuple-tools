//
//  HGCalEvent.cpp
//  xHGCalClustering
//
//  Created by Jeremi Niedziela on 07/06/2018.
//  Copyright © 2018 Jeremi Niedziela. All rights reserved.
//

#include "HGEvent.hpp"

HGEvent::HGEvent(TTree* _tree) : tree(_tree)
{
  genParticles = new HGGenParticles(_tree);
  recHits = new HGRecHits(_tree);
  simClusters = new HGSimClusters(_tree);
}

HGEvent::~HGEvent()
{
// tree should be deleted here
}
