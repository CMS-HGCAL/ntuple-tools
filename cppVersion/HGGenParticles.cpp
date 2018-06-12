//
//  HGGenParticle.cpp
//  xHGCalClustering
//
//  Created by Jeremi Niedziela on 07/06/2018.
//  Copyright © 2018 Jeremi Niedziela. All rights reserved.
//

#include "HGGenParticles.hpp"

HGGenParticles::HGGenParticles(TTree *_tree):reachedEE(nullptr)/*, tree(_tree)*/
{
  _tree->SetBranchAddress("genpart_reachedEE",&reachedEE);
  
}

HGGenParticles::~HGGenParticles()
{
  
}
