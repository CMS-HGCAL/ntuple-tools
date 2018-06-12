//
//  HGGenParticle.cpp
//
//  Created by Jeremi Niedziela on 07/06/2018.
//

#include "GenParticles.hpp"

GenParticles::GenParticles(TTree *_tree):reachedEE(nullptr)
{
  _tree->SetBranchAddress("genpart_reachedEE",&reachedEE);
}

GenParticles::~GenParticles()
{
  
}
