//
//  HGGenParticle.cpp
//
//  Created by Jeremi Niedziela on 07/06/2018.
//

#include "GenParticles.hpp"

GenParticles::GenParticles(TTree *_tree) :
reachedEE(nullptr)
{
  _tree->SetBranchAddress("genpart_reachedEE",&reachedEE);
  //... more branches can be added if needed
}

GenParticles::~GenParticles()
{
  if(reachedEE){ reachedEE->clear(); delete reachedEE;}
}

void GenParticles::Clean()
{
  if(reachedEE){reachedEE->clear(); delete reachedEE; reachedEE = nullptr;}
}
