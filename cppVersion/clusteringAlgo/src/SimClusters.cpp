//
//  HGSimClusters.cpp
//
//  Created by Jeremi Niedziela on 07/06/2018.
//

#include "SimClusters.hpp"

SimClusters::SimClusters(TTree *_tree):
eta(nullptr),
phi(nullptr),
energy(nullptr),
hits(nullptr)
{
  _tree->SetBranchAddress("simcluster_eta",&eta);
  _tree->SetBranchAddress("simcluster_phi",&phi);
  _tree->SetBranchAddress("simcluster_pt",&pt);
  _tree->SetBranchAddress("simcluster_energy",&energy);
  _tree->SetBranchAddress("simcluster_hits",&hits);
  // ... more branches can be add in the future if needed
}

SimClusters::~SimClusters()
{
  if(eta){ eta->clear(); delete eta;}
  if(phi){ phi->clear(); delete phi;}
  if(energy){ energy->clear(); delete energy;}
  if(hits){ hits->clear(); delete hits;}
}
