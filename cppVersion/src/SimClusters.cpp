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
//tree(_tree)
{
  _tree->SetBranchAddress("simcluster_eta",&eta);
  _tree->SetBranchAddress("simcluster_phi",&phi);
  _tree->SetBranchAddress("simcluster_energy",&energy);
  _tree->SetBranchAddress("simcluster_hits",&hits);
}

SimClusters::~SimClusters()
{
  
}
