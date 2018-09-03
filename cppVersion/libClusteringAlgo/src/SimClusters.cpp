//
//  HGSimClusters.cpp
//
//  Created by Jeremi Niedziela on 07/06/2018.
//

#include "SimClusters.hpp"

#include <iostream>

using namespace std;

SimClusters::SimClusters(TTree *_tree):
eta(nullptr),
phi(nullptr),
energy(nullptr),
pt(nullptr),
hits(nullptr),
layers(nullptr),
fractions(nullptr)
{
  _tree->SetBranchAddress("simcluster_eta",&eta);
  _tree->SetBranchAddress("simcluster_phi",&phi);
  _tree->SetBranchAddress("simcluster_pt",&pt);
  _tree->SetBranchAddress("simcluster_energy",&energy);
  _tree->SetBranchAddress("simcluster_hits",&hits);
  _tree->SetBranchAddress("simcluster_layers",&layers);
  _tree->SetBranchAddress("simcluster_fractions",&fractions);
  // ... more branches can be add in the future if needed
}

SimClusters::~SimClusters()
{
  if(eta){ eta->clear(); delete eta;}
  if(phi){ phi->clear(); delete phi;}
  if(energy){ energy->clear(); delete energy;}
  if(hits){ hits->clear(); delete hits;}
  if(layers){ layers->clear(); delete layers;}
  if(fractions){ fractions->clear(); delete fractions;}
}

void SimClusters::Print(int index)
{
  cout<<"Sim-cluster "<<index<<":";
  cout<<"\t p_T:"<<pt->at(index);
  cout<<"\t E:"<<energy->at(index);
  cout<<"\t eta:"<<eta->at(index);
  cout<<"\t phi:"<<phi->at(index);
  cout<<endl;
}

void SimClusters::Clean()
{
  if(eta){ eta->clear(); delete eta; eta = nullptr;}
  if(phi){ phi->clear(); delete phi; phi = nullptr;}
  if(energy){ energy->clear(); delete energy; energy = nullptr;}
  if(hits){ hits->clear(); delete hits; hits = nullptr;}
  if(layers){ layers->clear(); delete layers; layers = nullptr;}
  if(fractions){ fractions->clear(); delete fractions; fractions = nullptr;}
}

int SimClusters::GetNsimClustersInLayer(unsigned int iLayer)
{
  int nClusters =0;
  
  for(int iCluster=0;iCluster<layers->size();iCluster++){
    std::vector<unsigned int> thisCluster = layers->at(iCluster);
    
    for(int iHit=0;iHit<thisCluster.size();iHit++){
      if(thisCluster[iHit] == iLayer) nClusters++;
    }
    
  }
  return nClusters;
}
