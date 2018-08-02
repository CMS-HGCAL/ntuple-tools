//
//  Clusters2D.cpp
//
//  Created by Jeremi Niedziela on 04/07/2018.
//

#include "Clusters2D.hpp"

using namespace std;

Clusters2D::Clusters2D():
eta(new vector<float>()),
phi(new vector<float>()),
energy(new vector<float>()),
layer(new vector<int>())
{
  
}


Clusters2D::Clusters2D(TTree *_tree):
eta(nullptr),
phi(nullptr),
energy(nullptr),
layer(nullptr)
{
  _tree->SetBranchAddress("cluster2d_eta",&eta);
  _tree->SetBranchAddress("cluster2d_phi",&phi);
  _tree->SetBranchAddress("cluster2d_energy",&energy);
  _tree->SetBranchAddress("cluster2d_layer",&layer);
}

Clusters2D::~Clusters2D()
{
  if(eta){ eta->clear(); delete eta;}
  if(phi){ phi->clear(); delete phi;}
  if(energy){ energy->clear(); delete energy;}
  if(layer){ layer->clear(); delete layer;}
}



void Clusters2D::GetClustersInLayer(std::unique_ptr<Clusters2D> &clustersInLayer, int _layer)
{
  for(int i=0;i<N();i++){
    if(layer->at(i) == _layer){
      clustersInLayer->AddCluster(eta->at(i), phi->at(i), energy->at(i), layer->at(i));
    }
  }
}

void Clusters2D::AddCluster(float _eta, float _phi, float _energy, int _layer)
{
  eta->push_back(_eta);
  phi->push_back(_phi);
  energy->push_back(_energy);
  layer->push_back(_layer);
}
