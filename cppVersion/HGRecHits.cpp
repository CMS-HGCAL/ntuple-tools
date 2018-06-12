//
//  HGRecHits.cpp
//
//  Created by Jeremi Niedziela on 07/06/2018.
//

#include "HGRecHits.hpp"
#include <iostream>

using namespace std;

HGRecHits::HGRecHits():
eta(new vector<float>),
phi(new vector<float>),
energy(new vector<float>),
x(new vector<float>),
y(new vector<float>),
z(new vector<float>),
layer(new vector<int>),
detid(new vector<unsigned int>),
thickness(new vector<float>),
isHalf(new vector<bool>),
time(new vector<float>),
cluster2d(new vector<int>)
{
  
}

HGRecHits::HGRecHits(TTree *_tree):
eta(nullptr),
phi(nullptr),
energy(nullptr),
x(nullptr),
y(nullptr),
z(nullptr),
layer(nullptr),
detid(nullptr),
thickness(nullptr),
isHalf(nullptr),
time(nullptr),
cluster2d(nullptr)
//tree(_tree)
{
//  cout<<"tree:"<<_tree<<endl;
//  tree = _tree;
  _tree->SetBranchAddress("rechit_eta",&eta);
  _tree->SetBranchAddress("rechit_phi",&phi);
  _tree->SetBranchAddress("rechit_energy",&energy);
  _tree->SetBranchAddress("rechit_x",&x);
  _tree->SetBranchAddress("rechit_y",&y);
  _tree->SetBranchAddress("rechit_z",&z);
  _tree->SetBranchAddress("rechit_layer",&layer);
  _tree->SetBranchAddress("rechit_detid",&detid);
  _tree->SetBranchAddress("rechit_thickness",&thickness);
  _tree->SetBranchAddress("rechit_isHalf",&isHalf);
  _tree->SetBranchAddress("rechit_time",&time);
  _tree->SetBranchAddress("rechit_cluster2d",&cluster2d);
}

HGRecHits::~HGRecHits()
{
  
}

HGRecHit* HGRecHits::GetHit(int index)
{
  vector<float> &etaRef = *eta;
  vector<float> &phiRef = *phi;
  vector<float> &energyRef = *energy;
  vector<float> &xRef = *x;
  vector<float> &yRef = *y;
  vector<float> &zRef = *z;
  vector<int> &layerRef = *layer;
  vector<unsigned int> &detidRef = *detid;
  vector<float> &thicknessRef = *thickness;
  vector<bool> &isHalfRef = *isHalf;
  vector<float> &timeRef = *time;
  vector<int> &cluster2dRef = *cluster2d;
  
  return new HGRecHit(etaRef[index],phiRef[index],energyRef[index],xRef[index],yRef[index],zRef[index], layerRef[index],detidRef[index],thicknessRef[index],isHalfRef[index],timeRef[index],cluster2dRef[index]);
}

void HGRecHits::AddHit(HGRecHit *hit)
{
  eta->push_back(hit->eta);
  phi->push_back(hit->phi);
  energy->push_back(hit->energy);
  x->push_back(hit->x);
  y->push_back(hit->y);
  z->push_back(hit->z);
  layer->push_back(hit->layer);
  detid->push_back(hit->detid);
  thickness->push_back(hit->thickness);
  isHalf->push_back(hit->isHalf);
  time->push_back(hit->time);
  cluster2d->push_back(hit->cluster2d);
}

double HGRecHits::GetTotalEnergy()
{
  double totalEnergy=0;
  for(double e : (*energy)){
    totalEnergy+=e;
  }
  return totalEnergy;
}

double HGRecHits::GetXmin()
{
  return *min_element(x->begin(), x->end());
}
double HGRecHits::GetXmax()
{
  return *max_element(x->begin(), x->end());
}
double HGRecHits::GetYmin()
{
  return *min_element(y->begin(), y->end());
}
double HGRecHits::GetYmax()
{
  return *max_element(y->begin(), y->end());
}

//------------------------------------------------------------------------------------------------
// Single HGRecHit

HGRecHit::HGRecHit():
eta(0),
phi(0),
energy(0),
x(0),
y(0),
z(0),
layer(0),
detid(0),
thickness(0),
isHalf(false),
time(0),
cluster2d(0)
{
  
}

HGRecHit::HGRecHit(float _eta, float _phi, float _energy, float _x, float _y, float _z, int _layer, unsigned int _detid, float _thickness, bool _isHalf, float _time, int _cluster2d):
eta(_eta),
phi(_phi),
energy(_energy),
x(_x),
y(_y),
z(_z),
layer(_layer),
detid(_detid),
thickness(_thickness),
isHalf(_isHalf),
time(_time),
cluster2d(_cluster2d)
{
  
}

HGRecHit::~HGRecHit()
{
  
}







