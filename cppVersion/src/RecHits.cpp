//
//  RecHits.cpp
//
//  Created by Jeremi Niedziela on 07/06/2018.
//

#include "RecHits.hpp"
#include "Helpers.hpp"

#include <iostream>

using namespace std;

RecHits::RecHits():
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
  recHitCalib = new RecHitCalibration();
}

RecHits::RecHits(TTree *_tree):
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
{
  recHitCalib = new RecHitCalibration();
  
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

RecHits::~RecHits()
{
  if(recHitCalib) delete recHitCalib;
  
  if(eta){ eta->clear(); delete eta;}
  if(phi){ phi->clear(); delete phi;}
  if(energy){ energy->clear(); delete energy;}
  if(x){ x->clear(); delete x;}
  if(y){ y->clear(); delete y;}
  if(z){ z->clear(); delete z;}
  if(layer){ layer->clear(); delete layer;}
  if(detid){ detid->clear(); delete detid;}
  if(thickness){ thickness->clear(); delete thickness;}
  if(isHalf){ isHalf->clear(); delete isHalf;}
  if(time){ time->clear(); delete time;}
  if(cluster2d){ cluster2d->clear(); delete cluster2d;}
}

unique_ptr<RecHit> RecHits::GetHit(int index)
{
  unique_ptr<RecHit> hit(new RecHit(eta->at(index),phi->at(index),energy->at(index),x->at(index),y->at(index),z->at(index), layer->at(index),detid->at(index),thickness->at(index),isHalf->at(index),time->at(index),cluster2d->at(index)));
  
  return hit;
}

void RecHits::AddHit(unique_ptr<RecHit> &hit)
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

double RecHits::GetTotalEnergy()
{
  double totalEnergy=0;
  for(double e : (*energy)){
    totalEnergy+=e;
  }
  return totalEnergy;
}

double RecHits::GetXmin()
{
  return *min_element(x->begin(), x->end());
}
double RecHits::GetXmax()
{
  return *max_element(x->begin(), x->end());
}
double RecHits::GetYmin()
{
  return *min_element(y->begin(), y->end());
}
double RecHits::GetYmax()
{
  return *max_element(y->begin(), y->end());
}

unique_ptr<RecHits> RecHits::GetHitsAboveNoise(double ecut)
{
  unique_ptr<RecHits> hitsAboveNoise(new RecHits());
  unique_ptr<RecHit> hit;
  
  for(int iHit=0;iHit<N();iHit++){
    hit = GetHit(iHit);
    if(get<0>(hit->RecHitAboveThreshold(recHitCalib, ecut, true))){
      hitsAboveNoise->AddHit(hit);
    }
  }
  return hitsAboveNoise;
}

void RecHits::GetHitsPerSimCluster(vector<RecHits*> &hitsPerCluster, SimClusters *clusters, double energyMin)
{
  unique_ptr<RecHits> hitsAboveNoise = GetHitsAboveNoise(energyMin);
  vector<unsigned int> *hitsDetIDs = hitsAboveNoise->detid;
  
  for(int iCluster=0;iCluster<clusters->N();iCluster++){
    vector<unsigned int> hitsInClusterDetIDs = clusters->GetHits()->at(iCluster);
    vector<unsigned int> clusterToHitID;
    
    sort(hitsDetIDs->begin(), hitsDetIDs->end());
    sort(hitsInClusterDetIDs.begin(), hitsInClusterDetIDs.end());
    
    set_intersection(hitsDetIDs->begin(), hitsDetIDs->end(),
                     hitsInClusterDetIDs.begin(), hitsInClusterDetIDs.end(),
                     std::back_inserter(clusterToHitID));
    
    RecHits *hitsInThisCluster = new RecHits();
    
    for(unsigned int i : clusterToHitID){
      ptrdiff_t pos = distance(detid->begin(), find(detid->begin(), detid->end(), i));
      unique_ptr<RecHit> hit = GetHit((int)pos);
      hitsInThisCluster->AddHit(hit);
    }
    hitsPerCluster.push_back(hitsInThisCluster);
  }
}

void RecHits::GetRecHitsPerHexel(vector<RecHits*> &hitsClustered,vector<shared_ptr<Hexel>> &hexels, double energyMin)
{
  unique_ptr<RecHits> hitsAboveNoise = GetHitsAboveNoise(energyMin);
  
  vector<int> clusterIndices;
  vector<unsigned int> *hitDetIDs = hitsAboveNoise->detid; //ok
  vector<int> hexelDetIDs;
  
  for(auto hexel : hexels){
    if(find(clusterIndices.begin(), clusterIndices.end(), hexel->clusterIndex) == clusterIndices.end()){
      clusterIndices.push_back(hexel->clusterIndex);
    }
    hexelDetIDs.push_back(hexel->detid);//ok
  }
  
  vector<int> hexelToHitID;
  vector<int> new_vec_1 = vector<int>(hexelDetIDs);
  std::sort(std::begin(new_vec_1), std::end(new_vec_1));
  
  for (int i = 0; i < hitDetIDs->size(); ++i) {
    if (std::binary_search(std::begin(new_vec_1),
                           std::end(new_vec_1),
                           hitDetIDs->at(i))) {
      hexelToHitID.push_back(i);
    }
  }
  
  for(int clusterIndex=0;
      clusterIndex<= *max_element(clusterIndices.begin(),clusterIndices.end());
      clusterIndex++){
    hitsClustered.push_back(new RecHits());
  }
  
  for(int iHex=0;iHex<hexels.size();iHex++){
    int hitIndex = hexelToHitID[iHex];
    unique_ptr<RecHit> hit = GetHit(hitIndex);
    hitsClustered[hexels[iHex]->clusterIndex]->AddHit(hit);
  }
}

unique_ptr<RecHits> RecHits::GetHitsInLayer(int layer)
{
  unique_ptr<RecHits> hitsInLayer(new RecHits());
  
  for(int iHit=0;iHit<N();iHit++){
    unique_ptr<RecHit> hit = GetHit(iHit);
    if(hit->layer == layer){
      hitsInLayer->AddHit(hit);
    }
  }
  return hitsInLayer;
}

//------------------------------------------------------------------------------------------------
// Single RecHit

RecHit::RecHit():
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

RecHit::RecHit(float _eta, float _phi, float _energy, float _x, float _y, float _z, int _layer, unsigned int _detid, float _thickness, bool _isHalf, float _time, int _cluster2d):
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

RecHit::~RecHit()
{
  
}

unique_ptr<Hexel> RecHit::GetHexel()
{
  unique_ptr<Hexel> hexel(new Hexel(eta,phi,x,y,z,energy,thickness,time,detid,layer,cluster2d,isHalf));
  return hexel;
}

tuple<bool,double> RecHit::RecHitAboveThreshold(RecHitCalibration *calib, double ecut, double dependSensor)
{
  double sigmaNoise = 1.;
  
  if(dependSensor){
    int thickIndex = -1;
    
    if(layer <= lastLayerFH){  // EE + FH
      if(thickness > 99. and thickness < 101.)        thickIndex = 0;
      else if(thickness > 199. and thickness < 201.)  thickIndex = 1;
      else if(thickness > 299. and thickness < 301.)  thickIndex = 2;
      else cout<<"ERROR - silicon thickness has a nonsensical value"<<endl;
      // determine noise for each sensor/subdetector using RecHitCalibration library
    }
    sigmaNoise = 0.001 * calib->sigmaNoiseMeV(layer, thickIndex);  // returns threshold for EE, FH, BH (in case of BH thickIndex does not play a role)
  }
  bool aboveThreshold = energy >= ecut * sigmaNoise;  // this checks if rechit energy is above the threshold of ecut (times the sigma noise for the sensor, if that option is set)
  return make_tuple(aboveThreshold,sigmaNoise);
}





