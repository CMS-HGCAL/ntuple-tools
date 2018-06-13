//
//  RecHits.cpp
//
//  Created by Jeremi Niedziela on 07/06/2018.
//

#include "RecHits.hpp"

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
  
}

RecHit* RecHits::GetHit(int index)
{
  return new RecHit(eta->at(index),phi->at(index),energy->at(index),x->at(index),y->at(index),z->at(index), layer->at(index),detid->at(index),thickness->at(index),isHalf->at(index),time->at(index),cluster2d->at(index));
}

void RecHits::AddHit(RecHit *hit)
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

RecHits* RecHits::GetHitsAboveNoise(double ecut)
{
  double sigmaNoise = 1.;
  double thickIndex = -1;
  
  RecHits *hitsAboveNoise = new RecHits();
  RecHit *hit;
  
  for(int iHit=0;iHit<N();iHit++){
    hit = GetHit(iHit);
    int layer = hit->layer;
    float thickness = hit->thickness;
    
    if(layer <= 40){  // EE + FH
      if     (thickness > 99.  && thickness < 101.) thickIndex = 0;
      else if(thickness > 199. && thickness < 201.) thickIndex = 1;
      else if(thickness > 299. && thickness < 301.) thickIndex = 2;
      else cout<<"ERROR - silicon thickness has a nonsensical value"<<endl;
    }
    // determine noise for each sensor/subdetector using RecHitCalibration library
    sigmaNoise = 0.001 * recHitCalib->sigmaNoiseMeV(layer, thickIndex);  // returns threshold for EE, FH, BH (in case of BH thickIndex does not play a role)
    if(hit->energy >= ecut * sigmaNoise){
      // checks if rechit energy is above the threshold of ecut (times the sigma noise for the sensor, if that option is set)
      hitsAboveNoise->AddHit(hit);
    }
    delete hit;
  }
  return hitsAboveNoise;
}

// groups hits into array of clusters
void RecHits::GetHitsPerCluster(vector<RecHits*> &hitsPerCluster, SimClusters *clusters, double energyMin)
{
  RecHits *hitsAboveNoise = GetHitsAboveNoise(energyMin);
  vector<unsigned int> *hitsDetIDs = hitsAboveNoise->detid;
  
  for(int iCluster=0;iCluster<clusters->N();iCluster++){
    vector<unsigned int> hitsInClusterDetIDs = clusters->hits->at(iCluster);
    vector<unsigned int> clusterToHitID;
    
    sort(hitsDetIDs->begin(), hitsDetIDs->end());
    sort(hitsInClusterDetIDs.begin(), hitsInClusterDetIDs.end());
    
    set_intersection(hitsDetIDs->begin(), hitsDetIDs->end(),
                     hitsInClusterDetIDs.begin(), hitsInClusterDetIDs.end(),
                     std::back_inserter(clusterToHitID));
    
    RecHits *hitsInThisCluster = new RecHits();
    
    for(unsigned int i : clusterToHitID){
      ptrdiff_t pos = distance(detid->begin(), find(detid->begin(), detid->end(), i));
      hitsInThisCluster->AddHit(GetHit((int)pos));
    }
    hitsPerCluster.push_back(hitsInThisCluster);
  }
}

// groups hits associated with hexels into array of clusters
void RecHits::GetRecHitsPerHexel(vector<RecHits*> &hitsClustered,vector<Hexel*> hexels, double energyMin)
{
  RecHits *hitsAboveNoise = GetHitsAboveNoise(energyMin);
  
  vector<int> clusterIndices;
  vector<unsigned int> *hitDetIDs = hitsAboveNoise->detid; //ok
  vector<int> hexelDetIDs;
  
  for(Hexel *hexel : hexels){
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
    RecHit *hit = GetHit(hitIndex);
    hitsClustered[hexels[iHex]->clusterIndex]->AddHit(hit);
  }
}

RecHits* RecHits::GetHitsInLayer(int layer)
{
  RecHits *hitsInLayer = new RecHits();
  RecHit *hit = nullptr;
  
  for(int iHit=0;iHit<N();iHit++){
    hit = GetHit(iHit);
    if(!hit) continue;
    if(hit->layer == layer){
      hitsInLayer->AddHit(hit);
    }
    delete hit;
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

Hexel* RecHit::GetHexel()
{
  return new Hexel(eta,phi,x,y,z,energy,thickness,time,detid,layer,cluster2d,isHalf);
}






