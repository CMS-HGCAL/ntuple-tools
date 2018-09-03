//
//  RecHits.cpp
//
//  Created by Jeremi Niedziela on 07/06/2018.
//

#include "RecHits.hpp"
#include "Helpers.hpp"

#include <iostream>
#include <algorithm>

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
  recHitCalib = unique_ptr<RecHitCalibration>(new RecHitCalibration());
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
  recHitCalib = unique_ptr<RecHitCalibration>(new RecHitCalibration());
  
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


void RecHits::Print()
{
  cout<<"Rec hits collection:";
  cout<<"\tcenter eta:"<<GetCenterEta();
  cout<<"\tavg x:"<<(GetXmax()+GetXmin())/2.;
  cout<<"\tavg y:"<<(GetYmax()+GetYmin())/2.;
  cout<<"\tavg E:"<<GetTotalEnergy();
  cout<<endl;
  
  if(ConfigurationManager::Instance()->GetVerbosityLevel() > 1){
    for(int i=0;i<x->size();i++){
      cout<<"x:"<<x->at(i)<<"\ty:"<<y->at(i)<<"\tE:"<<energy->at(i)<<endl;
    }
  }
}

void RecHits::Clean()
{
  if(eta){ eta->clear(); delete eta; eta = new vector<float>;}
  if(phi){ phi->clear(); delete phi; phi = new vector<float>;}
  if(energy){ energy->clear(); delete energy; energy = new vector<float>;}
  if(x){ x->clear(); delete x; x = new vector<float>;}
  if(y){ y->clear(); delete y; y = new vector<float>;}
  if(z){ z->clear(); delete z; z = new vector<float>;}
  if(layer){ layer->clear(); delete layer; layer = new vector<int>;}
  if(detid){ detid->clear(); delete detid; detid = new vector<unsigned int>;}
  if(thickness){ thickness->clear(); delete thickness; thickness = new vector<float>;}
  if(isHalf){ isHalf->clear(); delete isHalf; isHalf = new vector<bool>;}
  if(time){ time->clear(); delete time; time = new vector<float>;}
  if(cluster2d){ cluster2d->clear(); delete cluster2d; cluster2d = new vector<int>;}
}

unique_ptr<RecHit> RecHits::GetHit(int index, double energyFraction)
{
  unique_ptr<RecHit> hit(new RecHit(eta->at(index),phi->at(index),energyFraction*energy->at(index),x->at(index),y->at(index),z->at(index), layer->at(index),detid->at(index),thickness->at(index),isHalf->at(index),time->at(index),cluster2d->at(index)));

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

void RecHits::AddHits(unique_ptr<RecHits> &hits)
{
  for(int i=0;i<hits->N();i++){
    unique_ptr<RecHit> hit = hits->GetHit(i);
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
  if(!x || x->size()==0) return 99999;
  return *min_element(x->begin(), x->end());
}
double RecHits::GetXmax()
{
  if(!x || x->size()==0) return 99999;
  return *max_element(x->begin(), x->end());
}
double RecHits::GetYmin()
{
  if(!y || y->size()==0) return 99999;
  return *min_element(y->begin(), y->end());
}
double RecHits::GetYmax()
{
  if(!y || y->size()==0) return 99999;
  return *max_element(y->begin(), y->end());
}

double RecHits::GetCenterEta()
{
  if(!eta || eta->size()==0) return 99999;
  
  double minEta = *min_element(eta->begin(), eta->end());
  double maxEta = *max_element(eta->begin(), eta->end());
  
  return (maxEta+minEta)/2.;
}

unique_ptr<RecHits> RecHits::GetHitsAboveNoise()
{
  unique_ptr<RecHits> hitsAboveNoise(new RecHits());
  unique_ptr<RecHit> hit;

  for(int iHit=0;iHit<N();iHit++){
    
    if(get<0>(RecHitAboveThreshold(iHit))){
      hit = GetHit(iHit);
      hitsAboveNoise->AddHit(hit);
    }
  }
  return hitsAboveNoise;
}

void RecHits::GetHitsPerSimCluster(vector<unique_ptr<RecHits>> &hitsPerCluster,
                                   shared_ptr<SimClusters> clusters)
{
  unique_ptr<RecHits> hitsAboveNoise = GetHitsAboveNoise();
  vector<unsigned int> *hitsDetIDs = hitsAboveNoise->detid;

  int nAssociatedHits = 0;
  
  for(int iCluster=0;iCluster<clusters->N();iCluster++){
    if(ConfigurationManager::Instance()->GetVerbosityLevel() > 0){
      clusters->Print(iCluster);
    }
    vector<float>        hitsInClusterFractions = clusters->GetFractions()->at(iCluster);
    vector<unsigned int> hitsInClusterDetIDs    = clusters->GetHits()->at(iCluster);
    vector<unsigned int> hitsInClusterDetIDsUnsorted = clusters->GetHits()->at(iCluster);
    vector<unsigned int> clusterToHitID;

    sort(hitsDetIDs->begin(), hitsDetIDs->end());
    sort(hitsInClusterDetIDs.begin(), hitsInClusterDetIDs.end());

    set_intersection(hitsDetIDs->begin(), hitsDetIDs->end(),
                     hitsInClusterDetIDs.begin(), hitsInClusterDetIDs.end(),
                     std::back_inserter(clusterToHitID));

    unique_ptr<RecHits> hitsInThisCluster = unique_ptr<RecHits>(new RecHits());

    for(unsigned int i : clusterToHitID){
      ptrdiff_t pos = distance(detid->begin(), find(detid->begin(), detid->end(), i));
      ptrdiff_t posCluster = distance(hitsInClusterDetIDsUnsorted.begin(), find(hitsInClusterDetIDsUnsorted.begin(), hitsInClusterDetIDsUnsorted.end(), i));
      double energyFraction = hitsInClusterFractions[(int)posCluster];
      unique_ptr<RecHit> hit = GetHit((int)pos,energyFraction);
      hitsInThisCluster->AddHit(hit);
    }
    nAssociatedHits += hitsInThisCluster->N();
    hitsPerCluster.push_back(move(hitsInThisCluster));
  }
  if(ConfigurationManager::Instance()->GetVerbosityLevel() > 0){
    cout<<"\nnum of rechits associated with sim-clusters:"<<nAssociatedHits<<"\n"<<endl;
  }
}

void RecHits::GetRecHitsPerHexel(vector<RecHits*> &hitsClustered,vector<shared_ptr<Hexel>> &hexels)
{
  unique_ptr<RecHits> hitsAboveNoise = GetHitsAboveNoise();

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

  for (uint i = 0; i < hitDetIDs->size(); ++i) {
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

  for(uint iHex=0;iHex<hexels.size();iHex++){
    int hitIndex = hexelToHitID[iHex];
    unique_ptr<RecHit> hit = GetHit(hitIndex);
    hitsClustered[hexels[iHex]->clusterIndex]->AddHit(hit);
  }
}

unique_ptr<RecHits> RecHits::GetHitsInLayer(int layerIndex)
{
  unique_ptr<RecHits> hitsInLayer(new RecHits());

  for(int iHit=0;iHit<N();iHit++){
    if(layer->at(iHit) == layerIndex){
      unique_ptr<RecHit> hit = GetHit(iHit);
      hitsInLayer->AddHit(hit);
    }
  }
  return hitsInLayer;
}

tuple<bool,double> RecHits::RecHitAboveThreshold(double iHit)
{
  ConfigurationManager *config = ConfigurationManager::Instance();
  
  double sigmaNoise = 1.;
  
  if(config->GetDependSensor()){
    int thickIndex = -1;
    
    if(layer->at(iHit) <= lastLayerFH){  // EE + FH
      if(thickness->at(iHit) > 99. && thickness->at(iHit) < 101.)        thickIndex = 0;
      else if(thickness->at(iHit) > 199. && thickness->at(iHit) < 201.)  thickIndex = 1;
      else if(thickness->at(iHit) > 299. && thickness->at(iHit) < 301.)  thickIndex = 2;
      else cout<<"ERROR - silicon thickness has a nonsensical value"<<endl;
      // determine noise for each sensor/subdetector using RecHitCalibration library
    }
    sigmaNoise = 0.001 * recHitCalib->sigmaNoiseMeV(layer->at(iHit), thickIndex);  // returns threshold for EE, FH, BH (in case of BH thickIndex does not play a role)
  }
  bool aboveThreshold = energy->at(iHit) >= config->GetEnergyMin() * sigmaNoise;  // this checks if rechit energy is above the threshold of ecut (times the sigma noise for the sensor, if that option is set)
  return make_tuple(aboveThreshold,sigmaNoise);
}

void RecHits::ShareCommonHits(unique_ptr<RecHits> &hits)
{
  // store DetIDs that are common for this RecHits and the provided ones
  std::vector<unsigned int> common;
  std::set_intersection(this->GetDetIDs()->begin(), this->GetDetIDs()->end(),
                        hits->GetDetIDs()->begin(), hits->GetDetIDs()->end(),
                        std::back_inserter(common));
  
  // Calculate energy of the core (without shared hits)
  double coreEnergyThis = this->GetTotalEnergy();
  double coreEnergyHits = hits->GetTotalEnergy();
  
  for(int i=0;i<common.size();i++){
    auto posThis = distance(this->GetDetIDs()->begin(),
                            find(this->GetDetIDs()->begin(),this->GetDetIDs()->end(), common[i]));
    
    auto posHits = distance(hits->GetDetIDs()->begin(),
                            find(hits->GetDetIDs()->begin(),hits->GetDetIDs()->end(), common[i]));
    
    coreEnergyThis += this->GetEnergyOfHit((int)posThis);
    coreEnergyHits += hits->GetEnergyOfHit((int)posHits);
  }
  
  double fractionThis = coreEnergyThis/(coreEnergyThis+coreEnergyHits);
  double fractionHits = 1-fractionThis;
  
  for(int i=0;i<common.size();i++){
    int posThis = (int)distance(this->GetDetIDs()->begin(),
                                find(this->GetDetIDs()->begin(),this->GetDetIDs()->end(), common[i]));
    
    int posHits = (int)distance(hits->GetDetIDs()->begin(),
                                find(hits->GetDetIDs()->begin(),hits->GetDetIDs()->end(), common[i]));
    
    this->SetEnergyOfHit(posThis, fractionThis * this->GetEnergyOfHit(posThis));
    hits->SetEnergyOfHit(posHits, fractionHits * hits->GetEnergyOfHit(posHits));
  }
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
  recHitCalibration = unique_ptr<RecHitCalibration>(new RecHitCalibration());
}

RecHit::~RecHit()
{

}

unique_ptr<Hexel> RecHit::GetHexel()
{
  unique_ptr<Hexel> hexel(new Hexel(eta,phi,x,y,z,energy,thickness,time,detid,layer,cluster2d,isHalf));
  return hexel;
}

tuple<bool,double> RecHit::RecHitAboveThreshold()
{
  ConfigurationManager *config = ConfigurationManager::Instance();
  
  double sigmaNoise = 1.;

  if(config->GetDependSensor()){
    int thickIndex = -1;

    if(layer <= lastLayerFH){  // EE + FH
      if(thickness > 99. && thickness < 101.)        thickIndex = 0;
      else if(thickness > 199. && thickness < 201.)  thickIndex = 1;
      else if(thickness > 299. && thickness < 301.)  thickIndex = 2;
      else cout<<"ERROR - silicon thickness has a nonsensical value"<<endl;
      // determine noise for each sensor/subdetector using RecHitCalibration library
    }
    sigmaNoise = 0.001 * recHitCalibration->sigmaNoiseMeV(layer, thickIndex);  // returns threshold for EE, FH, BH (in case of BH thickIndex does not play a role)
  }
  bool aboveThreshold = energy >= config->GetEnergyMin() * sigmaNoise;  // this checks if rechit energy is above the threshold of ecut (times the sigma noise for the sensor, if that option is set)
  return make_tuple(aboveThreshold,sigmaNoise);
}
