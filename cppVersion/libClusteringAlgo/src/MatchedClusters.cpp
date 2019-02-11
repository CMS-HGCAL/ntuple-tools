//
//  MatchedClusters.cpp
//
//  Created by Jeremi Niedziela on 27/08/2018.
//

#include "MatchedClusters.hpp"

#include <iostream>

using namespace std;

MatchedClusters::MatchedClusters()
{
  recClusters = make_unique<std::vector<shared_ptr<BasicCluster>>>();
  simClusters = make_unique<std::vector<shared_ptr<BasicCluster>>>();
  recHits = std::unique_ptr<RecHits>(new RecHits());
  simHits = std::unique_ptr<RecHits>(new RecHits());
}

MatchedClusters::~MatchedClusters()
{
  
}

void MatchedClusters::Print()
{
  cout<<"Matched cluster: "<<simClusters->size()<<" sim, "<<recClusters->size()<<" rec clusters in the collection."<<endl;
  
  cout<<"Merged sim properties:";
//  cout<<"\teta:"<<GetSimEta();
  cout<<"\tx:"<<GetSimX();
  cout<<"\ty:"<<GetSimY();
  cout<<"\tE:"<<GetSimEnergy();
  cout<<endl;
  
  for(int i=0;i<simClusters->size();i++){
    cout<<"\t\t\tSim "<<i<<":\tx:"<<simClusters->at(i)->GetX();
    cout<<"\ty:"<<simClusters->at(i)->GetY();
    cout<<"\tE:"<<simClusters->at(i)->GetEnergy();
    cout<<endl;
  }
  simHits->Print();

  cout<<"Merged rec properties:";
//  cout<<"\teta:"<<GetRecEta();
  cout<<"\tx:"<<GetRecX();
  cout<<"\ty:"<<GetRecY();
  cout<<"\tE:"<<GetRecEnergy();
  cout<<endl;
  
  for(int i=0;i<recClusters->size();i++){
    cout<<"\t\t\tRec "<<i<<":\tx:"<<recClusters->at(i)->GetX();
    cout<<"\ty:"<<recClusters->at(i)->GetY();
    cout<<"\tE:"<<recClusters->at(i)->GetEnergy();
    cout<<endl;
  }
  recHits->Print();
  
}

void MatchedClusters::AddRecCluster(int index, std::unique_ptr<RecHits> &hits)
{
  recClusters->push_back(GetBasicClusterFromRecHits(hits));
  recIndices.push_back(index);
  
  for(int iHit=0;iHit<hits->N();iHit++){
    recDetIDs.push_back(hits->GetDetIDofHit(iHit));
  }
  sort(recDetIDs.begin(),recDetIDs.end());
  recHits->AddHits(hits);
}

void MatchedClusters::AddSimCluster(int index,std::unique_ptr<RecHits> &hits)
{
  simClusters->push_back(GetBasicClusterFromRecHits(hits));
  simIndices.push_back(index);
  
  for(int iHit=0;iHit<hits->N();iHit++){
    simDetIDs.push_back(hits->GetDetIDofHit(iHit));
  }
  sort(simDetIDs.begin(),simDetIDs.end());
  simHits->AddHits(hits);
}

void MatchedClusters::Merge(shared_ptr<MatchedClusters> &clusters)
{
  for(auto basicCluster : *(clusters->simClusters)){simClusters->push_back(basicCluster);}
  for(auto index : clusters->simIndices){simIndices.push_back(index);}
  for(auto id : clusters->simDetIDs){simDetIDs.push_back(id);}
  simHits->AddHits(clusters->simHits);
  
  for(auto basicCluster : *(clusters->recClusters)){recClusters->push_back(basicCluster);}
  for(auto index : clusters->recIndices){recIndices.push_back(index);}
  for(auto id : clusters->recDetIDs){recDetIDs.push_back(id);}
  recHits->AddHits(clusters->recHits);
}

shared_ptr<BasicCluster> MatchedClusters::GetBasicClusterFromRecHits(std::unique_ptr<RecHits> &hits)
{
  double recEnergy = hits->GetTotalEnergy();
  int maxHit = hits->GetHighestEnergyHitIndex();
   
  double clusterX   = hits->GetXofHit(maxHit);
  double clusterY   = hits->GetYofHit(maxHit);
  double clusterEta = hits->GetCenterEta();
  
  // vector of squared distances and energies
  vector<tuple<float,float>> distances2andEnergies;
  
  for(int i=0;i<hits->N();i++){
    distances2andEnergies.push_back(make_tuple( pow(hits->GetXofHit(i)-clusterX,2)
                                               +pow(hits->GetYofHit(i)-clusterY,2),
                                               hits->GetEnergyOfHit(i)));
  }
  
  // sort hits by distance to cluster's center
  sort(begin(distances2andEnergies), end(distances2andEnergies),
          [](tuple<float, float> const &t1, tuple<float, float> const &t2) {return get<0>(t1) < get<0>(t2);}
       );
  
  double accumulatedEnergy = 0;
  double clusterR = 0;
  
  
  for(int i=0;i<distances2andEnergies.size();i++){
    accumulatedEnergy += get<1>(distances2andEnergies[i]);
    clusterR = get<0>(distances2andEnergies[i]);
    
    if(accumulatedEnergy > 0.68*recEnergy) break;
  }
  clusterR = sqrt(clusterR);
  
  if(clusterR==0) clusterR = 0.5; // this means that cluster has just one hit
  
  auto basicCluster = make_shared<BasicCluster>(recEnergy,clusterX,clusterY,0,clusterEta,clusterR);
  return basicCluster;
}

double MatchedClusters::GetSharedFractionWithRecHits(vector<unsigned int> &detIDs)
{
  // store DetIDs that are common for rec and sim clusters
  std::vector<unsigned int> common;
  std::set_intersection(detIDs.begin(), detIDs.end(),
                        recDetIDs.begin(), recDetIDs.end(),
                        std::back_inserter(common));
  
  // Sum up energy of shared hits
  double energySumShared=0;
  for(int i=0;i<common.size();i++){
    auto pos = distance(recDetIDs.begin(),
                        find(recDetIDs.begin(),recDetIDs.end(), common[i]));
    
    energySumShared += recHits->GetEnergyOfHit((int)pos);
  }
  
  // Return total shared energy
  return energySumShared;
}

double MatchedClusters::GetSharedFraction()
{
  // store DetIDs that are common for rec and sim clusters
  std::vector<unsigned int> common;
  std::set_intersection(simHits->GetDetIDs()->begin(), simHits->GetDetIDs()->end(),
                        recHits->GetDetIDs()->begin(), recHits->GetDetIDs()->end(),
                        std::back_inserter(common));
  
  // Sum up energy of shared hits
  double energySumShared=0;
  for(int i=0;i<common.size();i++){
    auto pos = distance(simHits->GetDetIDs()->begin(),
                        find(simHits->GetDetIDs()->begin(),simHits->GetDetIDs()->end(), common[i]));
    
    energySumShared += simHits->GetEnergyOfHit((int)pos);
  }
  
  // Return shared energy normalized to the total sim energy
  return energySumShared/GetSimEnergy();
}

shared_ptr<BasicCluster> MatchedClusters::GetRecClusterByIndex(int index){
  ptrdiff_t pos = distance(recIndices.begin(), find(recIndices.begin(), recIndices.end(), index));
  if(pos >= recClusters->size()) return nullptr;
  return recClusters->at(pos);
}

shared_ptr<BasicCluster> MatchedClusters::GetSimClusterByIndex(int index){
  ptrdiff_t pos = distance(simIndices.begin(), find(simIndices.begin(), simIndices.end(), index));
  if(pos >= simClusters->size()) return nullptr;
  return simClusters->at(pos);
}

double MatchedClusters::GetRecX(){
  if(recHits->N() <= 0) return 99999;
  return GetBasicClusterFromRecHits(recHits)->GetX();
}

double MatchedClusters::GetRecY(){
  if(recHits->N() <= 0) return 99999;
  return GetBasicClusterFromRecHits(recHits)->GetY();
}

double MatchedClusters::GetRecRadius(){
  if(recHits->N() <= 0) return 99999;
  return GetBasicClusterFromRecHits(recHits)->GetRadius();
}

double MatchedClusters::GetRecEta(){
  if(recHits->N() <= 0) return 99999;
  return GetBasicClusterFromRecHits(recHits)->GetEta();
}

double MatchedClusters::GetRecEnergy(){
  if(recHits->N() <= 0) return 99999;
  return GetBasicClusterFromRecHits(recHits)->GetEnergy();
}

double MatchedClusters::GetSimX(){
  if(simHits->N() <= 0) return 99999;
  return GetBasicClusterFromRecHits(simHits)->GetX();
}

double MatchedClusters::GetSimY(){
  if(simHits->N() <= 0) return 99999;
  return GetBasicClusterFromRecHits(simHits)->GetY();
}

double MatchedClusters::GetSimRadius(){
  if(simHits->N() <= 0) return 99999;
  return GetBasicClusterFromRecHits(simHits)->GetRadius();
}

double MatchedClusters::GetSimEta(){
  if(simHits->N() <= 0) return 99999;
  return GetBasicClusterFromRecHits(simHits)->GetEta();
}

double MatchedClusters::GetSimEnergy(){
  if(simHits->N() <= 0) return 99999;
  return GetBasicClusterFromRecHits(simHits)->GetEnergy();
}

bool MatchedClusters::ContainsSimCluster(int simClusterIndex){
  return (find(simIndices.begin(),simIndices.end(),simClusterIndex) != simIndices.end());
}

bool MatchedClusters::HasSimClusters(){
  return (simClusters->size()!=0);
}

bool MatchedClusters::HasRecClusters(){
  return (recClusters->size()!=0);
}



