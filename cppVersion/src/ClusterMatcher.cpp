//
//  ClusterMatcher.cpp
//
//  Created by Jeremi Niedziela on 05/07/2018.
//

#include "ClusterMatcher.hpp"

using namespace std;

ClusterMatcher::ClusterMatcher()
{
  
}

ClusterMatcher::~ClusterMatcher()
{
  
}

int ContainsSimCluster(vector<MatchedClusters*> &matched, int simClusterIndex){
  for(int i=0;i<matched.size();i++){
    if(find(matched[i]->simIndices.begin(),
            matched[i]->simIndices.begin(),simClusterIndex) != matched[i]->simIndices.end()){
      return i;
    }
  }
  return -1;
}

void ClusterMatcher::MatchClustersByDetID(vector<MatchedClusters*> &matched,
                                          vector<RecHits*> &recHitsPerCluster,
                                          vector<RecHits*> &simHitsPerCluster,
                                          int layer)
{
  vector<RecHits*> hitsMatchedToRecClusters;
  vector<vector<unsigned int>> recDetIDs;
  vector<vector<unsigned int>> simDetIDs;
  vector<double> Xs,Ys,Rs,Es;
  vector<double> Xsim,Ysim,Rsim,Esim;
  
  vector<unsigned long> assignedRecClusters;
  vector<int> simClustersParents;
  double maxDistance = ConfigurationManager::Instance()->GetMachingMaxDistance();
  vector<int> forRemoval;
  
  // loop over all rec clusters and add then to matched clusters vector
  for(uint recClusterIndex=0;recClusterIndex < recHitsPerCluster.size();recClusterIndex++){
    
    RecHits *recCluster = recHitsPerCluster[recClusterIndex];
    unique_ptr<RecHits> recHitsInLayerInCluster = recCluster->GetHitsInLayer(layer);
    
    if(recHitsInLayerInCluster->N()==0) continue;
    
    BasicCluster *basicCluster = GetBasicClusterFromRecHits(recHitsInLayerInCluster);
    
    vector<unsigned int> detIDsInRecCluster;
    vector<double> energiesInRecCluster;
    
    for(int iHit=0;iHit<recHitsInLayerInCluster->N();iHit++){
      detIDsInRecCluster.push_back(recHitsInLayerInCluster->GetDetIDofHit(iHit));
      double e = recHitsInLayerInCluster->GetEnergyOfHit(iHit);
      energiesInRecCluster.push_back(e);
    }
    sort(detIDsInRecCluster.begin(),detIDsInRecCluster.end());
    recDetIDs.push_back(detIDsInRecCluster);
    
    Xs.push_back(basicCluster->GetX());
    Ys.push_back(basicCluster->GetY());
    Rs.push_back(basicCluster->GetRadius());
    Es.push_back(basicCluster->GetEnergy());
    
    MatchedClusters *matchedCluster = new MatchedClusters();
    
    matchedCluster->recClusters->push_back(basicCluster);
    matchedCluster->recIndices.push_back(recClusterIndex);
    matchedCluster->AddRecDetIDs(detIDsInRecCluster);
    matchedCluster->AddRecEnergies(energiesInRecCluster);
    matchedCluster->recHits->AddHits(recHitsInLayerInCluster);
    matched.push_back(matchedCluster);
  }
  
  // for each sim cluster, find a rec cluster that shares the most hits. If distance between them is smaller than maxDistance, add this sim cluster to matching rec cluster
  for(uint simClusterIndex=0;simClusterIndex<simHitsPerCluster.size();simClusterIndex++){
    
    RecHits *simCluster = simHitsPerCluster[simClusterIndex];
    unique_ptr<RecHits> simHitsInLayerInCluster = simCluster->GetHitsInLayer(layer);
    
    if(simHitsInLayerInCluster->N()==0) continue;
    
    BasicCluster *basicCluster = GetBasicClusterFromRecHits(simHitsInLayerInCluster);
    
    vector<unsigned int> detIDsInSimCluster;
    vector<double> energiesInSimCluster;
    
    for(int iHit=0;iHit<simHitsInLayerInCluster->N();iHit++){
      detIDsInSimCluster.push_back(simHitsInLayerInCluster->GetDetIDofHit(iHit));
      energiesInSimCluster.push_back(simHitsInLayerInCluster->GetEnergyOfHit(iHit));
    }
    sort(detIDsInSimCluster.begin(),detIDsInSimCluster.end());
    simDetIDs.push_back(detIDsInSimCluster);
    
    Xsim.push_back(basicCluster->GetX());
    Ysim.push_back(basicCluster->GetY());
    Rsim.push_back(basicCluster->GetRadius());
    Esim.push_back(basicCluster->GetEnergy());
    
    int parentRecCluster = findMostDetIDsharingCluster(recDetIDs, detIDsInSimCluster);
    
    if(parentRecCluster < 0){
      if(ConfigurationManager::Instance()->GetVerbosityLevel() >= 1){
        cout<<"No rec cluster found for a sim cluster!!"<<endl;
      }
      continue;
    }
    
    double distance = sqrt(pow(Xs[parentRecCluster]-basicCluster->GetX(),2)+pow(Ys[parentRecCluster]-basicCluster->GetY(),2));
    
    if((distance > maxDistance) && maxDistance>=0) continue;
    
    matched[parentRecCluster]->simClusters->push_back(basicCluster);
    matched[parentRecCluster]->simIndices.push_back(simClusterIndex);
    matched[parentRecCluster]->AddSimDetIDs(detIDsInSimCluster);
    matched[parentRecCluster]->AddSimEnergies(energiesInSimCluster);
    matched[parentRecCluster]->simHits->AddHits(simHitsInLayerInCluster);
    assignedRecClusters.push_back(parentRecCluster);
  }
  
  // find rec clusters that were not matched with any sim clusters (they can be smaller clusters, which were not the "stronges" candidates, but still can be re-assigned to sim clusters in second interation)
  for(uint recClusterIndex=0;recClusterIndex < recHitsPerCluster.size();recClusterIndex++){
    if(find(assignedRecClusters.begin(),assignedRecClusters.end(),recClusterIndex)==assignedRecClusters.end()){
      
      
      RecHits *recCluster = recHitsPerCluster[recClusterIndex];
      unique_ptr<RecHits> recHitsInLayerInCluster = recCluster->GetHitsInLayer(layer);
      
      if(recHitsInLayerInCluster->N()==0) continue;
      
      BasicCluster *basicCluster = GetBasicClusterFromRecHits(recHitsInLayerInCluster);
      
      vector<unsigned int> detIDsInRecCluster;
      for(int iHit=0;iHit<recHitsInLayerInCluster->N();iHit++){
        detIDsInRecCluster.push_back(recHitsInLayerInCluster->GetDetIDofHit(iHit));
      }
      sort(detIDsInRecCluster.begin(),detIDsInRecCluster.end());
      
      int parentSimCluster = findMostDetIDsharingCluster(simDetIDs, detIDsInRecCluster);
      
      if(parentSimCluster < 0){
        if(ConfigurationManager::Instance()->GetVerbosityLevel() >= 1){
          cout<<"No sim cluster found for a rec cluster!!"<<endl;
        }
        continue;
      }
      
      double distance = sqrt(pow(Xsim[parentSimCluster]-basicCluster->GetX(),2)+pow(Ysim[parentSimCluster]-basicCluster->GetY(),2));
      
      if((distance > maxDistance) && maxDistance>=0) continue;
        
        
      int simAlreadyIn = ContainsSimCluster(matched, parentSimCluster);
        
      if(simAlreadyIn >= 0){
        matched[simAlreadyIn]->recClusters->push_back(basicCluster);
        matched[simAlreadyIn]->recIndices.push_back(recClusterIndex);
      }
    }
  }
  
  // remove matched clusters with no sim assigned to rec
//  for(int i=0;i<matched.size();i++){
//    if(matched[i]->simClusters->size() == 0){
//      matched.erase(matched.begin()+i);
//      i--;
//    }
//  }
}


void ClusterMatcher::MatchClustersClosest(vector<MatchedClusters*> &matched,
                                          vector<RecHits*> &recHitsPerCluster, vector<RecHits*> &simHitsPerCluster,
                                          int layer)
{
  vector<RecHits*> hitsMatchedToRecClusters;
  vector<double> Xs,Ys,Rs, Es;

  for(uint recClusterIndex=0;recClusterIndex < recHitsPerCluster.size();recClusterIndex++){
    
    RecHits *recCluster = recHitsPerCluster[recClusterIndex];
    unique_ptr<RecHits> recHitsInLayerInCluster = recCluster->GetHitsInLayer(layer);
    
    if(recHitsInLayerInCluster->N()==0) continue;
    
    BasicCluster *basicCluster = GetBasicClusterFromRecHits(recHitsInLayerInCluster);
    
    Xs.push_back(basicCluster->GetX());
    Ys.push_back(basicCluster->GetY());
    Rs.push_back(basicCluster->GetRadius());
    Es.push_back(basicCluster->GetEnergy());
    
    MatchedClusters *matchedCluster = new MatchedClusters();
    
    matchedCluster->recClusters->push_back(basicCluster);
    matched.push_back(matchedCluster);
  }
  
  double maxDistance = ConfigurationManager::Instance()->GetMachingMaxDistance();
  
  for(uint simClusterIndex=0;simClusterIndex<simHitsPerCluster.size();simClusterIndex++){
    
    RecHits *simCluster = simHitsPerCluster[simClusterIndex];
    unique_ptr<RecHits> simHitsInLayerInCluster = simCluster->GetHitsInLayer(layer);
    
    if(simHitsInLayerInCluster->N()==0) continue;
    
    BasicCluster *basicCluster = GetBasicClusterFromRecHits(simHitsInLayerInCluster);
    int parentRecCluster = findClosestCircle(Xs, Ys, Rs, basicCluster->GetX(), basicCluster->GetY());
    
    if(parentRecCluster < 0){
      if(ConfigurationManager::Instance()->GetVerbosityLevel() >= 1){
        cout<<"No rec cluster found for a sim cluster!!"<<endl;
      }
      continue;
    }
    double distance = sqrt(pow(Xs[parentRecCluster]-basicCluster->GetX(),2)+pow(Ys[parentRecCluster]-basicCluster->GetY(),2));
    
    if((distance <= maxDistance) || maxDistance < 0){
      matched[parentRecCluster]->simClusters->push_back(basicCluster);
    }
  }
}

void ClusterMatcher::MatchClustersAllToAll(vector<MatchedClusters*> &matched,
                                           vector<RecHits*> &recHitsPerCluster, vector<RecHits*> &simHitsPerCluster,
                                           int layer)
{
  for(uint recClusterIndex=0;recClusterIndex < recHitsPerCluster.size();recClusterIndex++){
    
    RecHits *recCluster = recHitsPerCluster[recClusterIndex];
    unique_ptr<RecHits> recHitsInLayerInCluster = recCluster->GetHitsInLayer(layer);
    
    if(recHitsInLayerInCluster->N()==0) continue;
    
    for(uint simClusterIndex=0;simClusterIndex < simHitsPerCluster.size();simClusterIndex++){
      
      RecHits *simCluster = simHitsPerCluster[simClusterIndex];
      unique_ptr<RecHits> simHitsInLayerInCluster = simCluster->GetHitsInLayer(layer);
      
      if(simHitsInLayerInCluster->N()==0) continue;
      
      MatchedClusters *matchedCluster = new MatchedClusters();
      matchedCluster->recClusters->push_back(GetBasicClusterFromRecHits(recHitsInLayerInCluster));
      
      BasicCluster *basicCluster = GetBasicClusterFromRecHits(simHitsInLayerInCluster);
      matchedCluster->simClusters->push_back(basicCluster);
      
      matched.push_back(matchedCluster);
    }
  }
}

BasicCluster* ClusterMatcher::GetBasicClusterFromRecHits(unique_ptr<RecHits> &hits)
{
  double recEnergy = hits->GetTotalEnergy();
  double xMax   = hits->GetXmax();
  double xMin   = hits->GetXmin();
  double yMax   = hits->GetYmax();
  double yMin   = hits->GetYmin();
  
  double clusterX = (xMax+xMin)/2.;
  double clusterY = (yMax+yMin)/2.;
  double clusterEta = hits->GetCenterEta();
  double clusterR = max((xMax-xMin)/2.,(yMax-yMin)/2.);
  
  BasicCluster *basicCluster = new BasicCluster(recEnergy,clusterX,clusterY,0,clusterEta,clusterR);
  return basicCluster;
}
