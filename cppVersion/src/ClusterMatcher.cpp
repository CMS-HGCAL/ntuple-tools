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

void ClusterMatcher::MatchClustersByDetID(vector<MatchedClusters*> &matched,
                                          vector<RecHits*> &recHitsPerCluster, vector<RecHits*> &simHitsPerCluster,
                                          int layer)
{
  vector<RecHits*> hitsMatchedToRecClusters;
  vector<vector<unsigned int>> recDetIDs;
  vector<double> Xs,Ys,Rs, Es;
  
  for(uint recClusterIndex=0;recClusterIndex < recHitsPerCluster.size();recClusterIndex++){
    
    RecHits *recCluster = recHitsPerCluster[recClusterIndex];
    unique_ptr<RecHits> recHitsInLayerInCluster = recCluster->GetHitsInLayer(layer);
    
    if(recHitsInLayerInCluster->N()==0) continue;
    
    BasicCluster *basicCluster = GetBasicClusterFromRecHits(recHitsInLayerInCluster);
    
    vector<unsigned int> detIDsInRecCluster;
    
    for(int iHit=0;iHit<recHitsInLayerInCluster->N();iHit++){
      detIDsInRecCluster.push_back(recHitsInLayerInCluster->GetDetIDofHit(iHit));
    }
    sort(detIDsInRecCluster.begin(),detIDsInRecCluster.end());
    recDetIDs.push_back(detIDsInRecCluster);
    
    Xs.push_back(basicCluster->GetX());
    Ys.push_back(basicCluster->GetY());
    Rs.push_back(basicCluster->GetRadius());
    Es.push_back(basicCluster->GetEnergy());
    
    MatchedClusters *matchedCluster = new MatchedClusters();
    
    matchedCluster->recCluster = basicCluster;
    matched.push_back(matchedCluster);
  }
   double maxDistance = ConfigurationManager::Instance()->GetMachingMaxDistance();
  
  for(uint simClusterIndex=0;simClusterIndex<simHitsPerCluster.size();simClusterIndex++){
    
    RecHits *simCluster = simHitsPerCluster[simClusterIndex];
    unique_ptr<RecHits> simHitsInLayerInCluster = simCluster->GetHitsInLayer(layer);
    
    if(simHitsInLayerInCluster->N()==0) continue;
    
    BasicCluster *basicCluster = GetBasicClusterFromRecHits(simHitsInLayerInCluster);
    
    vector<unsigned int> detIDsInSimCluster;
    for(int iHit=0;iHit<simHitsInLayerInCluster->N();iHit++){
      detIDsInSimCluster.push_back(simHitsInLayerInCluster->GetDetIDofHit(iHit));
    }
    sort(detIDsInSimCluster.begin(),detIDsInSimCluster.end());
    int parentRecCluster = findMostDetIDsharingCluster(recDetIDs, detIDsInSimCluster);
    
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
    
//    matched[parentRecCluster]->simClusters->push_back(basicCluster);
  }
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
    
    matchedCluster->recCluster = basicCluster;
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
      matchedCluster->recCluster = GetBasicClusterFromRecHits(recHitsInLayerInCluster);
      
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
