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

void ClusterMatcher::MatchClustersClosest(vector<MatchedClusters*> &matched,
                                          vector<RecHits*> &recHitsPerCluster, vector<RecHits*> &simHitsPerCluster,
                                          int layer)
{
  vector<int> alreadyAssociatedClusters;
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
    
    
    matched[parentRecCluster]->simClusters->push_back(basicCluster);
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
