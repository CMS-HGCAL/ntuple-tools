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
    if(matched[i]->ContainsSimCluster(simClusterIndex)) return i;
  }
  return -1;
}

void ClusterMatcher::MatchClustersByDetID(vector<MatchedClusters*> &matched,
                                          vector<RecHits*> &recHitsPerCluster,
                                          vector<RecHits*> &simHitsPerCluster,
                                          int layer)
{
  double maxDistance = ConfigurationManager::Instance()->GetMachingMaxDistance();
  int verbosityLevel = ConfigurationManager::Instance()->GetVerbosityLevel();
 
  if(verbosityLevel > 0){
    cout<<"\nMatching clusters in layer "<<layer<<endl<<endl;
  }
  
  // loop over all rec clusters and add them to matched clusters vector
  for(uint recClusterIndex=0;recClusterIndex < recHitsPerCluster.size();recClusterIndex++){
    
    // Get rec cluster in correct layer
    RecHits *recCluster = recHitsPerCluster[recClusterIndex];
    unique_ptr<RecHits> recHitsInLayerInCluster = recCluster->GetHitsInLayer(layer);
    if(recHitsInLayerInCluster->N()==0) continue;
    
    if(verbosityLevel > 0){
      cout<<"Rec hits in cluster:"<<endl;
      recHitsInLayerInCluster->Print();
    }
    
    // Create new matched cluster for each rec cluster
    MatchedClusters *matchedCluster = new MatchedClusters();
    matchedCluster->AddRecCluster(recClusterIndex,recHitsInLayerInCluster);
    matched.push_back(matchedCluster);
  }
  
  if(verbosityLevel > 0){
    cout<<"\n\nMatched clusters after step 1\n\n"<<endl;
    for(auto cluster : matched){
      cluster->Print();
      cout<<endl;
    }
  }
  
  // for each sim cluster, find a rec cluster that shares the most hits.
  // If distance between them is smaller than maxDistance, add this sim cluster to matching rec cluster
  for(uint simClusterIndex=0;simClusterIndex<simHitsPerCluster.size();simClusterIndex++){
    
    // Get sim cluster in correct layer
    RecHits *simCluster = simHitsPerCluster[simClusterIndex];
    unique_ptr<RecHits> simHitsInLayerInCluster = simCluster->GetHitsInLayer(layer);
    if(simHitsInLayerInCluster->N()==0) continue;
    
    if(verbosityLevel > 0){
      cout<<"Sim hits in cluster:"<<endl;
      simHitsInLayerInCluster->Print();
    }
    
    vector<unsigned int> detIDsInSimCluster = *(simHitsInLayerInCluster->GetDetIDs());
    
    // Find matched cluster that shares the most hits with this sim cluster
    double maxShared = 0;
    MatchedClusters *maxSharingMatchedCluster = nullptr;
    for(auto matchedCluster : matched){
      double shared = matchedCluster->GetSharedFractionWithRecHits(detIDsInSimCluster);
      if(shared > maxShared){
        maxShared = shared;
        maxSharingMatchedCluster = matchedCluster;
      }
    }
    
    MatchedClusters *matchedTmp = new MatchedClusters();
    matchedTmp->AddSimCluster(simClusterIndex, simHitsInLayerInCluster);
    if(!maxSharingMatchedCluster){
      if(ConfigurationManager::Instance()->GetVerbosityLevel() > 0){
        cout<<"\n========================================================"<<endl;
        cout<<"No rec cluster found for a sim cluster!!"<<endl;
        cout<<"Sim cluster:"<<endl;
        matchedTmp->Print();
        cout<<endl;
      }
      continue;
    }
   
    
    double distance = sqrt( pow(maxSharingMatchedCluster->GetRecX() - matchedTmp->GetSimX(),2)
                           +pow(maxSharingMatchedCluster->GetRecY() - matchedTmp->GetSimY(),2));
    
    if((distance > maxDistance) && maxDistance>=0) continue;
    
    maxSharingMatchedCluster->AddSimCluster(simClusterIndex,simHitsInLayerInCluster);
  }
  
  if(verbosityLevel > 0){
    cout<<"\n\nMatched clusters after step 2\n\n"<<endl;
    for(auto cluster : matched){
      cluster->Print();
      cout<<endl;
    }
  }
 
  // Iterate over rec clusters that didn't get any sim clusters assigned and try to match them to existing sets of sim-rec clusters
  for(MatchedClusters* cluster : matched){
    if(!cluster->HasSimClusters()){
      vector<unsigned int> recDetIDs = cluster->GetRecDetIDs();
      double maxShared = 0;
      MatchedClusters *maxSharingMatchedCluster = nullptr;
      
      for(MatchedClusters* otherCluster : matched){
        if(cluster == otherCluster) continue;
      
        double shared = otherCluster->GetSharedFractionWithRecHits(recDetIDs);
        if(shared > maxShared){
          maxShared = shared;
          maxSharingMatchedCluster = otherCluster;
        }
      }
      
      MatchedClusters *matchedTmp = new MatchedClusters();
      unique_ptr<RecHits> recHitsInThisCluster = unique_ptr<RecHits>(new RecHits());
      cluster->GetRecHits(recHitsInThisCluster);
      matchedTmp->AddSimCluster(cluster->GetFirstRecIndex(), recHitsInThisCluster);
      
      if(!maxSharingMatchedCluster){
        if(ConfigurationManager::Instance()->GetVerbosityLevel() > 0){
          cout<<"\n========================================================"<<endl;
          cout<<"No sim cluster found for a rec cluster!!"<<endl;
          cout<<"Rec cluster:"<<endl;
          matchedTmp->Print();
          cout<<endl;
        }
        continue;
      }
      
      double distance = sqrt( pow(maxSharingMatchedCluster->GetRecX() - matchedTmp->GetSimX(),2)
                             +pow(maxSharingMatchedCluster->GetRecY() - matchedTmp->GetSimY(),2));
      
      if((distance > maxDistance) && maxDistance>=0) continue;
      maxSharingMatchedCluster->Merge(cluster);
    }
  }
  
  if(verbosityLevel > 0){
    cout<<"\n\nMatched clusters after step 3\n\n"<<endl;
    for(auto cluster : matched){
      cluster->Print();
      cout<<endl;
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
  vector<double> Xs,Ys,Rs, Es;
  double maxDistance = ConfigurationManager::Instance()->GetMachingMaxDistance();
  
  
  for(uint recClusterIndex=0;recClusterIndex < recHitsPerCluster.size();recClusterIndex++){
    
    // Get rec cluster in correct layer
    RecHits *recCluster = recHitsPerCluster[recClusterIndex];
    unique_ptr<RecHits> recHitsInLayerInCluster = recCluster->GetHitsInLayer(layer);
    if(recHitsInLayerInCluster->N()==0) continue;
    
    // Build new mached cluster for each rec cluster
    MatchedClusters *matchedCluster = new MatchedClusters();
    matchedCluster->AddRecCluster(recClusterIndex, recHitsInLayerInCluster);
    matched.push_back(matchedCluster);
    
    // Save coordinates of this rec cluster
    BasicCluster *basicCluster = matchedCluster->GetRecClusterByIndex(recClusterIndex);
    Xs.push_back(basicCluster->GetX());
    Ys.push_back(basicCluster->GetY());
    Rs.push_back(basicCluster->GetRadius());
    Es.push_back(basicCluster->GetEnergy());
  }
  
  
  for(uint simClusterIndex=0;simClusterIndex<simHitsPerCluster.size();simClusterIndex++){
    
    // Get sim cluster in correct layer
    RecHits *simCluster = simHitsPerCluster[simClusterIndex];
    unique_ptr<RecHits> simHitsInLayerInCluster = simCluster->GetHitsInLayer(layer);
    if(simHitsInLayerInCluster->N()==0) continue;
    
    // Get coordinates of this sim cluster
    unique_ptr<MatchedClusters> matchedTmp = unique_ptr<MatchedClusters>(new MatchedClusters());
    matchedTmp->AddSimCluster(simClusterIndex, simHitsInLayerInCluster);
    BasicCluster *basicCluster = matchedTmp->GetSimClusterByIndex(simClusterIndex);
    
    // Check which of the rec clusters is the closest to this sim cluster
    int parentRecCluster = findClosestCircle(Xs, Ys, Rs, basicCluster->GetX(), basicCluster->GetY());
    
    if(parentRecCluster < 0){
      if(ConfigurationManager::Instance()->GetVerbosityLevel() >= 1){
        cout<<"Matching clusters by distance -- No rec cluster found for a sim cluster!!"<<endl;
      }
      continue;
    }
    
    // Calculate the distance to the closest rec cluster
    double distance = sqrt( pow(Xs[parentRecCluster]-basicCluster->GetX(),2)
                           +pow(Ys[parentRecCluster]-basicCluster->GetY(),2));
    
    // If distance to the closest rec cluster is below the limit, add this sim cluster to the matched clusters
    if((distance <= maxDistance) || maxDistance < 0){
      matched[parentRecCluster]->AddSimCluster(simClusterIndex, simHitsInLayerInCluster);
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
      matchedCluster->AddRecCluster(recClusterIndex, recHitsInLayerInCluster);
      matchedCluster->AddSimCluster(simClusterIndex, simHitsInLayerInCluster);
      
      matched.push_back(matchedCluster);
    }
  }
}

