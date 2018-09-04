//
//  ClusterMatcher.cpp
//
//  Created by Jeremi Niedziela on 05/07/2018.
//

#include "ClusterMatcher.hpp"

#include <TCanvas.h>
#include <TGraph.h>
#include <TAxis.h>

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
                                          vector<unique_ptr<RecHits>> &recHitsPerCluster,
                                          vector<unique_ptr<RecHits>> &simHitsPerCluster,
                                          bool draw)
{
  double maxDistance = ConfigurationManager::Instance()->GetMachingMaxDistance();
  int verbosityLevel = ConfigurationManager::Instance()->GetVerbosityLevel();
 
  // loop over all rec clusters and add them to matched clusters vector
  for(uint recClusterIndex=0;recClusterIndex < recHitsPerCluster.size();recClusterIndex++){
    if(verbosityLevel > 1){
      cout<<"Rec hits in cluster:"<<endl;
      recHitsPerCluster[recClusterIndex]->Print();
    }
    // Create new matched cluster for each rec cluster
    MatchedClusters *matchedCluster = new MatchedClusters();
    matchedCluster->AddRecCluster(recClusterIndex,recHitsPerCluster[recClusterIndex]);
    matched.push_back(matchedCluster);
  }
  
  if(verbosityLevel > 1){
    cout<<"\n\nMatched clusters after step 1\n\n"<<endl;
    for(auto cluster : matched){
      cluster->Print();
      cout<<endl;
    }
  }
  
  // for each sim cluster, find a rec cluster that shares the most hits.
  // If distance between them is smaller than maxDistance, add this sim cluster to matching rec cluster
  for(uint simClusterIndex=0;simClusterIndex<simHitsPerCluster.size();simClusterIndex++){
    if(verbosityLevel > 1){
      cout<<"Sim hits in cluster:"<<endl;
      simHitsPerCluster[simClusterIndex]->Print();
    }
    
    vector<unsigned int> detIDsInSimCluster = *(simHitsPerCluster[simClusterIndex]->GetDetIDs());
    
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
    matchedTmp->AddSimCluster(simClusterIndex, simHitsPerCluster[simClusterIndex]);
    if(!maxSharingMatchedCluster){
      if(ConfigurationManager::Instance()->GetVerbosityLevel() > 1){
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
    
    maxSharingMatchedCluster->AddSimCluster(simClusterIndex,simHitsPerCluster[simClusterIndex]);
  }
  
  if(verbosityLevel > 1){
    cout<<"\n\nMatched clusters after step 2\n\n"<<endl;
    for(auto cluster : matched){
      cluster->Print();
      cout<<endl;
    }
  }
 
  // Iterate over rec clusters that didn't get any sim clusters assigned and try to match them to existing sets of sim-rec clusters
  vector<int> matchedToErase;
  
  for(int iRec=0;iRec<matched.size();iRec++){
    MatchedClusters* failedRecCluster = matched[iRec];
    
    if(!failedRecCluster->HasSimClusters()){
      
      double maxShared = 0;
      MatchedClusters *maxSharingMatchedCluster = nullptr;
      
      for(MatchedClusters* otherCluster : matched){
        if(failedRecCluster == otherCluster) continue;
      
        vector<unsigned int> simDetIDs = otherCluster->GetSimDetIDs();
        
        double shared = failedRecCluster->GetSharedFractionWithRecHits(simDetIDs);
        if(shared > maxShared){
          maxShared = shared;
          maxSharingMatchedCluster = otherCluster;
        }
      }
      
      if(!maxSharingMatchedCluster){
        if(ConfigurationManager::Instance()->GetVerbosityLevel() > 1){
          cout<<"\n========================================================"<<endl;
          cout<<"No sim cluster found for a rec cluster!!"<<endl;
          cout<<"Rec cluster:"<<endl;
          failedRecCluster->Print();
          cout<<endl;
        }
        continue;
      }
      
      double distance = sqrt( pow(maxSharingMatchedCluster->GetSimX() - failedRecCluster->GetRecX(),2)
                             +pow(maxSharingMatchedCluster->GetSimY() - failedRecCluster->GetRecY(),2));
      
      if((distance > maxDistance) && maxDistance>=0) continue;
      maxSharingMatchedCluster->Merge(failedRecCluster);
      matchedToErase.push_back(iRec);
    }
  }
  
  if(verbosityLevel > 1){
    cout<<"\n\nMatched clusters after step 3\n\n"<<endl;
    for(auto cluster : matched){
      cluster->Print();
    }
    cout<<endl;
  }
  
  vector<MatchedClusters*> newMached;
  
  for(int i=0;i<matched.size();i++){
    if(find(matchedToErase.begin(),matchedToErase.end(),i) == matchedToErase.end()){
      newMached.push_back(matched[i]);
    }
  }
  
  matched = newMached;
  
  if(verbosityLevel > 1){
    cout<<"\n\nMatched clusters after step 4\n\n"<<endl;
    for(auto cluster : matched){
      cluster->Print();
    }
    cout<<endl;
  }
  
  if(draw) DrawMatched(matched);
}


void ClusterMatcher::MatchClustersClosest(vector<MatchedClusters*> &matched,
                                          vector<unique_ptr<RecHits>> &recHitsPerCluster,
                                          vector<unique_ptr<RecHits>> &simHitsPerCluster)
{
  vector<double> Xs,Ys,Rs, Es;
  double maxDistance = ConfigurationManager::Instance()->GetMachingMaxDistance();
  
  
  for(uint recClusterIndex=0;recClusterIndex < recHitsPerCluster.size();recClusterIndex++){
    // Build new mached cluster for each rec cluster
    MatchedClusters *matchedCluster = new MatchedClusters();
    matchedCluster->AddRecCluster(recClusterIndex, recHitsPerCluster[recClusterIndex]);
    matched.push_back(matchedCluster);
    
    // Save coordinates of this rec cluster
    BasicCluster *basicCluster = matchedCluster->GetRecClusterByIndex(recClusterIndex);
    Xs.push_back(basicCluster->GetX());
    Ys.push_back(basicCluster->GetY());
    Rs.push_back(basicCluster->GetRadius());
    Es.push_back(basicCluster->GetEnergy());
  }
  
  
  for(uint simClusterIndex=0;simClusterIndex<simHitsPerCluster.size();simClusterIndex++){
    // Get coordinates of this sim cluster
    unique_ptr<MatchedClusters> matchedTmp = unique_ptr<MatchedClusters>(new MatchedClusters());
    matchedTmp->AddSimCluster(simClusterIndex, simHitsPerCluster[simClusterIndex]);
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
      matched[parentRecCluster]->AddSimCluster(simClusterIndex, simHitsPerCluster[simClusterIndex]);
    }
  }
}

void ClusterMatcher::MatchClustersAllToAll(vector<MatchedClusters*> &matched,
                                           vector<unique_ptr<RecHits>> &recHitsPerCluster,
                                           vector<unique_ptr<RecHits>> &simHitsPerCluster)
{
  for(uint recClusterIndex=0;recClusterIndex < recHitsPerCluster.size();recClusterIndex++){
    for(uint simClusterIndex=0;simClusterIndex < simHitsPerCluster.size();simClusterIndex++){

      MatchedClusters *matchedCluster = new MatchedClusters();
      matchedCluster->AddRecCluster(recClusterIndex, recHitsPerCluster[recClusterIndex]);
      matchedCluster->AddSimCluster(simClusterIndex, simHitsPerCluster[simClusterIndex]);
      matched.push_back(matchedCluster);
    }
  }
}

void ClusterMatcher::DrawMatched(std::vector<MatchedClusters*> &matched)
{
  vector<TGraph*> recGraph;
  TGraph *recHitsGraph = new TGraph();
  TGraph *simGraph = new TGraph();
  TGraph *simGraphDouble = new TGraph();
  
  int iPointRec=0;
  int iPointSim=0;
  int iPointSimDouble=0;
  int iPointRecHits=0;
  int iCluster=0;
  double scale = 35;
  
  for(auto cluster : matched){
    recGraph.push_back(new TGraph());
    
    std::unique_ptr<RecHits> recHits = std::unique_ptr<RecHits>(new RecHits());
    std::unique_ptr<RecHits> simHits = std::unique_ptr<RecHits>(new RecHits());
    
    cluster->GetRecHits(recHits);
    cluster->GetSimHits(simHits);
    
    recGraph[iCluster]->SetPoint(iPointRec++, cluster->GetRecX(), cluster->GetRecY());
    
    recGraph[iCluster]->SetMarkerSize(scale*cluster->GetRecRadius());
    recGraph[iCluster]->SetMarkerColor(kRed);
    recGraph[iCluster]->SetMarkerStyle(24);
    
    for(int i=0;i<recHits->N();i++){
      recHitsGraph->SetPoint(iPointRecHits++, recHits->GetX()->at(i), recHits->GetY()->at(i));
    }
    for(int i=0;i<simHits->N();i++){
      double x = simHits->GetX()->at(i);
      double y = simHits->GetY()->at(i);
      bool alreadyIn = false;
      
      for(int n=0;n<simGraph->GetN();n++){
        double xx,yy;
        simGraph->GetPoint(n, xx, yy);
        if(xx == x && yy == y){
          alreadyIn = true;
        }
        
      }
      if(alreadyIn) simGraphDouble->SetPoint(iPointSimDouble++, x, y);
      else          simGraph->SetPoint(iPointSim++, x, y);
    }
  }

  recHitsGraph->SetMarkerSize(0.15*scale);
  recHitsGraph->SetMarkerStyle(25);
  recHitsGraph->SetMarkerColor(kRed);
  
  simGraph->SetMarkerSize(0.05*scale);
  simGraph->SetMarkerStyle(20);
  simGraph->SetMarkerColor(kGreen+2);
  
  simGraphDouble->SetMarkerSize(0.10*scale);
  simGraphDouble->SetMarkerStyle(24);
  simGraphDouble->SetMarkerColor(kGreen+2);
  
  TCanvas *c1 = new TCanvas("c1","c1",1800,1000);
  c1->cd();
  
  simGraph->Draw("AP");
//  simGraph->GetXaxis()->SetLimits(-56, -48);
//  simGraph->GetYaxis()->SetRangeUser(44,60);
  
  simGraphDouble->Draw("Psame");
  recHitsGraph->Draw("Psame");
  
  for(int i=0;i<recGraph.size();i++){
    recGraph[i]->Draw("Psame");
  }
  
  c1->Update();
}


