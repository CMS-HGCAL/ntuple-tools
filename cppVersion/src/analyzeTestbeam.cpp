//
//  analyzeTestbeam.cpp
//
//  Created by Jeremi Niedziela on 30/11/2018.
//

#include "Event.hpp"
#include "RecHitCalibration.hpp"
#include "ImagingAlgo.hpp"
#include "Helpers.hpp"
#include "ConfigurationManager.hpp"
#include "ClusterMatcher.hpp"
#include "MatchedClusters.hpp"

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TApplication.h>

#include <string>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <utility>
#include <memory>
#include <fstream>

using namespace std;

enum EMonitor1D{
  kResolution,
  kSeparation,
  kContainment,
  kNmonitors1D
};

enum EMonitor2D{
  kResolutionVsEta,
  kNmonitors2D
};

// Limits of 1D monitoring histograms {nBinsX, minX, maxX}
double monitorLimits1D[kNmonitors1D][3] = {
  {1000,-5,5},    // resolution
  {1000,0,10},    // separation
  {1000,-1,1},    // containment
};

// Limits of 2D monitoring histograms {nBinsX, minX, maxX, nBinsY, minY, maxY}
double monitorLimits2D[kNmonitors2D][6] = {
  {100,1.5,3.2,100,-1.5,1.0}, // resolution vs. eta
};

const char* monitorNames1D[kNmonitors1D] = {
  "resolution",
  "separation",
  "containment",
};

const char* monitorNames2D[kNmonitors2D] = {
  "resolutionVsEta",
};

int main(int argc, char* argv[])
{
  if(argc != 2){
    cout<<"Usage: analyzeTestbeam path_to_config"<<endl;
    exit(0);
  }
  ConfigurationManager *config = ConfigurationManager::Instance(argv[1]);
  cout<<endl;config->Print();
  
  gROOT->ProcessLine(".L loader.C+");
  TApplication theApp("App", &argc, argv);
  
  std::system(("mkdir -p "+config->GetOutputPath()).c_str());
  
  ImagingAlgo *algo = new ImagingAlgo();
  
  TH1D *monitors1D[kNmonitors1D];
  TH2D *monitors2D[kNmonitors2D];
  
  for(int i=0;i<kNmonitors1D;i++){
    monitors1D[i] = new TH1D(monitorNames1D[i],monitorNames1D[i],
                             monitorLimits1D[i][0],monitorLimits1D[i][1],monitorLimits1D[i][2]);
  }
  
  for(int i=0;i<kNmonitors2D;i++){
    monitors2D[i] = new TH2D(monitorNames2D[i],monitorNames2D[i],
                             monitorLimits2D[i][0],monitorLimits2D[i][1],monitorLimits2D[i][2],
                             monitorLimits2D[i][3],monitorLimits2D[i][4],monitorLimits2D[i][5]);
  }
  
  
  for(int nTupleIter=config->GetMinNtuple();nTupleIter<=config->GetMaxNtuple();nTupleIter++){
    cout<<"\nCurrent ntup: "<<nTupleIter<<endl;
    
    TFile *inFile = TFile::Open(Form("%s%i.root",config->GetInputPath().c_str(),nTupleIter));
    if(!inFile) continue;
    TTree *tree = (TTree*)inFile->Get("rechitntupler/hits");
    if(!tree) continue;
    
    long long nEvents = tree->GetEntries();
    cout<<"n entries:"<<nEvents<<endl;
    
    unique_ptr<Event> hgCalEvent(new Event(tree));
    
    // start event loop
    for(int iEvent=0;iEvent<nEvents;iEvent++){
      if(iEvent>config->GetMaxEventsPerTuple()) break;
      
      hgCalEvent->GoToEvent(iEvent);
      if(!hgCalEvent->IsTestBeam()){
        cout<<"The event doesn't come from the testbeam of the ntuple structure changed significantly."<<endl;
        cout<<"I don't know how to process this event..."<<endl;
        continue;
      }
      cout<<"Current event:"<<iEvent<<endl;
      
      // re-run clustering with HGCalAlgo
      shared_ptr<RecHits> recHitsRaw = hgCalEvent->GetRecHits();
      vector<shared_ptr<Hexel>> recClusters;
      algo->getRecClusters(recClusters, recHitsRaw);
      
      vector<vector<vector<unique_ptr<Hexel>>>> clusters;
      vector<shared_ptr<BasicCluster>> rec3Dclusters;
      algo->makeClusters(clusters, recHitsRaw);
      algo->make3DClusters(rec3Dclusters, clusters);
      
      vector<RecHits*> recHitsPerClusterArray;
      if(recClusters.size()>0){
        recHitsRaw->GetRecHitsPerHexel(recHitsPerClusterArray, recClusters);
      }
      
      if(config->GetVerbosityLevel() > 0){
        cout<<"\nReconstructed hits grouped by clusters:"<<endl;
        for(RecHits *hits : recHitsPerClusterArray){
          hits->Print();
        }
      }
      
      
      for(int layer=config->GetMinLayer();layer<config->GetMaxLayer();layer++){
        
        // Take rec clusters in this layer, check if there are any
        vector<unique_ptr<RecHits>> recHitsInClusterInLayer;
        
        for(RecHits *cluster : recHitsPerClusterArray){
          unique_ptr<RecHits> clusterInLayer = cluster->GetHitsInLayer(layer);
          if(clusterInLayer->N() == 0) continue;
          recHitsInClusterInLayer.push_back(move(clusterInLayer));
        }

        cout<<"Layer "<<layer<<endl;
        for(auto &cluster : recHitsInClusterInLayer){
          cluster->Print();
        }
      }
      
      for(auto cluster : rec3Dclusters){
        cluster->Print();
      }
      
      recHitsPerClusterArray.clear();
    }
    inFile->Close();
    delete inFile;
  }
  cout<<endl<<endl;
  
  // Create output file
  string outpath = config->GetOutputPath();
  
  // Save monitors
  for(int i=0;i<kNmonitors1D;i++){
    monitors1D[i]->SaveAs((outpath+"/"+monitorNames1D[i]+".root").c_str());
  }
  for(int i=0;i<kNmonitors2D;i++){
    monitors2D[i]->SaveAs((outpath+"/"+monitorNames2D[i]+".root").c_str());
  }
  
  delete algo;
  theApp.Run();
  
  return 0;
}

