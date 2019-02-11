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
#include <TH3D.h>
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
  kTotalEnergyMIP,
  kNmonitors1D
};

enum EMonitor2D{
  kPositionXY,
  kNmonitors2D
};

enum EMonitor3D{
  kPositionXYZ,
  kTracksImpactXYZ,
  kNmonitors3D
};

// Limits of 1D monitoring histograms {nBinsX, minX, maxX}
double monitorLimits1D[kNmonitors1D][3] = {
  {1000,0,10000},    // total energy
};

// Limits of 2D monitoring histograms {nBinsX, minX, maxX, nBinsY, minY, maxY}
double monitorLimits2D[kNmonitors2D][6] = {
  {100,-50,50,100,-50,50}, // position XY
};

// Limits of 3D monitoring histograms {nBinsX, minX, maxX, nBinsY, minY, maxY, nBinsZ, minZ, maxZ}
double monitorLimits3D[kNmonitors3D][9] = {
  {100,-50,50,100,-50,50,100,-50,50}, // position XYZ
  {100,-50,50,100,-50,50,100,-50,50}, // tracks impact points XYZ
};

const char* monitorNames1D[kNmonitors1D] = {
  "totalEnergyMIP",
};

const char* monitorNames2D[kNmonitors2D] = {
  "positionXY",
};

const char* monitorNames3D[kNmonitors3D] = {
  "positionXYZ",
  "tracksImpactXYZ",
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
  TH3D *monitors3D[kNmonitors3D];
  
  for(int i=0;i<kNmonitors1D;i++){
    monitors1D[i] = new TH1D(monitorNames1D[i],monitorNames1D[i],
                             monitorLimits1D[i][0],monitorLimits1D[i][1],monitorLimits1D[i][2]);
  }
  
  for(int i=0;i<kNmonitors2D;i++){
    monitors2D[i] = new TH2D(monitorNames2D[i],monitorNames2D[i],
                             monitorLimits2D[i][0],monitorLimits2D[i][1],monitorLimits2D[i][2],
                             monitorLimits2D[i][3],monitorLimits2D[i][4],monitorLimits2D[i][5]);
  }
  
  for(int i=0;i<kNmonitors3D;i++){
    monitors3D[i] = new TH3D(monitorNames3D[i],monitorNames3D[i],
                             monitorLimits3D[i][0],monitorLimits3D[i][1],monitorLimits3D[i][2],
                             monitorLimits3D[i][3],monitorLimits3D[i][4],monitorLimits3D[i][5],
                             monitorLimits3D[i][6],monitorLimits3D[i][7],monitorLimits3D[i][8]);
  }
  
  TTree *outputTree = new TTree("tree","tree");
  double slopeX, slopeY, offsetX, offsetY;
  
  outputTree->Branch("slopeX", &slopeX);
  outputTree->Branch("slopeY", &slopeY);
  outputTree->Branch("offsetX", &offsetX);
  outputTree->Branch("offsetY", &offsetY);
  
  for(int nTupleIter=config->GetMinNtuple();nTupleIter<=config->GetMaxNtuple();nTupleIter++){
    cout<<"\nCurrent ntup: "<<nTupleIter<<endl;
    
    TFile *inFile = TFile::Open(Form("%s%i.root",config->GetInputPath().c_str(),nTupleIter));
    if(!inFile) continue;
    TTree *tree = (TTree*)inFile->Get("rechitntupler/hits");
    if(!tree) continue;
    TTree *testBeamTree = (TTree*)inFile->Get("trackimpactntupler/impactPoints");
    if(!testBeamTree) continue;
    
    long long nEvents = tree->GetEntries();
    cout<<"n entries:"<<nEvents<<endl;
    
    unique_ptr<Event> hgCalEvent(new Event(tree, testBeamTree));
    
    // start event loop
    for(int iEvent=0;iEvent<nEvents;iEvent++){
//      if(iEvent != 2) continue;
      
      if(iEvent>config->GetMaxEventsPerTuple()) break;
      
      hgCalEvent->GoToEvent(iEvent);
      
      if(!hgCalEvent->IsTestBeam()){
        cout<<"The event doesn't come from the testbeam of the ntuple structure changed significantly."<<endl;
        cout<<"I don't know how to process this event..."<<endl;
        continue;
      }
      if(iEvent%100 == 0){
        cout<<"Current event:"<<iEvent<<endl;
      }
      if(config->GetVerbosityLevel() > 0){
        hgCalEvent->Print();
      }
      
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
      
      double avgX = 0;
      double avgY = 0;
      int iClusters = 0;
      for(int layer=config->GetMinLayer();layer<config->GetMaxLayer();layer++){
        vector<unique_ptr<RecHits>> recHitsInClusterInLayer;
        
        for(RecHits *cluster : recHitsPerClusterArray){
          unique_ptr<RecHits> clusterInLayer = cluster->GetHitsInLayer(layer);
          if(clusterInLayer->N() == 0) continue;
          recHitsInClusterInLayer.push_back(move(clusterInLayer));
        }
        for(unique_ptr<RecHits> &cluster : recHitsInClusterInLayer){
          iClusters++;
          avgX += (cluster->GetXmax() - cluster->GetXmin())/2.;
          avgY += (cluster->GetYmax() - cluster->GetYmin())/2.;
        }
      }
      avgX /= iClusters;
      avgY /= iClusters;
      
//      cout<<"avg x:"<<avgX<<"\ty:"<<avgY<<endl;
      
      double totalEnergyMIP = 0;
      for(auto cluster : rec3Dclusters){
        if(config->GetVerbosityLevel() > 0){
          cluster->Print();
        }
        totalEnergyMIP += cluster->GetEnergy();
        monitors2D[kPositionXY]->Fill(cluster->GetX(),cluster->GetY());
        monitors3D[kPositionXYZ]->Fill(cluster->GetX(),cluster->GetY(),cluster->GetZ());
      }
      
      shared_ptr<TestbeamTrack> track = hgCalEvent->GetTestbeamTrack();
      
      for(int iLayer=0;iLayer<28;iLayer++){
        monitors3D[kTracksImpactXYZ]->Fill(track->GetX(iLayer),
                                           track->GetY(iLayer),
                                           iLayer);
      }
      monitors1D[kTotalEnergyMIP]->Fill(totalEnergyMIP);
      
      slopeX  = track->GetSlopeX();
      offsetX = track->GetOffsetX();
      slopeY  = track->GetSlopeY();
      offsetY = track->GetOffsetY();
      
      outputTree->Fill();
      
      recHitsPerClusterArray.clear();
      
    }
    inFile->Close();
    delete inFile;
  }
  
  // Save monitors
  
  string outpath = config->GetOutputPath();
  TFile *outFile = new TFile((outpath+"/results.root").c_str(),"recreate");
  outFile->cd();
  for(int i=0;i<kNmonitors1D;i++){monitors1D[i]->Write();}
  for(int i=0;i<kNmonitors2D;i++){monitors2D[i]->Write();}
  for(int i=0;i<kNmonitors3D;i++){monitors3D[i]->Write();}
  outputTree->Write();
  outFile->Close();
  
  delete algo;
//  theApp.Run();
  
  return 0;
}

