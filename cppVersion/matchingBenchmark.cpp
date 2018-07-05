#include "Event.hpp"
#include "RecHitCalibration.hpp"
#include "ImagingAlgo.hpp"
#include "Helpers.hpp"
#include "ConfigurationManager.hpp"

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>

#include <string>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <utility>
#include <memory>

using namespace std;

int main(int argc, char* argv[])
{
  if(argc != 2){
    cout<<"Usage: matchingBenchmark path_to_config"<<endl;
    exit(0);
  }
  string configPath(argv[1]);
  ConfigurationManager *config = ConfigurationManager::Instance(configPath);
  
  gROOT->ProcessLine(".L loader.C+");

  std::system(("mkdir -p "+config->GetOutputPath()).c_str());

  ImagingAlgo *algo = new ImagingAlgo();
  
  for(int nTupleIter=config->GetMinNtuple();nTupleIter<=config->GetMaxNtuple();nTupleIter++){
    cout<<"\nCurrent ntup: "<<nTupleIter<<endl;

    TFile *inFile = TFile::Open(Form("%s%i.root",config->GetInputPath().c_str(),nTupleIter));
    TTree *tree = (TTree*)inFile->Get("ana/hgc");
    long long nEvents = tree->GetEntries();
    
    cout<<"\n\nLoading ntuple...";
    unique_ptr<Event> hgCalEvent(new Event(tree));
    cout<<"n entries:"<<nEvents<<endl;
    // start event loop
    for(int iEvent=0;iEvent<nEvents;iEvent++){
      if(iEvent>config->GetMaxEventsPerTuple()) break;
      
      hgCalEvent->GoToEvent(iEvent);

      // check if particles reached EE
      bool skipEvent = false;
      for(auto reachedEE : *(hgCalEvent->GetGenParticles()->GetReachedEE())){
        if(reachedEE==0){
          skipEvent = true;
          break;
        }
      }
      if(skipEvent) continue;

      string eventDir = config->GetOutputPath()+"/ntup"+to_string(nTupleIter)+"/event"+to_string(iEvent);
      std::system(("mkdir -p "+eventDir).c_str());

      cout<<"\nCurrent event:"<<iEvent<<endl;

      shared_ptr<RecHits> recHitsRaw = hgCalEvent->GetRecHits();
      shared_ptr<SimClusters> simClusters = hgCalEvent->GetSimClusters();

      // get simulated hits associated with a cluster
      vector<RecHits*> simHitsPerClusterArray;
      recHitsRaw->GetHitsPerSimCluster(simHitsPerClusterArray, simClusters, config->GetEnergyMin());
      
      // re-run clustering with HGCalAlgo, save to file
      std::vector<shared_ptr<Hexel>> recClusters;
      algo->getRecClusters(recClusters, recHitsRaw);
//      cout<<"num of rechits clustered with imaging algo:"<<recClusters.size()<<endl;
      
      // recClusters -> array of hexel objects
      vector<RecHits*> recHitsPerClusterArray;
      recHitsRaw->GetRecHitsPerHexel(recHitsPerClusterArray, recClusters, config->GetEnergyMin());
      
      // perform final analysis, fill in histograms and save to files
      TH2D *energyComparisonNoMatchingHist = new TH2D("no matching","no matching",100,0,100,100,0,100);
      TH2D *energyComparisonClosestHist = new TH2D("closest rec cluster","closest rec cluster",100,0,100,100,0,100);

      for(int layer=config->GetMinLayer();layer<config->GetMaxLayer();layer++){
        vector<int> alreadyAssociatedClusters;
        vector<RecHits*> hitsMatchedToRecClusters;
        vector<double> Xs,Ys,Rs, Es, Eassoc;
        bool skip = false;
        
        for(uint recClusterIndex=0;recClusterIndex<recHitsPerClusterArray.size();recClusterIndex++){

          RecHits *recCluster = recHitsPerClusterArray[recClusterIndex];
          unique_ptr<RecHits> recHitsInLayerInCluster = recCluster->GetHitsInLayer(layer);
          skip = false;
          if(recHitsInLayerInCluster->N()==0){
            skip = true;
            continue;
          }
          double recEnergy = recHitsInLayerInCluster->GetTotalEnergy();
          double xMaxRec   = recHitsInLayerInCluster->GetXmax();
          double xMinRec   = recHitsInLayerInCluster->GetXmin();
          double yMaxRec   = recHitsInLayerInCluster->GetYmax();
          double yMinRec   = recHitsInLayerInCluster->GetYmin();

          double recClusterX = xMinRec+(xMaxRec-xMinRec)/2.;
          double recClusterY = yMinRec+(yMaxRec-yMinRec)/2.;
          double recClusterR = max((xMaxRec-xMinRec)/2.,(yMaxRec-yMinRec)/2.);

//          double assocSimEnergy = 0;

          Xs.push_back(recClusterX);
          Ys.push_back(recClusterY);
          Rs.push_back(recClusterR);
          Es.push_back(recEnergy);
          Eassoc.push_back(0.0);
          
          for(uint simClusterIndex=0;simClusterIndex<simHitsPerClusterArray.size();simClusterIndex++){
            RecHits *simCluster = simHitsPerClusterArray[simClusterIndex];
            unique_ptr<RecHits> simHitsInLayerInCluster = simCluster->GetHitsInLayer(layer);
            
            if(simHitsInLayerInCluster->N()==0) continue;
            double simEnergy = simHitsInLayerInCluster->GetTotalEnergy();
            
            if(recEnergy*simEnergy != 0){
              energyComparisonNoMatchingHist->Fill(recEnergy,simEnergy);
            }
          }
          
        }
        if(skip) continue;
        for(uint simClusterIndex=0;simClusterIndex<simHitsPerClusterArray.size();simClusterIndex++){
          RecHits *simCluster = simHitsPerClusterArray[simClusterIndex];
          unique_ptr<RecHits> simHitsInLayerInCluster = simCluster->GetHitsInLayer(layer);
          
          if(simHitsInLayerInCluster->N()==0) continue;
          
          double simEnergy = simHitsInLayerInCluster->GetTotalEnergy();
          double xMaxSim   = simHitsInLayerInCluster->GetXmax();
          double xMinSim   = simHitsInLayerInCluster->GetXmin();
          double yMaxSim   = simHitsInLayerInCluster->GetYmax();
          double yMinSim   = simHitsInLayerInCluster->GetYmin();
          
          double simClusterX = xMinSim+(xMaxSim-xMinSim)/2.;
          double simClusterY = yMinSim+(yMaxSim-yMinSim)/2.;
          // double simClusterR = max((xMaxSim-xMinSim)/2.,(yMaxSim-yMinSim)/2.);
          
          int parentRecCluster = findClosestCircle(Xs, Ys, Rs, simClusterX, simClusterY);
          
          if(Es[parentRecCluster]*simEnergy != 0){
            energyComparisonNoMatchingHist->Fill(Es[parentRecCluster],simEnergy);
          }
          
          if(parentRecCluster < 0){
            cout<<"No rec cluster found for a sim cluster!!"<<endl;
            continue;
          }
          
          
          Eassoc[parentRecCluster] += simEnergy;
        }
        
        for(uint recClusterIndex=0;recClusterIndex<recHitsPerClusterArray.size();recClusterIndex++){
          
          RecHits *recCluster = recHitsPerClusterArray[recClusterIndex];
          unique_ptr<RecHits> recHitsInLayerInCluster = recCluster->GetHitsInLayer(layer);
          
          if(recHitsInLayerInCluster->N()==0) continue;
          
          double recEnergy = recHitsInLayerInCluster->GetTotalEnergy();
          if(recEnergy*Eassoc[recClusterIndex] != 0){
            energyComparisonClosestHist->Fill(recEnergy,Eassoc[recClusterIndex]);
          }
          
        }
      }
      
      
      energyComparisonNoMatchingHist->SaveAs(Form("%s/energyComparisonNoMatchingHist.root",eventDir.c_str()));
      energyComparisonClosestHist->SaveAs(Form("%s/energyComparisonClosestHist.root",eventDir.c_str()));
      
      recHitsPerClusterArray.clear();
      simHitsPerClusterArray.clear();

    }
    delete tree;
    inFile->Close();
    delete inFile;
  }
  delete algo;

  return 0;
}
