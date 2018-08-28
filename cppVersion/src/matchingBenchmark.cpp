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
  ClusterMatcher *matcher = new ClusterMatcher();
  
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
      if(config->GetReachedEEonly()){
        bool skipEvent = false;
        for(auto reachedEE : *(hgCalEvent->GetGenParticles()->GetReachedEE())){
          if(reachedEE==0){
            skipEvent = true;
            break;
          }
        }
        if(skipEvent) continue;
      }
      string eventDir = config->GetOutputPath()+"/ntup"+to_string(nTupleIter)+"/event"+to_string(iEvent);
      std::system(("mkdir -p "+eventDir).c_str());

      cout<<"\nCurrent event:"<<iEvent<<endl;

      shared_ptr<RecHits> recHitsRaw = hgCalEvent->GetRecHits();
      shared_ptr<SimClusters> simClusters = hgCalEvent->GetSimClusters();

      // get simulated hits associated with a cluster
      vector<RecHits*> simHitsPerClusterArray;
      recHitsRaw->GetHitsPerSimCluster(simHitsPerClusterArray, simClusters);
      
      // re-run clustering with HGCalAlgo, save to file
      std::vector<shared_ptr<Hexel>> recClusters;
      algo->getRecClusters(recClusters, recHitsRaw);
//      cout<<"num of rechits clustered with imaging algo:"<<recClusters.size()<<endl;
      
      // recClusters -> array of hexel objects
      vector<RecHits*> recHitsPerClusterArray;
      recHitsRaw->GetRecHitsPerHexel(recHitsPerClusterArray, recClusters);
      
      // perform final analysis, fill in histograms and save to files
      TH2D *energyComparisonNoMatchingHist = new TH2D("no matching","no matching",500,0,100,500,0,100);
      TH2D *energyComparisonClosestHist = new TH2D("closest rec cluster","closest rec cluster",500,0,100,500,0,100);

      
      int nZeroSim = 0;
      for(int layer=config->GetMinLayer();layer<config->GetMaxLayer();layer++){
        
        vector<MatchedClusters*> unmatchedClusters;
        vector<MatchedClusters*> matchedClusters;
        
        matcher->MatchClustersByDetID(matchedClusters,recHitsPerClusterArray,simHitsPerClusterArray,layer);
        matcher->MatchClustersAllToAll(unmatchedClusters,recHitsPerClusterArray,simHitsPerClusterArray,layer);
        
        for(MatchedClusters *clusters : unmatchedClusters){
            energyComparisonNoMatchingHist->Fill(clusters->GetRecEnergy(),
                                                 clusters->GetSimEnergy());
        }
        
        for(MatchedClusters *clusters : matchedClusters){
           if(!clusters->HasSimClusters()) continue;
            energyComparisonClosestHist->Fill(clusters->GetRecEnergy(),
                                              clusters->GetSimEnergy());
          
          if(clusters->GetSimEnergy() < 0.0001) nZeroSim++;
          
//          cout<<"E sim:"<<clusters->GetTotalSimEnergy()<<"\tErec:"<<clusters->recCluster->GetEnergy()<<endl;
        }
      }
      
      cout<<"\nN zero sim:"<<nZeroSim<<"\n\n"<<endl;
      
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
  delete matcher;
  return 0;
}
