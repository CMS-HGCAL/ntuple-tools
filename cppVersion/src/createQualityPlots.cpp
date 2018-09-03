//
//  createQualityPlots.cpp
//
//  Created by Jeremi Niedziela on 28/06/2018.
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

// specify here an event to stop at and plot matched clusters (for debugging)
int plotEvent = -1;
int plotLayer = -1;

enum EMonitor1D{
  kResolution,
  kSeparation,
  kContainment,
  kDeltaNclusters,
  kFakeVsLayer,
  kNmonitors1D
};

enum EMonitor2D{
  kResolutionVsEta,
  kErecVsEsimUnmatched,
  kErecVsEsimDetIdMatching,
  kNrecVsNsim,
  kNmonitors2D
};

// Limits of 1D monitoring histograms {nBinsX, minX, maxX}
double monitorLimits1D[kNmonitors1D][3] = {
  {1000,-5,5},    // resolution
  {1000,0,10},    // separation
  {1000,-1,1},    // containment
  {2000,-10,10},  // ΔN_clusters
  {100,0,100},    // fake rec clusters vs. layer
};

// Limits of 2D monitoring histograms {nBinsX, minX, maxX, nBinsY, minY, maxY}
double monitorLimits2D[kNmonitors2D][6] = {
  {100,1.5,3.2,100,-1.5,1.0}, // resolution vs. eta
  {500,0,100,500,0,100},       // Erec vs. Esim (unmatched)
  {500,0,100,500,0,100},      // Erec vs. Esim (matching by detID)
  {20,1,20,20,1,20},          // Nrec vs. Nsim
};

const char* monitorNames1D[kNmonitors1D] = {
  "resolution",
  "separation",
  "containment",
  "deltaNclusters",
  "fakeVsLayer",
};

const char* monitorNames2D[kNmonitors2D] = {
  "resolutionVsEta",
  "ErecVsEsimUnmatched",
  "ErecVsEsimDetIdMatching",
  "NrecVsNsim",
};

int main(int argc, char* argv[])
{
  if((argc != 2) && (argc != 22)){
    cout<<"Usage: createQualityPlots path_to_config"<<endl;
    exit(0);
  }
  
  ConfigurationManager *config = nullptr;
  
  if(argc == 2){
    string configPath(argv[1]);
    config = ConfigurationManager::Instance(configPath);
  }
  if(argc == 22){
    config = ConfigurationManager::Instance(atoi(argv[1]), // depend sensor
                                            argv[2], // input path
                                            argv[3], // output path
                                            atof(argv[4]), // deltac EE
                                            atof(argv[5]), // deltac FH
                                            atof(argv[6]), // deltac BH
                                            atof(argv[7]), // energy threshold
                                            atof(argv[8]), // crit dist EE
                                            atof(argv[9]),// crit dist FH
                                            atof(argv[10]),// crit dist BH
                                            atof(argv[11]),// kappa
                                            atoi(argv[12]),// verbosity
                                            atoi(argv[13]),// min n tuple
                                            atoi(argv[14]),// max n tuple
                                            atoi(argv[15]),// min layer
                                            atoi(argv[16]),// max layer
                                            atoi(argv[17]),// events per tuple
                                            argv[18],// energy density function
                                            atoi(argv[19]),// reached EE only
                                            atof(argv[20]),// matching distance
                                            argv[21] // score output path
                                            );
    
  }
  
  cout<<endl;config->Print();
  
  
  gROOT->ProcessLine(".L loader.C+");
  TApplication theApp("App", &argc, argv);
  
  std::system(("mkdir -p "+config->GetOutputPath()).c_str());
  
  ImagingAlgo *algo = new ImagingAlgo();
  ClusterMatcher *matcher = new ClusterMatcher();
  
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
  
  int failure1=0, failure2=0, failure3=0;
  int nTotalLayerEvents = 0;
  int nTotalMatchedClusters = 0;
  
  for(int nTupleIter=config->GetMinNtuple();nTupleIter<=config->GetMaxNtuple();nTupleIter++){
    cout<<"\nCurrent ntup: "<<nTupleIter<<endl;
    
    TFile *inFile = TFile::Open(Form("%s%i.root",config->GetInputPath().c_str(),nTupleIter));
    if(!inFile) continue;
    TTree *tree = (TTree*)inFile->Get("ana/hgc");
    if(!tree) continue;
    long long nEvents = tree->GetEntries();
    cout<<"n entries:"<<nEvents<<endl;
    
    unique_ptr<Event> hgCalEvent(new Event(tree));
//    gDebug = 2;
    
    // start event loop
    for(int iEvent=0;iEvent<nEvents;iEvent++){
      if(iEvent>config->GetMaxEventsPerTuple()) break;
      
      hgCalEvent->GoToEvent(iEvent);
      
      cout<<"Current event:"<<iEvent<<endl;
      
      auto genParticles = hgCalEvent->GetGenParticles();
      
      if(config->GetVerbosityLevel() > 0){
        for(int iGen=0;iGen<genParticles->N();iGen++){
          genParticles->Print(iGen);
        }
        cout<<endl;
      }
      
      // check if particles reached EE
      if(config->GetReachedEEonly()){
        bool skipEvent = false;
        for(auto reachedEE : *(genParticles->GetReachedEE())){
          if(reachedEE==0){
            skipEvent = true;
            break;
          }
        }
        if(skipEvent){
          if(config->GetVerbosityLevel()>0){
            cout<<"\n\nSkipping event because there were particles that converted before reaching EE\n\n"<<endl;
          }
          continue;
        }
        if(config->GetVerbosityLevel()>0){
          cout<<"\n\nEvent OK - all particles reached EE\n\n"<<endl;
        }
        
      }
      
      shared_ptr<RecHits> recHitsRaw = hgCalEvent->GetRecHits();
      shared_ptr<SimClusters> simClusters = hgCalEvent->GetSimClusters();
      
      // get simulated hits associated with a cluster
      vector<unique_ptr<RecHits>> simHitsPerClusterArray;
      recHitsRaw->GetHitsPerSimCluster(simHitsPerClusterArray, simClusters);
      
      if(config->GetVerbosityLevel() > 0){
        cout<<"Simulated hits grouped by clusters:"<<endl;
        for(auto &hits : simHitsPerClusterArray){
          hits->Print();
        }
      }
      
      // re-run clustering with HGCalAlgo
      std::vector<shared_ptr<Hexel>> recClusters;
      algo->getRecClusters(recClusters, recHitsRaw);
      
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
        
        // Take sim clusters in this layer, check if there are any
        vector<unique_ptr<RecHits>> simHitsInClusterInLayer;
        for(auto &cluster : simHitsPerClusterArray){
          unique_ptr<RecHits> clusterInLayer = cluster->GetHitsInLayer(layer);
          if(clusterInLayer->N() == 0) continue;
          simHitsInClusterInLayer.push_back(move(clusterInLayer));
        }
        if(simHitsInClusterInLayer.size() == 0){
          // this is OK, somethimes there may be nothing simulated in given layer, for instance when particles stop before the end of the calorimeter
          if(config->GetVerbosityLevel() > 1){
            cout<<"No sim clusters in layer:"<<layer<<endl;
          }
          continue;
        }
        
        // Some sim clusters will share a fraction of hits. Try to re-assign shared hits energy.
        // This will not work very well if single hit is share by 3 or more clusters...
//        for(int i=0;i<simHitsInClusterInLayer.size();i++){
//          for(int j=i+1;j<simHitsInClusterInLayer.size();j++){
//            simHitsInClusterInLayer[i]->ShareCommonHits(simHitsInClusterInLayer[j]);
//          }
//        }
        
        // If this event-layer makes sense at all (there are some sim clusters in it), count it for the denominator of failures fraction
        nTotalLayerEvents++;
        
        // Take rec clusters in this layer, check if there are any
        vector<unique_ptr<RecHits>> recHitsInClusterInLayer;
        for(RecHits *cluster : recHitsPerClusterArray){
          unique_ptr<RecHits> clusterInLayer = cluster->GetHitsInLayer(layer);
          if(clusterInLayer->N() == 0) continue;
          recHitsInClusterInLayer.push_back(move(clusterInLayer));
        }
        if(recHitsInClusterInLayer.size() == 0){
          // This is not good, one should minimize possibility for that to happen. This means that there were no rec clusters found despite there were some sim clusters in this layer
          if(config->GetVerbosityLevel() > 1){
            cout<<"No rec clusters in layer:"<<layer<<" despite there are "<<simHitsInClusterInLayer.size()<<" sim clusters in this layer!!"<<endl;
          }
          failure1++;
          continue;
        }
        
        monitors2D[kNrecVsNsim]->Fill(simHitsPerClusterArray.size(),recHitsPerClusterArray.size());
        monitors1D[kDeltaNclusters]->Fill(((int)simHitsPerClusterArray.size()-(int)recHitsPerClusterArray.size())/(double)simHitsPerClusterArray.size());
        
        // Match rec clusters with sim clusters by det ID, check if there is at least one such pair
        vector<MatchedClusters*> matchedClusters;
        vector<MatchedClusters*> unmatchedClusters;
        
        bool draw=false;
        if(iEvent==plotEvent && layer==plotLayer) draw = true;
        
        matcher->MatchClustersByDetID(matchedClusters,recHitsInClusterInLayer,simHitsInClusterInLayer,draw);
        matcher->MatchClustersAllToAll(unmatchedClusters,recHitsInClusterInLayer,simHitsInClusterInLayer);
        
        if(iEvent==plotEvent && layer==plotLayer) goto finish; // psss, don't tell anyone...
        
        if(matchedClusters.size() == 0){
          // This is bad, because there are some sim and some rec clusters in this layer, but not even a single pair sim-rec was found
          if(config->GetVerbosityLevel() > 1){
            cout<<"Zero matched clusters created in layer:"<<layer<<endl;
            cout<<"Sim clusters:"<<endl;
            for(auto &cluster : simHitsInClusterInLayer){cluster->Print();}
            cout<<"Rec clusters:"<<endl;
            for(auto &cluster : recHitsInClusterInLayer){cluster->Print();}
          }
          failure2++;
          continue;
        }
        
        for(MatchedClusters *clusters : matchedClusters){
          nTotalMatchedClusters++;
          
          if(!clusters->HasSimClusters()){
            // This may happen when fake rec cluster was found and there are no simulated clusters around. One should minimize probability of that to happen
            if(config->GetVerbosityLevel() > 1){
              cout<<"\n\nNo sim clusters in matched clusters in layer:"<<layer<<endl;
            }
            monitors1D[kFakeVsLayer]->Fill(layer);
            failure3++;
            continue;
          }
          if(!clusters->HasRecClusters()){
            // Should never happen, as we check earlier if there are rec clusters in this layer and matched clusters are created from rec clusters
            cout<<"ERROR -- No rec clusters in matched clusters in layer:"<<layer<<endl;
            continue;
          }
          if(clusters->GetRecRadius() == 0){
            // Should never happen, as we set radius of 0.5 cm for single-hit clusters
            cout<<"ERROR -- Rec cluster with size zero"<<endl;
            continue;
          }
          
          double recEnergy = clusters->GetRecEnergy();
          double recEta    = clusters->GetRecEta();
          double simEnergy = clusters->GetSimEnergy();
        
          monitors2D[kResolutionVsEta]->Fill(fabs(recEta),(recEnergy-simEnergy)/simEnergy);
          monitors1D[kResolution]->Fill((recEnergy-simEnergy)/simEnergy);
          monitors1D[kContainment]->Fill(clusters->GetSharedFraction());
          monitors2D[kErecVsEsimDetIdMatching]->Fill(clusters->GetRecEnergy(),clusters->GetSimEnergy());
        }
        
        for(MatchedClusters *clusters : unmatchedClusters){
          monitors2D[kErecVsEsimUnmatched]->Fill(clusters->GetRecEnergy(),clusters->GetSimEnergy());
        }
        
        for(uint i=0;i<matchedClusters.size();i++){
          if(!matchedClusters[i]->HasRecClusters()) continue;
          for(uint j=(i+1);j<matchedClusters.size();j++){
            if(!matchedClusters[j]->HasRecClusters()) continue;
            
            double distance = sqrt(pow(matchedClusters[i]->GetRecX()-matchedClusters[j]->GetRecX(),2)+
                                   pow(matchedClusters[i]->GetRecY()-matchedClusters[j]->GetRecY(),2));
            
            double sigma1 = matchedClusters[i]->GetRecRadius();
            double sigma2 = matchedClusters[j]->GetRecRadius();
            
            if(sigma1+sigma2 == 0){
              // This should basically never happen, for single-hit clusters we set radius of 0.5 cm
              cout<<"sum zero size"<<endl;
              continue;
            }
            
            if(distance/(sigma1+sigma2) > 10) monitors1D[kSeparation]->Fill(9.9);
            else                              monitors1D[kSeparation]->Fill(distance/(sigma1+sigma2));
          }
        }
      }
      
      recHitsPerClusterArray.clear();
      simHitsPerClusterArray.clear();
    }
    inFile->Close();
    delete inFile;
  }
  cout<<endl<<endl;
  
finish:
  
  // Create output file
  string outpath = config->GetOutputPath();
  ofstream outputFile;
  cout<<"writing output to:"<<config->GetScoreOutputPath()<<endl;
  outputFile.open(config->GetScoreOutputPath());
  
  // Save 1D monitors and put mean values in a text file
  for(int i=0;i<kNmonitors1D;i++){
    monitors1D[i]->SaveAs((outpath+"/"+monitorNames1D[i]+".root").c_str());
    double mean, sigma;
    
    if(monitors1D[i] && monitors1D[i]->GetEntries()>0 && monitors1D[i]->GetStdDev() != 0){
      mean = monitors1D[i]->GetMean();
      sigma= monitors1D[i]->GetStdDev();
    }
    else{
      mean = sigma = 999999;
    }
    outputFile<<mean<<endl;
    outputFile<<sigma<<endl;
    cout<<"Average "<<monitorNames1D[i]<<" per event:"<<mean<<endl;
    cout<<monitorNames1D[i]<<" sigma:"<<sigma<<endl;
  }
  
  // Save 2D monitors
  for(int i=0;i<kNmonitors2D;i++){
      monitors2D[i]->SaveAs((outpath+"/"+monitorNames2D[i]+".root").c_str());
  }
  
  // Save fraction of fails in the text file
  outputFile<<failure1/(double)(nTotalLayerEvents)<<endl;
  cout<<"\% of event-layers where algo failed to find sim clusters:"<<failure1/(double)(nTotalLayerEvents)<<endl;
  
  outputFile<<failure2/(double)(nTotalLayerEvents)<<endl;
  cout<<"\% of event-layers were rec and sim clusters couldn't be matched:"<<failure2/(double)(nTotalLayerEvents)<<endl;
  
  outputFile<<failure3/(double)(nTotalMatchedClusters)<<endl;
  cout<<"\% of fake rec clusters (those that don't match any sim cluster):"<<failure3/(double)(nTotalMatchedClusters)<<endl;
  
  outputFile.close();
  
  // copy file with results also to the location where histograms are stored
  system(("cp "+config->GetScoreOutputPath()+" "+outpath+"/output.txt").c_str());
  
  delete algo;
  delete matcher;
  
  theApp.Run();
  return 0;
}

