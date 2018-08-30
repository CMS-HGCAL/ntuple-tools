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

#include <string>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <utility>
#include <memory>
#include <fstream>

using namespace std;

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
  
  std::system(("mkdir -p "+config->GetOutputPath()).c_str());
  
  ImagingAlgo *algo = new ImagingAlgo();
  ClusterMatcher *matcher = new ClusterMatcher();
  
  TH1D *deltaE = new TH1D("resolution","resolution",1000,-5,5);
  TH1D *separation = new TH1D("separation","separation",1000,0,10);
  TH1D *containment = new TH1D("containment","containment",200,-1,1);
  TH1D *deltaN = new TH1D("numberClusters","numberClusters",2000,-10,10);
  TH1D *clusterRadius = new TH1D("clusterRadius","clusterRadius",3000,-1,300);
  
  TH2D *ErecEsimVsEta = new TH2D("ErecEsim vs. eta","ErecEsim vs. eta",100,1.5,3.2,100,0,2.5);
  TH2D *sigmaEvsEtaEsim = new TH2D("sigma(E)Esim vs. eta","sigma(E)Esim vs. eta",100,1.5,3.2,100,-1.5,1.0);
  TH2D *NrecNsim = new TH2D("NrecNsim","NrecNsim",20,1,20,20,1,20);
  TH2D *energyComparisonNoMatchingHist = new TH2D("no matching","no matching",500,0,100,500,0,100);
  TH2D *energyComparisonClosestHist = new TH2D("closest rec cluster","closest rec cluster",500,0,100,500,0,100);
  
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
      vector<RecHits*> simHitsPerClusterArray;
      recHitsRaw->GetHitsPerSimCluster(simHitsPerClusterArray, simClusters);
      
      if(config->GetVerbosityLevel() > 0){
        cout<<"Simulated hits grouped by clusters:"<<endl;
        for(RecHits *hits : simHitsPerClusterArray){
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
        for(RecHits *cluster : simHitsPerClusterArray){
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
        
        NrecNsim->Fill(simHitsPerClusterArray.size(),recHitsPerClusterArray.size());
        deltaN->Fill(((int)simHitsPerClusterArray.size()-(int)recHitsPerClusterArray.size())/(double)simHitsPerClusterArray.size());
        
        // Match rec clusters with sim clusters by det ID, check if there is at least one such pair
        vector<MatchedClusters*> matchedClusters;
        vector<MatchedClusters*> unmatchedClusters;
        
        matcher->MatchClustersByDetID(matchedClusters,recHitsInClusterInLayer,simHitsInClusterInLayer);
        matcher->MatchClustersAllToAll(unmatchedClusters,recHitsPerClusterArray,simHitsPerClusterArray,layer);
        
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
            failure3++;
            continue;
          }
          if(!clusters->HasRecClusters()){
            // Should never happen, as we check earlier if there are rec clusters in this layer and matched clusters are created from rec clusters
            cout<<"ERROR -- No rec clusters in matched clusters in layer:"<<layer<<endl;
            continue;
          }
          clusterRadius->Fill(clusters->GetRecRadius());
          
          if(clusters->GetRecRadius() == 0){
            // Should never happen, as we set radius of 0.5 cm for single-hit clusters
            cout<<"ERROR -- Rec cluster with size zero"<<endl;
            continue;
          }
          
          double recEnergy = clusters->GetRecEnergy();
          double recEta    = clusters->GetRecEta();
          double simEnergy = clusters->GetSimEnergy();
        
          ErecEsimVsEta->Fill(fabs(recEta),recEnergy/simEnergy);
          sigmaEvsEtaEsim->Fill(fabs(recEta),(recEnergy-simEnergy)/simEnergy);
          deltaE->Fill((recEnergy-simEnergy)/simEnergy);
          containment->Fill(clusters->GetSharedFraction());
          energyComparisonClosestHist->Fill(clusters->GetRecEnergy(),clusters->GetSimEnergy());
        }
        
        for(MatchedClusters *clusters : unmatchedClusters){
          energyComparisonNoMatchingHist->Fill(clusters->GetRecEnergy(),clusters->GetSimEnergy());
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
            
            if(distance/(sigma1+sigma2) > 10) separation->Fill(9.9);
            else                              separation->Fill(distance/(sigma1+sigma2));
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
  
  string outpath = config->GetOutputPath();
  
  deltaE->SaveAs((outpath+"/resolution.root").c_str());
  separation->SaveAs((outpath+"/separation.root").c_str());
  containment->SaveAs((outpath+"/containment.root").c_str());
  deltaN->SaveAs((outpath+"/deltaN.root").c_str());
  clusterRadius->SaveAs((outpath+"/clusterRadius.root").c_str());
  
  ErecEsimVsEta->SaveAs((outpath+"/ErecEsimVsEta.root").c_str());
  sigmaEvsEtaEsim->SaveAs((outpath+"/simgaEVsEtaEsim.root").c_str());
  NrecNsim->SaveAs((outpath+"/NrecNsim.root").c_str());
  energyComparisonNoMatchingHist->SaveAs((outpath+"/energyComparisonNoMatchingHist.root").c_str());
  energyComparisonClosestHist->SaveAs((outpath+"/energyComparisonClosestHist.root").c_str());
  
  ofstream outputFile;
  cout<<"writing output to:"<<config->GetScoreOutputPath()<<endl;
  outputFile.open(config->GetScoreOutputPath());
  
  if(deltaE && deltaE->GetEntries()>0 && deltaE->GetStdDev() != 0){
    outputFile<<deltaE->GetMean()<<endl;
    outputFile<<deltaE->GetStdDev()<<endl;
    cout<<"Average resolution per event:"<<deltaE->GetMean()<<endl;
    cout<<"Resolution sigma:"<<deltaE->GetStdDev()<<endl;
  }
  else{
    outputFile<<999999<<endl;
    outputFile<<999999<<endl;
    cout<<"Average resolution per event:"<<999999<<endl;
    cout<<"Resolution sigma:"<<999999<<endl;
  }
  if(separation && separation->GetEntries()>0 && separation->GetStdDev() != 0){
    outputFile<<separation->GetMean()<<endl;
    outputFile<<separation->GetStdDev()<<endl;
    cout<<"Average separation per event:"<<separation->GetMean()<<endl;
    cout<<"Separation sigma:"<<separation->GetStdDev()<<endl;
  }
  else{
    outputFile<<999999<<endl;
    outputFile<<999999<<endl;
    cout<<"Average separation per event:"<<999999<<endl;
    cout<<"Separation sigma:"<<999999<<endl;
  }
  if(containment && containment->GetEntries()>0 && containment->GetStdDev() != 0){
    outputFile<<containment->GetMean()<<endl;
    outputFile<<containment->GetStdDev()<<endl;
    cout<<"Average containment per event:"<<containment->GetMean()<<endl;
    cout<<"Containment sigma:"<<containment->GetStdDev()<<endl;
  }
  else{
    outputFile<<999999<<endl;
    outputFile<<999999<<endl;
    cout<<"Average containment per event:"<<999999<<endl;
    cout<<"Containment sigma:"<<999999<<endl;
  }
  if(deltaN && deltaN->GetEntries()>0 && deltaN->GetStdDev() != 0){
    outputFile<<deltaN->GetMean()<<endl;
    outputFile<<deltaN->GetStdDev()<<endl;
    cout<<"Average difference in N clusters (sim-rec) per event:"<<deltaN->GetMean()<<endl;
    cout<<"N clusters (sim-rec) sigma:"<<deltaN->GetStdDev()<<endl;
  }
  else{
    outputFile<<999999<<endl;
    outputFile<<999999<<endl;
    cout<<"Average difference in N clusters (sim-rec) per event:"<<999999<<endl;
    cout<<"N clusters (sim-rec) sigma:"<<999999<<endl;
  }
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
  return 0;
}

