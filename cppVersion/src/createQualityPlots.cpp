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
  TH1D *separation = new TH1D("separation","separation",500,0,5);
  TH1D *containment = new TH1D("containment","containment",200,-1,1);
  TH1D *deltaN = new TH1D("numberClusters","numberClusters",2000,-10,10);
  int emptyMatchedClusters = 0;
  int zeroSizeClusters = 0;
  int noMatchedClusters = 0;
  int nTotalMatchedClusters = 0;
  int nTotalLayerEvents = 0;
  
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
      
      cout<<"\nCurrent event:"<<iEvent<<"\n\n"<<endl;
      
      auto genParticles = hgCalEvent->GetGenParticles();
      
      for(int iGen=0;iGen<genParticles->N();iGen++){
        if(config->GetVerbosityLevel() > 0){
          genParticles->Print(iGen);
        }
      }
      cout<<endl;
      
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
      
      string eventDir = config->GetOutputPath()+"/ntup"+to_string(nTupleIter)+"/event"+to_string(iEvent);
      std::system(("mkdir -p "+eventDir).c_str());
      
      
      
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
      
      if(recClusters.size() == 0){
        cout<<"ERROR - algorithm couldn't find any rec clusters!!!"<<endl;
        continue;
      }
      
      vector<RecHits*> recHitsPerClusterArray;
      recHitsRaw->GetRecHitsPerHexel(recHitsPerClusterArray, recClusters);
    
      if(config->GetVerbosityLevel() > 0){
        cout<<"\nReconstructed hits grouped by clusters:"<<endl;
        for(RecHits *hits : recHitsPerClusterArray){
          hits->Print();
        }
      }
      
      
      // perform final analysis, fill in histograms and save to files
      TH2D *ErecEsimVsEta = new TH2D("ErecEsim vs. eta","ErecEsim vs. eta",100,1.5,3.2,100,0,2.5);
      TH2D *sigmaEvsEta = new TH2D("sigma(E) vs. eta","sigma(E) vs. eta",100,1.5,3.2,100,-10,10);
      TH2D *sigmaEvsEtaEsim = new TH2D("sigma(E)Esim vs. eta","sigma(E)Esim vs. eta",100,1.5,3.2,100,-1.5,1.0);
      TH1D *twoSeparation = new TH1D("two clusters separation","two clusters separation",1000,0,10);
      TH1D *twoSeparationJer = new TH1D("two clusters separation","two clusters separation",1000,0,10);
      TH2D *NrecNsim = new TH2D("NrecNsim","NrecNsim",20,1,20,20,1,20);
      
      double totalSimEnergy = 0;
      double totalRecEnergy = 0;
      
      
      double total3DsimEnergy = 0;
      
      for(int i=0;i<simClusters->N();i++){
        total3DsimEnergy += simClusters->GetEnergy(i);
      }
      
      for(int layer=config->GetMinLayer();layer<config->GetMaxLayer();layer++){
        
        if(simClusters->GetNsimClustersInLayer(layer) == 0){
          // this means that there were no 2d clusters simulated in this layer, so there's nothing to look for there (in the future one should check number of fake clusters reconstructed even in layers where there was nothing simulated)
          continue;
        }
        
        vector<MatchedClusters*> matchedClusters;
        matcher->MatchClustersByDetID(matchedClusters,recHitsPerClusterArray,simHitsPerClusterArray,layer);
        
        nTotalMatchedClusters+=matchedClusters.size();
        nTotalLayerEvents++;
        
        if(matchedClusters.size()==0){
          noMatchedClusters++;
          continue;
        }
        
        for(MatchedClusters *clusters : matchedClusters){
          if(clusters->GetRecRadius() == 0){
            zeroSizeClusters++;
            // maybe skip those?
          }
          
          if(!clusters->HasSimClusters() ||
             !clusters->HasRecClusters()){
            emptyMatchedClusters++;
            continue;
          }
          double recEnergy = clusters->GetRecEnergy();
          double recEta    = clusters->GetRecEta();
          double simEnergy = clusters->GetSimEnergy();
          
          totalRecEnergy += recEnergy;
          totalSimEnergy += simEnergy;

          ErecEsimVsEta->Fill(fabs(recEta),recEnergy/simEnergy);
          sigmaEvsEta->Fill(fabs(recEta),(recEnergy-simEnergy)/recEnergy);
          sigmaEvsEtaEsim->Fill(fabs(recEta),(recEnergy-simEnergy)/simEnergy);
          deltaE->Fill((recEnergy-simEnergy)/simEnergy);
          containment->Fill(clusters->GetSharedFraction());
        }
        NrecNsim->Fill(simHitsPerClusterArray.size(),recHitsPerClusterArray.size());
        deltaN->Fill(((int)simHitsPerClusterArray.size()-(int)recHitsPerClusterArray.size())/(double)simHitsPerClusterArray.size());
        
        for(uint i=0;i<matchedClusters.size();i++){
          for(uint j=(i+1);j<matchedClusters.size();j++){
            
            double distance = sqrt(pow(matchedClusters[i]->GetRecX()-matchedClusters[j]->GetRecX(),2)+
                                   pow(matchedClusters[i]->GetRecY()-matchedClusters[j]->GetRecY(),2));
            
            double sigma1 = matchedClusters[i]->GetRecRadius();
            double sigma2 = matchedClusters[j]->GetRecRadius();
            
            if(sigma1+sigma2 == 0){
              continue;
            }
            
            twoSeparation->Fill(distance/sqrt(sigma1*sigma1+sigma2*sigma2));
            twoSeparationJer->Fill(distance/(sigma1+sigma2));
            separation->Fill(distance/(sigma1+sigma2));
          }
        }
      }
      
      ErecEsimVsEta->SaveAs(Form("%s/ErecEsimVsEta.root",eventDir.c_str()));
      sigmaEvsEta->SaveAs(Form("%s/simgaEVsEta.root",eventDir.c_str()));
      sigmaEvsEtaEsim->SaveAs(Form("%s/simgaEVsEtaEsim.root",eventDir.c_str()));
      twoSeparation->SaveAs(Form("%s/twoSeparation.root",eventDir.c_str()));
      twoSeparationJer->SaveAs(Form("%s/twoSeparationJer.root",eventDir.c_str()));
      NrecNsim->SaveAs(Form("%s/NrecNsim.root",eventDir.c_str()));
      
      recHitsPerClusterArray.clear();
      simHitsPerClusterArray.clear();
    }
    inFile->Close();
    delete inFile;
  }
  cout<<endl<<endl;
  
  deltaE->SaveAs(Form("%s/resolution.root",config->GetOutputPath().c_str()));
  separation->SaveAs(Form("%s/separation.root",config->GetOutputPath().c_str()));
  containment->SaveAs(Form("%s/containment.root",config->GetOutputPath().c_str()));
  deltaN->SaveAs(Form("%s/deltaN.root",config->GetOutputPath().c_str()));
  
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
  outputFile<<emptyMatchedClusters/(double)(nTotalMatchedClusters)<<endl;
  cout<<"\% of event-layers with empty (sim or rec) clusters in matched clusters:"<<emptyMatchedClusters/(double)(nTotalMatchedClusters)<<endl;
  outputFile<<zeroSizeClusters/(double)(nTotalMatchedClusters)<<endl;
  cout<<"\% of event-layers with zero size matched clusters:"<<zeroSizeClusters/(double)(nTotalMatchedClusters)<<endl;
  outputFile<<noMatchedClusters/(double)(nTotalMatchedClusters)<<endl;
  cout<<"N events with no matched clusters:"<<noMatchedClusters/(double)(nTotalLayerEvents)<<endl;
  
  outputFile.close();
  
  
  delete algo;
  delete matcher;
  return 0;
}

