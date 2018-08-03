//
//  Chromosome.cpp
//
//  Created by Jeremi Niedziela on 02/08/2018.
//

#include "Chromosome.hpp"
#include "ConfigurationManager.hpp"
#include "ImagingAlgo.hpp"
#include "ClusterMatcher.hpp"
#include "Event.hpp"
#include "Helpers.hpp"

#include <TMath.h>
#include <iostream>
#include <TROOT.h>
#include <TH1D.h>
#include <TFile.h>

using namespace std;

Chromosome::Chromosome()
{
  uniqueID = reinterpret_cast<uint64_t>(this);
  configPath = "tmp/config_"+to_string(uniqueID)+".md";
  clusteringOutputPath = "tmp/output_"+to_string(uniqueID)+".txt";
  executionTime = 9999999;
  score = -999999;
} 

Chromosome::~Chromosome()
{
  
}

Chromosome* Chromosome::GetRandom()
{
  Chromosome *result = new Chromosome();
  
  result->SetCriticalDistanceEE(RandFloat(0.0, 100.0));
  result->SetCriticalDistanceFH(RandFloat(0.0, 100.0));
  result->SetCriticalDistanceBH(RandFloat(0.0, 100.0));
  result->SetDependSensor(RandBool());
  result->SetReachedEE(RandBool());
  result->SetKernel(RandInt(0, 2));
  result->SetDeltacEE(RandFloat(0.0, 100.0));
  result->SetDeltacFH(RandFloat(0.0, 100.0));
  result->SetDeltacBH(RandFloat(0.0, 100.0));
  result->SetKappa(RandFloat(0.0, 10000.0));
  result->SetEnergyMin(RandFloat(0.0, 100.0));
  result->SetMatchingDistance(RandFloat(0.0, 100.0));
  result->SetMinClusters(RandInt(0, 10));
  
  return result;
}

void Chromosome::SaveToBitChromosome()
{
  int currentShift = 0;
  bitChromosome[0] = 0;
  bitChromosome[1] = 0;
  bitChromosome[2] = 0;
  
  // load values to 1st chromosome
  ShiftIntoChromosome(criticalDistanceEE,currentShift,0);
  ShiftIntoChromosome(criticalDistanceFH,currentShift,0);
  ShiftIntoChromosome(criticalDistanceBH,currentShift,0);
  ShiftIntoChromosome(dependSensor,currentShift,0);
  ShiftIntoChromosome(reachedEE,currentShift,0);
  
  // load values to 2st chromosome
  currentShift = 0;
  ShiftIntoChromosome(kernel,currentShift,1);
  ShiftIntoChromosome(deltacEE,currentShift,1);
  ShiftIntoChromosome(deltacFH,currentShift,1);
  ShiftIntoChromosome(deltacBH,currentShift,1);
  
  // load values to 3st chromosome
  currentShift = 0;
  ShiftIntoChromosome(kappa,currentShift,2);
  ShiftIntoChromosome(energyMin,currentShift,2);
  ShiftIntoChromosome(matchingDistance,currentShift,2);
  ShiftIntoChromosome(minClusters,currentShift,2);
  
}

void Chromosome::ReadFromBitChromosome()
{
  int currentShift = 0;
  
  // read 1st chromosome
  SetValueFromChromosome(criticalDistanceEE, currentShift, 0);
  SetValueFromChromosome(criticalDistanceFH, currentShift, 0);
  SetValueFromChromosome(criticalDistanceBH, currentShift, 0);
  SetValueFromChromosome(dependSensor, currentShift, 0);
  SetValueFromChromosome(reachedEE, currentShift, 0);
  
  // read 1st chromosome
  currentShift = 0;
  SetValueFromChromosome(kernel, currentShift, 1);
  SetValueFromChromosome(deltacEE, currentShift, 1);
  SetValueFromChromosome(deltacFH, currentShift, 1);
  SetValueFromChromosome(deltacBH, currentShift, 1);
  
  // read 2st chromosome
  currentShift = 0;
  SetValueFromChromosome(kappa, currentShift, 2);
  SetValueFromChromosome(energyMin, currentShift, 2);
  SetValueFromChromosome(matchingDistance, currentShift, 2);
  SetValueFromChromosome(minClusters, currentShift, 2);
}

void Chromosome::Print()
{
  cout<<"================================================="<<endl;
  cout<<"Chromosome "<<uniqueID<<endl;
  
  cout<<"Critical distance:"<<endl;
  cout<<"\tEE:"<<criticalDistanceEE/1000.<<endl;
  cout<<"\tFH:"<<criticalDistanceFH/1000.<<endl;
  cout<<"\tBH:"<<criticalDistanceBH/1000.<<endl;
  
  cout<<"Depend on sensor:"<<dependSensor<<endl;
  cout<<"Reached EE only:"<<reachedEE<<endl;
  
  cout<<"Kernel: ";
  if(kernel==0)       cout<<"step"<<endl;
  else if(kernel==1)  cout<<"gaus"<<endl;
  else                cout<<"exp"<<endl;
  
  cout<<"Critical #delta:"<<endl;
  cout<<"\tEE:"<<deltacEE/1000.<<endl;
  cout<<"\tFH:"<<deltacFH/1000.<<endl;
  cout<<"\tBH:"<<deltacBH/1000.<<endl;
  
  cout<<"kappa:"<<kappa/1000.<<endl;
  cout<<"energy threshold:"<<energyMin/100000.<<endl;
  cout<<"max matching distance:"<<matchingDistance/1000.<<endl;
  cout<<"min clusters:"<<minClusters<<endl;
  
  clusteringOutput.Print();
  cout<<"execution time:"<<executionTime<<endl;
  cout<<"score:"<<score<<endl;
  cout<<"================================================="<<endl;
}


template<class T>
void Chromosome::ShiftIntoChromosome(T value, int &shift, int chromoIndex)
{
  uint64_t mask = (uint64_t)value << shift;
  bitChromosome[chromoIndex] |= mask;
  shift += BitSize(value);
}

template<class T>
void Chromosome::SetValueFromChromosome(T &value, int &shift, int chromoIndex)
{
  uint64_t mask = 0;
  for(int i=0;i<BitSize(value);i++){mask |= 1ull << (i+shift);}
  value = (bitChromosome[chromoIndex] & mask) >> shift;
  shift += BitSize(value);
}


void Chromosome::StoreInConfig()
{
  system(("cp baseConfig.md "+configPath).c_str());
  
  UpdateParamValue(configPath, "depend_sensor",GetDependSensor());
  if(GetKernel() == 0)      UpdateParamValue(configPath, "energy_density_function","step");
  else if(GetKernel() == 1) UpdateParamValue(configPath, "energy_density_function","gaus");
  else                      UpdateParamValue(configPath, "energy_density_function","exp");
  
  UpdateParamValue(configPath, "critial_distance_EE",GetCriticalDistanceEE());
  UpdateParamValue(configPath, "critial_distance_FH",GetCriticalDistanceFH());
  UpdateParamValue(configPath, "critial_distance_BH",GetCriticalDistanceBH());
  UpdateParamValue(configPath, "deltac_EE",GetDeltacEE());
  UpdateParamValue(configPath, "deltac_FH",GetDeltacFH());
  UpdateParamValue(configPath, "deltac_BH",GetDeltacBH());
  UpdateParamValue(configPath, "kappa",GetKappa());
  UpdateParamValue(configPath, "energy_min",GetEnergyMin());
  UpdateParamValue(configPath, "min_clusters",GetMinClusters());
  UpdateParamValue(configPath, "reachedEE_only",GetReachedEE());
  UpdateParamValue(configPath, "matching_max_distance",GetMatchingDistance());
  UpdateParamValue(configPath, "score_output_path",clusteringOutputPath);
}

void Chromosome::RunClustering()
{
  cout<<"Running clusterization"<<endl;
  
  auto start = now();
//  Clusterize(configPath);
  system(("./createQualityPlots "+configPath+" > /dev/null 2>&1").c_str());
//  system(("./createQualityPlots "+configPath).c_str());
  auto end = now();
  executionTime = duration(start,end);
  clusteringOutput = ReadOutput(clusteringOutputPath);
  
  system(("rm "+configPath).c_str());
  system(("rm "+clusteringOutputPath).c_str());
  
  cout<<"Done. Execution time:"<<executionTime<<endl;
}

void Chromosome::CalculateScore()
{
  StoreInConfig();
  RunClustering();
  
  score =   clusteringOutput.containmentMean              // we want high containment
          - 0.1*fabs(clusteringOutput.containmentSigma)   // with small spread (but not that important)
          - fabs(clusteringOutput.resolutionMean)         // with good resolution
          - 0.1*fabs(clusteringOutput.resolutionSigma)    // also with small spread (with lower weight)
          - clusteringOutput.separationMean               // small separation factor
          - 0.1*fabs(clusteringOutput.separationSigma);   // with small spread
  
  if(executionTime > 60){
    score -= (executionTime-30); // add additional penalty for super long execution
  }
  
}

Chromosome* Chromosome::ProduceChildWith(Chromosome *partner)
{
  Chromosome *child = new Chromosome();
  
  // combine chromosomes of parents in a random way
  for(int i=0;i<3;i++){
    int crossingPoint = RandInt(0, 63);
    uint64_t newBitChromosome = 0;
    
    uint64_t maskA = 0;
    for(int j=0;j<BitSize(maskA)-crossingPoint;j++){maskA |= 1ull << (j+crossingPoint);}
    
    uint64_t maskB = 0;
    for(int j=(int)BitSize(maskB)-crossingPoint;j<BitSize(maskB);j++){maskB |= 1ull << j;}
    
    newBitChromosome = bitChromosome[0] & maskA;
    newBitChromosome |= partner->GetBitChromosome(0) & maskB;
    
    child->SetBitChromosome(i, newBitChromosome);
  }
  
  
  
  child->ReadFromBitChromosome();
  return child;
}


void Chromosome::Clusterize(string configPath)
{
  ConfigurationManager *config = ConfigurationManager::Instance(configPath);
  gROOT->ProcessLine(".L loader.C+");
  
  ImagingAlgo *algo = new ImagingAlgo();
  ClusterMatcher *matcher = new ClusterMatcher();
  
  TH1D *deltaE = new TH1D("Erec-Esim/Esim","Erec-Esim/Esim",100,-1.0,2.0);
  TH1D *separation = new TH1D("separation","separation",100,0,10);
  TH1D *containment = new TH1D("separation","separation",100,0,10);
  
  for(int nTupleIter=config->GetMinNtuple();nTupleIter<=config->GetMaxNtuple();nTupleIter++){

    TFile *inFile = TFile::Open(Form("%s%i.root",config->GetInputPath().c_str(),nTupleIter));
    if(!inFile) continue;
    TTree *tree = (TTree*)inFile->Get("ana/hgc");
    if(!tree) continue;
    long long nEvents = tree->GetEntries();
    cout<<"n entries:"<<nEvents<<endl;
    
    unique_ptr<Event> hgCalEvent(new Event(tree));
    
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
      shared_ptr<RecHits> recHitsRaw = hgCalEvent->GetRecHits();
      shared_ptr<SimClusters> simClusters = hgCalEvent->GetSimClusters();
      
      // get simulated hits associated with a cluster
      vector<RecHits*> simHitsPerClusterArray;
      recHitsRaw->GetHitsPerSimCluster(simHitsPerClusterArray, simClusters);
      
      // re-run clustering with HGCalAlgo
      std::vector<shared_ptr<Hexel>> recClusters;
      algo->getRecClusters(recClusters, recHitsRaw);
      
      vector<RecHits*> recHitsPerClusterArray;
      recHitsRaw->GetRecHitsPerHexel(recHitsPerClusterArray, recClusters);
      
      
      // perform final analysis, fill in histograms and save to files
      for(int layer=config->GetMinLayer();layer<config->GetMaxLayer();layer++){
        
        vector<MatchedClusters*> matchedClusters;
        matcher->MatchClustersByDetID(matchedClusters,recHitsPerClusterArray,simHitsPerClusterArray,layer);
        if(matchedClusters.size()==0) continue;
        
        for(MatchedClusters *clusters : matchedClusters){
          if(clusters->simClusters->size() == 0) continue;
          
          double recEnergy = clusters->GetTotalRecEnergy();
          double simEnergy = clusters->GetTotalSimEnergy();
          
          deltaE->Fill((recEnergy-simEnergy)/simEnergy);
          containment->Fill(clusters->GetSharedFraction());
        }
        
        for(uint i=0;i<matchedClusters.size();i++){
          for(uint j=(i+1);j<matchedClusters.size();j++){
            
            BasicCluster *recCluster1 = matchedClusters[i]->GetMergedRecCluster();
            BasicCluster *recCluster2 = matchedClusters[j]->GetMergedRecCluster();
            
            double distance = sqrt(pow(recCluster1->GetX()-recCluster2->GetX(),2)+
                                   pow(recCluster1->GetY()-recCluster2->GetY(),2));
            
            double sigma1 = recCluster1->GetRadius();
            double sigma2 = recCluster2->GetRadius();
            
            if(sigma1+sigma2 == 0) continue;
            
            separation->Fill(distance/(sigma1+sigma2));
          }
        }
      }
      recHitsPerClusterArray.clear();
      simHitsPerClusterArray.clear();
    }
    
    ofstream outputFile;
    outputFile.open(config->GetScoreOutputPath());
    outputFile<<deltaE->GetMean()<<endl;
    outputFile<<deltaE->GetStdDev()<<endl;
    outputFile<<separation->GetMean()<<endl;
    outputFile<<separation->GetStdDev()<<endl;
    outputFile<<containment->GetMean()<<endl;
    outputFile<<containment->GetStdDev()<<endl;
    outputFile.close();
    
    delete tree;
    inFile->Close();
    delete inFile;
    
  }
  delete algo;
  delete matcher;
}



