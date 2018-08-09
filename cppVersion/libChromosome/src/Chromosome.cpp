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
#include "GeneticHelpers.hpp"

#include <TMath.h>
#include <iostream>
#include <TROOT.h>
#include <TH1D.h>
#include <TFile.h>

using namespace std;

#define mutationChance 0.002
#define critialExecutionTime 30.0
#define executionTimeout 50

Chromosome::Chromosome() :
criticalDistanceEE(0.0),
criticalDistanceFH(0.0),
criticalDistanceBH(0.0),
dependSensor(false),
reachedEE(false),
kernel(0),
deltacEE(0.0),
deltacFH(0.0),
deltacBH(0.0),
kappa(0.0),
energyMin(0.0),
matchingDistance(0.0),
minClusters(0)
{
  int nChromosomes = 6;
  for(int i=0;i<nChromosomes;i++){bitChromosome.push_back(0);}
  
  clusteringOutput.resolutionMean = 999999;
  clusteringOutput.resolutionSigma = 999999;
  clusteringOutput.separationMean = 999999;
  clusteringOutput.separationSigma = 999999;
  clusteringOutput.containmentMean = 999999;
  clusteringOutput.containmentSigma = 999999;
    
  uniqueID = reinterpret_cast<uint64_t>(this);
  configPath = "tmp/config_"+to_string(uniqueID)+".md";
  clusteringOutputPath = "tmp/output_"+to_string(uniqueID)+".txt";
  executionTime = 99999;
  score = -99999;
} 

Chromosome::~Chromosome()
{
  
}

Chromosome* Chromosome::GetRandom()
{
  Chromosome *result = new Chromosome();
  
  result->SetCriticalDistanceEE(RandFloat(0.0, criticalDistanceEEmax));
  result->SetCriticalDistanceFH(RandFloat(0.0, criticalDistanceFHmax));
  result->SetCriticalDistanceBH(RandFloat(0.0, criticalDistanceBHmax));
  result->SetDependSensor(RandBool());
  result->SetReachedEE(RandBool());
  result->SetKernel(RandInt(0, 2));
  result->SetDeltacEE(RandFloat(0.0, deltacEEmax));
  result->SetDeltacFH(RandFloat(0.0, deltacFHmax));
  result->SetDeltacBH(RandFloat(0.0, deltacBHmax));
  result->SetKappa(RandFloat(0.0, kappaMax));
  result->SetEnergyMin(RandFloat(0.0, energyThresholdMax));
  result->SetMatchingDistance(RandFloat(0.0, matchingDistanceMax));
  result->SetMinClusters(RandInt(0, 10));
  
  return result;
}

void Chromosome::SaveToBitChromosome()
{
  for(int i=0;i<bitChromosome.size();i++){bitChromosome[i] = 0;}
  
  // load values to 1st chromosome
  int currentShift = 0;
  int chromoIndex = 0;
  ShiftIntoChromosome(criticalDistanceEE,currentShift,chromoIndex);
  ShiftIntoChromosome(criticalDistanceFH,currentShift,chromoIndex);
  
  currentShift = 0;
  chromoIndex++;
  ShiftIntoChromosome(criticalDistanceBH,currentShift,chromoIndex);
  ShiftIntoChromosome(dependSensor,currentShift,chromoIndex);
  ShiftIntoChromosome(reachedEE,currentShift,chromoIndex);
  
  
  currentShift = 0;
  chromoIndex++;
  ShiftIntoChromosome(kernel,currentShift,chromoIndex);
  ShiftIntoChromosome(deltacEE,currentShift,chromoIndex);
  
  currentShift = 0;
  chromoIndex++;
  ShiftIntoChromosome(deltacFH,currentShift,chromoIndex);
  ShiftIntoChromosome(deltacBH,currentShift,chromoIndex);
  
  currentShift = 0;
  chromoIndex++;
  ShiftIntoChromosome(kappa,currentShift,chromoIndex);
  ShiftIntoChromosome(energyMin,currentShift,chromoIndex);
  
  currentShift = 0;
  chromoIndex++;
  ShiftIntoChromosome(matchingDistance,currentShift,chromoIndex);
  ShiftIntoChromosome(minClusters,currentShift,chromoIndex);
}

void Chromosome::ReadFromBitChromosome()
{
  // read 1st chromosome
  int currentShift = 0;
  int chromoIndex = 0;
  SetValueFromChromosome(criticalDistanceEE, currentShift, chromoIndex);
  SetValueFromChromosome(criticalDistanceFH, currentShift, chromoIndex);
  
  currentShift = 0;
  chromoIndex++;
  SetValueFromChromosome(criticalDistanceBH, currentShift, chromoIndex);
  SetValueFromChromosome(dependSensor, currentShift, chromoIndex);
  SetValueFromChromosome(reachedEE, currentShift, chromoIndex);

  currentShift = 0;
  chromoIndex++;
  SetValueFromChromosome(kernel, currentShift, chromoIndex);
  SetValueFromChromosome(deltacEE, currentShift, chromoIndex);
  
  currentShift = 0;
  chromoIndex++;
  SetValueFromChromosome(deltacFH, currentShift, chromoIndex);
  SetValueFromChromosome(deltacBH, currentShift, chromoIndex);
  
  currentShift = 0;
  chromoIndex++;
  SetValueFromChromosome(kappa, currentShift, chromoIndex);
  SetValueFromChromosome(energyMin, currentShift, chromoIndex);
  
  currentShift = 0;
  chromoIndex++;
  SetValueFromChromosome(matchingDistance, currentShift, chromoIndex);
  SetValueFromChromosome(minClusters, currentShift, chromoIndex);
}

void Chromosome::Print()
{
  cout<<"================================================="<<endl;
  cout<<"Chromosome "<<uniqueID<<endl;
  
  cout<<"Critical distance:";
  cout<<"\tEE:"<<criticalDistanceEE/1000.;
  cout<<"\tFH:"<<criticalDistanceFH/1000.;
  cout<<"\tBH:"<<criticalDistanceBH/1000.<<endl;
  
  cout<<"Depend on sensor:"<<dependSensor<<endl;
  cout<<"Reached EE only:"<<reachedEE<<endl;
  
  cout<<"Kernel: ";
  if(kernel==0)       cout<<"step"<<endl;
  else if(kernel==1)  cout<<"gaus"<<endl;
  else                cout<<"exp"<<endl;
  
  cout<<"Critical #delta:";
  cout<<"\tEE:"<<deltacEE/1000.;
  cout<<"\tFH:"<<deltacFH/1000.;
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
//  cout<<"Running clusterization"<<endl;
  
  auto start = now();
//  Clusterize(configPath);
  system(("./createQualityPlots "+configPath+" > /dev/null 2>&1").c_str());
  
//  system(("./execWithTimeout.sh -t "+to_string(executionTimeout)+" ./createQualityPlots "+configPath+" > /dev/null 2>&1").c_str());
  
//  system(("./createQualityPlots "+configPath).c_str());
  auto end = now();
  executionTime = duration(start,end);
  clusteringOutput = ReadOutput(clusteringOutputPath);
  
  system(("rm "+configPath).c_str());
  system(("rm "+clusteringOutputPath).c_str());
  
//  cout<<"Done. Execution time:"<<executionTime<<endl;
}

void Chromosome::CalculateScore()
{
  StoreInConfig();
  RunClustering();
  
  double distance =     fabs(clusteringOutput.containmentMean-1)
                      +      clusteringOutput.containmentSigma
                      + fabs(clusteringOutput.resolutionMean)
                      +      clusteringOutput.resolutionSigma
                      +      clusteringOutput.separationMean
                      +      clusteringOutput.separationSigma;
  
  score = 1./distance;
  
  if(clusteringOutput.resolutionMean > 1000){ // this means that clustering failed completely
    score = 0;
  }
  
  // the time measurement we have now is not realiable, can't be used to punish population members...
//  if(executionTime > critialExecutionTime){
//    score -= pow((executionTime-critialExecutionTime),2); // add additional penalty for super long execution
//  }
}

Chromosome* Chromosome::ProduceChildWith(Chromosome *partner)
{
  Chromosome *child = new Chromosome();
  
  // combine chromosomes of parents in a random way
  // this is a single-point crossover
  for(int i=0;i<bitChromosome.size();i++){
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
  
  // perform child genes mutation
  for(int i=0;i<bitChromosome.size();i++){
    uint64_t bits = child->GetBitChromosome(i);
    
    for(int iBit=0;iBit<BitSize(bits);iBit++){
      double r = RandDouble(0, 1);
      if(r < mutationChance){ // 10% chance for mutation
        ReverseBit(bits, iBit);
      }
    }
    child->SetBitChromosome(i, bits);
  }
  
  child->ReadFromBitChromosome();
  return child;
}


void Chromosome::Clusterize(string configPath)
{
//  ConfigurationManager *config = ConfigurationManager::Instance(configPath);
  
  ImagingAlgo *algo = new ImagingAlgo(configPath);
  ClusterMatcher *matcher = new ClusterMatcher();
  
  TH1D *deltaE = new TH1D("Erec-Esim/Esim","Erec-Esim/Esim",100,-1.0,2.0);
  TH1D *separation = new TH1D("separation","separation",100,0,10);
  TH1D *containment = new TH1D("separation","separation",100,0,10);
  
  int minNtuple, maxNtuple, eventsPerTuple, minLayer, maxLayer;
  string inputPath, scoreOutputPath;
  bool reachedEEonly;
  
  GetParamFomeConfig(configPath, "min_Ntuple", minNtuple);
  GetParamFomeConfig(configPath, "max_Ntuple", maxNtuple);
  GetParamFomeConfig(configPath, "input_path", inputPath);
  GetParamFomeConfig(configPath, "analyze_events_per_tuple", eventsPerTuple);
  GetParamFomeConfig(configPath, "reachedEE_only", reachedEEonly);
  GetParamFomeConfig(configPath, "min_layer", minLayer);
  GetParamFomeConfig(configPath, "max_layer", maxLayer);
  GetParamFomeConfig(configPath, "score_output_path", scoreOutputPath);
  
  for(int nTupleIter=minNtuple; nTupleIter<=maxNtuple; nTupleIter++){

    TFile *inFile = TFile::Open(Form("%s%i.root",inputPath.c_str(),nTupleIter));
    if(!inFile) continue;
    TTree *tree = (TTree*)inFile->Get("ana/hgc");
    if(!tree) continue;
    long long nEvents = tree->GetEntries();
    cout<<"n entries:"<<nEvents<<endl;
    
    unique_ptr<Event> hgCalEvent(new Event(tree));
    
    // start event loop
    for(int iEvent=0;iEvent<nEvents;iEvent++){
      if(iEvent>eventsPerTuple) break;
      
      hgCalEvent->GoToEvent(iEvent);
      
      // check if particles reached EE
      if(reachedEEonly){
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
      for(int layer=minLayer;layer<maxLayer;layer++){
        
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
    outputFile.open(scoreOutputPath);
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



