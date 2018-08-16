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
#include <TROOT.h>
#include <TH1D.h>
#include <TFile.h>

#include <iostream>
#include <stdio.h>
#include <unistd.h>

using namespace std;

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
  normalizedScore = -999999;
  
  mutationChance = 0.002;
  severityFactor = 1.0;
} 

Chromosome::~Chromosome()
{
  
}

Chromosome* Chromosome::GetRandom()
{
  Chromosome *result = new Chromosome();
  
  result->SetCriticalDistanceEE(RandFloat(criticalDistanceEEmin, criticalDistanceEEmax));
  result->SetCriticalDistanceFH(RandFloat(criticalDistanceFHmin, criticalDistanceFHmax));
  result->SetCriticalDistanceBH(RandFloat(criticalDistanceBHmin, criticalDistanceBHmax));
  result->SetDependSensor(RandBool());
  result->SetReachedEE(RandBool());
  result->SetKernel(RandInt(kernelMin, kernelMax));
  result->SetDeltacEE(RandFloat(deltacEEmin, deltacEEmax));
  result->SetDeltacFH(RandFloat(deltacFHmin, deltacFHmax));
  result->SetDeltacBH(RandFloat(deltacBHmin, deltacBHmax));
  result->SetKappa(RandFloat(kappaMin, kappaMax));
  result->SetEnergyMin( RandFloat(result->GetDependSensor() ? energyThresholdMin : energyThresholdMinNoSensor,
                                  result->GetDependSensor() ? energyThresholdMax : energyThresholdMaxNoSensor));
  result->SetMatchingDistance(RandFloat(matchingDistanceMin, matchingDistanceMax));
  result->SetMinClusters(RandInt(minClustersMin, minClustersMax));
  
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
  ShiftIntoChromosome(criticalDistanceBH,currentShift,chromoIndex);
  ShiftIntoChromosome(dependSensor,currentShift,chromoIndex);
  ShiftIntoChromosome(reachedEE,currentShift,chromoIndex);
  
  
  currentShift = 0;
  chromoIndex++;
  ShiftIntoChromosome(kernel,currentShift,chromoIndex);
  ShiftIntoChromosome(deltacEE,currentShift,chromoIndex);
  ShiftIntoChromosome(deltacFH,currentShift,chromoIndex);
  ShiftIntoChromosome(deltacBH,currentShift,chromoIndex);
  
  currentShift = 0;
  chromoIndex++;
  ShiftIntoChromosome(kappa,currentShift,chromoIndex);
  ShiftIntoChromosome(energyMin,currentShift,chromoIndex);
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
  SetValueFromChromosome(criticalDistanceBH, currentShift, chromoIndex);
  SetValueFromChromosome(dependSensor, currentShift, chromoIndex);
  SetValueFromChromosome(reachedEE, currentShift, chromoIndex);

  currentShift = 0;
  chromoIndex++;
  SetValueFromChromosome(kernel, currentShift, chromoIndex);
  SetValueFromChromosome(deltacEE, currentShift, chromoIndex);
  SetValueFromChromosome(deltacFH, currentShift, chromoIndex);
  SetValueFromChromosome(deltacBH, currentShift, chromoIndex);
  
  currentShift = 0;
  chromoIndex++;
  SetValueFromChromosome(kappa, currentShift, chromoIndex);
  SetValueFromChromosome(energyMin, currentShift, chromoIndex);
  SetValueFromChromosome(matchingDistance, currentShift, chromoIndex);
  SetValueFromChromosome(minClusters, currentShift, chromoIndex);
}

void Chromosome::Print()
{
  cout<<"================================================="<<endl;
  cout<<"Chromosome "<<GetUniqueID()<<endl;
  
  cout<<"Critical distance:";
  cout<<"\tEE:"<<GetCriticalDistanceEE();
  cout<<"\tFH:"<<GetCriticalDistanceFH();
  cout<<"\tBH:"<<GetCriticalDistanceBH()<<endl;
  
  cout<<"Depend on sensor:"<<GetDependSensor()<<endl;
  cout<<"Reached EE only:"<<GetReachedEE()<<endl;
  
  cout<<"Kernel: ";
  if(GetKernel()==0)        cout<<"step"<<endl;
  else if(GetKernel()==1)   cout<<"gaus"<<endl;
  else                      cout<<"exp"<<endl;
  
  cout<<"Critical #delta:";
  cout<<"\tEE:"<<GetDeltacEE();
  cout<<"\tFH:"<<GetDeltacFH();
  cout<<"\tBH:"<<GetDeltacBH()<<endl;
  
  cout<<"kappa:"<<GetKappa()<<endl;
  cout<<"energy threshold:"<<GetEnergyMin()<<endl;
  cout<<"max matching distance:"<<GetMatchingDistance()<<endl;
  cout<<"min clusters:"<<GetMinClusters()<<endl;
  
  clusteringOutput.Print();
  cout<<"execution time:"<<executionTime<<endl;
  cout<<"score:"<<score<<endl;
  cout<<"normalized score:"<<normalizedScore<<endl;
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

void Chromosome::CalculateScore()
{
  clusteringOutput = ReadOutput(clusteringOutputPath);
  
  system(("rm -f "+configPath).c_str());
  system(("rm -f "+clusteringOutputPath).c_str());
  
  double distance =     fabs(clusteringOutput.containmentMean-1)
                      +      clusteringOutput.containmentSigma
                      + fabs(clusteringOutput.resolutionMean)
                      +      clusteringOutput.resolutionSigma
                      +      clusteringOutput.separationMean
                      +      clusteringOutput.separationSigma;
  
  score = severityFactor/distance;
  
  if(clusteringOutput.resolutionMean > 1000){ // this means that clustering failed completely
    score = 0;
  }
  if(score < 1E-5){ // just round down to zero if score it extremaly poor
    score = 0;
  }
  
  if(score==0){
    cout<<"This chromosome failed completely:"<<endl;
    Print();
  }
}

uint64_t Chromosome::SinglePointCrossover(uint64_t a, uint64_t b)
{
  int crossingPoint = RandInt(0, 63);
  uint64_t newBitChromosome = 0;
  
  uint64_t maskA = 0;
  for(int j=0;j<BitSize(maskA)-crossingPoint;j++){maskA |= 1ull << (j+crossingPoint);}
  
  uint64_t maskB = 0;
  for(int j=(int)BitSize(maskB)-crossingPoint;j<BitSize(maskB);j++){maskB |= 1ull << j;}
  
  newBitChromosome = a & maskA;
  newBitChromosome |= b & maskB;
  
  return newBitChromosome;
}

Chromosome* Chromosome::ProduceChildWith(Chromosome *partner)
{
  Chromosome *child = new Chromosome();
  
  // combine chromosomes of parents in a random way
  
  // this is a single-point crossover
  for(int i=0;i<bitChromosome.size();i++){
    uint64_t newBitChromosome = SinglePointCrossover(GetBitChromosome(i), partner->GetBitChromosome(i));
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
  // populate fields
  child->ReadFromBitChromosome();
  
  // make sure that after crossing and mutation all parameters are within limits
  int wasOutside = child->BackToLimits();
  if(wasOutside > 0){
//    cout<<"produced child was outside of limits ("<<wasOutside<<")"<<endl;
  }
  
  // update the bits after updating values!!
  child->SaveToBitChromosome();
  
  // set other parameters
  child->SetMutationChance(mutationChance);
  child->SetSeverityFactor(severityFactor);
  
  return child;
}

int Chromosome::BackToLimits()
{
  int wasOutside = 0;
  
  if(GetCriticalDistanceEE() < criticalDistanceEEmin || GetCriticalDistanceEE() > criticalDistanceEEmax){
    SetCriticalDistanceEE(RandFloat(criticalDistanceEEmin, criticalDistanceEEmax));
    wasOutside++;
    cout<<"a"<<endl;
  }
  
  if(GetCriticalDistanceFH() < criticalDistanceFHmin || GetCriticalDistanceFH() > criticalDistanceFHmax){
    SetCriticalDistanceFH(RandFloat(criticalDistanceFHmin, criticalDistanceFHmax));
    wasOutside++;
    cout<<"b"<<endl;
  }
  if(GetCriticalDistanceBH() < criticalDistanceBHmin || GetCriticalDistanceBH() > criticalDistanceBHmax){
    SetCriticalDistanceBH(RandFloat(criticalDistanceBHmin, criticalDistanceBHmax));
    wasOutside++;
    cout<<"c"<<endl;
  }
  
  if(GetKernel() < kernelMin || GetKernel() > kernelMax){
    SetKernel(RandInt(kernelMin, kernelMax));
    wasOutside++;
    cout<<"d"<<endl;
  }
    
  if(GetDeltacEE() < deltacEEmin || GetDeltacEE() > deltacEEmax){
    SetDeltacEE(RandFloat(deltacEEmin, deltacEEmax));
    wasOutside++;
    cout<<"e"<<endl;
  }
  if(GetDeltacFH() < deltacFHmin || GetDeltacFH() > deltacFHmax){
    SetDeltacFH(RandFloat(deltacFHmin, deltacFHmax));
    wasOutside++;
    cout<<"f"<<endl;
  }
  if(GetDeltacBH() < deltacBHmin || GetDeltacBH() > deltacBHmax){
    SetDeltacBH(RandFloat(deltacBHmin, deltacBHmax));
    wasOutside++;
    cout<<"g"<<endl;
  }
    
  if(GetKappa() < kappaMin || GetKappa() > kappaMax){
    SetKappa(RandFloat(kappaMin, kappaMax));
    wasOutside++;
    cout<<"h"<<endl;
  }
  
  if(GetDependSensor()){
    if(GetEnergyMin() < energyThresholdMin || GetEnergyMin() > energyThresholdMax){
      SetEnergyMin(RandFloat(energyThresholdMin, energyThresholdMax));
      wasOutside++;
      cout<<"i"<<endl;
    }
  }
  else{
    if(GetEnergyMin() < energyThresholdMinNoSensor || GetEnergyMin() > energyThresholdMaxNoSensor){
      SetEnergyMin(RandFloat(energyThresholdMinNoSensor, energyThresholdMaxNoSensor));
      wasOutside++;
      cout<<"j";
    }
  }
    
  if(GetMatchingDistance() < matchingDistanceMin || GetMatchingDistance() > matchingDistanceMax){
    SetMatchingDistance(RandFloat(matchingDistanceMin, matchingDistanceMax));
    wasOutside++;
    cout<<"k"<<endl;
  }
  if(GetMinClusters() < minClustersMin || GetMinClusters() > minClustersMax){
    SetMinClusters(RandInt(minClustersMin, minClustersMax));
    wasOutside++;
    cout<<"l"<<endl;
  }
  return wasOutside;
}
