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
dependSensor(true),
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
 
  crossover = kSinglePoint;
  
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
  result->SetReachedEE(RandBool());
  result->SetKernel(RandInt(kernelMin, kernelMax));
  result->SetDeltacEE(RandFloat(deltacEEmin, deltacEEmax));
  result->SetDeltacFH(RandFloat(deltacFHmin, deltacFHmax));
  result->SetDeltacBH(RandFloat(deltacBHmin, deltacBHmax));
  result->SetKappa(RandFloat(kappaMin, kappaMax));
  result->SetEnergyMin(RandFloat(energyThresholdMin, energyThresholdMax));
  result->SetMatchingDistance(RandFloat(matchingDistanceMin, matchingDistanceMax));
  result->SetMinClusters(RandInt(minClustersMin, minClustersMax));
  
  result->SetDependSensor(true); // it has to depend on the sensor thickness
  
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


void Chromosome::StoreInConfig(string path)
{
  string currentConfigPath = path=="" ? configPath : path;
  
  system(("cp baseConfig.md "+currentConfigPath).c_str());
  
  UpdateParamValue(currentConfigPath, "depend_sensor",GetDependSensor());
  if(GetKernel() == 0)      UpdateParamValue(currentConfigPath, "energy_density_function","step");
  else if(GetKernel() == 1) UpdateParamValue(currentConfigPath, "energy_density_function","gaus");
  else                      UpdateParamValue(currentConfigPath, "energy_density_function","exp");
  
  UpdateParamValue(currentConfigPath, "critial_distance_EE",GetCriticalDistanceEE());
  UpdateParamValue(currentConfigPath, "critial_distance_FH",GetCriticalDistanceFH());
  UpdateParamValue(currentConfigPath, "critial_distance_BH",GetCriticalDistanceBH());
  UpdateParamValue(currentConfigPath, "deltac_EE",GetDeltacEE());
  UpdateParamValue(currentConfigPath, "deltac_FH",GetDeltacFH());
  UpdateParamValue(currentConfigPath, "deltac_BH",GetDeltacBH());
  UpdateParamValue(currentConfigPath, "kappa",GetKappa());
  UpdateParamValue(currentConfigPath, "energy_min",GetEnergyMin());
  UpdateParamValue(currentConfigPath, "min_clusters",GetMinClusters());
  UpdateParamValue(currentConfigPath, "reachedEE_only",GetReachedEE());
  UpdateParamValue(currentConfigPath, "matching_max_distance",GetMatchingDistance());
  UpdateParamValue(currentConfigPath, "score_output_path",clusteringOutputPath);
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
//                      + fabs(clusteringOutput.deltaNclustersMean)
//                      +      clusteringOutput.deltaNclustersSigma;
  
  score = severityFactor/distance;
  
  if(   clusteringOutput.resolutionMean     > 1000
     || clusteringOutput.separationMean     > 1000
     || clusteringOutput.containmentMean    > 1000
     || clusteringOutput.deltaNclustersMean > 1000
    )
  { // this means that clustering failed completely
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

vector<uint64_t> Chromosome::SinglePointCrossover(uint64_t a, uint64_t b)
{
  int crossingPoint = RandInt(0, 63);
  vector<uint64_t> newBitChromosomes = {0,0};
  
  uint64_t maskA = 0;
  for(int j=0;j<BitSize(maskA)-crossingPoint;j++){maskA |= 1ull << (j+crossingPoint);}
  
  uint64_t maskB = 0;
  for(int j=(int)BitSize(maskB)-crossingPoint;j<BitSize(maskB);j++){maskB |= 1ull << j;}
  
  newBitChromosomes[0] = a & maskA;
  newBitChromosomes[0] |= b & maskB;
  
  newBitChromosomes[1] = b & maskA;
  newBitChromosomes[1] |= a & maskB;
  
  return newBitChromosomes;
}

vector<Chromosome*> Chromosome::ProduceChildWith(Chromosome *partner)
{
  vector<Chromosome*> children = {new Chromosome(),new Chromosome()};
  // combine chromosomes of parents in a random way
  
  if(crossover == kMultiPoint){ // single-point crossover in each chromosome
    for(int i=0;i<bitChromosome.size();i++){
      vector<uint64_t> newBitChromosomes = SinglePointCrossover(GetBitChromosome(i), partner->GetBitChromosome(i));
      
      children[0]->SetBitChromosome(i, newBitChromosomes[0]);
      children[1]->SetBitChromosome(i, newBitChromosomes[1]);
    }
  }
  else if(crossover == kSinglePoint){ // true single-point crossover
    int crossingChromo = RandInt(0, (int)bitChromosome.size()-1);
    
    for(int i=0;i<bitChromosome.size();i++){
    
      if(i < crossingChromo){
        children[0]->SetBitChromosome(i, this->GetBitChromosome(i));
        children[1]->SetBitChromosome(i, partner->GetBitChromosome(i));
      }
      else if(i == crossingChromo){
        vector<uint64_t> newBitChromosomes = SinglePointCrossover(this->GetBitChromosome(i),
                                                                  partner->GetBitChromosome(i));
        children[0]->SetBitChromosome(i, newBitChromosomes[0]);
        children[1]->SetBitChromosome(i, newBitChromosomes[1]);
      }
      else{
        children[0]->SetBitChromosome(i, partner->GetBitChromosome(i));
        children[1]->SetBitChromosome(i, this->GetBitChromosome(i));
      }
    }
  }
  else if(crossover == kUniform){
    
  }
  
  // perform child genes mutation
  for(int iChild=0;iChild<2;iChild++){
    for(int i=0;i<bitChromosome.size();i++){
      uint64_t bits = children[iChild]->GetBitChromosome(i);
      
      for(int iBit=0;iBit<BitSize(bits);iBit++){
        double r = RandDouble(0, 1);
        if(r < mutationChance){ // 10% chance for mutation
          ReverseBit(bits, iBit);
        }
      }
      children[iChild]->SetBitChromosome(i, bits);
    }
    // populate fields
    children[iChild]->ReadFromBitChromosome();
    
    // make sure that after crossing and mutation all parameters are within limits
    int wasOutside = children[iChild]->BackToLimits();
    if(wasOutside > 0){
      cout<<"ERROR -- produced child was outside of limits ("<<wasOutside<<" times)"<<endl;
      cout<<"with the current implementation it should never happen!!"<<endl;
    }
    
    // update the bits after updating values!!
    children[iChild]->SaveToBitChromosome();
    
    // set other parameters
    children[iChild]->SetMutationChance(mutationChance);
    children[iChild]->SetSeverityFactor(severityFactor);
    children[iChild]->SetCrossover(crossover);
  }
  
  return children;
}

int Chromosome::BackToLimits()
{
  int wasOutside = 0;
  
  if(GetCriticalDistanceEE() < criticalDistanceEEmin || GetCriticalDistanceEE() > criticalDistanceEEmax){
    SetCriticalDistanceEE(RandFloat(criticalDistanceEEmin, criticalDistanceEEmax));
    wasOutside++;
  }
  
  if(GetCriticalDistanceFH() < criticalDistanceFHmin || GetCriticalDistanceFH() > criticalDistanceFHmax){
    SetCriticalDistanceFH(RandFloat(criticalDistanceFHmin, criticalDistanceFHmax));
    wasOutside++;
  }
  if(GetCriticalDistanceBH() < criticalDistanceBHmin || GetCriticalDistanceBH() > criticalDistanceBHmax){
    SetCriticalDistanceBH(RandFloat(criticalDistanceBHmin, criticalDistanceBHmax));
    wasOutside++;
  }
  
  if(GetKernel() < kernelMin || GetKernel() > kernelMax){
    SetKernel(RandInt(kernelMin, kernelMax));
    wasOutside++;
  }
    
  if(GetDeltacEE() < deltacEEmin || GetDeltacEE() > deltacEEmax){
    SetDeltacEE(RandFloat(deltacEEmin, deltacEEmax));
    wasOutside++;
  }
  if(GetDeltacFH() < deltacFHmin || GetDeltacFH() > deltacFHmax){
    SetDeltacFH(RandFloat(deltacFHmin, deltacFHmax));
    wasOutside++;
  }
  if(GetDeltacBH() < deltacBHmin || GetDeltacBH() > deltacBHmax){
    SetDeltacBH(RandFloat(deltacBHmin, deltacBHmax));
    wasOutside++;
  }
    
  if(GetKappa() < kappaMin || GetKappa() > kappaMax){
    SetKappa(RandFloat(kappaMin, kappaMax));
    wasOutside++;
  }
  
  if(GetEnergyMin() < energyThresholdMin || GetEnergyMin() > energyThresholdMax){
    SetEnergyMin(RandFloat(energyThresholdMin, energyThresholdMax));
    wasOutside++;
  }
    
  if(GetMatchingDistance() < matchingDistanceMin || GetMatchingDistance() > matchingDistanceMax){
    SetMatchingDistance(RandFloat(matchingDistanceMin, matchingDistanceMax));
    wasOutside++;
  }
  if(GetMinClusters() < minClustersMin || GetMinClusters() > minClustersMax){
    SetMinClusters(RandInt(minClustersMin, minClustersMax));
    wasOutside++;
  }
  return wasOutside;
}
