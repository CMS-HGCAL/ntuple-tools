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
executionTime(99999),
score(-99999),
normalizedScore(-99999),
mutationChance(0.002),
severityFactor(1.0),
crossover(kSinglePoint),
inputDataPath("")
{
  int nChromosomes = 3;
  for(int i=0;i<nChromosomes;i++){bitChromosome.push_back(0);}
  
  clusteringOutput = ClusteringOutput();

  for(int i=0;i<kNparams;i++){params.push_back(0.0);}
  
  uniqueID = reinterpret_cast<uint64_t>(this);
  configPath = "tmp/config_"+to_string(uniqueID)+".md";
  clusteringOutputPath = "tmp/output_"+to_string(uniqueID)+".txt";
}

Chromosome::~Chromosome()
{
  
}

Chromosome* Chromosome::GetRandom()
{
  Chromosome *result = new Chromosome();
  
  for(int i=0;i<kNparams;i++){
    result->SetParam((EParam)i,RandDouble(paramMin[i],paramMax[i]));
  }
  
  return result;
}

void Chromosome::SaveToBitChromosome()
{
  for(int i=0;i<bitChromosome.size();i++){bitChromosome[i] = 0;}
  
  int currentShift = 0;
  int chromoIndex = -1;
  
  for(int iPar=0;iPar<kNparams;iPar++){
    if(iPar%4 == 0){
      currentShift = 0;
      chromoIndex++;
    }
    ShiftIntoChromosome(params[iPar],currentShift,chromoIndex);
  }
}

void Chromosome::ReadFromBitChromosome()
{
  int currentShift = 0;
  int chromoIndex = -1;
  
  for(int iPar=0;iPar<kNparams;iPar++){
    if(iPar%4 == 0){
      currentShift = 0;
      chromoIndex++;
    }
    SetValueFromChromosome(params[iPar], currentShift, chromoIndex);
  }
}

void Chromosome::Print()
{
  cout<<"================================================="<<endl;
  cout<<"Chromosome "<<GetUniqueID()<<endl;
  
  cout<<"Critical distance:";
  cout<<"\tEE:"<<GetParam(kCriticalDistanceEE);
  cout<<"\tFH:"<<GetParam(kCriticalDistanceFH);
  cout<<"\tBH:"<<GetParam(kCriticalDistanceBH)<<endl;
  
  cout<<"Kernel: ";
  if(round(GetParam(kKernel)) == 0)        cout<<"step"<<endl;
  else if(round(GetParam(kKernel)) == 1)   cout<<"gaus"<<endl;
  else                                     cout<<"exp"<<endl;
  
  cout<<"Critical #delta:";
  cout<<"\tEE:"<<GetParam(kDeltacEE);
  cout<<"\tFH:"<<GetParam(kDeltacFH);
  cout<<"\tBH:"<<GetParam(kDeltacBH)<<endl;
  
  cout<<"kappa:"<<GetParam(kKappa)<<endl;
  cout<<"energy threshold:"<<GetParam(kEnergyThreshold)<<endl;
  cout<<"max matching distance:"<<GetParam(kMatchingDistance)<<endl;
  
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
  
  UpdateParamValue(currentConfigPath, "input_path",inputDataPath);
  UpdateParamValue(currentConfigPath, "min_layer",minLayer);
  UpdateParamValue(currentConfigPath, "max_layer",maxLayer);
  
  int kernelIndex = round(GetParam(kKernel));
  
  if(kernelIndex == 0)      UpdateParamValue(currentConfigPath, "energy_density_function","step");
  else if(kernelIndex == 1) UpdateParamValue(currentConfigPath, "energy_density_function","gaus");
  else                      UpdateParamValue(currentConfigPath, "energy_density_function","exp");
  
  UpdateParamValue(currentConfigPath, "critical_distance_EE",GetParam(kCriticalDistanceEE));
  UpdateParamValue(currentConfigPath, "critical_distance_FH",GetParam(kCriticalDistanceFH));
  UpdateParamValue(currentConfigPath, "critical_distance_BH",GetParam(kCriticalDistanceBH));
  UpdateParamValue(currentConfigPath, "deltac_EE",GetParam(kDeltacEE));
  UpdateParamValue(currentConfigPath, "deltac_FH",GetParam(kDeltacFH));
  UpdateParamValue(currentConfigPath, "deltac_BH",GetParam(kDeltacBH));
  UpdateParamValue(currentConfigPath, "kappa",GetParam(kKappa));
  UpdateParamValue(currentConfigPath, "energy_min",GetParam(kEnergyThreshold));
  UpdateParamValue(currentConfigPath, "matching_max_distance",GetParam(kMatchingDistance));
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
                      +      clusteringOutput.separationSigma
                      + fabs(clusteringOutput.deltaNclustersMean)
                      +      clusteringOutput.deltaNclustersSigma
                      +      clusteringOutput.nRecoFailed
                      +      clusteringOutput.nCantMatchRecSim
                      +      clusteringOutput.nFakeRec;
  
  
  
  score = severityFactor/distance;
  
  if(   clusteringOutput.resolutionMean     > 1000
     || clusteringOutput.separationMean     > 1000
     || clusteringOutput.containmentMean    > 1000
     || clusteringOutput.deltaNclustersMean > 1000
     || clusteringOutput.nRecoFailed        > 1000
     || clusteringOutput.nCantMatchRecSim   > 1000
     || clusteringOutput.nFakeRec           > 1000
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

vector<uint64_t> Chromosome::SinglePointCrossover(uint64_t a, uint64_t b, bool fixed)
{
  int crossingPoint;
  if(fixed){
    int crossingParam = RandInt(0, 3); // this is the index of parameter to cross after
    crossingPoint = crossingParam * 16; // crossing point will be 0, 16, 32 or 48, preserving parameter content
  }
  else{
    crossingPoint = RandInt(0, 63);
  }
  
  vector<uint64_t> newBitChromosomes = {0,0};
  
  uint64_t maskA = 0;
  for(int j=crossingPoint;j<BitSize(maskA);j++){maskA |= (1ull << j);}
  
  uint64_t maskB = ~maskA;
  
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
  else if(crossover == kSinglePoint || crossover == kFixedSinglePoint){ // true single-point crossover (can be fixed)
    int crossingChromo = RandInt(0, (int)bitChromosome.size()-1);
    
    for(int i=0;i<bitChromosome.size();i++){
    
      if(i < crossingChromo){
        children[0]->SetBitChromosome(i, this->GetBitChromosome(i));
        children[1]->SetBitChromosome(i, partner->GetBitChromosome(i));
      }
      else if(i == crossingChromo){
        vector<uint64_t> newBitChromosomes = SinglePointCrossover(this->GetBitChromosome(i),
                                                                  partner->GetBitChromosome(i),
                                                                  (crossover==kFixedSinglePoint));
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
    for(int i=0;i<bitChromosome.size();i++){
      uint64_t dadBits = this->GetBitChromosome(i);
      uint64_t momBits = partner->GetBitChromosome(i);
      
      for(int j=0;j<64;j++){
        bool cross = RandBool();
        if(cross){
          ReverseBit(dadBits, j);
          ReverseBit(momBits, j);
        }
      }
      children[0]->SetBitChromosome(i, dadBits);
      children[1]->SetBitChromosome(i, momBits);
    }
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
    
    int wasOutside = children[iChild]->BackToLimits();
    if(wasOutside) cout<<"Was outside:"<<wasOutside<<endl;
    // update the bits after updating values!!
    children[iChild]->SaveToBitChromosome();
    
    // set other parameters
    children[iChild]->SetMutationChance(mutationChance);
    children[iChild]->SetSeverityFactor(severityFactor);
    children[iChild]->SetCrossover(crossover);
    children[iChild]->SetInputDataPath(inputDataPath);
    children[iChild]->SetMinLayer(minLayer);
    children[iChild]->SetMaxLayer(maxLayer);
  }
  
  return children;
}


int Chromosome::BackToLimits()
{
  int wasOutside = 0;
  
  for(int i=0;i<kNparams;i++){
    if(GetParam((EParam)i) <= paramMin[i] || GetParam((EParam)i) >= paramMax[i]){
      SetParam((EParam)i, RandDouble(paramMin[i], paramMax[i]));
      wasOutside++;
    }
  }
  return wasOutside;
}
