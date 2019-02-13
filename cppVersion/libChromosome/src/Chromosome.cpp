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
normalizedScore(-99999)
{
  int nChromosomes = 3;
  for(int i=0;i<nChromosomes;i++){bitChromosome.push_back(0);}
  
  
  clusteringOutput = ClusteringOutput();

  for(int i=0;i<kNparams;i++){
    params.push_back(0.0);
    isParamFixed.push_back(false);
  }
  
  uniqueID = reinterpret_cast<uint64_t>(this);
  configPath = "tmp/config_"+to_string(uniqueID)+".md";
  clusteringOutputPath = "tmp/output_"+to_string(uniqueID)+".txt";
}

Chromosome::~Chromosome()
{
  
}

void Chromosome::SetParam(EParam par, double val)
{
  params[par] = static_cast<uint16_t>(std::numeric_limits<uint16_t>::max()/(paramMax[par]-paramMin[par])*(val-paramMin[par]));
}

void Chromosome::FixParam(EParam par, double val)
{
  SetParam(par, val);
  isParamFixed[par] = true;
}

double Chromosome::GetParam(EParam par)
{
  return paramMin[par] + (double)params[par]*(paramMax[par]-paramMin[par])/std::numeric_limits<uint16_t>::max();
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
