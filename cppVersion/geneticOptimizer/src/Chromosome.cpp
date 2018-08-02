//
//  Chromosome.cpp
//
//  Created by Jeremi Niedziela on 02/08/2018.
//

#include "Chromosome.hpp"

#include <iostream>

using namespace std;

Chromosome::Chromosome()
{
  uniqueID = reinterpret_cast<uint64_t>(this);
}

Chromosome::Chromosome(Chromosome &c)
{
  for(int i=0;i<3;i++){
    bitChromosome[i] = c.bitChromosome[i];
  }
  uniqueID = reinterpret_cast<uint64_t>(this);
}

Chromosome::~Chromosome()
{
  
}

Chromosome Chromosome::GetRandom()
{
  Chromosome result;
  
  result.SetCriticalDistanceEE(RandFloat(0.0, 100.0));
  result.SetCriticalDistanceFH(RandFloat(0.0, 100.0));
  result.SetCriticalDistanceBH(RandFloat(0.0, 100.0));
  result.SetDependSensor(RandBool());
  result.SetReachedEE(RandBool());
  result.SetKernel(RandInt(0, 2));
  result.SetDeltacEE(RandFloat(0.0, 100.0));
  result.SetDeltacFH(RandFloat(0.0, 100.0));
  result.SetDeltacBH(RandFloat(0.0, 100.0));
  result.SetKappa(RandFloat(0.0, 10000.0));
  result.SetEnergyMin(RandFloat(0.0, 100.0));
  result.SetMatchingDistance(RandFloat(0.0, 100.0));
  result.SetMinClusters(RandInt(0, 100));
  
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
  cout<<"Chromosome "<<uniqueID<<endl;
  
  cout<<"Critical distance:"<<endl;
  cout<<"\tEE:"<<criticalDistanceEE/1000.<<endl;
  cout<<"\tFH:"<<criticalDistanceFH/1000.<<endl;
  cout<<"\tBH:"<<criticalDistanceBH/1000.<<endl;
  
  cout<<"Depend on sensor:"<<dependSensor<<endl;
  cout<<"Reached EE only:"<<reachedEE<<endl;
  
  cout<<"Kernel: ";
  if(kernel==0) cout<<"step"<<endl;
  if(kernel==1) cout<<"gaus"<<endl;
  if(kernel==2) cout<<"exp"<<endl;
  
  cout<<"Critical #delta:"<<endl;
  cout<<"\tEE:"<<deltacEE/1000.<<endl;
  cout<<"\tFH:"<<deltacFH/1000.<<endl;
  cout<<"\tBH:"<<deltacBH/1000.<<endl;
  
  cout<<"kappa:"<<kappa/1000.<<endl;
  cout<<"energy threshold:"<<energyMin/100000.<<endl;
  cout<<"max matching distance:"<<matchingDistance/1000.<<endl;
  cout<<"min clusters:"<<minClusters<<endl;
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
  string configPath = "tmp/config_"+to_string(uniqueID)+".md";
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
  
}





