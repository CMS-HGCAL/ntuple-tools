//
//  ChromosomeProcessor.cpp
//
//  Created by Jeremi Niedziela on 13/02/2019.
//

#include "ChromosomeProcessor.hpp"

using namespace std;

ChromosomeProcessor::ChromosomeProcessor(double _mutationChance, double _severityFactor, ECrossover _crossover) :
mutationChance(_mutationChance),
severityFactor(_severityFactor),
crossover(_crossover)
{
  
}

ChromosomeProcessor::~ChromosomeProcessor()
{
  
}

shared_ptr<Chromosome> ChromosomeProcessor::GetRandomChromosome()
{
  auto result = make_shared<Chromosome>();
  
  for(int i=0;i<kNparams;i++){
    result->SetParam((EParam)i,RandDouble(paramMin[i],paramMax[i]));
  }
  return result;
}

void ChromosomeProcessor::StoreChromosomeInConfig(const shared_ptr<Chromosome> chromo, string path)
{
  string currentConfigPath = path=="" ? chromo->configPath : path;
  
  system(("cp baseConfig.md "+currentConfigPath).c_str());
  
  int kernelIndex = round(chromo->GetParam(kKernel));
  
  if(kernelIndex == 0)      UpdateParamValue(currentConfigPath, "energy_density_function","step");
  else if(kernelIndex == 1) UpdateParamValue(currentConfigPath, "energy_density_function","gaus");
  else                      UpdateParamValue(currentConfigPath, "energy_density_function","exp");
  
  UpdateParamValue(currentConfigPath, "critical_distance_EE",chromo->GetParam(kCriticalDistanceEE));
  UpdateParamValue(currentConfigPath, "critical_distance_FH",chromo->GetParam(kCriticalDistanceFH));
  UpdateParamValue(currentConfigPath, "critical_distance_BH",chromo->GetParam(kCriticalDistanceBH));
  UpdateParamValue(currentConfigPath, "deltac_EE",chromo->GetParam(kDeltacEE));
  UpdateParamValue(currentConfigPath, "deltac_FH",chromo->GetParam(kDeltacFH));
  UpdateParamValue(currentConfigPath, "deltac_BH",chromo->GetParam(kDeltacBH));
  UpdateParamValue(currentConfigPath, "kappa",chromo->GetParam(kKappa));
  UpdateParamValue(currentConfigPath, "energy_min",chromo->GetParam(kEnergyThreshold));
  UpdateParamValue(currentConfigPath, "matching_max_distance",chromo->GetParam(kMatchingDistance));
  UpdateParamValue(currentConfigPath, "score_output_path",chromo->clusteringOutputPath);
}

void ChromosomeProcessor::CalculateScore(shared_ptr<Chromosome> chromo)
{
  chromo->clusteringOutput = ReadOutput(chromo->GetClusteringOutputPath());
  
  system(("rm -f "+chromo->configPath).c_str());
  system(("rm -f "+chromo->clusteringOutputPath).c_str());
  
  ClusteringOutput output = chromo->clusteringOutput;
  
  double distance =     fabs(output.containmentMean-1)
                      +      output.containmentSigma
                      + fabs(output.resolutionMean)
                      +      output.resolutionSigma
                      +      output.separationMean
                      +      output.separationSigma
                      + fabs(output.deltaNclustersMean)
                      +      output.deltaNclustersSigma
                      +      output.nRecoFailed
                      +      output.nCantMatchRecSim
                      +      output.nFakeRec;
  
  
  
  chromo->score = severityFactor/distance;
  
  if(   output.resolutionMean     > 1000
     || output.separationMean     > 1000
     || output.containmentMean    > 1000
     || output.deltaNclustersMean > 1000
     || output.nRecoFailed        > 1000
     || output.nCantMatchRecSim   > 1000
     || output.nFakeRec           > 1000
     )
  { // this means that clustering failed completely
    chromo->score = 0;
  }
  if(chromo->score < 1E-5){ // just round down to zero if score it extremaly poor
    chromo->score = 0;
  }
  
  if(chromo->score==0){
    cout<<"This chromosome failed completely:"<<endl;
    chromo->Print();
  }
}

pair<shared_ptr<Chromosome>, shared_ptr<Chromosome>>
ChromosomeProcessor::CrossChromosomes(const shared_ptr<Chromosome> mom, const shared_ptr<Chromosome> dad)
{
  pair<shared_ptr<Chromosome>, shared_ptr<Chromosome>> children(make_shared<Chromosome>(),
                                                                make_shared<Chromosome>());
  
  // Perform crossover accoring to the strategy specified by the user
  int crossingChromo = RandInt(0, (int)mom->bitChromosome.size()-1);
  pair<uint64_t,uint64_t>  newBitChromosomes;
  
  for(int i=0;i<mom->bitChromosome.size();i++){
    if(crossover == kMultiPoint){
      newBitChromosomes = SinglePointCrossover(mom->bitChromosome[i], dad->bitChromosome[i]);
    }
    else if(crossover == kSinglePoint || crossover == kFixedSinglePoint){
      if(i == crossingChromo){
        newBitChromosomes = SinglePointCrossover(mom->bitChromosome[i], dad->bitChromosome[i],
                                                 (crossover==kFixedSinglePoint));
      }
      else if(i < crossingChromo)   newBitChromosomes = make_pair(mom->bitChromosome[i], dad->bitChromosome[i]);
      else                          newBitChromosomes = make_pair(dad->bitChromosome[i], mom->bitChromosome[i]);
    }
    else if(crossover == kUniform){
      newBitChromosomes = make_pair(mom->bitChromosome[i], dad->bitChromosome[i]);
      
      for(int j=0;j<64;j++){
        if(RandBool()){
          ReverseBit(newBitChromosomes.first, j);
          ReverseBit(newBitChromosomes.second, j);
        }
      }
    }
    children.first->bitChromosome[i]  = newBitChromosomes.first;
    children.second->bitChromosome[i] = newBitChromosomes.second;
  }
  
  // Perform child genes mutation
  for(int i=0;i<mom->bitChromosome.size();i++){
    pair<uint64_t, uint64_t> bits(children.first->bitChromosome[i], children.second->bitChromosome[i]);
    
    for(int iBit=0;iBit<BitSize(bits.first);iBit++){
      if(RandDouble(0, 1) < mutationChance)  ReverseBit(bits.first, iBit);
      if(RandDouble(0, 1) < mutationChance)  ReverseBit(bits.second, iBit);
    }
    children.first->bitChromosome[i]  = bits.first;
    children.second->bitChromosome[i] = bits.second;
  }
  
  // Populate fields
  children.first->ReadFromBitChromosome();
  children.second->ReadFromBitChromosome();
  
  // Make sure that new parameters are within limits
  if(children.first->BackToLimits() || children.second->BackToLimits()){
    cout<<"Param was outside (should never happen!)"<<endl;
  }
  
  // Fix parameters that should be fixed
  for(int iPar=0;iPar<kNparams;iPar++){
    if(mom->isParamFixed[iPar]){
      children.first->FixParam((EParam)iPar,mom->GetParam((EParam)iPar));
      children.second->FixParam((EParam)iPar,mom->GetParam((EParam)iPar));
    }
  }
  
  // Update the bits after updating values!!
  children.first->SaveToBitChromosome();
  children.second->SaveToBitChromosome();
  
  return children;
}

pair<uint64_t,uint64_t> ChromosomeProcessor::SinglePointCrossover(uint64_t a, uint64_t b, bool fixed)
{
  int crossingPoint;
  if(fixed){
    int crossingParam = RandInt(0, 3);  // this is the index of parameter to cross after
    crossingPoint = crossingParam * 16; // crossing point will be 0, 16, 32 or 48, preserving parameter content
  }
  else{
    crossingPoint = RandInt(0, 63);
  }
  
  pair<uint64_t,uint64_t> newBitChromosomes(0,0);
  
  uint64_t maskA = 0;
  for(int j=crossingPoint;j<BitSize(maskA);j++){maskA |= (1ull << j);}
  
  uint64_t maskB = ~maskA;
  
  newBitChromosomes.first  = a & maskA;
  newBitChromosomes.first |= b & maskB;
  
  newBitChromosomes.second  = b & maskA;
  newBitChromosomes.second |= a & maskB;
  
  return newBitChromosomes;
}
