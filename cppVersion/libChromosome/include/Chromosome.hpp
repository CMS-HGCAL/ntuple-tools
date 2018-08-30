//
//  Chromosome.hpp
//
//  Created by Jeremi Niedziela on 02/08/2018.
//

#ifndef Chromosome_hpp
#define Chromosome_hpp

#include "GeneticHelpers.hpp"

class Chromosome
{
public:
  enum ECrossover{
    kUniform,       ///< each bit has a chance to be exchaned between parents
    kSinglePoint,   ///< one chromosome gets crossed, the rest stays the same or is exchaned intact
    kFixedSinglePoint,///< one chromosome gets crossed at point that doesn't modify any of the parameters, the rest stays the same or is exchaned intact
    kMultiPoint,     ///< each chromosome is crossed at a random point
    kNcrossover
  };
  
  inline static std::string crossoverName[kNcrossover] = {
    "Uniform",
    "Single point",
    "Fixed single point",
    "Multi point"
  };
  
  /// Default constructor
  Chromosome();
  
  /// Default destructor
  ~Chromosome();
  
  static Chromosome* GetRandom();
  
  inline void SetParam(EParam par, double val){
    params[par] = static_cast<uint16_t>(std::numeric_limits<uint16_t>::max()/(paramMax[par]-paramMin[par])*(val-paramMin[par]));
  }
  
  inline void SetBitChromosome(int i, uint64_t bits){bitChromosome[i] = bits;}
  
  inline void SetScore(double val){score = val;}
  inline void SetNormalizedScore(double val){normalizedScore = val;}
  inline void SetExecutionTime(double val){executionTime = val;}
  inline void SetMutationChance(double val){mutationChance = val;}
  inline void SetSeverityFactor(double val){severityFactor = val;}
  inline void SetCrossover(ECrossover val){crossover = val;}
  inline void SetInputDataPath(std::string val){inputDataPath = val;}
  
  
  // Getters
  inline double GetParam(EParam par){
    return paramMin[par] + (double)params[par]*(paramMax[par]-paramMin[par])/std::numeric_limits<uint16_t>::max();
  }
  inline uint64_t GetBitChromosome(int i){return bitChromosome[i];}
  
  inline double  GetScore(){return score;}
  inline double  GetNormalizedScore(){return normalizedScore;}
  inline uint64_t GetUniqueID(){return uniqueID;}
  inline std::string GetConfigPath(){return configPath;}
  inline std::string GetClusteringOutputPath(){return clusteringOutputPath;}
  
  void SaveToBitChromosome();
  void ReadFromBitChromosome();
  
  void StoreInConfig(std::string path="");
  
  void Print();
  void CalculateScore();
  
  std::vector<Chromosome*> ProduceChildWith(Chromosome *partner);
private:
  std::vector<uint16_t> params;
  
  std::vector<uint64_t> bitChromosome;
  
  // other variables not stored in bit chromosome
  uint64_t uniqueID;
  std::string configPath;
  std::string clusteringOutputPath;
  std::string inputDataPath;
  double mutationChance;
  double severityFactor;  // larger the value, more easily population members will die
  ECrossover crossover;
  
  ClusteringOutput clusteringOutput;
  double executionTime;
  double score;
  double normalizedScore;
  
  template<class T>
  void ShiftIntoChromosome(T value, int &shift, int chromoIndex);
  
  template<class T>
  void SetValueFromChromosome(T &value, int &shift, int chromoIndex);
  
  std::vector<uint64_t> SinglePointCrossover(uint64_t a, uint64_t b, bool fixed = false);
  
  int BackToLimits();
};

#endif /* Chromosome_hpp */
