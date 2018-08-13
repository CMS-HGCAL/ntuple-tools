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
  /// Default constructor
  Chromosome();
  
  /// Default destructor
  ~Chromosome();
  
  static Chromosome* GetRandom();
  
  // Setters
  inline void SetCriticalDistanceEE(float val){criticalDistanceEE = (uint32_t)(val*1000);}
  inline void SetCriticalDistanceFH(float val){criticalDistanceFH = (uint32_t)(val*1000);}
  inline void SetCriticalDistanceBH(float val){criticalDistanceBH = (uint32_t)(val*1000);}
  inline void SetDependSensor(bool val){dependSensor = val;}
  inline void SetReachedEE(bool val){reachedEE = val;}
  inline void SetKernel(int val){kernel = (uint8_t)val;} // 0 - step, 1 - gaus, 2 - exp
  inline void SetDeltacEE(float val){deltacEE = (uint32_t)(val*1000);}
  inline void SetDeltacFH(float val){deltacFH = (uint32_t)(val*1000);}
  inline void SetDeltacBH(float val){deltacBH = (uint32_t)(val*1000);}
  inline void SetKappa(float val){kappa = (uint32_t)(val*1000);}
  inline void SetEnergyMin(float val){energyMin = (uint32_t)(val*100000);}
  inline void SetMatchingDistance(float val){matchingDistance = (uint32_t)(val*1000);}
  inline void SetMinClusters(int val){minClusters = (uint8_t)val;}
  
  inline void SetScore(double val){score = val;}
  inline void SetNormalizedScore(double val){normalizedScore = val;}
  
  inline void SetBitChromosome(int i, uint64_t bits){bitChromosome[i] = bits;}
  
  inline void SetExecutionTime(double val){executionTime = val;}
  
  
  // Getters
  inline float  GetCriticalDistanceEE(){return criticalDistanceEE/1000.;}
  inline float  GetCriticalDistanceFH(){return criticalDistanceFH/1000.;}
  inline float  GetCriticalDistanceBH(){return criticalDistanceBH/1000.;}
  inline bool   GetDependSensor(){return dependSensor;}
  inline bool   GetReachedEE(){return reachedEE;}
  inline int    GetKernel(){return kernel;} // 0 - step, 1 - gaus, 2 - exp
  inline float  GetDeltacEE(){return deltacEE/1000.;}
  inline float  GetDeltacFH(){return deltacFH/1000.;}
  inline float  GetDeltacBH(){return deltacBH/1000.;}
  inline float  GetKappa(){return kappa/1000.;}
  inline float  GetEnergyMin(){return energyMin/100000.;}
  inline float  GetMatchingDistance(){return matchingDistance/1000.;}
  inline int    GetMinClusters(){return minClusters;}
  
  inline double  GetScore(){return score;}
  inline double  GetNormalizedScore(){return normalizedScore;}
  inline uint64_t GetBitChromosome(int i){return bitChromosome[i];}
  
  inline std::string GetConfigPath(){return configPath;}
  inline uint64_t GetUniqueID(){return uniqueID;}
  
  void SaveToBitChromosome();
  void ReadFromBitChromosome();
  
  void StoreInConfig();
  
  void Print();
  void CalculateScore();
  
  Chromosome* ProduceChildWith(Chromosome *partner);
  void BackToLimits();
private:
  // 1st chromosome
  uint32_t criticalDistanceEE;
  uint32_t criticalDistanceFH;
  
  // 2nd chromosome
  uint32_t criticalDistanceBH;
  bool dependSensor;          // 2 bits
  bool reachedEE;             // 2 bits
  
  // 3rd chromosome
  uint8_t kernel;
  uint32_t deltacEE;
  
  // 4th chromosome
  uint32_t deltacFH;
  uint32_t deltacBH;
  
  // 5th chromosome
  uint32_t kappa;
  uint32_t energyMin;
  
  // 6th chromosome
  uint32_t matchingDistance;
  uint8_t  minClusters;
  
  std::vector<uint64_t> bitChromosome; // 64 bit
  
  uint64_t uniqueID;
  std::string configPath;
  std::string clusteringOutputPath;
  
  ClusteringOutput clusteringOutput;
  double executionTime;
  double score;
  double normalizedScore;
  
  template<class T>
  void ShiftIntoChromosome(T value, int &shift, int chromoIndex);
  
  template<class T>
  void SetValueFromChromosome(T &value, int &shift, int chromoIndex);
  
  void Clusterize(std::string configPath);
};

#endif /* Chromosome_hpp */
