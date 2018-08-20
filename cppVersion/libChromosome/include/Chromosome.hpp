//
//  Chromosome.hpp
//
//  Created by Jeremi Niedziela on 02/08/2018.
//

#ifndef Chromosome_hpp
#define Chromosome_hpp

#include "GeneticHelpers.hpp"

template <typename T,typename K>
class AlgoParam;

class Chromosome
{
public:
  enum ECrossover{
    kUniform,       ///< each bit has a chance to be exchaned between parents
    kSinglePoint,   ///< one chromosome gets crossed, the rest stays the same or is exchaned intact
    kMultiPoint     ///< each chromosome is crossed at a random point
  };
  
  enum EParam{
    kCriticalDistanceEE,
    kCriticalDistanceFH,
    kCriticalDiscanceBH
  };
  
  /// Default constructor
  Chromosome();
  
  /// Default destructor
  ~Chromosome();
  
  static Chromosome* GetRandom();
  
  // Setters
  inline void SetCriticalDistanceEE(double val){
    criticalDistanceEE = DoubleInRangeToUint(val, criticalDistanceEEmin, criticalDistanceEEmax, criticalDistanceEE);
  }
  inline void SetCriticalDistanceFH(double val){
    criticalDistanceFH = DoubleInRangeToUint(val, criticalDistanceFHmin, criticalDistanceFHmax, criticalDistanceFH);
  }
  inline void SetCriticalDistanceBH(double val){
    criticalDistanceBH = DoubleInRangeToUint(val, criticalDistanceBHmin, criticalDistanceBHmax, criticalDistanceBH);
  }
  inline void SetDependSensor(bool val){
    dependSensor = DoubleInRangeToUint(val, 0, 1, dependSensor);
  }
  inline void SetReachedEE(bool val){
    reachedEE = DoubleInRangeToUint(val, 0, 1, reachedEE);
  }

  inline void SetDeltacEE(double val){
    deltacEE = DoubleInRangeToUint(val, deltacEEmin, deltacEEmax, deltacEE);
  }
  inline void SetDeltacFH(double val){
    deltacFH = DoubleInRangeToUint(val, deltacFHmin, deltacFHmax, deltacFH);
  }
  inline void SetDeltacBH(double val){
    deltacBH = DoubleInRangeToUint(val, deltacBHmin, deltacBHmax, deltacBH);
  }
  inline void SetKappa(double val){
    kappa = DoubleInRangeToUint(val, kappaMin, kappaMax, kappa);
  }
  inline void SetEnergyMin(double val){
    energyMin = DoubleInRangeToUint(val, energyThresholdMin, energyThresholdMax, energyMin);
  }
  inline void SetMatchingDistance(double val){
    matchingDistance = DoubleInRangeToUint(val, matchingDistanceMin, matchingDistanceMax, matchingDistance);
  }
  inline void SetKernel(int val){// 0 - step, 1 - gaus, 2 - exp
    kernel = DoubleInRangeToUint(val, kernelMin, kernelMax, kernel);
  }
  inline void SetMinClusters(int val){
    minClusters = DoubleInRangeToUint(val, minClustersMin, minClustersMax, minClusters);
  }
  
  inline void SetBitChromosome(int i, uint64_t bits){bitChromosome[i] = bits;}
  
  inline void SetScore(double val){score = val;}
  inline void SetNormalizedScore(double val){normalizedScore = val;}
  inline void SetExecutionTime(double val){executionTime = val;}
  inline void SetMutationChance(double val){mutationChance = val;}
  inline void SetSeverityFactor(double val){severityFactor = val;}
  inline void SetCrossover(ECrossover val){crossover = val;}
  
  // Getters
  inline double GetCriticalDistanceEE(){
    return UintToDoubleInRange(criticalDistanceEE, criticalDistanceEEmin, criticalDistanceEEmax);
  }
  inline double  GetCriticalDistanceFH(){
    return UintToDoubleInRange(criticalDistanceFH, criticalDistanceFHmin, criticalDistanceFHmax);
  }
  inline double  GetCriticalDistanceBH(){
    return UintToDoubleInRange(criticalDistanceBH, criticalDistanceBHmin, criticalDistanceBHmax);
  }
  inline bool   GetDependSensor(){
    double val = UintToDoubleInRange(dependSensor, 0, 1);
    return static_cast<bool>(val);
  }
  inline bool   GetReachedEE(){
    double val = UintToDoubleInRange(reachedEE, 0, 1);
    return static_cast<bool>(val);
  }
  inline double  GetDeltacEE(){
    return UintToDoubleInRange(deltacEE, deltacEEmin, deltacEEmax);
  }
  inline double  GetDeltacFH(){
    return UintToDoubleInRange(deltacFH, deltacFHmin, deltacFHmax);
  }
  inline double  GetDeltacBH(){
    return UintToDoubleInRange(deltacBH, deltacBHmin, deltacBHmax);
  }
  inline double  GetKappa(){
    return UintToDoubleInRange(kappa, kappaMin, kappaMax);
  }
  inline double  GetEnergyMin(){
    return UintToDoubleInRange(energyMin, energyThresholdMin, energyThresholdMax);
  }
  inline double  GetMatchingDistance(){
    return UintToDoubleInRange(matchingDistance, matchingDistanceMin, matchingDistanceMax);
  }
  inline int GetKernel(){// 0 - step, 1 - gaus, 2 - exp
    double val = UintToDoubleInRange(kernel, kernelMin, kernelMax);
    return round(val);
  }
  inline int GetMinClusters(){
    double val = UintToDoubleInRange(minClusters, minClustersMin, minClustersMax);
    return round(val);
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
  // 1st chromosome
  uint16_t criticalDistanceEE;
  uint16_t criticalDistanceFH;
  uint16_t criticalDistanceBH;
  uint16_t reachedEE;
  
  // 3rd chromosome
  uint8_t kernel;
  uint16_t deltacEE;
  uint16_t deltacFH;
  uint16_t deltacBH;
  
  // 5th chromosome
  uint16_t kappa;
  uint16_t energyMin;
  uint16_t matchingDistance;
  uint16_t  minClusters;
  
  std::vector<uint64_t> bitChromosome; // 64 bit
  
  
  // other variables not stored in bit chromosome
  uint64_t uniqueID;
  std::string configPath;
  std::string clusteringOutputPath;
  bool dependSensor;
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
  
  template<class T>
  T DoubleInRangeToUint(double val, double min, double max, T output);
  
  template<class T>
  double UintToDoubleInRange(T input, double min, double max);
  
  std::vector<uint64_t> SinglePointCrossover(uint64_t a, uint64_t b);
  int BackToLimits();
};

template <class T,class K>
class AlgoParam{
public:
  AlgoParam(){};
  ~AlgoParam(){};
  
  T bitValue;
  
  K realValue;
  K min;
  K max;
};


#endif /* Chromosome_hpp */
