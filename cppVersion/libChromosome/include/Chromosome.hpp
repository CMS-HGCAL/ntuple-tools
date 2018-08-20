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
    kMultiPoint     ///< each chromosome is crossed at a random point
  };
  
  /// Default constructor
  Chromosome();
  
  /// Default destructor
  ~Chromosome();
  
  static Chromosome* GetRandom();
  
  // Setters
  inline void SetCriticalDistanceEE(double val){
    criticalDistanceEE = (uint16_t)(std::numeric_limits<uint16_t>::max()/(criticalDistanceEEmax-criticalDistanceEEmin)*(val-criticalDistanceEEmin));
  }
  inline void SetCriticalDistanceFH(double val){
    criticalDistanceFH = (uint16_t)(std::numeric_limits<uint16_t>::max()/(criticalDistanceFHmax-criticalDistanceFHmin)*(val-criticalDistanceFHmin));
  }
  inline void SetCriticalDistanceBH(double val){
    criticalDistanceBH = (uint16_t)(std::numeric_limits<uint16_t>::max()/(criticalDistanceBHmax-criticalDistanceBHmin)*(val-criticalDistanceBHmin));
  }
  inline void SetDependSensor(bool val){dependSensor = val;}
  inline void SetReachedEE(bool val){reachedEE = val;}
  

  inline void SetDeltacEE(double val){
    deltacEE = (uint16_t)(std::numeric_limits<uint16_t>::max()/(deltacEEmax-deltacEEmin)*(val-deltacEEmin));
  }
  inline void SetDeltacFH(double val){
    deltacFH = (uint16_t)(std::numeric_limits<uint16_t>::max()/(deltacFHmax-deltacFHmin)*(val-deltacFHmin));
  }
  inline void SetDeltacBH(double val){
    deltacBH = (uint16_t)(std::numeric_limits<uint16_t>::max()/(deltacBHmax-deltacBHmin)*(val-deltacBHmin));
  }
  inline void SetKappa(double val){
    kappa = (uint16_t)(std::numeric_limits<uint16_t>::max()/(kappaMax-kappaMin)*(val-kappaMin));
  }
  inline void SetEnergyMin(double val){
    energyMin = (uint16_t)(std::numeric_limits<uint16_t>::max()/(energyThresholdMax-energyThresholdMin)*(val-energyThresholdMin));
  }
  inline void SetMatchingDistance(double val){
    matchingDistance = (uint16_t)(std::numeric_limits<uint16_t>::max()/(matchingDistanceMax-matchingDistanceMin)*(val-matchingDistanceMin));
  }
  
  inline void SetKernel(int val){// 0 - step, 1 - gaus, 2 - exp
    kernel = (uint8_t)(std::numeric_limits<uint8_t>::max()/(kernelMax-kernelMin)*((double)val-kernelMin));
  }
  inline void SetMinClusters(int val){
    minClusters = (uint8_t)(std::numeric_limits<uint8_t>::max()/(minClustersMax-minClustersMin)*((double)val-minClustersMin));
  }
  
  inline void SetBitChromosome(int i, uint64_t bits){bitChromosome[i] = bits;}
  
  inline void SetScore(double val){score = val;}
  inline void SetNormalizedScore(double val){normalizedScore = val;}
  inline void SetExecutionTime(double val){executionTime = val;}
  inline void SetMutationChance(double val){mutationChance = val;}
  inline void SetSeverityFactor(double val){severityFactor = val;}
  inline void SetCrossover(ECrossover val){crossover = val;}
  
  // Getters
  inline double  GetCriticalDistanceEE(){
    return criticalDistanceEEmin + criticalDistanceEE*(criticalDistanceEEmax-criticalDistanceEEmin)/std::numeric_limits<uint16_t>::max();
  }
  inline double  GetCriticalDistanceFH(){
    return criticalDistanceFHmin + criticalDistanceFH*(criticalDistanceFHmax-criticalDistanceFHmin)/std::numeric_limits<uint16_t>::max();
  }
  inline double  GetCriticalDistanceBH(){
    return criticalDistanceBHmin + criticalDistanceBH*(criticalDistanceBHmax-criticalDistanceBHmin)/std::numeric_limits<uint16_t>::max();
  }
  inline bool   GetDependSensor(){return dependSensor;}
  inline bool   GetReachedEE(){return reachedEE;}

  inline double  GetDeltacEE(){
    return deltacEEmin + deltacEE*(deltacEEmax-deltacEEmin)/std::numeric_limits<uint16_t>::max();
  }
  inline double  GetDeltacFH(){
    return deltacFHmin + deltacFH*(deltacFHmax-deltacFHmin)/std::numeric_limits<uint16_t>::max();
  }
  inline double  GetDeltacBH(){
    return deltacBHmin + deltacBH*(deltacBHmax-deltacBHmin)/std::numeric_limits<uint16_t>::max();
  }
  inline double  GetKappa(){
    return kappaMin + kappa*(kappaMax-kappaMin)/std::numeric_limits<uint16_t>::max();
  }
  inline double  GetEnergyMin(){
    return energyThresholdMin + energyMin*(energyThresholdMax-energyThresholdMin)/std::numeric_limits<uint16_t>::max();
  }
  inline double  GetMatchingDistance(){
    return matchingDistanceMin + matchingDistance*(matchingDistanceMax-matchingDistanceMin)/std::numeric_limits<uint16_t>::max();
  }
  inline int GetKernel(){// 0 - step, 1 - gaus, 2 - exp
    return round(kernelMin + kernel*(kernelMax-kernelMin)/std::numeric_limits<uint8_t>::max());
  }
  inline int GetMinClusters(){
    return round(minClustersMin + minClusters*(minClustersMax-minClustersMin)/std::numeric_limits<uint8_t>::max());
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
  bool reachedEE;             // 2 bits
  
  // 3rd chromosome
  uint8_t kernel;
  uint16_t deltacEE;
  uint16_t deltacFH;
  uint16_t deltacBH;
  
  // 5th chromosome
  uint16_t kappa;
  uint16_t energyMin;
  uint16_t matchingDistance;
  uint8_t  minClusters;
  
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
  
  void Clusterize(std::string configPath);
  
  std::vector<uint64_t> SinglePointCrossover(uint64_t a, uint64_t b);
  int BackToLimits();
};

#endif /* Chromosome_hpp */
