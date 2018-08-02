//
//  Chromosome.hpp
//
//  Created by Jeremi Niedziela on 02/08/2018.
//

#ifndef Chromosome_hpp
#define Chromosome_hpp

#include "Helpers.hpp"

class Chromosome
{
public:
  /// Default constructor
  Chromosome();
  
  /// Copy contructor. Copies only the bit chromosome.
  /// Call ReadFromBitChromosome() to populate other fields.
  /// uniqueID of the copy will be different.
  Chromosome(Chromosome &c);
  
  /// Default destructor
  ~Chromosome();
  
  static Chromosome GetRandom();
  
  // Setters
  inline void SetCriticalDistanceEE(float val){criticalDistanceEE = (uint16_t)(val*1000);}
  inline void SetCriticalDistanceFH(float val){criticalDistanceFH = (uint16_t)(val*1000);}
  inline void SetCriticalDistanceBH(float val){criticalDistanceBH = (uint16_t)(val*1000);}
  inline void SetDependSensor(bool val){dependSensor = val;}
  inline void SetReachedEE(bool val){reachedEE = val;}
  inline void SetKernel(int val){kernel = (uint16_t)val;} // 0 - step, 1 - gaus, 2 - exp
  inline void SetDeltacEE(float val){deltacEE = (uint16_t)(val*1000);}
  inline void SetDeltacFH(float val){deltacFH = (uint16_t)(val*1000);}
  inline void SetDeltacBH(float val){deltacBH = (uint16_t)(val*1000);}
  inline void SetKappa(float val){kappa = (uint16_t)(val*1000);}
  inline void SetEnergyMin(float val){energyMin = (uint16_t)(val*100000);}
  inline void SetMatchingDistance(float val){matchingDistance = (uint16_t)(val*1000);}
  inline void SetMinClusters(float val){minClusters = (uint16_t)val;}
  
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
  inline float  GetMinClusters(){return minClusters;}
  
  void SaveToBitChromosome();
  void ReadFromBitChromosome();
  
  void Print();
  void StoreInConfig();
private:
  // 52 bits (1st choromosome)
  uint16_t criticalDistanceEE;
  uint16_t criticalDistanceFH;
  uint16_t criticalDistanceBH;
  bool dependSensor;          // 2 bits
  bool reachedEE;             // 2 bits
  
  // 64 bits (2nd chromosome)
  uint16_t kernel;
  uint16_t deltacEE;
  uint16_t deltacFH;
  uint16_t deltacBH;
  
  // 64 bits (3d chromosome)
  uint16_t kappa;
  uint16_t energyMin;
  uint16_t matchingDistance;
  uint16_t minClusters;
  
  uint64_t bitChromosome[3]; // 64 bit
  
  uint64_t uniqueID;
  
  template<class T>
  void ShiftIntoChromosome(T value, int &shift, int chromoIndex);
  
  template<class T>
  void SetValueFromChromosome(T &value, int &shift, int chromoIndex);
  
  
};

#endif /* Chromosome_hpp */
