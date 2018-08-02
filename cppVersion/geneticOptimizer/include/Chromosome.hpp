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
  Chromosome();
  ~Chromosome();
  
  // 52 bits (1st choromosome)
  uint16_t criticalDistanceEE;
  uint16_t criticalDistanceFH;
  uint16_t criticalDistanceBH;
  bool dependSensor;          // 2 bits
  bool reachedEE;             // 2 bits
  
  // 56 bits (2nd chromosome)
  uint16_t kernel;
  uint16_t deltacEE;
  uint16_t deltacFH;
  uint16_t deltacBH;
  
  // 56 bits (3d chromosome)
  uint16_t kappa;
  uint16_t energyMin;
  uint16_t matchingDistance;
  uint16_t minClusters;
  
  uint64_t bitChromosome[3]; // 64 bit
  
  template<class T>
  void ShiftIntoChromosome(T value, int &shift, int chromoIndex);
  
  template<class T>
  void SetValueFromChromosome(T &value, int &shift, int chromoIndex);
  
  void SetupBitChromosome();
  void GetFromBitChromosome();
  
  void Print();
};

#endif /* Chromosome_hpp */
