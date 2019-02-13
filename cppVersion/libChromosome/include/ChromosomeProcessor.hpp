//
//  ChromosomeProcessor.hpp
//
//  Created by Jeremi Niedziela on 13/02/2019.
//

#ifndef ChromosomeProcessor_hpp
#define ChromosomeProcessor_hpp

#include "GeneticHelpers.hpp"
#include "Chromosome.hpp"

class ChromosomeProcessor {
public:
  /// Default constructor
  ChromosomeProcessor(double _mutationChance,
                      double _severityFactor,
                      ECrossover _crossover);
  
  /// Default destructor
  ~ChromosomeProcessor();
  
  std::shared_ptr<Chromosome> GetRandomChromosome();
  
  void CalculateScore(std::shared_ptr<Chromosome> chromo);
  
  std::pair<std::shared_ptr<Chromosome>, std::shared_ptr<Chromosome>>
  CrossChromosomes(const std::shared_ptr<Chromosome> mom,
                   const std::shared_ptr<Chromosome> dad);
  
  void StoreChromosomeInConfig(const std::shared_ptr<Chromosome> chromo,
                               std::string path="");
  
private:
  
  std::pair<uint64_t,uint64_t> SinglePointCrossover(uint64_t a, uint64_t b, bool fixed=false);
  
  double mutationChance;
  double severityFactor;   // larger the value, more easily population members will die
  ECrossover crossover;
};

#endif /* ChromosomeProcessor_hpp */
