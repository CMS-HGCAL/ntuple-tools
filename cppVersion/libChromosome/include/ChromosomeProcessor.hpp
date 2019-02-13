//
//  ChromosomeProcessor.hpp
//
//  Created by Jeremi Niedziela on 13/02/2019.
//

#ifndef ChromosomeProcessor_hpp
#define ChromosomeProcessor_hpp

#include "GeneticHelpers.hpp"
#include "Chromosome.hpp"

/// Performs operations on chromosomes, such as crossover and mutation, but can also assign
/// score, or store chromosome as an Imaging Algo config file.
class ChromosomeProcessor {
public:
  /// Default constructor
  /// \_mutationChance Probability that each bit will flip during mutation step
  /// \_severityFactor Determines how easily will weak creatures die
  /// \_crossover      Crossover strategy (see ECrossover definition)
  ChromosomeProcessor(double _mutationChance,
                      double _severityFactor,
                      ECrossover _crossover);
  
  /// Default destructor
  ~ChromosomeProcessor();
  
  /// Returns a chromosome with random parameters within limits defined in GeneticHelpers.
  /// Parameters are NOT automatically stored in the bits, user has to call SaveToBitChromosome()
  /// on the resulting objects.
  std::shared_ptr<Chromosome> GetRandomChromosome();
  
  /// Sets score of the chromosome, which is determined based on the Imaging Algo output.
  /// The output and config files for this chromosome will be removed from the disk.
  void CalculateScore(std::shared_ptr<Chromosome> chromo);
  
  /// Performs crossover and mutation of two chromosomes: mom and dad.
  /// Will automatically store params in the bits of child chromosomes.
  /// \return Returns a pair of new chromosomes resulting from crossing of mom and dad.
  std::pair<std::shared_ptr<Chromosome>, std::shared_ptr<Chromosome>>
  CrossChromosomes(const std::shared_ptr<Chromosome> mom,
                   const std::shared_ptr<Chromosome> dad);
  
  /// Saves parameters of chromosome in form of an Imagin Algo config file.
  /// If path not explicitely specified, will use path of the chromosome.
  void StoreChromosomeInConfig(const std::shared_ptr<Chromosome> chromo,
                               std::string path="");
  
private:
  
  /// Swaps bits between a and b at randomly chosen position.
  /// If fixed==true, bits belonging to one parameter will never be cut and thus a and b
  /// will only exchange some parameters, without modifying their values.
  std::pair<uint64_t,uint64_t> SinglePointCrossover(uint64_t a, uint64_t b, bool fixed=false);
  
  double mutationChance;  ///< Probability that each bit will flip during mutation step
  double severityFactor;  ///< Determines how easily will weak creatures die
  ECrossover crossover;   ///< Crossover strategy (see ECrossover definition)
};

#endif /* ChromosomeProcessor_hpp */
