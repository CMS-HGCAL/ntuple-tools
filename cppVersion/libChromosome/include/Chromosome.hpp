//
//  Chromosome.hpp
//
//  Created by Jeremi Niedziela on 02/08/2018.
//

#ifndef Chromosome_hpp
#define Chromosome_hpp

#include "GeneticHelpers.hpp"

/// Representation of Imaging Algo parameters that can be used for optimization with use
/// of a genetic algorithm. Allows to save and read double parameters in a vector of 64bit
/// unsigned integers.
class Chromosome
{
public:  
  /// Default constructor
  Chromosome();
  
  /// Default destructor
  ~Chromosome();
  
  /// Prints basic information about this creature
  void Print();
  
  /// Sets value of parameter. It will NOT be automatically added to bitChromosome - user has to
  /// call SaveToBitChromosome() manually.
  /// \par Type of parameter
  /// \val Value of the parameter
  void SetParam(EParam par, double val);
  
  /// Fixes value of parameter. It will NOT be automatically added to bitChromosome - user has to
  /// call SaveToBitChromosome() manually.
  /// \par Type of parameter
  /// \val Value of the parameter
  void FixParam(EParam par, double val);
  
  /// Returns value of the parameter. This should be the only way of accessing parameters!
  /// Don't try to access it directly, even within this class.
  double GetParam(EParam par);
  
  /// Saves all parameters in a vector of 64bit integers that will be used for crossover and mutation.
  void  SaveToBitChromosome();
  
  /// Restores all parameters from a vector of 64bit integers used for crossover and mutation.
  void  ReadFromBitChromosome();
  
  // Trivial setters
  inline void SetNormalizedScore(double val){normalizedScore = val;}
  inline void SetExecutionTime(double val){executionTime = val;}
  inline void SetClusteringOutput(ClusteringOutput val){clusteringOutput = val;}
  
  // Trivial getters
  inline double       GetScore(){return score;}
  inline double       GetNormalizedScore(){return normalizedScore;}
  inline uint64_t     GetUniqueID(){return uniqueID;}
  inline std::string  GetConfigPath(){return configPath;}
  inline std::string  GetClusteringOutputPath(){return clusteringOutputPath;}
  
  
private:
  std::vector<uint16_t> params;       ///< Vector of chromosome parameters (see EParam)
  std::vector<bool>     isParamFixed; ///< Tells if a parameter should be fixed
  
  std::vector<uint64_t> bitChromosome;///< Bits representing all parameters (used for crossover
                                      ///  and mutation
  
  // other variables not stored in bit chromosome
  uint64_t    uniqueID;             ///< Unique identifier of this creature
  std::string configPath;           ///< Path to file in which this creature will be stored
  std::string clusteringOutputPath; ///< Path to file to which results of clustering with params
                                    ///  of this creature should be stored
  
  ClusteringOutput clusteringOutput;
  double executionTime;
  double score;
  double normalizedScore;
  
  template<class T>
  void ShiftIntoChromosome(T value, int &shift, int chromoIndex);
  
  template<class T>
  void SetValueFromChromosome(T &value, int &shift, int chromoIndex);
  
  int BackToLimits();
  
  friend class ChromosomeProcessor;
};

#endif /* Chromosome_hpp */
