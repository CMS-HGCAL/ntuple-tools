//
//  ConfigurationManager.hpp
//
//  Created by Jeremi Niedziela on 02/07/2018.
//

#ifndef ConfigurationManager_h
#define ConfigurationManager_h

#include "Helpers.hpp"

#include <TEnv.h>

#include <string>

class ConfigurationManager {
public:  
  static ConfigurationManager* Instance(std::string _configPath="");

  /// Returns the path to the input files with ntuples
  inline std::string GetInputPath(){return inputPath;}
  
  /// Returns a path under which the output files should be stored
  inline std::string GetOutputPath(){return outputPath;}
  
  /// Tells if sensor dependance should be taken into account
  inline bool GetDependSensor(){return dependSensor;}
  
  /// Returns radius of circle in which to look for hexels (cartesian coordiantes in cm, per detector: EE, FH, BH)
  double GetDeltac(EDet det);
  
  /// Returns requested min number of 2D clusters + 1
  inline int GetMinClusters(){return minClusters;}
  
  /// Returns cut on energy (minimum)
  inline double GetEnergyMin(){return minEnergy;}
  
  /// Returns maximum distance to include hits in the energy density calculation.
  double GetCriticalDistance(EDet det);
  
  /// Returns kappa parameter defining minimum energy density to consider a hit as a cluster seed
  /// Critical energy density rho_c = max(rho)/kappa
  inline double GetKappa(){return kappa;}
  
  /// Verbosity level of the algo (0 - no output, 1 - basic output, 2 - debug)
  inline int GetVerbosityLevel(){return verbosityLevel;}
  
  /// Start running algo from this Ntuple
  inline int GetMinNtuple(){return minNtuple;}
  
  /// Last Ntuple to run the algo on
  inline int GetMaxNtuple(){return maxNtuple;}
  
  /// Start algorithm in this layer
  inline int GetMinLayer(){return minLayer;}
  
  /// Finish alsorithm in this layer
  inline int GetMaxLayer(){return maxLayer;}

  /// Stop analyzing each ntuple after that many events
  inline int GetMaxEventsPerTuple(){return eventsPerTuple;}

  /// Returns a name of the function to be used for the energy density calculation
  inline std::string GetEnergyDensityFunction(){return energyDensityFunction;}
  
  /// Tells if events with at least one particle converting before EE should be skipped or not
  inline bool GetReachedEEonly(){return reachedEEonly;}
  
  /// Tells how far sim cluster can be from a rec cluster to be matched with it
  inline double GetMachingMaxDistance(){return matchingMaxDistance;}
  
  /// Where to store resulting algo characteristics (resolution, separation, containment)
  inline std::string GetScoreOutputPath(){return scoreOutputPath;}
  
private:
  ConfigurationManager(std::string _configPath);
  ~ConfigurationManager(){};
  
  static ConfigurationManager *instance;
  
  TEnv *settings;           ///< Object storing current configuration
  std::string configPath;   ///< Current path to the configuration
  
  bool dependSensor;
  std::string inputPath;
  std::string outputPath;
  double deltac_EE;
  double deltac_FH;
  double deltac_BH;
  double minEnergy;
  int minClusters;
  double criticalDistance_EE;
  double criticalDistance_FH;
  double criticalDistance_BH;
  double kappa;
  int verbosityLevel;
  int minNtuple;
  int maxNtuple;
  int minLayer;
  int maxLayer;
  int eventsPerTuple;
  std::string energyDensityFunction;
  bool reachedEEonly;
  double matchingMaxDistance;
  std::string scoreOutputPath;
};

#endif /* ConfigurationManager_h */
