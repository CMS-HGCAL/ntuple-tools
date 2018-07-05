//
//  ConfigurationManager.hpp
//
//  Created by Jeremi Niedziela on 02/07/2018.
//

#ifndef ConfigurationManager_h
#define ConfigurationManager_h

#include <TEnv.h>

#include <string>

class ConfigurationManager {
public:
  static ConfigurationManager* Instance(std::string _configPath="");

  /// Returns the path to the input files with ntuples
  std::string GetInputPath();
  
  /// Returns a path under which the output files should be stored
  std::string GetOutputPath();
  
  /// Tells if sensor dependance should be taken into account
  bool GetDependSensor();
  
  /// Returns radius of circle in which to look for hexels (cartesian coordiantes in cm, per detector: EE, FH, BH)
  std::vector<double> GetDeltac();
  
  /// Returns requested min number of 2D clusters + 1
  int GetMinClusters();
  
  /// Returns cut on energy (minimum)
  double GetEnergyMin();
  
  double GetKappa();
  
  /// Verbosity level of the algo (0 - no output, 1 - basic output, 2 - debug)
  int GetVerbosityLevel();
  
  /// Start running algo from this Ntuple
  int GetMinNtuple();
  
  /// Last Ntuple to run the algo on
  int GetMaxNtuple();
  
  /// Start algorithm in this layer
  int GetMinLayer();
  
  /// Finish alsorithm in this layer
  int GetMaxLayer();

  /// Stop analyzing each ntuple after that many events
  int GetMaxEventsPerTuple();

private:
  ConfigurationManager(std::string _configPath);
  ~ConfigurationManager(){};
  
  static ConfigurationManager *instance;
  
  TEnv *settings;           ///< Object storing current configuration
  std::string configPath;   ///< Current path to the configuration
};

#endif /* ConfigurationManager_h */