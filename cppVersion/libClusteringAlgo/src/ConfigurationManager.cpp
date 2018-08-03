//
//  ConfigurationManager.cpp
//
//  Created by Jeremi Niedziela on 02/07/2018.
//

#include "ConfigurationManager.hpp"

#include <TString.h>

#include <iostream>

using namespace std;

ConfigurationManager* ConfigurationManager::instance  = nullptr;

ConfigurationManager* ConfigurationManager::Instance(string _configPath)
{
  if(instance){
    if(_configPath != ""){
      cout<<"WARNING - Configuration Manager was already created, but new config path was specified later on. Carefull - the new path will be ignored!!"<<endl;
    }
  }
  else{
    instance = new ConfigurationManager(_configPath);
  }
  return instance;
}

ConfigurationManager::ConfigurationManager(string _configPath) :
configPath(_configPath)
{
  settings = new TEnv();
  
  if(settings->ReadFile(configPath.c_str(), kEnvUser) < 0){
    cout<<"ERROR - could not load config file:"<<configPath<<endl;
    exit(0);
  }
  
  TString inpath = settings->GetValue("input_path","unknown");
  
  if(inpath == "unknown"){
    cout<<"WARNING -- could not read input data path from config file:"<<configPath<<endl;
    inputPath = "";
  }
  inputPath = inpath.Data();
  
  TString outpath = settings->GetValue("output_path","unknown");
  
  if(outpath == "unknown"){
    cout<<"WARNING -- could not read output data path from config file:"<<configPath<<endl;
    outputPath = "";
  }
  outputPath = outpath.Data();
  
  dependSensor = settings->GetValue("depend_sensor",true);
  
  deltac_EE = settings->GetValue("deltac_EE",2.0);
  deltac_FH = settings->GetValue("deltac_FH",2.0);
  deltac_BH = settings->GetValue("deltac_BH",2.0);
  
  minEnergy = settings->GetValue("energy_min",3.0);
  minClusters = settings->GetValue("min_clusters",3);
  
  criticalDistance_EE = settings->GetValue("critical_distance_EE",2.0);
  criticalDistance_FH = settings->GetValue("critical_distance_FH",2.0);
  criticalDistance_BH = settings->GetValue("critical_distance_BH",5.0);
  
  kappa = settings->GetValue("kappa",9.0);
  
  verbosityLevel = settings->GetValue("verbosity_level",0);
  
  minNtuple = settings->GetValue("min_Ntuple",0);
  maxNtuple = settings->GetValue("max_Ntuple",0);
  minLayer = settings->GetValue("min_layer",0);
  maxLayer = settings->GetValue("max_layer",0);
  
  eventsPerTuple = settings->GetValue("analyze_events_per_tuple",99999);
  energyDensityFunction = settings->GetValue("energy_density_function","step");
  reachedEEonly = settings->GetValue("reachedEE_only",1);
  
  matchingMaxDistance = settings->GetValue("matching_max_distance",-1);
  
  scoreOutputPath = settings->GetValue("score_output_path","./output.txt");
}

double ConfigurationManager::GetDeltac(EDet det)
{
  if(det == kEE) return deltac_EE;
  if(det == kFH) return deltac_FH;
  return deltac_BH;
}

double ConfigurationManager::GetCriticalDistance(EDet det)
{
  if(det == kEE) return criticalDistance_EE;
  if(det == kFH) return criticalDistance_FH;
  return criticalDistance_BH;
}
