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
}

string ConfigurationManager::GetInputPath()
{
  TString path = settings->GetValue("input_path","unknown");
  
  if(path == "unknown"){
    cout<<"WARNING -- could not read input data path from config file:"<<configPath<<endl;
    return "";
  }
  return path.Data();
}

string ConfigurationManager::GetOutputPath()
{
  TString path = settings->GetValue("output_path","unknown");
  
  if(path == "unknown"){
    cout<<"WARNING -- could not read output data path from config file:"<<configPath<<endl;
    return "";
  }
  return path.Data();
}

bool ConfigurationManager::GetDependSensor()
{
  return settings->GetValue("depend_sensor",true);
}

double ConfigurationManager::GetDeltac(EDet det)
{
  
  if(det == kEE) return settings->GetValue("deltac_EE",2.0);
  if(det == kFH) return settings->GetValue("deltac_FH",2.0);
  return settings->GetValue("deltac_BH",2.0);
}

int ConfigurationManager::GetMinClusters()
{
  return settings->GetValue("min_clusters",3);
}

double ConfigurationManager::GetEnergyMin()
{
  return settings->GetValue("energy_min",3.0);
}

double ConfigurationManager::GetCriticalDistance(EDet det)
{
  if(det == kEE) return settings->GetValue("critical_distance_EE",2.0);
  if(det == kFH) return settings->GetValue("critical_distance_FH",2.0);
  return settings->GetValue("critical_distance_BH",5.0);
}

double ConfigurationManager::GetKappa()
{
  return settings->GetValue("kappa",9.0);
}

int ConfigurationManager::GetVerbosityLevel()
{
  return settings->GetValue("verbosity_level",0);
}

int ConfigurationManager::GetMinNtuple()
{
  return settings->GetValue("min_Ntuple",0);
}

int ConfigurationManager::GetMaxNtuple()
{
  return settings->GetValue("max_Ntuple",0);
}

int ConfigurationManager::GetMinLayer()
{
  return settings->GetValue("min_layer",0);
}

int ConfigurationManager::GetMaxLayer()
{
  return settings->GetValue("max_layer",0);
}

int ConfigurationManager::GetMaxEventsPerTuple()
{
  return settings->GetValue("analyze_events_per_tuple",99999);
}






