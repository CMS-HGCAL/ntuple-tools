//
//  ConfigurationManager.cpp
//
//  Created by Jeremi Niedziela on 02/07/2018.
//

#include "ConfigurationManager.hpp"

#include <TString.h>

#include <iostream>

using namespace std;

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

vector<double> ConfigurationManager::GetDeltac()
{
  vector<double> deltac;
  deltac.push_back(settings->GetValue("deltac_EE",2.0));
  deltac.push_back(settings->GetValue("deltac_FH",2.0));
  deltac.push_back(settings->GetValue("deltac_BH",2.0));
  return deltac;
}

int ConfigurationManager::GetMinClusters()
{
  return settings->GetValue("min_clusters",3);
}

double ConfigurationManager::GetEnergyMin()
{
  if(GetDependSensor()){
    return settings->GetValue("energy_min_with_sensor_dependance",3.0);
  }
  return settings->GetValue("energy_min_no_sensor_dependance",0.060);
}

double ConfigurationManager::GetKappa()
{
  if(GetDependSensor()){
    return settings->GetValue("kappa_with_sensor_dependance",9.0);
  }
  return settings->GetValue("kappa_no_sensor_dependance",10.0);
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






