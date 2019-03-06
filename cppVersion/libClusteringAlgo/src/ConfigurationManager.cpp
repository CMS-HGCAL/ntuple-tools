//
//  ConfigurationManager.cpp
//
//  Created by Jeremi Niedziela on 02/07/2018.
//

#include "ConfigurationManager.hpp"

#include <TString.h>

#include <iostream>

using namespace std;

shared_ptr<ConfigurationManager> ConfigurationManager::instance  = nullptr;

shared_ptr<ConfigurationManager> ConfigurationManager::Instance(string _configPath)
{
  if(instance){
    if(_configPath != ""){
      cout<<"WARNING - Configuration Manager was already created, but new config path was specified later on. Carefull - the new path will be ignored!!"<<endl;
    }
  }
  else{
    instance = make_shared<ConfigurationManager>(_configPath);
  }
  return instance;
}

shared_ptr<ConfigurationManager> ConfigurationManager::Instance(bool _dependSensor,
                                                                std::string _inputPath,
                                                                std::string _outputPath,
                                                                double _deltac_EE,
                                                                double _deltac_FH,
                                                                double _deltac_BH,
                                                                double _minEnergy,
                                                                double _criticalDistance_EE,
                                                                double _criticalDistance_FH,
                                                                double _criticalDistance_BH,
                                                                double _assignmentDistance_EE,
                                                                double _assignmentDistance_FH,
                                                                double _assignmentDistance_BH,
                                                                double _kappa,
                                                                int _verbosityLevel,
                                                                int _minNtuple,
                                                                int _maxNtuple,
                                                                int _minLayer,
                                                                int _maxLayer,
                                                                int _eventsPerTuple,
                                                                std::string _energyDensityFunction,
                                                                bool _reachedEEonly,
                                                                double _matchingMaxDistance,
                                                                std::string _scoreOutputPath,
                                                                bool _doHalo)
{
  if(instance){
    cout<<"WARNING - Configuration Manager was already created!!"<<endl;
  }
  else{
    instance = make_shared<ConfigurationManager>(_dependSensor,
                                                 _inputPath,
                                                 _outputPath,
                                                 _deltac_EE,
                                                 _deltac_FH,
                                                 _deltac_BH,
                                                 _minEnergy,
                                                 _criticalDistance_EE,
                                                 _criticalDistance_FH,
                                                 _criticalDistance_BH,
                                                 _assignmentDistance_EE,
                                                 _assignmentDistance_FH,
                                                 _assignmentDistance_BH,
                                                 _kappa,
                                                 _verbosityLevel,
                                                 _minNtuple,
                                                 _maxNtuple,
                                                 _minLayer,
                                                 _maxLayer,
                                                 _eventsPerTuple,
                                                 _energyDensityFunction,
                                                 _reachedEEonly,
                                                 _matchingMaxDistance,
                                                 _scoreOutputPath,
                                                 _doHalo
                                                 );
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
  
  criticalDistance_EE = settings->GetValue("critical_distance_EE",2.0);
  criticalDistance_FH = settings->GetValue("critical_distance_FH",2.0);
  criticalDistance_BH = settings->GetValue("critical_distance_BH",5.0);
  
  assignmentDistance_EE = settings->GetValue("assignment_distance_EE",2.0);
  assignmentDistance_FH = settings->GetValue("assignment_distance_FH",2.0);
  assignmentDistance_BH = settings->GetValue("assignment_distance_BH",2.0);
  
  kappa = settings->GetValue("kappa",9.0);
  
  verbosityLevel = settings->GetValue("verbosity_level",0);
  
  minNtuple = settings->GetValue("min_Ntuple",0);
  maxNtuple = settings->GetValue("max_Ntuple",0);
  minLayer = settings->GetValue("min_layer",0);
  maxLayer = settings->GetValue("max_layer",0);
  
  eventsPerTuple = settings->GetValue("analyze_events_per_tuple",99999);
  energyDensityFunction = settings->GetValue("energy_density_function","step");
  reachedEEonly = settings->GetValue("reachedEE_only",1);
  
  matchingMaxDistance = settings->GetValue("matching_max_distance",-1.0);
  
  scoreOutputPath = settings->GetValue("score_output_path","./output.txt");
  doHalo = settings->GetValue("do_halo",0);
}

ConfigurationManager::ConfigurationManager(bool _dependSensor,
                                           std::string _inputPath,
                                           std::string _outputPath,
                                           double _deltac_EE,
                                           double _deltac_FH,
                                           double _deltac_BH,
                                           double _minEnergy,
                                           double _criticalDistance_EE,
                                           double _criticalDistance_FH,
                                           double _criticalDistance_BH,
                                           double _assignmentDistance_EE,
                                           double _assignmentDistance_FH,
                                           double _assignmentDistance_BH,
                                           double _kappa,
                                           int _verbosityLevel,
                                           int _minNtuple,
                                           int _maxNtuple,
                                           int _minLayer,
                                           int _maxLayer,
                                           int _eventsPerTuple,
                                           std::string _energyDensityFunction,
                                           bool _reachedEEonly,
                                           double _matchingMaxDistance,
                                           std::string _scoreOutputPath,
                                           bool _doHalo)
{
  dependSensor = _dependSensor;
//  cout<<"depend sensor:"<<dependSensor<<endl;
  inputPath = _inputPath;
//  cout<<"input:"<<inputPath<<endl;
  outputPath = _outputPath;
//  cout<<"output:"<<outputPath<<endl;
  deltac_EE = _deltac_EE;
  deltac_FH = _deltac_FH;
  deltac_BH = _deltac_BH;
  minEnergy = _minEnergy;
  criticalDistance_EE = _criticalDistance_EE;
  criticalDistance_FH = _criticalDistance_FH;
  criticalDistance_BH = _criticalDistance_BH;
  assignmentDistance_EE = _assignmentDistance_EE;
  assignmentDistance_FH = _assignmentDistance_FH;
  assignmentDistance_BH = _assignmentDistance_BH;
  kappa = _kappa;
  verbosityLevel = _verbosityLevel;
  minNtuple = _minNtuple;
  maxNtuple = _maxNtuple;
  minLayer = _minLayer;
  maxLayer = _maxLayer;
  eventsPerTuple = _eventsPerTuple;
  energyDensityFunction = _energyDensityFunction;
  reachedEEonly = _reachedEEonly;
  matchingMaxDistance = _matchingMaxDistance;
  scoreOutputPath = _scoreOutputPath;
  doHalo = _doHalo;
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

double ConfigurationManager::GetAssignmentDistance(EDet det)
{
  if(det == kEE) return assignmentDistance_EE;
  if(det == kFH) return assignmentDistance_FH;
  return assignmentDistance_BH;
}

void ConfigurationManager::Print()
{
  cout<<"Current configuration:"<<endl;
  cout<<"Depend sensor:"<<dependSensor<<endl;
  cout<<"Reached EE only:"<<reachedEEonly<<endl;
  cout<<"Critical distance: "<<criticalDistance_EE<<"\t"<<criticalDistance_FH<<"\t"<<criticalDistance_BH<<endl;
  cout<<"Critical delta: "<<deltac_EE<<"\t"<<deltac_FH<<"\t"<<deltac_BH<<endl;
  cout<<"Kappa:"<<kappa<<endl;
  cout<<"Energy threshold:"<<minEnergy<<endl;
  cout<<"Kernel:"<<energyDensityFunction<<endl;
  cout<<"Matching distance:"<<matchingMaxDistance<<endl;
  cout<<"Data path:"<<inputPath<<endl;
}
