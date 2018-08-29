//
//  Helpers.hpp
//
//  Created by Jeremi Niedziela on 13/06/2018.
//

#ifndef Helpers_h
#define Helpers_h

#include <cstdint>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <bitset>

#include <TMath.h>

enum EParam{
  kCriticalDistanceEE,
  kCriticalDistanceFH,
  kCriticalDistanceBH,
  kDeltacEE,
  kDeltacFH,
  kDeltacBH,
  kKappa,
  kEnergyThreshold,
  kMatchingDistance,
  kKernel,
  kNparams
};

const double paramMin[kNparams] = {
  0.0,  ///< min critical distance EE
  2.0,  ///< min critical distance FH (all dies below 2.0)
  5.0,  ///< min critical distance BH (all dies below 5.0)
  0.01, ///< min critical delta EE (too small - algorithm cannot find clusters)
  0.01, ///< min critical delta FH
  0.01, ///< min critical delta BH
  1.0,  ///< min kappa  // below 1.0 algorithm doesn't work
  2.0,  ///< min energy threshold
  0.0,  ///< min matching distance
 -0.49  ///< min kernel index
};

const double paramMax[kNparams] = {
  20.0, ///< max critical distance EE (almost all dies above 20)
  30.0, ///< max critical distance FH
  50.0, ///< max critical distance BH
  30.0, ///< max critical delta EE (above ~30 problems start to occur)
  30.0, ///< max critical delta FH
  40.0, ///< max critical delta BH
  150.0,///< max kappa (above 150 almost all creatures fail completely)
  10.0, ///< max energy threshold
  30.0, ///< max matching distance
  2.49  ///< max kernel index
};

const double paramStart[kNparams] = {
  2.00, ///< initial critical distance EE
  2.00, ///< initial critical distance FH
  2.00, ///< initial critical distance BH
  2.00, ///< initial critical delta EE
  2.00, ///< initial critical delta FH
  5.00, ///< initial critical delta BH
  9.00, ///< initial kappa
  3.00, ///< initial energy threshold
  5.00, ///< initial matching distance
  0.0   ///< initial kernel index (0 - step, 1 - gaus, 2 - exp)
};
  
inline const char* paramTitle[kNparams] = {
  "critDistEE",
  "critDistFH",
  "critDistBH",
  "deltaEE",
  "deltaFH",
  "deltaBH",
  "kappa",
  "eMin",
  "matchingDist",
  "kernel"
};

enum EDet {
  kEE,  ///< electromagneric endcap (silicon)
  kFH,  ///< front hadronic endcap (silicon)
  kBH   ///< back hadronic endcap (plastic)
};

inline double RandDouble(double min, double max)
{
  return min + static_cast<double>(rand()) /( static_cast<double>(RAND_MAX/(max-min)));
}

inline float RandFloat(float min, float max)
{
  return min + static_cast<float>(rand()) /( static_cast<float>(RAND_MAX/(max-min)));
}

inline int RandInt(int min, int max)
{
  return min + (rand() % static_cast<int>(max - min + 1));
}

inline bool RandBool()
{
  return rand() % 2;
}

inline void ReverseBit(uint64_t &bits, int pos)
{
  uint64_t mask = pow(2,pos);
  if((bits&mask)==0) bits = bits|mask;
  else bits = bits&(~mask);
}

inline void PrintBits(uint64_t bits)
{
  std::bitset<64> x(bits);
  std::cout<<x<<std::endl;
}

template<class T>
T BitSize(T&)
{
  return std::numeric_limits<T>::digits;
}

template<class T>
inline void UpdateParamValue(std::string configPath, std::string keyToReplace, T newValue){
  std::string tmpName = "tmp/tmp_"+std::to_string(RandInt(0, 100000))+".md";
  
//  std::cout<<"tmp name:"<<tmpName<<std::endl;
//  std::cout<<"config path:"<<configPath<<std::endl;
  
  std::ifstream is_file(configPath);
  std::ofstream outputFile;
  outputFile.open(tmpName);
  
  std::string line;
  while(getline(is_file, line)){
    std::istringstream is_line(line);
    std::string key;
    
    if( std::getline(is_line, key, ':')){
      std::string value;
      
      if(std::getline(is_line, value)){
        if(key==keyToReplace) outputFile<<key<<":\t"<<newValue<<std::endl;
        else                  outputFile<<line<<std::endl;
      }
      else outputFile<<line<<std::endl;
    }
    else  outputFile<<line<<std::endl;
  }
  outputFile.close();
//  std::cout<<"Executing: "<<("mv "+tmpName+" "+configPath)<<std::endl;
  system(("mv "+tmpName+" "+configPath).c_str());
}

template<class T>
void GetParamFromConfig(std::string configPath, std::string keyToFind, T &returnValue){
  std::ifstream is_file(configPath);
  
  std::string line;
  while(getline(is_file, line)){
    std::istringstream is_line(line);
    std::string key;
    if( std::getline(is_line, key, ':')){
      std::string value;
      if(std::getline(is_line, value)){
        if(key==keyToFind){
          std::istringstream is_value(value);
          returnValue = decltype(returnValue)(is_value);
          return;
        }
      }
    }
  }
  return;
}

struct ClusteringOutput {
  ClusteringOutput(){
    resolutionMean = 99999;
    resolutionSigma = 99999;
    separationMean = 99999;
    separationSigma = 99999;
    containmentMean = 99999;
    containmentSigma = 99999;
    deltaNclustersMean = 99999;
    deltaNclustersSigma = 99999;
    nRecoFailed = 99999;
    nCantMatchRecSim = 99999;
    nFakeRec = 99999;
  }
  
  void Print(){
    std::cout<<"Resolution:"<<resolutionMean<<" +/- "<<resolutionSigma<<std::endl;
    std::cout<<"Separation:"<<separationMean<<" +/- "<<separationSigma<<std::endl;
    std::cout<<"Containment:"<<containmentMean<<" +/- "<<containmentSigma<<std::endl;
    std::cout<<"Delta N clusters:"<<deltaNclustersMean<<" +/- "<<deltaNclustersSigma<<std::endl;
    std::cout<<"\% of event-layers where algo failed to find sim clusters:"<<nRecoFailed<<std::endl;
    std::cout<<"\% of event-layers were rec and sim clusters couldn't be matched:"<<nCantMatchRecSim<<std::endl;
    std::cout<<"\% of fake rec clusters (those that don't match any sim cluster):"<<nFakeRec<<std::endl;
  }
  
  double resolutionMean;
  double resolutionSigma;
  double separationMean;
  double separationSigma;
  double containmentMean;
  double containmentSigma;
  double deltaNclustersMean;
  double deltaNclustersSigma;
  double nRecoFailed;
  double nCantMatchRecSim;
  double nFakeRec;
  
};

inline ClusteringOutput ReadOutput(std::string fileName){
  std::ifstream is_file(fileName);
  std::string line;
  int iter=0;
  ClusteringOutput output = ClusteringOutput();
  
  while(getline(is_file, line)){
    std::istringstream is_line(line);
    if(iter==0) is_line >> output.resolutionMean;
    if(iter==1) is_line >> output.resolutionSigma;
    if(iter==2) is_line >> output.separationMean;
    if(iter==3) is_line >> output.separationSigma;
    if(iter==4) is_line >> output.containmentMean;
    if(iter==5) is_line >> output.containmentSigma;
    if(iter==6) is_line >> output.deltaNclustersMean;
    if(iter==7) is_line >> output.deltaNclustersSigma;
    if(iter==8) is_line >> output.nRecoFailed;
    if(iter==9) is_line >> output.nCantMatchRecSim;
    if(iter==10)is_line >> output.nFakeRec;
    iter++;
  }
  return output;
}

/// Calculate duration between two events
/// \param t0 Start time
/// \param t1 End time
/// \return Difference between events t0 and t1 in seconds
template<class T>
double duration(T t0,T t1)
{
  auto elapsed_secs = t1-t0;
  typedef std::chrono::duration<float> float_seconds;
  auto secs = std::chrono::duration_cast<float_seconds>(elapsed_secs);
  return secs.count();
}

/// Returns current time
inline std::chrono::time_point<std::chrono::steady_clock> now()
{
  return std::chrono::steady_clock::now();
//  return std::chrono::system_clock::now();
}

  
#endif /* Helpers_h */
