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

inline uint64_t reverse_bit (uint64_t bits, int pos)
{
  if(bits&(2^pos)==0) return (bits|(2^pos));
  return (bits&(~2^pos));
}


inline void PrintBits(uint64_t bits)
{
  std::bitset<80> x(bits);
  std::cout<<x<<std::endl;
}

template<class T>
T BitSize(T&)
{
  return std::numeric_limits<T>::digits;
}

template<class T>
inline void UpdateParamValue(std::string configPath, std::string keyToReplace, T newValue){
  std::string tmpName = "tmp/tmp_"+std::to_string(RandInt(0, 10000000))+".md";
  
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
void GetParamFomeConfig(std::string configPath, std::string keyToFind, T &returnValue){
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

//inline double GetParamFomeConfig(std::string configPath, std::string keyToFind){
//  std::ifstream is_file(configPath);
//
//  std::string line;
//  while(getline(is_file, line)){
//    std::istringstream is_line(line);
//    std::string key;
//    if( std::getline(is_line, key, ':')){
//      std::string value;
//      if(std::getline(is_line, value)){
//        if(key==keyToFind){
//          return stod(value);
//        }
//      }
//    }
//  }
//  return -999999999;
//}

struct ClusteringOutput {
  ClusteringOutput(){
    resolutionMean = 99999;
    resolutionSigma = 99999;
    separationMean = 99999;
    separationSigma = 99999;
    containmentMean = 99999;
    containmentSigma = -99999;
  }
  
  void Print(){
    std::cout<<"Resolution:"<<resolutionMean<<" +/- "<<resolutionSigma<<std::endl;
    std::cout<<"Separation:"<<separationMean<<" +/- "<<separationSigma<<std::endl;
    std::cout<<"Containment:"<<containmentMean<<" +/- "<<containmentSigma<<std::endl;
  }
  
  double resolutionMean;
  double resolutionSigma;
  double separationMean;
  double separationSigma;
  double containmentMean;
  double containmentSigma;
};

inline ClusteringOutput ReadOutput(std::string fileName){
  std::ifstream is_file(fileName);
  std::string line;
  int iter=0;
  ClusteringOutput output;
  
  while(getline(is_file, line)){
    std::istringstream is_line(line);
    if(iter==0) is_line >> output.resolutionMean;
    if(iter==1) is_line >> output.resolutionSigma;
    if(iter==2) is_line >> output.separationMean;
    if(iter==3) is_line >> output.separationSigma;
    if(iter==4) is_line >> output.containmentMean;
    if(iter==5) is_line >> output.containmentSigma;
    iter++;
  }
  return output;
}

/// Calculate duration between two events
/// \param t0 Start time
/// \param t1 End time
/// \return Difference between events t0 and t1 in seconds
inline double duration(std::chrono::time_point<std::chrono::system_clock> t0,
                       std::chrono::time_point<std::chrono::system_clock> t1)
{
  auto elapsed_secs = t1-t0;
  typedef std::chrono::duration<float> float_seconds;
  auto secs = std::chrono::duration_cast<float_seconds>(elapsed_secs);
  return secs.count();
}

/// Returns current time
inline std::chrono::time_point<std::chrono::system_clock> now()
{
  return std::chrono::system_clock::now();
}

  
#endif /* Helpers_h */
