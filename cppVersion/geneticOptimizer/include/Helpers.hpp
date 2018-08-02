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

inline float RandFloat(float min, float max)
{
  return min + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(max-min)));
}

inline int RandInt(int min, int max)
{
  return min + (rand() % static_cast<int>(max - min + 1));
}

inline bool RandBool()
{
  return rand() % 2;
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
  std::ifstream is_file(configPath);
  std::ofstream outputFile;
  std::string tmpName = "tmp_"+std::to_string(RandInt(0, 10000000))+".md";
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
  system(("mv "+tmpName+" "+configPath).c_str());
}


inline double GetParamFomeConfig(std::string configPath, std::string keyToFind){
  std::ifstream is_file(configPath);
  
  std::string line;
  while(getline(is_file, line)){
    std::istringstream is_line(line);
    std::string key;
    if( std::getline(is_line, key, ':')){
      std::string value;
      if(std::getline(is_line, value)){
        if(key==keyToFind){
          return stod(value);
        }
      }
    }
  }
  return -999999999;
}

struct ClusteringOutput {
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

  
#endif /* Helpers_h */
