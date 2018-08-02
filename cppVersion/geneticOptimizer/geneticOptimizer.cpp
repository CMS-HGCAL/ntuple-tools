#include "include/Helpers.hpp"

using namespace std;

const string configPath = "../configs/autoGenConfig.md";
const string outputPath = "autoGenOutput.txt";

struct Chromosome {
  // 52 bits (1st choromosome)
  uint16_t criticalDistanceEE;
  uint16_t criticalDistanceFH;
  uint16_t criticalDistanceBH;
  bool dependSensor;          // 2 bits
  bool reachedEE;             // 2 bits

  // 56 bits (2nd chromosome)
  uint16_t kernel;
  uint16_t deltacEE;
  uint16_t deltacFH;
  uint16_t deltacBH;
  
  // 56 bits (3d chromosome)
  uint16_t kappa;
  uint16_t energyMin;
  uint16_t matchingDistance;
  uint16_t minClusters;
  
  uint64_t bitChromosome[3]; // 64 bit
  
  template<class T>
  void ShiftIntoChromosome(T value, int &shift, int chromoIndex){
    uint64_t mask = (uint64_t)value << shift;
    bitChromosome[chromoIndex] |= mask;
    shift += BitSize(value);
  }
  
  template<class T>
  void SetValueFromChromosome(T &value, int &shift, int chromoIndex){
    uint64_t mask = 0;
    for(int i=0;i<BitSize(value);i++){mask |= 1ull << (i+shift);}
    value = (bitChromosome[chromoIndex] & mask) >> shift;
    shift += BitSize(value);
  }
  
  inline void SetupBitChromosome(){
    int currentShift = 0;
    bitChromosome[0] = 0;
    bitChromosome[1] = 0;
    bitChromosome[2] = 0;
    
    // load values to 1st chromosome
    ShiftIntoChromosome(criticalDistanceEE,currentShift,0);
    ShiftIntoChromosome(criticalDistanceFH,currentShift,0);
    ShiftIntoChromosome(criticalDistanceBH,currentShift,0);
    ShiftIntoChromosome(dependSensor,currentShift,0);
    ShiftIntoChromosome(reachedEE,currentShift,0);
    
    // load values to 2st chromosome
    currentShift = 0;
    ShiftIntoChromosome(kernel,currentShift,1);
    ShiftIntoChromosome(deltacEE,currentShift,1);
    ShiftIntoChromosome(deltacFH,currentShift,1);
    ShiftIntoChromosome(deltacBH,currentShift,1);
    
    // load values to 3st chromosome
    currentShift = 0;
    ShiftIntoChromosome(kappa,currentShift,2);
    ShiftIntoChromosome(energyMin,currentShift,2);
    ShiftIntoChromosome(matchingDistance,currentShift,2);
    ShiftIntoChromosome(minClusters,currentShift,2);

  }
  
  inline void GetFromBitChromosome(){
    int currentShift = 0;
    
    // read 1st chromosome
    SetValueFromChromosome(criticalDistanceEE, currentShift, 0);
    SetValueFromChromosome(criticalDistanceFH, currentShift, 0);
    SetValueFromChromosome(criticalDistanceBH, currentShift, 0);
    SetValueFromChromosome(dependSensor, currentShift, 0);
    SetValueFromChromosome(reachedEE, currentShift, 0);
    
    // read 1st chromosome
    currentShift = 0;
    SetValueFromChromosome(kernel, currentShift, 1);
    SetValueFromChromosome(deltacEE, currentShift, 1);
    SetValueFromChromosome(deltacFH, currentShift, 1);
    SetValueFromChromosome(deltacBH, currentShift, 1);
    
    // read 2st chromosome
    currentShift = 0;
    SetValueFromChromosome(kappa, currentShift, 2);
    SetValueFromChromosome(energyMin, currentShift, 2);
    SetValueFromChromosome(matchingDistance, currentShift, 2);
    SetValueFromChromosome(minClusters, currentShift, 2);
  }
  
  void Print(){
    cout<<"Critical distance:"<<endl;
    cout<<"\tEE:"<<criticalDistanceEE/1000.<<endl;
    cout<<"\tFH:"<<criticalDistanceFH/1000.<<endl;
    cout<<"\tBH:"<<criticalDistanceBH/1000.<<endl;
    
    cout<<"Depend on sensor:"<<dependSensor<<endl;
    cout<<"Reached EE only:"<<reachedEE<<endl;
    
    cout<<"Kernel: ";
    if(kernel==0) cout<<"step"<<endl;
    if(kernel==1) cout<<"gaus"<<endl;
    if(kernel==2) cout<<"exp"<<endl;
    
    cout<<"Critical #delta:"<<endl;
    cout<<"\tEE:"<<deltacEE/1000.<<endl;
    cout<<"\tFH:"<<deltacFH/1000.<<endl;
    cout<<"\tBH:"<<deltacBH/1000.<<endl;
    
    cout<<"kappa:"<<kappa/1000.<<endl;
    cout<<"energy threshold:"<<energyMin/100000.<<endl;
    cout<<"max matching distance:"<<matchingDistance/1000.<<endl;
    cout<<"min clusters:"<<minClusters<<endl;
  }
};

void myfuncf(Int_t&, Double_t*, Double_t &f, Double_t *par, Int_t)
{
  double chi2 = 0;
  
  UpdateParamValue(configPath, "depend_sensor",par[0] < 0.5 ? 0 : 1);
  if(par[1]<0.5)      UpdateParamValue(configPath, "energy_density_function","step");
  else if(par[1]<1.5) UpdateParamValue(configPath, "energy_density_function","gaus");
  else                UpdateParamValue(configPath, "energy_density_function","exp");
  
  UpdateParamValue(configPath, "critial_distance_EE",par[2]);
  UpdateParamValue(configPath, "critial_distance_FH",par[3]);
  UpdateParamValue(configPath, "critial_distance_BH",par[4]);
  UpdateParamValue(configPath, "deltac_EE",par[5]);
  UpdateParamValue(configPath, "deltac_FH",par[6]);
  UpdateParamValue(configPath, "deltac_BH",par[7]);
  UpdateParamValue(configPath, "kappa",par[8]);
  UpdateParamValue(configPath, "energy_min",par[9]);
  UpdateParamValue(configPath, "min_clusters",par[10]);
  UpdateParamValue(configPath, "reachedEE_only",par[11] < 0.5 ? 0 : 1);
  UpdateParamValue(configPath, "matching_max_distance",par[12]);
  
  cout<<"Running clusterization"<<endl;
  system(("./createQualityPlots "+configPath+" > /dev/null 2>&1").c_str());
//  system(("./createQualityPlots "+configPath).c_str());
  cout<<"Clusterization output:"<<endl;
  
  ClusteringOutput output = ReadOutput(outputPath);
  output.Print();
  
  chi2 = fabs(output.resolutionMean) + output.separationMean + 1/output.containmentMean;
  cout<<"\n\nchi2:"<<chi2<<"\n\n"<<endl;
  
  f = chi2;
  return;
}


int main()
{
//  UpdateParamValue(configPath, "analyze_events_per_tuple",10);
//  GetParamFomeConfig(configPath, "depend_sensor"), 1, 0, 1); // 0 - no dependance, 1 - depend on sensor
  
  Chromosome ch1;
  ch1.dependSensor = 1;
  ch1.reachedEE = 0;
  ch1.criticalDistanceEE = 2.53*1000;
  ch1.criticalDistanceFH = 4.53*1000;
  ch1.criticalDistanceBH = 6.53*1000;
  ch1.kernel = 2;
  ch1.deltacEE = 123.412*1000;
  ch1.deltacFH = 3.524621*1000;
  ch1.deltacBH = 4123.431*1000;
  ch1.energyMin = 0.00012*100000;
  ch1.minClusters = 3;
  ch1.matchingDistance = 4.32*1000;
  ch1.SetupBitChromosome();
  
  Chromosome ch2;
  ch2.bitChromosome[0] = ch1.bitChromosome[0];
  ch2.bitChromosome[1] = ch1.bitChromosome[1];
  ch2.bitChromosome[2] = ch1.bitChromosome[2];
  ch2.GetFromBitChromosome();
  
  
  cout<<"Chromosome 1:"<<endl;
  ch1.Print();
  
  cout<<"\n\nChromosome 2:"<<endl;
  ch2.Print();
  
  
  return 0;
}
