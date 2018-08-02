#include "include/Helpers.hpp"

#include "Chromosome.hpp"

#include <TMath.h>

using namespace std;

const string configPath = "../configs/autoGenConfig.md";
const string outputPath = "autoGenOutput.txt";


void myfuncf(int&, double*, double &f, double *par, int)
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
