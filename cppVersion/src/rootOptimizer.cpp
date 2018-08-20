#include "GeneticHelpers.hpp"
#include "ConfigurationManager.hpp"

#include <TFitter.h>

using namespace std;

const string configPath = "bestGenetic.md";
const string outputPath = "autoGenOutput.txt";

int nEventsPerTest = 30;   ///< On how many events each population member will be tested

int minNtuple = 1;
int maxNtuple = 1;

string dataPath = "../../data/MultiParticleInConeGunProducer_PDGid22_nPart1_Pt6p57_Eta2p2_InConeDR0p10_PDGid22_predragm_cmssw1020pre1_20180730/NTUP/partGun_PDGid22_x96_Pt6.57To6.57_NTUP_ ";

string histsOutputPath = "../clusteringResultsCXX/rootOptimizer/";

void myfuncf(Int_t&, Double_t*, Double_t &f, Double_t *par, Int_t)
{
  double chi2 = 0;
  
  string kernel;
  if(round(par[1]) == 0)      kernel = "step";
  else if(round(par[1]) == 1) kernel = "gaus";
  else                        kernel = "exp";
  
  string command = "./createQualityPlots "
  +to_string(par[0] < 0.5 ? 0 : 1)+" "
  +dataPath+" "
  +histsOutputPath+" "
  +to_string(par[5])+" "
  +to_string(par[6])+" "
  +to_string(par[7])+" "
  +to_string(par[9])+" "
  +to_string(par[10])+" "
  +to_string(par[2])+" "
  +to_string(par[3])+" "
  +to_string(par[4])+" "
  +to_string(par[8])+" "
  +"0 " // verbosity
  +to_string(minNtuple)+" " // min n tuple
  +to_string(maxNtuple)+" " // max n tuple
  +"0 " // min layer
  +"52 " // max layer
  +to_string(nEventsPerTest)+" "
  +kernel+" "
  +to_string(par[11] < 0.5 ? 0 : 1)+" "
  +to_string(par[12])+" "
  +outputPath+" "
  +" > /dev/null 2>&1";
  
  cout<<"Running clusterization"<<endl;
  system(command.c_str());
  cout<<"Clusterization output:"<<endl;
  
  ClusteringOutput output = ReadOutput(outputPath);
  output.Print();
  
  chi2 =   fabs(output.containmentMean-1)
          +      output.containmentSigma
          + fabs(output.resolutionMean)
          +      output.resolutionSigma
          +      output.separationMean
          +      output.separationSigma;
  cout<<"\n\nchi2:"<<chi2<<"\n\n"<<endl;
  cout<<"Score (GA equivalent):"<<10.0/chi2<<endl;
  
  f = chi2;
  return;
}


int main()
{
  TVirtualFitter::SetDefaultFitter("Minuit");
  const int nPar = 13;
  TFitter *fitter = new TFitter(nPar);
  fitter->SetFCN(myfuncf);
  
  // set fitter params
  fitter->SetParameter(0, "depend_sensor",          true, 1.0, 0, 1); // 0 - no dependance, 1 - depend on sensor
  fitter->FixParameter(0);
  fitter->SetParameter(1, "energy_density_function",kernelStart, 1.0, 0, 2); // 0 - step, 1 - gaus, 2 - exp
  fitter->SetParameter(2, "critial_distance_EE",    criticalDistanceEEstart,0.1,
                                                    criticalDistanceEEmin, criticalDistanceEEmax);
  fitter->SetParameter(3, "critial_distance_FH",    criticalDistanceFHstart, 0.1,
                                                    criticalDistanceFHmin, criticalDistanceFHmax);
  fitter->SetParameter(4, "critial_distance_BH",    criticalDistanceBHstart, 0.1,
                                                    criticalDistanceBHmin, criticalDistanceBHmax);
  fitter->SetParameter(5, "deltac_EE",              deltacEEstart, 0.1, deltacEEmin, deltacEEmax);
  fitter->SetParameter(6, "deltac_FH",              deltacFHstart, 0.1, deltacFHmin, deltacFHmax);
  fitter->SetParameter(7, "deltac_BH",	            deltacBHstart, 0.1, deltacBHmin, deltacBHmax);
  
  fitter->SetParameter(8, "kappa",                  kappaStart, 0.1, kappaMin, kappaMax);
  fitter->SetParameter(9, "energy_min",             energyThresholdStart, 0.01, energyThresholdMin, energyThresholdMax);
  fitter->SetParameter(10,"min_clusters",           minClustersStart, 1.0, minClustersMin, minClustersMax);
  fitter->SetParameter(11,"reachedEE_only",         true, 1, 0, 1);
  fitter->SetParameter(12,"matching_max_distance",  matchingDistanceStart, 0.1,
                                                    matchingDistanceMin, matchingDistanceMax);
  
  double args = 0; // put to 0 for results only, or to -1 for no garbage
  fitter->ExecuteCommand( "SET PRINTOUT"  , &args, 1);
  //  fitter->ExecuteCommand( "SET NOWARNINGS", &args, 0);
  fitter->ExecuteCommand( "SET PRINT"     , &args, 1);
  //  double fitterError[1] = {5.0};
  //  fitter->ExecuteCommand( "SET ERR", fitterError, 1);
  double strategyLevel[1] = {2};
  fitter->ExecuteCommand( "SET STR", strategyLevel, 1);
  double arglist[1] = {0};
  fitter->ExecuteCommand("MIGRAD", arglist, 0);
  
  return 0;
}
