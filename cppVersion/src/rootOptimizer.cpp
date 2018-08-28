#include "GeneticHelpers.hpp"
#include "ConfigurationManager.hpp"

#include <TFitter.h>

using namespace std;

const string configPath = "bestGenetic.md";
const string outputPath = "autoGenOutput.txt";

int nEventsPerTest = 100;   ///< On how many events each population member will be tested

double severityFactor = 10.0;

bool dependSensor = true;
bool reachedEE = true;

int minNtuple = 1;
int maxNtuple = 1;

string dataPath = "../../data/MultiParticleInConeGunProducer_PDGid22_nPart1_Pt6p57_Eta2p2_InConeDR0p10_PDGid22_predragm_cmssw1020pre1_20180730/NTUP/partGun_PDGid22_x96_Pt6.57To6.57_NTUP_ ";

string histsOutputPath = "../clusteringResultsCXX/rootOptimizer/";

void myfuncf(Int_t&, Double_t*, Double_t &f, Double_t *par, Int_t)
{
  double chi2 = 0;
  
  string kernel;
  if(round(par[kKernel]) == 0)       kernel = "step";
  else if(round(par[kKernel]) == 1)  kernel = "gaus";
  else                               kernel = "exp";
  
  string command = "./createQualityPlots "
  +to_string(dependSensor)+" "
  +dataPath+" "
  +histsOutputPath+" "
  +to_string(par[kDeltacEE])+" "
  +to_string(par[kDeltacFH])+" "
  +to_string(par[kDeltacBH])+" "
  +to_string(par[kEnergyThreshold])+" "
  +to_string(par[kCriticalDistanceEE])+" "
  +to_string(par[kCriticalDistanceFH])+" "
  +to_string(par[kCriticalDistanceBH])+" "
  +to_string(par[kKappa])+" "
  +"0 " // verbosity
  +to_string(minNtuple)+" " // min n tuple
  +to_string(maxNtuple)+" " // max n tuple
  +"0 " // min layer
  +"52 " // max layer
  +to_string(nEventsPerTest)+" "
  +kernel+" "
  +to_string(reachedEE)+" "
  +to_string(par[kMatchingDistance])+" "
  +outputPath+" "
  +" > /dev/null 2>&1";
  
  cout<<"Running clusterization:"<<endl;
  cout<<command<<endl;
  system(command.c_str());
  cout<<"Clusterization output:"<<endl;
  
  ClusteringOutput output = ReadOutput(outputPath);
  output.Print();
  
  chi2 =    fabs(output.containmentMean-1)
          +      output.containmentSigma
          + fabs(output.resolutionMean)
          +      output.resolutionSigma
          +      output.separationMean
          +      output.separationSigma
          + fabs(output.deltaNclustersMean)
          +      output.deltaNclustersSigma
          +      output.nEmptyMatched
          +      output.nNoMached
          +      output.nZeroSize;
  
  cout<<"\n\nchi2:"<<chi2<<"\n\n"<<endl;
  cout<<"Score (GA equivalent):"<<severityFactor/chi2<<endl;
  
  f = chi2;
  return;
}


int main()
{
  TVirtualFitter::SetDefaultFitter("Minuit");
  const int nPar = 13;
  TFitter *fitter = new TFitter(nPar);
  fitter->SetFCN(myfuncf);
  
  for(int iPar=0;iPar<kNparams;iPar++){
    cout<<"Setting parameter:"<<endl;
    cout<<iPar<<"\t"<<paramTitle[iPar]<<"\t"<<paramStart[iPar]<<"\t"<<0.1<<"\t"<<paramMin[iPar]<<"\t"<<paramMax[iPar]<<endl;
    fitter->SetParameter(iPar, paramTitle[iPar], paramStart[iPar], 0.1, paramMin[iPar], paramMax[iPar]);
  }
  
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
