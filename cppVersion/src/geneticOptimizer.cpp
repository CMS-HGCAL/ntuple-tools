#include "GeneticHelpers.hpp"

#include "Chromosome.hpp"

#include <TMath.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TROOT.h>
#include <TFile.h>

#include <future>
#include <thread>
#include <assert.h>
#include <chrono>
#include <utility>
#include <atomic>
#include <mutex>
#include <cstring>
#include <pthread.h>
#include <random>
#include <unistd.h>
#include <signal.h>

#include <stdio.h>
#include <dirent.h>

using namespace std;

string baseResultsPath;

string baseResultsSearchPath = "geneticResults/qcd/";
string baseResultsDirName = "results_";

int populationSize = 10;  ///< Size of the population, will stay the same for all generations. Make it an even number, otherwise there may be some complications.
int maxBatchSize = 20;  ///< execute this number of jobs simultaneously
int nGenerations = 100;     ///< Number of iterations
int nEventsPerTest = 5;   ///< On how many events per ntuple each population member will be tested

int processTimeout = 200; ///< this is a timeout for the test of whole population in given generation, give it at least 2-3 seconds per member per event (processTimeout ~ 2*maxBatchSize*nEventsPerTest)

double mutationChance = 0.003;
double severityFactor = 10.0; // larger the value, more easily population members will die (and the more good solutions will be promoted)

bool dependSensor = true;
bool reachedEE = false;

Chromosome::ECrossover crossoverStrategy = Chromosome::kFixedSinglePoint;

int minNtuple = 1;
int maxNtuple = 1;

int minLayer = 1;
int maxLayer = 53;

bool fixEnergyThreshold = true;
double energyThreshold = 3.0;
bool fixKernelFunction = true;
int kernelFunction = 0; // 0 - step, 1 - gauss, 2 - exp
bool fixMatchingDistance = true;
double matchingDistance = 10.0;

// two-photons
//string dataPath = "../../data/MultiParticleInConeGunProducer_PDGid22_nPart1_Pt6p57_Eta2p2_InConeDR0p10_PDGid22_predragm_cmssw1020pre1_20180730/NTUP/partGun_PDGid22_x96_Pt6.57To6.57_NTUP_";

// two-pions
//string dataPath = "../../data/MultiParticleInConeGunProducer_SinglePion_Pt80_Eta2_InConePion_DeltaR0p1_clange_20171102/NTUP/partGun_PDGid211_x120_Pt80.0To80.0_NTUP_";

// QCD event
string dataPath = "../../data/eventQCD_wf24031p0_Pt_80_120_14TeV_2023D28_noPU_jniedzie_20190208/NTUP/eventQCD_x1320_1.0To35.0_NTUP_";

string outputPath = "../clusteringResultsCXX/geneticOptimizerQCD/";

mt19937 randGenerator;

TGraph *scoresMean;
TFile *outfile;

TH2D *paramHists[kNparams];
TH2D *paramHistsWgt[kNparams];

TH1D *paramHistsFailed[kNparams];
TH1D *paramHistsPassed[kNparams];

atomic<bool> allKidsFinished;

int scheduleClustering(shared_ptr<Chromosome> chromo)
{
  int kernelIndex = round(chromo->GetParam(kKernel));
  string kernel;
  if(kernelIndex == 0) kernel = "step";
  if(kernelIndex == 1) kernel = "gaus";
  if(kernelIndex == 2) kernel = "exp";
  
  string command = "./createQualityPlots "
  +to_string(dependSensor)+" "
  +dataPath+" "
  +outputPath+" "
  +to_string(chromo->GetParam(kDeltacEE))+" "
  +to_string(chromo->GetParam(kDeltacFH))+" "
  +to_string(chromo->GetParam(kDeltacBH))+" "
  +to_string(chromo->GetParam(kEnergyThreshold))+" "
  +to_string(chromo->GetParam(kCriticalDistanceEE))+" "
  +to_string(chromo->GetParam(kCriticalDistanceFH))+" "
  +to_string(chromo->GetParam(kCriticalDistanceBH))+" "
  +to_string(chromo->GetParam(kKappa))+" "
  +"0 " // verbosity
  +to_string(minNtuple)+" "
  +to_string(maxNtuple)+" "
  +to_string(minLayer)+" "
  +to_string(maxLayer)+" "
  +to_string(nEventsPerTest)+" "
  +kernel+" "
  +to_string(reachedEE)+" "
  +to_string(chromo->GetParam(kMatchingDistance))+" "
  +chromo->GetClusteringOutputPath()+" "
  +" > /dev/null 2>&1";
  
  cout<<".";
  int pid = ::fork();
  
  if(pid==0){
    system(command.c_str());
    exit(0);
  }
  else return pid;
}

void waitGently(vector<int> childPid)
{
  for(int i=0;i<childPid.size();i++){
    int status;
    waitpid(childPid[i],&status,0);
//    cout<<"Child with pid "<<childPid[i]<<" finished"<<endl;
  }
  allKidsFinished = true;
}

void killChildrenAfterTimeout(vector<int> childPid, int timeout)
{
  int timeElapsed=0;
  while(timeElapsed < timeout && !allKidsFinished){
    sleep(1);
    timeElapsed++;
    cout<<"Time elapsed:"<<timeElapsed<<" (killing after "<<processTimeout<<" s.)"<<endl;
  }
  if(!allKidsFinished){
    cout<<"\nKilling all children\n\n"<<endl;
    for(int pid : childPid){
      kill(pid,SIGKILL);
    }
  }
  else{
    cout<<"\nIt seems that all child processes already finished\n\n"<<endl;
  }
}

int GetWeightedRandom(discrete_distribution<double> dist)
{
  return dist(randGenerator);
}

void TestPopulation(vector<shared_ptr<Chromosome>> population,
                    TH1D *hist, discrete_distribution<double> &dist, int generation)
{
  std::vector<int> childPid;
  
  for(int iBatch=0;iBatch<ceil(populationSize/(double)maxBatchSize);iBatch++){
    cout<<"Batch "<<iBatch<<"/"<<ceil(populationSize/(double)maxBatchSize)-1<<endl;
    
    for(int i=0;i<maxBatchSize;i++){
      if(iBatch*maxBatchSize+i >= populationSize) break;
      int pid = scheduleClustering(population[iBatch*maxBatchSize+i]);
      if(pid > 0) childPid.push_back(pid);
    }

    std::cout<<"\n\nall forks created\n\n"<<std::endl;
  
    allKidsFinished = false;
  
    thread *waitThread = new thread(waitGently,childPid);
    thread *killThread = new thread(killChildrenAfterTimeout,childPid,processTimeout);
  
    waitThread->join();
    killThread->join();
  }
    
  vector<double> scores;
  double minScore=99999, maxScore=-99999;
  int bestGuyIndex = -1;
  
  for(int i=0;i<populationSize;i++){
    population[i]->CalculateScore();
    double score = population[i]->GetScore();
    if(score < minScore && score > 1E-5) minScore = score; // make sure to remove veeery bad results
    if(score > maxScore){
      maxScore = score;
      bestGuyIndex = i;
    }
    if(score > 1E-5) hist->Fill(score);
    scores.push_back(score);
  }
  
  // re-assing normalized points to population members
  vector<double> scoresNormalized;
  double normScore;

  for(int i=0;i<populationSize;i++){
    normScore = (scores[i]-minScore)/(maxScore-minScore);
    if(normScore < 1E-5) normScore = 0;
    
    cout<<"ID: "<<population[i]->GetUniqueID()<<"\tscore:\t"<<scores[i]<<"\tnormalized:\t"<<normScore<<endl;
    
    population[i]->SetNormalizedScore(normScore);
    scoresNormalized.push_back(normScore);
  }

  cout<<"\n\n================================================="<<endl;
  cout<<"The best guy in this generation was..."<<endl;
  population[bestGuyIndex]->Print();
  population[bestGuyIndex]->StoreInConfig(baseResultsPath+"bestGenetic"+to_string(generation)+".md");
  
  dist = discrete_distribution<double>(scoresNormalized.begin(), scoresNormalized.end());
}

void SaveHists()
{
  outfile = new TFile((baseResultsPath+"geneticHists.root").c_str(),"update");
  outfile->cd();
  
  scoresMean->Write("scores",TObject::kOverwrite);
  
  for(int i=0;i<kNparams;i++){
    paramHists[i]->Write(paramTitle[i],TObject::kOverwrite);
    paramHistsWgt[i]->Write(Form("%swgt",paramTitle[i]),TObject::kOverwrite);
    paramHistsFailed[i]->Write(Form("%s_failed",paramTitle[i]),TObject::kOverwrite);
    paramHistsPassed[i]->Write(Form("%s_passed",paramTitle[i]),TObject::kOverwrite);
  }
  
  outfile->Close();
  delete outfile;
}

void SaveConfigurationToFile(){
  ofstream outputFile;
  outputFile.open(baseResultsPath+"geneticConfig.txt");
  
  outputFile<<"Configuration of the genetic optimizer:"<<endl;
  outputFile<<"Number of creatures in each population:\t"<<populationSize<<endl;
  outputFile<<"Number of generations:\t"<<nGenerations<<endl;
  outputFile<<"Test N events per tuples:\t"<<nEventsPerTest<<endl;
  outputFile<<"N tuple min:\t"<<minNtuple<<endl;
  outputFile<<"N tuple max:\t"<<maxNtuple<<endl;
  outputFile<<"Mutation probability:\t"<<mutationChance<<endl;
  outputFile<<"Severity_factor:\t"<<severityFactor<<endl;
  outputFile<<"Crossover_strategy:\t"<<Chromosome::crossoverName[crossoverStrategy]<<endl;
  outputFile<<"Timeout for each generation:\t"<<processTimeout<<" (s)"<<endl;
  outputFile<<"Sensor dependance:\t"<<dependSensor<<endl;
  outputFile<<"Reached EE only:\t"<<reachedEE<<endl;
  outputFile<<"\nInput data path:\n"<<dataPath<<"\n"<<endl;
  
  outputFile.close();
}

void SetBaseResultsPath()
{
  DIR *dir = opendir(baseResultsSearchPath.c_str());
  
  if(!dir){
    cout<<"could not open directory:"<<baseResultsSearchPath<<endl;
    cout<<"trying to create it..."<<endl;
    system(("mkdir -p "+baseResultsSearchPath).c_str());
    
    dir = opendir(baseResultsSearchPath.c_str());
    
    if(!dir){
      cout<<"Failed! Please create this direcoty manually and try again."<<endl;
      exit(0);
    }
    else{
      cout<<"Successfully created output directory"<<endl;
    }
  }
  
  struct dirent *ent;
  int maxIndex = -1;
  
  while( (ent = readdir(dir)) ){
    string fileName = ent->d_name;
    
    size_t pos = fileName.find(baseResultsDirName);
    if(pos != std::string::npos){
      string indexString = fileName.substr(fileName.find(baseResultsDirName) + 8);
      int index = stoi(indexString);
      if(index > maxIndex) maxIndex=index;
    }
  }
  closedir(dir);
  
  baseResultsPath = baseResultsSearchPath+baseResultsDirName+to_string(maxIndex+1)+"/";
}

int main(int argc, char* argv[])
{
  if(argc > 1 && argc!=9){
    cout<<"Usage:"<<endl;
    cout<<"./geneticOptimizer"<<endl;
    cout<<"or with custom parameters:"<<endl;
    cout<<"./geneticOptimizer populationSize nGenerations nEventsPerTest processTimeout mutationChance severityFactor dataPath outputPath"<<endl;
    exit(0);
  }
  if(argc == 9){
    populationSize = atoi(argv[1]);
    nGenerations = atoi(argv[2]);
    nEventsPerTest = atoi(argv[3]);
    processTimeout = atoi(argv[4]);
    mutationChance = atof(argv[5]);
    severityFactor = atof(argv[6]);
    dataPath = argv[7];
    outputPath = argv[8];
  }
  
  gROOT->ProcessLine(".L loader.C+");
  TApplication theApp("App", &argc, argv);

  SetBaseResultsPath();
  system(("mkdir -p "+baseResultsPath).c_str());
  SaveConfigurationToFile();
  
  TH1D *scoresDist[nGenerations];
  scoresMean = new TGraph();
  
  for(int i=0;i<kNparams;i++){
    paramHists[i] = new TH2D(paramTitle[i],paramTitle[i],nGenerations,0,nGenerations,100,paramMin[i],paramMax[i]);
    paramHistsWgt[i] = new TH2D(Form("%sWgt",paramTitle[i]),Form("%sWgt",paramTitle[i]),nGenerations,0,nGenerations,100,paramMin[i],paramMax[i]);
    
    paramHistsFailed[i] = new TH1D(Form("%s_failed",paramTitle[i]),
                                   Form("%s_failed",paramTitle[i]),100,paramMin[i],paramMax[i]);
    
    paramHistsPassed[i] = new TH1D(Form("%s_passed",paramTitle[i]),
                                   Form("%s_passed",paramTitle[i]),100,paramMin[i],paramMax[i]);
    
  }
  
  // make sure to recreate output file
  outfile = new TFile((baseResultsPath+"geneticHists.root").c_str(),"recreate");
  outfile->Close();
  delete outfile;
  
  vector<shared_ptr<Chromosome>> population;
  discrete_distribution<double> scores;
  
  for(int iGeneration=0;iGeneration<nGenerations;iGeneration++){
    cout<<"\n\nGeneration "<<iGeneration<<"\n\n"<<endl;
    
    if(iGeneration==0){
      // draw initial population
      for(int i=0;i<populationSize;i++){
        auto chromo = Chromosome::GetRandom();
        if(fixEnergyThreshold)  chromo->FixParam(kEnergyThreshold, energyThreshold);
        if(fixKernelFunction)   chromo->FixParam(kKernel, kernelFunction);
        if(fixMatchingDistance) chromo->FixParam(kMatchingDistance, matchingDistance);
        
        chromo->SaveToBitChromosome();
        chromo->SetMutationChance(mutationChance);
        chromo->SetSeverityFactor(severityFactor);
        chromo->SetCrossover(crossoverStrategy);
        chromo->SetInputDataPath(dataPath);
        chromo->SetMinLayer(minLayer);
        chromo->SetMaxLayer(maxLayer);
        population.push_back(chromo);
      }
    }
    else{
      // draw new population
      vector<shared_ptr<Chromosome>> oldPopulation;
      for(auto &chromo : population){
        oldPopulation.push_back(chromo);
      }
      population.clear();
      
      for(int i=0;i<populationSize/2;i++){
        shared_ptr<Chromosome> mom, dad;
        
        int iter=0;
        do {
          if(iter>10) cout<<"mom==dad for "<<iter<<" times. Probably something went wrong..."<<endl;
          if(iter>20) exit(0);
          iter++;
          
          mom = oldPopulation[(int)scores(randGenerator)];
          dad = oldPopulation[(int)scores(randGenerator)];
        } while (mom==dad);

        auto children = dad->ProduceChildWith(mom);
        population.push_back(children[0]);
        population.push_back(children[1]);
      }
    }
    
    scoresDist[iGeneration] = new TH1D(Form("score dist gen[%i]",iGeneration),
                                       Form("score dist gen[%i]",iGeneration),
                                       200,-20,20);
    
    TestPopulation(population, scoresDist[iGeneration], scores, iGeneration);
    
    for(int i=0;i<populationSize;i++){
      for(int iPar=0;iPar<kNparams;iPar++){
        double val = population[i]->GetParam((EParam)iPar);
        double score = population[i]->GetScore();
        double scoreNormalized = population[i]->GetNormalizedScore();
        
        if(iPar == kKernel) val = round(val);
        
        paramHists[iPar]->Fill(iGeneration,val);
        paramHistsWgt[iPar]->Fill(iGeneration,val,score);
        
        if(scoreNormalized > 0.1) paramHistsPassed[iPar]->Fill(val);
        else                      paramHistsFailed[iPar]->Fill(val);
      }
    }
    scoresMean->SetPoint(iGeneration, iGeneration, scoresDist[iGeneration]->GetMean());
    
    // Save updated histograms after each iteration
    SaveHists();
    
//    if(iGeneration == 0 || iGeneration == nGenerations-1){
//      cout<<"\n\nPopulation:\n\n"<<endl;
//      for(int i=0;i<populationSize;i++){
//        population[i]->Print();
//      }
//    }
  }
  
  theApp.Run();
  return 0;
}
