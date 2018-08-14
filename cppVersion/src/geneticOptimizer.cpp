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

using namespace std;

int populationSize = 100;  ///< Size of the population, will stay the same for all generations
int nGenerations = 100;     ///< Number of iterations
int nEventsPerTest = 10;   ///< On how many events each population member will be tested

int processTimeout = 300; ///< this is a timeout for the test of whole population in given generation, give it at least 2-3 seconds per member per event ( processTimeout ~ 2*populationSize

double mutationChance = 0.002;
double severityFactor = 1.0;

string dataPath = "../../data/MultiParticleInConeGunProducer_PDGid22_nPart1_Pt6p57_Eta2p2_InConeDR0p10_PDGid22_predragm_cmssw1020pre1_20180730/NTUP/partGun_PDGid22_x96_Pt6.57To6.57_NTUP_ ";

string outputPath = "../clusteringResultsCXX/geneticOptimizer/";

mt19937 randGenerator;

TGraph *scoresMean;
TFile *outfile;
TH2D *criticalDistanceEE, *criticalDistanceFH, *criticalDistanceBH, *dependSensor, *reachedEE, *kernel, *deltacEE, *deltacFH, *deltacBH, *kappa, *energyMin, *matchingDistance, *minClusters;

TH2D *criticalDistanceEEwgt, *criticalDistanceFHwgt, *criticalDistanceBHwgt, *dependSensorwgt, *reachedEEwgt, *kernelwgt, *deltacEEwgt, *deltacFHwgt, *deltacBHwgt, *kappawgt, *energyMinwgt, *matchingDistancewgt, *minClusterswgt;


atomic<bool> allKidsFinished;

int scheduleClustering(Chromosome *chromo)
{
  string kernel;
  if(chromo->GetKernel() == 0) kernel = "step";
  if(chromo->GetKernel() == 1) kernel = "gaus";
  if(chromo->GetKernel() == 2) kernel = "exp";
  
  string command = "./createQualityPlots "
  +to_string(chromo->GetDependSensor())+" "
  +dataPath+" "
  +outputPath+" "
  +to_string(chromo->GetDeltacEE())+" "
  +to_string(chromo->GetDeltacFH())+" "
  +to_string(chromo->GetDeltacBH())+" "
  +to_string(chromo->GetEnergyMin())+" "
  +to_string(chromo->GetMinClusters())+" "
  +to_string(chromo->GetCriticalDistanceEE())+" "
  +to_string(chromo->GetCriticalDistanceFH())+" "
  +to_string(chromo->GetCriticalDistanceBH())+" "
  +to_string(chromo->GetKappa())+" "
  +"0 " // verbosity
  +"1 " // min n tuple
  +"1 " // max n tuple
  +"0 " // min layer
  +"52 " // max layer
  +to_string(nEventsPerTest)+" "
  +kernel+" "
  +to_string(chromo->GetReachedEE())+" "
  +to_string(chromo->GetMatchingDistance())+" "
  +chromo->GetClusteringOutputPath()+" "
  +" > /dev/null 2>&1";
  
  cout<<"Forking"<<endl;
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
    cout<<"Child with pid "<<childPid[i]<<" finished"<<endl;
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

void TestPopulation(vector<Chromosome*> population, TH1D *hist, discrete_distribution<double> &dist)
{
   std::vector<int> childPid;
  
  for(int i=0;i<population.size();i++){
    int pid = scheduleClustering(population[i]);
    if(pid > 0) childPid.push_back(pid);
  }

  std::cout<<"\n\nall forks created\n\n"<<std::endl;
  
  allKidsFinished = false;
  
  thread *waitThread = new thread(waitGently,childPid);
  thread *killThread = new thread(killChildrenAfterTimeout,childPid,processTimeout);
  
  waitThread->join();
  killThread->join();
  
  vector<double> scores;
  double minScore=99999, maxScore=-99999;
  
  for(int i=0;i<populationSize;i++){
    population[i]->CalculateScore();
    double score = population[i]->GetScore();
    if(score < minScore) minScore = score; // make sure to remove veeery bad results
    if(score > maxScore) maxScore = score;
    hist->Fill(score);
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

  dist = discrete_distribution<double>(scoresNormalized.begin(), scoresNormalized.end());
}

void SaveHists()
{
  outfile = new TFile("geneticHists.root","update");
  outfile->cd();
  
  scoresMean->Write("scores",TObject::kOverwrite);
  criticalDistanceEE->Write("critDistEE",TObject::kOverwrite);
  criticalDistanceFH->Write("critDistFH",TObject::kOverwrite);
  criticalDistanceBH->Write("critDistBH",TObject::kOverwrite);
  dependSensor->Write("sensor",TObject::kOverwrite);
  reachedEE->Write("reachedEE",TObject::kOverwrite);
  kernel->Write("kernel",TObject::kOverwrite);
  deltacEE->Write("deltaEE",TObject::kOverwrite);
  deltacFH->Write("deltaFH",TObject::kOverwrite);
  deltacBH->Write("deltaBH",TObject::kOverwrite);
  kappa->Write("kappa",TObject::kOverwrite);
  energyMin->Write("eMin",TObject::kOverwrite);
  matchingDistance->Write("matchingDist",TObject::kOverwrite);
  minClusters->Write("minCluster",TObject::kOverwrite);

  criticalDistanceEEwgt->Write("critDistEEwgt",TObject::kOverwrite);
  criticalDistanceFHwgt->Write("critDistFHwgt",TObject::kOverwrite);
  criticalDistanceBHwgt->Write("critDistBHwgt",TObject::kOverwrite);
  dependSensorwgt->Write("sensorwgt",TObject::kOverwrite);
  reachedEEwgt->Write("reachedEEwgt",TObject::kOverwrite);
  kernelwgt->Write("kernelwgt",TObject::kOverwrite);
  deltacEEwgt->Write("deltaEEwgt",TObject::kOverwrite);
  deltacFHwgt->Write("deltaFHwgt",TObject::kOverwrite);
  deltacBHwgt->Write("deltaBHwgt",TObject::kOverwrite);
  kappawgt->Write("kappawgt",TObject::kOverwrite);
  energyMinwgt->Write("eMinwgt",TObject::kOverwrite);
  matchingDistancewgt->Write("matchingDistwgt",TObject::kOverwrite);
  minClusterswgt->Write("minClusterwgt",TObject::kOverwrite);
  
  outfile->Close();
  delete outfile;
}

int main(int argc, char* argv[])
{
  if(argc > 1 && argc!=9){
    cout<<"Usage:"<<endl;
    cout<<"./geneticOptimizer"<<endl;
    cout<<"or with custom parameters:"<<endl;
    cout<<"./ geneticOptimizer populationSize nGenerations nEventsPerTest processTimeout mutationChance severityFactor dataPath outputPath"<<endl;
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
  
  // Set number of events for each member to be tested on
  UpdateParamValue("baseConfig.md", "analyze_events_per_tuple",nEventsPerTest);
  
  srand((unsigned int)time(0));
  randGenerator.seed((unsigned int)time(0));
  
  TH1D *scoresDist[nGenerations];
  scoresMean = new TGraph();
  criticalDistanceEE = new TH2D("critical distance EE","critical distance EE",nGenerations,0,nGenerations,100,0.0,criticalDistanceEEmax);
  criticalDistanceFH = new TH2D("critical distance FH","critical distance FH",nGenerations,0,nGenerations,100,0.0,criticalDistanceFHmax);
  criticalDistanceBH = new TH2D("critical distance BH","critical distance BH",nGenerations,0,nGenerations,100,0.0,criticalDistanceBHmax);
  dependSensor = new TH2D("depend sensor","depend sensor",nGenerations,0,nGenerations,2,0,2);
  reachedEE = new TH2D("reached EE","reached EE",nGenerations,0,nGenerations,2, 0, 2);
  kernel = new TH2D("kernel","kernel",nGenerations,0,nGenerations,3, 0, 3);
  deltacEE = new TH2D("delta c EE","delta c EE",nGenerations,0,nGenerations,100, 0.0,deltacEEmax);
  deltacFH = new TH2D("delta c FH","delta c FH",nGenerations,0,nGenerations,100, 0.0,deltacFHmax);
  deltacBH = new TH2D("delta c BH","delta c BH",nGenerations,0,nGenerations,100, 0.0,deltacBHmax);
  kappa = new TH2D("kappa","kappa",nGenerations,0,nGenerations,10000, 0.0,kappaMax);
  energyMin = new TH2D("energy threshold","energy threshold",nGenerations,0,nGenerations,10000,0.0,energyThresholdMax);
  matchingDistance = new TH2D("matching distance","matching distance",nGenerations,0,nGenerations,100, 0.0,matchingDistanceMax);
  minClusters = new TH2D("min clusters","min clusters",nGenerations,0,nGenerations,10, 0,10);
  
  criticalDistanceEEwgt = new TH2D("critical distance EE wgt","critical distance EE wgt",nGenerations,0,nGenerations,100,0.0,criticalDistanceEEmax);
  criticalDistanceFHwgt = new TH2D("critical distance FH wgt","critical distance FH wgt",nGenerations,0,nGenerations,100,0.0,criticalDistanceFHmax);
  criticalDistanceBHwgt = new TH2D("critical distance BH wgt","critical distance BH wgt",nGenerations,0,nGenerations,100,0.0,criticalDistanceBHmax);
  dependSensorwgt = new TH2D("depend sensor wgt","depend sensor wgt",nGenerations,0,nGenerations,2,0,2);
  reachedEEwgt = new TH2D("reached EE wgt","reached EE wgt",nGenerations,0,nGenerations,2, 0, 2);
  kernelwgt = new TH2D("kernel wgt","kernel wgt",nGenerations,0,nGenerations,3, 0, 3);
  deltacEEwgt = new TH2D("delta c EE wgt","delta c EE wgt",nGenerations,0,nGenerations,100, 0.0,deltacEEmax);
  deltacFHwgt = new TH2D("delta c FH wgt","delta c FH wgt",nGenerations,0,nGenerations,100, 0.0,deltacFHmax);
  deltacBHwgt = new TH2D("delta c BH wgt","delta c BH wgt",nGenerations,0,nGenerations,100, 0.0,deltacBHmax);
  kappawgt = new TH2D("kappa wgt","kappa wgt",nGenerations,0,nGenerations,10000, 0.0,kappaMax);
  energyMinwgt = new TH2D("energy threshold wgt","energy threshold wgt",nGenerations,0,nGenerations,10000,0.0,energyThresholdMax);
  matchingDistancewgt = new TH2D("matching distance wgt","matching distance wgt",nGenerations,0,nGenerations,100, 0.0,matchingDistanceMax);
  minClusterswgt = new TH2D("min clusters wgt","min clusters wgt",nGenerations,0,nGenerations,10, 0,10);
  
  // make sure to recreate output file
  outfile = new TFile("geneticHists.root","recreate");
  outfile->Close();
  delete outfile;
  
  vector<Chromosome*> population;
  discrete_distribution<double> scores;
  
  for(int iGeneration=0;iGeneration<nGenerations;iGeneration++){
    cout<<"\n\nGeneration "<<iGeneration<<"\n\n"<<endl;
    
    if(iGeneration==0){
      // draw initial population
      for(int i=0;i<populationSize;i++){
        Chromosome *chromo;
        chromo = Chromosome::GetRandom();
        chromo->SaveToBitChromosome();
        chromo->SetMutationChance(mutationChance);
        chromo->SetSeverityFactor(severityFactor);
        population.push_back(chromo);
      }
    }
    else{
      // draw new population
      population.clear();
      
      for(int i=0;i<populationSize;i++){
        Chromosome *mom, *dad;
        
        int iter=0;
        do {
          if(iter>10) cout<<"mom==dad for "<<iter<<" times. Probably something went wrong..."<<endl;
          if(iter>20) exit(0);
          iter++;
          
          mom = population[(int)scores(randGenerator)];
          dad = population[(int)scores(randGenerator)];
        } while (mom==dad);

        Chromosome *child = dad->ProduceChildWith(mom);
        population.push_back(child);
      }
    }
    
    scoresDist[iGeneration] = new TH1D(Form("score dist gen[%i]",iGeneration),
                                       Form("score dist gen[%i]",iGeneration),
                                       200,-20,20);
    
    TestPopulation(population, scoresDist[iGeneration], scores);
    
    for(int i=0;i<populationSize;i++){
      criticalDistanceEE->Fill(iGeneration,population[i]->GetCriticalDistanceEE());
      criticalDistanceFH->Fill(iGeneration,population[i]->GetCriticalDistanceFH());
      criticalDistanceBH->Fill(iGeneration,population[i]->GetCriticalDistanceBH());
      dependSensor->Fill(iGeneration,population[i]->GetDependSensor());
      reachedEE->Fill(iGeneration,population[i]->GetReachedEE());
      kernel->Fill(iGeneration,population[i]->GetKernel());
      deltacEE->Fill(iGeneration,population[i]->GetDeltacEE());
      deltacFH->Fill(iGeneration,population[i]->GetDeltacFH());
      deltacBH->Fill(iGeneration,population[i]->GetDeltacBH());
      kappa->Fill(iGeneration,population[i]->GetKappa());
      energyMin->Fill(iGeneration,population[i]->GetEnergyMin());
      matchingDistance->Fill(iGeneration,population[i]->GetMatchingDistance());
      minClusters->Fill(iGeneration,population[i]->GetMinClusters());
      
      criticalDistanceEEwgt->Fill(iGeneration,population[i]->GetCriticalDistanceEE(),population[i]->GetScore());
      criticalDistanceFHwgt->Fill(iGeneration,population[i]->GetCriticalDistanceFH(),population[i]->GetScore());
      criticalDistanceBHwgt->Fill(iGeneration,population[i]->GetCriticalDistanceBH(),population[i]->GetScore());
      dependSensorwgt->Fill(iGeneration,population[i]->GetDependSensor(),population[i]->GetScore());
      reachedEEwgt->Fill(iGeneration,population[i]->GetReachedEE(),population[i]->GetScore());
      kernelwgt->Fill(iGeneration,population[i]->GetKernel(),population[i]->GetScore());
      deltacEEwgt->Fill(iGeneration,population[i]->GetDeltacEE(),population[i]->GetScore());
      deltacFHwgt->Fill(iGeneration,population[i]->GetDeltacFH(),population[i]->GetScore());
      deltacBHwgt->Fill(iGeneration,population[i]->GetDeltacBH(),population[i]->GetScore());
      kappawgt->Fill(iGeneration,population[i]->GetKappa(),population[i]->GetScore());
      energyMinwgt->Fill(iGeneration,population[i]->GetEnergyMin(),population[i]->GetScore());
      matchingDistancewgt->Fill(iGeneration,population[i]->GetMatchingDistance(),population[i]->GetScore());
      minClusterswgt->Fill(iGeneration,population[i]->GetMinClusters(),population[i]->GetScore());
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
  
  cout<<"plotting"<<endl;
  TCanvas *c1 = new TCanvas("c1","c1",2880,1800);
  c1->Divide(4,4);
  
  c1->cd(1);
  scoresMean->SetMarkerSize(1.0);
  scoresMean->SetMarkerStyle(20);
  scoresMean->SetMarkerColor(kGreen+2);
  scoresMean->Draw("AP");

  c1->cd(2);
  criticalDistanceEE->Draw("colz");
  c1->cd(3);
  criticalDistanceFH->Draw("colz");
  c1->cd(4);
  criticalDistanceBH->Draw("colz");
  c1->cd(5);
  dependSensor->Draw("colz");
  c1->cd(6);
  reachedEE->Draw("colz");
  c1->cd(7);
  kernel->Draw("colz");
  c1->cd(8);
  deltacEE->Draw("colz");
  c1->cd(9);
  deltacFH->Draw("colz");
  c1->cd(10);
  deltacBH->Draw("colz");
  c1->cd(11);
  kappa->Draw("colz");
  c1->cd(12);
  energyMin->Draw("colz");
  c1->cd(13);
  matchingDistance->Draw("colz");
  c1->cd(14);
  minClusters->Draw("colz");
  c1->cd(15);
  
  c1->Update();
  
  theApp.Run();
  return 0;
}
