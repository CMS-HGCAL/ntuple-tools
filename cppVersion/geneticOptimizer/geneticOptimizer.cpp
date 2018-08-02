#include "include/Helpers.hpp"

#include "Chromosome.hpp"

#include <TMath.h>
#include <TH1D.h>
#include <TApplication.h>
#include <TCanvas.h>

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

using namespace std;

const string configPath = "../configs/autoGenConfig.md";
const string outputPath = "autoGenOutput.txt";

const int initPopulationSize = 10;
mt19937 gen;

void threadFunction(Chromosome *chromo)
{
  chromo->CalculateScore();
}

int GetWeightedRandom(discrete_distribution<double> dist)
{
  return dist(gen);
}

int main(int argc, char* argv[])
{
  TApplication theApp("App", &argc, argv);
  
  srand((unsigned int)time(0));
  gen.seed((unsigned int)time(0));
  
  Chromosome* initPopulation[initPopulationSize];
  
  // draw initial population
  for(int i=0;i<initPopulationSize;i++){
    Chromosome *chromo;
    chromo = Chromosome::GetRandom();
    chromo->SaveToBitChromosome();
    chromo->Print();
    initPopulation[i] = chromo;
  }
  TH1D *scoresDist = new TH1D("score dist","score dist",200,-20,20);
  
  // Starting Threads & move the future object in lambda function by reference
  thread *threads[initPopulationSize];

//  thread::native_handle_type threadHandles[initPopulationSize];
  
  for(int i=0;i<initPopulationSize;i++){
    threads[i] = new thread(threadFunction, initPopulation[i]);
//    threadHandles[i] = threads[i]->native_handle();
//    threads[i]->detach();
  }
  
  //Wait for 10 sec and then ask thread to join
//  this_thread::sleep_for(std::chrono::seconds(10));
//  cout<<"\n\nTime passed\n\n"<<endl;
  
  vector<double> scores;
  vector<double> scoresNormalized;
  double minScore=99999, maxScore=-99999;
  
  for(int i=0;i<initPopulationSize;i++){
//    pthread_cancel(threadHandles[i]);
//    pthread_kill(threadHandles[i], 1);
//    if(threads[i]->joinable()) threads[i]->join();
    
    threads[i]->join();
    
    cout<<"filling scores"<<endl;
    double score = initPopulation[i]->GetScore();
    if(score < minScore && score > -100) minScore = score; // make sure to remove veeery bad results
    if(score > maxScore) maxScore = score;
    scoresDist->Fill(score);
    scores.push_back(score);
  }

  // re-assing normalized points to population members
  
  cout<<"min score:"<<minScore<<endl;
  cout<<"max score:"<<maxScore<<endl;
  
  double normScore;
  
  for(int i=0;i<initPopulationSize;i++){
    if(scores[i] > -100) normScore = (scores[i]-minScore)/(maxScore-minScore);
    else normScore = 0; // kill extremally weak members of population
      
    cout<<"Original score:"<<scores[i]<<endl;
    cout<<"Normalized score:"<<normScore<<endl;
    initPopulation[i]->SetScore(normScore);
    scoresNormalized.push_back(normScore);
  }

  cout<<"plotting"<<endl;
  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  c1->cd();
  scoresDist->Draw();
  
  // draw new population
  Chromosome *newPopulation[initPopulationSize];
  discrete_distribution<double> dist(scoresNormalized.begin(), scoresNormalized.end());
  
  for(int i=0;i<initPopulationSize;i++){
    Chromosome *chromo1 = initPopulation[GetWeightedRandom(dist)];
    Chromosome *chromo2 = initPopulation[GetWeightedRandom(dist)];
    
    cout<<"Pair "<<i<<endl;
    chromo1->Print();
    chromo2->Print();
  }
  
  
  c1->Update();
  theApp.Run();
  return 0;
}
