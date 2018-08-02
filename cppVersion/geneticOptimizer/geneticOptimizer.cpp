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

using namespace std;

const string configPath = "../configs/autoGenConfig.md";
const string outputPath = "autoGenOutput.txt";

const int initPopulationSize = 3;

void threadFunction(Chromosome *chromo)
{
  chromo->CalculateScore();
}

int main(int argc, char* argv[])
{
  TApplication theApp("App", &argc, argv);
  
  srand((unsigned int)time(0));
  
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
  
  // calculate scores for each member of initial population
  vector<float> initPopulationScores;
  float minScore=99999999, maxScore=-99999999;

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
  
  for(int i=0;i<initPopulationSize;i++){
//    pthread_cancel(threadHandles[i]);
//    pthread_kill(threadHandles[i], 1);
//    if(threads[i]->joinable()) threads[i]->join();
    
    threads[i]->join();
    
    float score = initPopulation[i]->GetScore();
    initPopulationScores.push_back(score);
    if(score < minScore) minScore = score;
    if(score > maxScore) maxScore = score;
    scoresDist->Fill(score);
  }

  // re-assing normalized points to population members
  for(int i=0;i<initPopulationSize;i++){
    initPopulation[i]->SetScore((initPopulationScores[i]-minScore)/(maxScore-minScore));
  }
    
  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  c1->cd();
  scoresDist->Draw();
  
  c1->Update();
  theApp.Run();
  return 0;
}
