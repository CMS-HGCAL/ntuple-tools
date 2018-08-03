#include "GeneticHelpers.hpp"

#include "Chromosome.hpp"

#include <TMath.h>
#include <TH1D.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TGraph.h>

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

const string configPath = "baseConfig.md";
const string outputPath = "autoGenOutput.txt";

const int populationSize = 10;  ///< Size of the population, will stay the same for all generations
const int nGenerations = 2;     ///< Number of iterations
const int nEventsPerTest = 1;   ///< On how many events each population member will be tested

mt19937 randGenerator;

void threadFunction(Chromosome *chromo)
{
  chromo->CalculateScore();
}

int GetWeightedRandom(discrete_distribution<double> dist)
{
  return dist(randGenerator);
}

void TestPopulation(Chromosome *population[populationSize], TH1D *hist, discrete_distribution<double> &dist)
{
  // Starting Threads & move the future object in lambda function by reference
  thread *threads[populationSize];
  
  //  thread::native_handle_type threadHandles[populationSize];
  
  for(int i=0;i<populationSize;i++){
    threads[i] = new thread(threadFunction, population[i]);
    //    threadHandles[i] = threads[i]->native_handle();
    //    threads[i]->detach();
  }
  
  //Wait for 10 sec and then ask thread to join
  //  this_thread::sleep_for(std::chrono::seconds(10));
  //  cout<<"\n\nTime passed\n\n"<<endl;
  
  vector<double> scores;
  vector<double> scoresNormalized;
  double minScore=99999, maxScore=-99999;
  
  for(int i=0;i<populationSize;i++){
    //    pthread_cancel(threadHandles[i]);
    //    pthread_kill(threadHandles[i], 1);
    //    if(threads[i]->joinable()) threads[i]->join();
    
    threads[i]->join();
    double score = population[i]->GetScore();
    if(score < minScore && score > -100) minScore = score; // make sure to remove veeery bad results
    if(score > maxScore) maxScore = score;
    hist->Fill(score);
    scores.push_back(score);
  }
  
  // re-assing normalized points to population members
  double normScore;

  for(int i=0;i<populationSize;i++){
    if(scores[i] > -100) normScore = (scores[i]-minScore)/(maxScore-minScore);
    else normScore = 0; // kill extremally weak members of population
    
    population[i]->SetScore(normScore);
    scoresNormalized.push_back(normScore);
  }

  dist = discrete_distribution<double>(scoresNormalized.begin(), scoresNormalized.end());
}



int main(int argc, char* argv[])
{
  TApplication theApp("App", &argc, argv);
  
  // Set number of events for each member to be tested on
  UpdateParamValue("baseConfig.md", "analyze_events_per_tuple",nEventsPerTest);
  
  srand((unsigned int)time(0));
  randGenerator.seed((unsigned int)time(0));
  
  TH1D *scoresDist[nGenerations];
  
  Chromosome* population[populationSize];
  discrete_distribution<double> scores;
  
  for(int iGeneration=0;iGeneration<nGenerations;iGeneration++){
    cout<<"\n\nGeneration "<<iGeneration<<"\n\n"<<endl;
    
    if(iGeneration==0){
      // draw initial population
      for(int i=0;i<populationSize;i++){
        Chromosome *chromo;
        chromo = Chromosome::GetRandom();
        chromo->SaveToBitChromosome();
        population[i] = chromo;
      }
    }
    else{
      // draw new population
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
        delete population[i];
        population[i] = child;
      }
    }
    
    scoresDist[iGeneration] = new TH1D(Form("score dist gen[%i]",iGeneration),
                                       Form("score dist gen[%i]",iGeneration),
                                       200,-20,20);
    
    TestPopulation(population, scoresDist[iGeneration], scores);
    
    if(iGeneration == 0 || iGeneration == nGenerations-1){
      cout<<"\n\nPopulation:\n\n"<<endl;
      for(int i=0;i<populationSize;i++){
        population[i]->Print();
      }
    }
  }
  
  cout<<"plotting"<<endl;
  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  TGraph *scoresMean = new TGraph();
  
  c1->Divide(1,2);
  c1->cd(1);
  for(int iGeneration=0;iGeneration<nGenerations;iGeneration++){
    scoresDist[iGeneration]->Draw(iGeneration==0 ? "" : "same");
    scoresDist[iGeneration]->SetLineColor(iGeneration);
    
    scoresMean->SetPoint(iGeneration, iGeneration, scoresDist[iGeneration]->GetMean());
  }
  c1->cd(2);
  scoresMean->SetMarkerSize(1.0);
  scoresMean->SetMarkerStyle(20);
  scoresMean->SetMarkerColor(kGreen+2);
  scoresMean->Draw("AP");
  
  c1->Update();
  theApp.Run();
  return 0;
}
