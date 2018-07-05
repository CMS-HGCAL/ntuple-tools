#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TRandom.h>
#include <TROOT.h>
#include <TList.h>
#include <TSystemDirectory.h>

#include <iostream>

using namespace std;

string baseDir = "clusteringResultsCXX/twoPions_Pt80_Eta2_DeltaR0p4";

vector<TH2D*> getHistsWithName(const char* histFileName, const char* histName)
{
  vector<TH2D*> result;
  

  TSystemDirectory mainDir(baseDir.c_str(),baseDir.c_str());
  TList *ntupleDirsList = mainDir.GetListOfFiles();
  
  ntupleDirsList->Print();
  
  TSystemDirectory *ntupleDir;
  TIter next(ntupleDirsList);
  
  while((ntupleDir=(TSystemDirectory*)next())){
    TString name(ntupleDir->GetName());
    if(name.Contains("ntup")){
      TSystemDirectory eventsDir(Form("%s/%s",baseDir.c_str(),name.Data()),
                                 Form("%s/%s",baseDir.c_str(),name.Data()));
      TList *eventsDirsList = eventsDir.GetListOfFiles();
      
      TSystemDirectory *eventDir;
      TIter nextEvent(eventsDirsList);
      while((eventDir=(TSystemDirectory*)nextEvent())){
        TString eventName(eventDir->GetName());
        if(eventName.Contains("event")){
          cout<<Form("%s/%s/%s/%s.root",baseDir.c_str(),name.Data(),eventName.Data(),histFileName)<<endl;
          TFile *inFile = TFile::Open(Form("%s/%s/%s/%s.root",baseDir.c_str(),name.Data(),eventName.Data(),histFileName));
          if(!inFile) continue;
          TH2D *hist = (TH2D*)inFile->Get(histName);
          if(!hist){
            cout<<"no hist: "<<histName<<"!!"<<endl;
            continue;
          }
          result.push_back(hist);
        }
      }
    }
  }
  return result;
}

void plotMatchingBenchmark()
{
  TCanvas *canvas = new TCanvas("canvas","canvas",1200,800);
  canvas->Divide(2,2);
  TF1 *fun = new TF1("fun","x",0,100);
  
  vector<TH2D*> inputEnergyComparisonHists = getHistsWithName("energyComparisonNoMatchingHist","no matching");
  
  if(inputEnergyComparisonHists.size() > 0){
    TH2D *mergedEnergyComparisonHist = new TH2D(*inputEnergyComparisonHists[0]);

    for(int iter=1;iter<inputEnergyComparisonHists.size();iter++){
      mergedEnergyComparisonHist->Add(inputEnergyComparisonHists[iter]);
    }
    
    canvas->cd(1);
    mergedEnergyComparisonHist->Draw("colz");
    
    mergedEnergyComparisonHist->GetXaxis()->SetTitle("E_{rec} (GeV)");
    mergedEnergyComparisonHist->GetYaxis()->SetTitle("E_{sim} (Gev)");
    mergedEnergyComparisonHist->GetZaxis()->SetRangeUser(0,500);
    
    fun->Draw("same");
  }

  vector<TH2D*> inputEnergyComparisonOverlapHists = getHistsWithName("energyComparisonClosestHist","closest rec cluster");
  
  if(inputEnergyComparisonOverlapHists.size() >0){
    TH2D *mergedEnergyComparisonOverlapHist = new TH2D(*inputEnergyComparisonOverlapHists[0]);
    
    for(int iter=1;iter<inputEnergyComparisonOverlapHists.size();iter++){
      mergedEnergyComparisonOverlapHist->Add(inputEnergyComparisonOverlapHists[iter]);
    }
    canvas->cd(2);
    mergedEnergyComparisonOverlapHist->Draw("colz");
    
    mergedEnergyComparisonOverlapHist->GetXaxis()->SetTitle("E_{rec} (GeV)");
    mergedEnergyComparisonOverlapHist->GetYaxis()->SetTitle("E_{sim} (GeV)");
    mergedEnergyComparisonOverlapHist->GetZaxis()->SetRangeUser(0,30);
    
    fun->Draw("same");
  }
  
}

