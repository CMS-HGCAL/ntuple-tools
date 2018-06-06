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

const int nLayers = 40;
const int nClusters = 40;

vector<TH2D*> getHistsWithName(const char* histFileName, const char* histName)
{
  vector<TH2D*> result;
  string baseDir = "clusteringResults";
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
          cout<<Form("%s/%s/%s/energyComparisonHist.root",baseDir.c_str(),name.Data(),eventName.Data())<<endl;
          TFile *inFile = TFile::Open(Form("%s/%s/%s/%s.root",baseDir.c_str(),name.Data(),eventName.Data(),histFileName));
          if(!inFile) continue;
          TH2D *hist = (TH2D*)inFile->Get(histName);
          
          result.push_back(hist);
        }
      }
    }
  }
  return result;
}

void analyzeClustering()
{
  vector<TH2D*> inputEnergyComparisonHists = getHistsWithName("energyComparisonHist","energy comparison");
  TH2D *mergedEnergyComparisonHist = new TH2D(*inputEnergyComparisonHists[0]);
  
  vector<TH2D*> inputEnergyComparisonOverlapHists = getHistsWithName("energyComparisonOverlapHist","energy comparison overlap.");
  TH2D *mergedEnergyComparisonOverlapHist = new TH2D(*inputEnergyComparisonOverlapHists[0]);
  
  for(int iter=1;iter<inputEnergyComparisonHists.size();iter++){
    mergedEnergyComparisonHist->Add(inputEnergyComparisonHists[iter]);
  }
  
  for(int iter=1;iter<inputEnergyComparisonOverlapHists.size();iter++){
    mergedEnergyComparisonOverlapHist->Add(inputEnergyComparisonOverlapHists[iter]);
  }
  
  TCanvas *canvas = new TCanvas("canvas","canvas",1200,800);
  canvas->Divide(2,1);
  
  canvas->cd(1);
  mergedEnergyComparisonHist->Draw("colz");
  mergedEnergyComparisonHist->GetXaxis()->SetTitle("E_{rec}");
  mergedEnergyComparisonHist->GetYaxis()->SetTitle("E_{sim}");
  mergedEnergyComparisonHist->GetZaxis()->SetRangeUser(0,100);
  
  TF1 *fun = new TF1("fun","x",0,100);
  fun->Draw("same");
  
  canvas->cd(2);
  mergedEnergyComparisonOverlapHist->Draw("colz");
  mergedEnergyComparisonOverlapHist->GetXaxis()->SetTitle("E_{rec}");
  mergedEnergyComparisonOverlapHist->GetYaxis()->SetTitle("E_{sim}");
  mergedEnergyComparisonOverlapHist->GetZaxis()->SetRangeUser(0,100);
  
  fun->Draw("same");
  
  
  
  
  
}

