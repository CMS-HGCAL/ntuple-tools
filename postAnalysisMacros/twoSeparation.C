//
//  twoSeparation.c
//
//  Created by Jeremi Niedziela on 05/07/2018.
//


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

const int nDeltaRbins = 7;

string baseDir = "clusteringResultsCXX/twoPions_Pt80_Eta2_DeltaR";
string pathSuffix[nDeltaRbins] = {"0p1","0p15","0p2","0p25","0p3","0p35","0p4"};

const double deltaR[nDeltaRbins] = {0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4};
const int colors[nDeltaRbins] = {kRed, kGreen+1, kOrange+1, kBlue, kBlack, kCyan, kMagenta};

vector<TH1D*> get1DHistsWithName(const char* histFileName, const char* histName, const char* histDir)
{
  vector<TH1D*> result;
  
  TSystemDirectory mainDir(histDir,histDir);
  TList *ntupleDirsList = mainDir.GetListOfFiles();
  
  ntupleDirsList->Print();
  
  TSystemDirectory *ntupleDir;
  TIter next(ntupleDirsList);
  
  while((ntupleDir=(TSystemDirectory*)next())){
    TString name(ntupleDir->GetName());
    if(name.Contains("ntup")){
      TSystemDirectory eventsDir(Form("%s/%s",histDir,name.Data()),
                                 Form("%s/%s",histDir,name.Data()));
      TList *eventsDirsList = eventsDir.GetListOfFiles();
      
      TSystemDirectory *eventDir;
      TIter nextEvent(eventsDirsList);
      while((eventDir=(TSystemDirectory*)nextEvent())){
        TString eventName(eventDir->GetName());
        if(eventName.Contains("event")){
          cout<<Form("%s/%s/%s/%s.root",histDir,name.Data(),eventName.Data(),histFileName)<<endl;
          TFile *inFile = TFile::Open(Form("%s/%s/%s/%s.root",histDir,name.Data(),eventName.Data(),histFileName));
          if(!inFile) continue;
          TH1D *hist = (TH1D*)inFile->Get(histName);
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

void twoSeparation()
{
  TCanvas *canvas = new TCanvas("canvas","canvas",1200,800);
    canvas->Divide(2,1);
  canvas->cd(1);
  
  TF1 *fitFun = new TF1("fitFun","[0]/x+[1]",0.1,100);
  fitFun->SetParameter(0,1);
  fitFun->SetParameter(1,0);
  
  TGraph *distWidthVsDeltaR = new TGraph();
  
  for(int iDelta=0;iDelta<nDeltaRbins;iDelta++){
    string histDir = baseDir + pathSuffix[iDelta];
    
    vector<TH1D*> inputSeparation = get1DHistsWithName("twoSeparation","two clusters separation",histDir.c_str());
    if(inputSeparation.size() > 0){
      TH1D *mergedSeparationHists = new TH1D(*inputSeparation[0]);
      
      for(int iter=1;iter<inputSeparation.size();iter++){
        mergedSeparationHists->Add(inputSeparation[iter]);
      }
      
      mergedSeparationHists->SetLineColor(colors[iDelta]);
      mergedSeparationHists->Draw(iDelta == 0 ? "" : "same");
      mergedSeparationHists->GetXaxis()->SetTitle("sep");
      mergedSeparationHists->GetXaxis()->SetRangeUser(0,10);
      mergedSeparationHists->Fit(fitFun);
      
      distWidthVsDeltaR->SetPoint(iDelta,deltaR[iDelta],fitFun->GetParameter(0));
    }
  }
  
  canvas->cd(2);
  
  distWidthVsDeltaR->SetMarkerStyle(20);
  distWidthVsDeltaR->SetMarkerSize(1.0);
  distWidthVsDeltaR->SetMarkerColor(kGreen+1);
  distWidthVsDeltaR->Draw("AP");
  distWidthVsDeltaR->GetXaxis()->SetTitle("#Delta R");
  distWidthVsDeltaR->GetYaxis()->SetTitle("Fit offset");
}

