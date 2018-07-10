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
#include <TGraphErrors.h>
#include <TRandom.h>
#include <TROOT.h>
#include <TList.h>
#include <TSystemDirectory.h>
#include <TLegend.h>

#include <iostream>

using namespace std;

const int nDeltaRbins = 7;

string baseDir = "clusteringResultsCXX/twoPions_Pt80_Eta2_DeltaR";
string pathSuffix[nDeltaRbins] = {"0p1","0p15","0p2","0p25","0p3","0p35","0p4"};
bool runDir[nDeltaRbins]       = {  1  ,  1   ,  1  ,  1   ,  1  ,  1   ,  1  };

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
  double legendW = 0.15;
  double legendH = 0.5;
  double legendX = 0.75;
  double legendY = 0.25;
  TLegend *leg = new TLegend(legendX,legendY,legendX+legendW,legendY+legendH);
  leg->SetHeader("#Delta R(#eta,#phi):");
  
  TCanvas *canvas = new TCanvas("canvas","canvas",2800,1800);
  canvas->Divide(2,2);
  canvas->cd(1);
  
//  TF1 *fitFun = new TF1("fitFun","[0]/x+[1]",0.1,100);
//  TF1 *fitFun = new TF1("fitFun","[0]/(sqrt(2*TMath::Pi())*[1])*exp(-pow(x-[2],2)/(2*[1]*[1]))+[3]/x+[4]",0.1,100);
  TF1 *fitFun = new TF1("fitFun","[0]/(sqrt(2*TMath::Pi())*[1])*exp(-pow(x-[2],2)/(2*[1]*[1]))+[4]",0.1,100);
  fitFun->SetParameter(0,260); // scaling of gaussian
  fitFun->SetParameter(1,0.5);   // sigma
  fitFun->SetParameter(2,1);   // miu
//  fitFun->SetParameter(3,30); // scaling of 1/x
  fitFun->SetParameter(4,20);  // offset
  
  TGraphErrors *distWidthVsDeltaR = new TGraphErrors();
  TGraphErrors *distMiuVsDeltaR = new TGraphErrors();
  
  for(int iDelta=0;iDelta<nDeltaRbins;iDelta++){
    if(!runDir[iDelta]) continue;
    string histDir = baseDir + pathSuffix[iDelta];
    
    vector<TH1D*> inputSeparation = get1DHistsWithName("twoSeparationJer","two clusters separation",histDir.c_str());
    if(inputSeparation.size() > 0){
      TH1D *mergedSeparationHists = new TH1D(*inputSeparation[0]);
      
      for(int iter=1;iter<inputSeparation.size();iter++){
        mergedSeparationHists->Add(inputSeparation[iter]);
      }
      mergedSeparationHists->SetTitle("Two-cluster separation");
      mergedSeparationHists->SetLineColor(colors[iDelta]);
      mergedSeparationHists->Rebin(5);
      mergedSeparationHists->Draw(iDelta == 0 ? "" : "same");
      mergedSeparationHists->GetXaxis()->SetTitle("Separation");
      mergedSeparationHists->GetXaxis()->SetRangeUser(0,4);
      mergedSeparationHists->GetYaxis()->SetRangeUser(0,600);
      mergedSeparationHists->GetXaxis()->SetTitleSize(0.06);
      mergedSeparationHists->GetXaxis()->SetTitleOffset(0.75);
      
      mergedSeparationHists->Fit(fitFun,"","",0.03,100);
      fitFun->SetLineColor(colors[iDelta]);
      fitFun->DrawCopy("same");
      distWidthVsDeltaR->SetPoint(iDelta,deltaR[iDelta],fitFun->GetParameter(1));
      distWidthVsDeltaR->SetPointError(iDelta,0,fitFun->GetParError(1));
      distMiuVsDeltaR->SetPoint(iDelta,deltaR[iDelta],fitFun->GetParameter(2));
      distMiuVsDeltaR->SetPointError(iDelta,0,fitFun->GetParError(2));
      
      leg->AddEntry(mergedSeparationHists,Form("%.2f",deltaR[iDelta]),"pl");
    }
  }
  leg->Draw();
  
  canvas->cd(2);
  
  distWidthVsDeltaR->SetMarkerStyle(20);
  distWidthVsDeltaR->SetMarkerSize(1.0);
  distWidthVsDeltaR->SetMarkerColor(kGreen+2);
  distWidthVsDeltaR->Draw("APE");
  distWidthVsDeltaR->SetTitle("Gaussian width vs. #Delta R(#eta,#phi)");
  distWidthVsDeltaR->GetXaxis()->SetTitle("#Delta R(#eta,#phi)");
  distWidthVsDeltaR->GetYaxis()->SetTitle("#sigma");
  
  distWidthVsDeltaR->GetXaxis()->SetTitleSize(0.06);
  distWidthVsDeltaR->GetXaxis()->SetTitleOffset(0.75);
  distWidthVsDeltaR->GetYaxis()->SetTitleSize(0.06);
  distWidthVsDeltaR->GetYaxis()->SetTitleOffset(0.75);
  
  
  canvas->cd(3);
  
  distMiuVsDeltaR->SetMarkerStyle(20);
  distMiuVsDeltaR->SetMarkerSize(1.0);
  distMiuVsDeltaR->SetMarkerColor(kGreen+2);
  distMiuVsDeltaR->Draw("APE");
  distMiuVsDeltaR->SetTitle("Gaussian position vs. #Delta R(#eta,#phi)");
  distMiuVsDeltaR->GetXaxis()->SetTitle("#Delta R(#eta,#phi)");
  distMiuVsDeltaR->GetYaxis()->SetTitle("#mu");
  
  distMiuVsDeltaR->GetXaxis()->SetTitleSize(0.06);
  distMiuVsDeltaR->GetXaxis()->SetTitleOffset(0.75);
  distMiuVsDeltaR->GetYaxis()->SetTitleSize(0.06);
  distMiuVsDeltaR->GetYaxis()->SetTitleOffset(0.75);
}

