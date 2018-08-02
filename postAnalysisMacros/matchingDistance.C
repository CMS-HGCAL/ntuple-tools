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

const int nDistBins = 7;

string baseDir = "clusteringResultsCXX/twoPions_jniedzie_max_";
string pathSuffix[nDistBins] = {"1000","10","5","4","3","2","1"};
bool runDir[nDistBins]       = {  1   , 1  , 1 , 1 , 1 , 1 , 0 };

const double maxDist[nDistBins] = {1000, 10, 5, 4, 3, 2, 1};
const int colors[nDistBins] = {kRed, kGreen+1, kOrange+1, kBlue, kBlack, kCyan, kMagenta};

vector<TObject*> getHistsWithName(const char* histFileName, const char* histName, const char* histDir)
{
  vector<TObject*> result;
  
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
          TObject *hist = (TObject*)inFile->Get(histName);
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

void matchingDistance()
{
  double legendW = 0.15;
  double legendH = 0.5;
  double legendX = 0.75;
  double legendY = 0.25;
  TLegend *leg = new TLegend(legendX,legendY,legendX+legendW,legendY+legendH);
  leg->SetHeader("Max distance:");
  
  TCanvas *canvas = new TCanvas("canvas","canvas",2800,1800);
  canvas->Divide(1,1);
  
  
  for(int iDist=0;iDist<nDistBins;iDist++){
    if(!runDir[iDist]) continue;
    
    canvas->cd(1);
    string histDir = baseDir + pathSuffix[iDist];
    
    vector<TObject*> inputSigmaE = getHistsWithName("simgaEVsEtaEsim","sigma(E)Esim vs. eta",histDir.c_str());
    if(inputSigmaE.size() > 0){
      TH2D *mergedSigmaE = new TH2D(*(TH2D*)inputSigmaE[0]);
      
      for(int iter=1;iter<inputSigmaE.size();iter++){
        mergedSigmaE->Add((TH2D*)inputSigmaE[iter]);
      }
      
      mergedSigmaE->SetLineColor(colors[iDist]);
      mergedSigmaE->GetXaxis()->SetTitle("|#eta|");
      mergedSigmaE->GetYaxis()->SetTitle("#Delta E/E_{sim}");
      
      TH1D *sigmaE = mergedSigmaE->ProjectionY();
      sigmaE->Rebin(2);
      sigmaE->SetLineColor(colors[iDist]);
      sigmaE->SetTitle("#DeltaE/E_{sim} distribution");
      sigmaE->SetName("#DeltaE/E_{sim} distribution");
      
      sigmaE->GetXaxis()->SetRangeUser(-1,0.4);
      sigmaE->GetXaxis()->SetTitleSize(0.06);
      sigmaE->GetXaxis()->SetTitleOffset(0.75);
      
      sigmaE->GetYaxis()->SetRangeUser(0,450);
      
      sigmaE->DrawCopy(iDist == 0 ? "" : "same");
      
      leg->AddEntry(sigmaE,Form("%.0f cm",maxDist[iDist]),"pl");
      continue;
      
      canvas->cd(2);
      string histDir = baseDir + pathSuffix[iDist];
      
      vector<TObject*> inputSeparation = getHistsWithName("twoSeparation","two clusters separation",histDir.c_str());
      if(inputSeparation.size() > 0){
        TH1D *mergedSeparation = new TH1D(*(TH1D*)inputSeparation[0]);
        
        for(int iter=1;iter<inputSeparation.size();iter++){
          mergedSeparation->Add((TH1D*)inputSeparation[iter]);
        }
        
        mergedSeparation->SetLineColor(colors[iDist]);
        mergedSeparation->GetXaxis()->SetTitle("Separation");
        //        mergedSeparation->Rebin(2);
        mergedSeparation->GetYaxis()->SetRangeUser(0,500);
        
        for(int i=1;i<=mergedSeparation->GetNbinsX();i++){mergedSeparation->SetBinContent(i,mergedSeparation->GetBinContent(i)+40*iDist);}
        mergedSeparation->DrawCopy(iDist == 0 ? "" : "same");
        
      }
      
      canvas->cd(3);
      
      vector<TObject*> inputResponse = getHistsWithName("ErecEsimVsEta","ErecEsim vs. eta",histDir.c_str());
      if(inputResponse.size() > 0){
        TH2D *mergedResponse = new TH2D(*(TH2D*)inputResponse[0]);
        
        for(int iter=1;iter<inputResponse.size();iter++){
          mergedResponse->Add((TH2D*)inputResponse[iter]);
        }
        
        mergedResponse->SetLineColor(colors[iDist]);
        mergedResponse->GetXaxis()->SetTitle("|#eta|");
        mergedResponse->GetYaxis()->SetTitle("E_{rec}/E_{sim}");
        
        TH1D *response = mergedResponse->ProjectionY();
//        response->Rebin(2);
        response->SetLineColor(colors[iDist]);
        response->GetXaxis()->SetRangeUser(0,1.4);
        response->GetYaxis()->SetRangeUser(0,350);
        response->DrawCopy(iDist == 0 ? "" : "same");
      }
      
      
    }
  }
  canvas->cd(1);
  leg->Draw();
  canvas->cd(2);
  leg->Draw();
  canvas->cd(3);
  leg->Draw();

}

