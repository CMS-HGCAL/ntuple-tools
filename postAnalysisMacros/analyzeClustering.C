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
#include <TObject.h>

#include <iostream>

using namespace std;

const int nLayers = 40;
const int nClusters = 40;

// include particles decaying before EE:
//string baseDir = "clusteringResults/noReachedEE_centerWithinRecCluster/";

// play with rec cluster radius:
//string baseDir = "clusteringResults/reachedEE_centerWithin_80pc_RecCluster/";
//string baseDir = "clusteringResults/reachedEE_centerWithin_90pc_RecCluster/";
//string baseDir = "clusteringResults/reachedEE_centerWithinRecCluster/";
//string baseDir = "clusteringResults/reachedEE_centerWithin_120pc_RecCluster/";
//string baseDir = "clusteringResults/reachedEE_centerWithin_130pc_RecCluster/";

// play with energy threshold:
//string baseDir = "clusteringResults/reachedEE_centerWithinRecCluster_noice5s/";
//string baseDir = "clusteringResults/reachedEE_centerWithinRecCluster_noice1s/";

// play with overlapping circles:
//string baseDir = "clusteringResults/reachedEE_centersOverlap/";


//string baseDir = "clusteringResults/reachedEE_centerWithinRecCluster_inverted/";
//string baseDir = "clusteringResults/reachedEE_centerWithin_80pc_RecCluster_inverted/";

// from C++ version
//string baseDir = "clusteringResultsCXX";
//string baseDir = "clusteringResultsCXX/fixedSamples";
//string baseDir = "clusteringResultsCXX/oldSamples";
//string baseDir = "clusteringResultsCXX/twoPions_Pt80_Eta2_DeltaR0p4/";
string baseDir = "clusteringResultsCXX/twoPions_jniedzie";
//string baseDir = "clusteringResultsCXX/geneticOptimizer";
//string baseDir = "clusteringResultsCXX/singleGamma";

vector<TObject*> getHistsWithName(const char* histFileName, const char* histName)
{
  vector<TObject*> result;
  

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
          TObject *hist = inFile->Get(histName);
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

void analyzeClustering()
{
  TCanvas *canvas = new TCanvas("canvas","canvas",2880,1800);
  canvas->Divide(3,3);
  TF1 *fun = new TF1("fun","x",0,100);
  
  
  vector<TObject*> inputEnergyComparisonHists = getHistsWithName("energyComparisonNoMatchingHist","no matching");
  
  if(inputEnergyComparisonHists.size() > 0){
    TH2D *mergedEnergyComparisonHist = new TH2D(*(TH2D*)inputEnergyComparisonHists[0]);

    for(int iter=1;iter<inputEnergyComparisonHists.size();iter++){
      mergedEnergyComparisonHist->Add((TH2D*)inputEnergyComparisonHists[iter]);
    }
    
    canvas->cd(1);
    mergedEnergyComparisonHist->Draw("colz");
    mergedEnergyComparisonHist->SetTitle("No matching");
    
    mergedEnergyComparisonHist->GetXaxis()->SetTitle("E_{rec} (GeV)");
    mergedEnergyComparisonHist->GetYaxis()->SetTitle("E_{sim} (Gev)");
    mergedEnergyComparisonHist->GetZaxis()->SetRangeUser(0,500);
    
    mergedEnergyComparisonHist->GetXaxis()->SetTitleSize(0.06);
    mergedEnergyComparisonHist->GetXaxis()->SetTitleOffset(0.75);
    mergedEnergyComparisonHist->GetYaxis()->SetTitleSize(0.06);
    mergedEnergyComparisonHist->GetYaxis()->SetTitleOffset(0.75);
    
    fun->Draw("same");
  }

  vector<TObject*> inputEnergyComparisonOverlapHists = getHistsWithName("energyComparisonClosestHist","closest rec cluster");
  
  if(inputEnergyComparisonOverlapHists.size() >0){
    TH2D *mergedEnergyComparisonOverlapHist = new TH2D(*(TH2D*)inputEnergyComparisonOverlapHists[0]);
    
    for(int iter=1;iter<inputEnergyComparisonOverlapHists.size();iter++){
      mergedEnergyComparisonOverlapHist->Add((TH2D*)inputEnergyComparisonOverlapHists[iter]);
    }
    canvas->cd(2);
    mergedEnergyComparisonOverlapHist->Draw("colz");
    mergedEnergyComparisonOverlapHist->SetTitle("Matching by detID");
    
    mergedEnergyComparisonOverlapHist->GetXaxis()->SetTitle("E_{rec} (GeV)");
    mergedEnergyComparisonOverlapHist->GetYaxis()->SetTitle("E_{sim} (GeV)");
    
    mergedEnergyComparisonOverlapHist->GetXaxis()->SetTitleSize(0.06);
    mergedEnergyComparisonOverlapHist->GetXaxis()->SetTitleOffset(0.75);
    mergedEnergyComparisonOverlapHist->GetYaxis()->SetTitleSize(0.06);
    mergedEnergyComparisonOverlapHist->GetYaxis()->SetTitleOffset(0.75);
    
    mergedEnergyComparisonOverlapHist->GetZaxis()->SetRangeUser(0,30);
    
    fun->Draw("same");
  }
  
  vector<TObject*> inputErecEsimVsEta = getHistsWithName("ErecEsimVsEta","ErecEsim vs. eta");
  if(inputErecEsimVsEta.size() > 0){
    
    TH2D *mergedErecEsimVsEtaHists = new TH2D(*(TH2D*)inputErecEsimVsEta[0]);
    
    for(int iter=1;iter<inputErecEsimVsEta.size();iter++){
      mergedErecEsimVsEtaHists->Add((TH2D*)inputErecEsimVsEta[iter]);
    }
    
    canvas->cd(3);
    mergedErecEsimVsEtaHists->Draw("colz");
    
    mergedErecEsimVsEtaHists->GetXaxis()->SetTitle("|#eta|");
    mergedErecEsimVsEtaHists->GetYaxis()->SetTitle("E_{rec}/E_{sim}");
    
    canvas->cd(4);
    TH1D *sigmaEvsEsim = mergedErecEsimVsEtaHists->ProjectionY();
    sigmaEvsEsim->Draw();
    sigmaEvsEsim->GetXaxis()->SetRangeUser(0,2);
  }
  
  vector<TObject*> inputsigmaEsimVsEta = getHistsWithName("simgaEvsEta","sigma(E) vs. eta");
  if(inputsigmaEsimVsEta.size() > 0){
    
    TH2D *mergedSigmaEsimVsEtaHists = new TH2D(*(TH2D*)inputsigmaEsimVsEta[0]);
    
    for(int iter=1;iter<inputsigmaEsimVsEta.size();iter++){
      mergedSigmaEsimVsEtaHists->Add((TH2D*)inputsigmaEsimVsEta[iter]);
    }
    
    canvas->cd(5);
    mergedSigmaEsimVsEtaHists->Draw("colz");
    mergedSigmaEsimVsEtaHists->GetXaxis()->SetTitle("|#eta|");
    mergedSigmaEsimVsEtaHists->GetYaxis()->SetTitle("(E_{rec}-E_{sim})/E_{rec}");
  }
  
  vector<TObject*> inputsigmaEsimVsEtaEsim = getHistsWithName("simgaEvsEtaEsim","sigma(E)Esim vs. eta");
  if(inputsigmaEsimVsEtaEsim.size() > 0){
    
    TH2D *mergedSigmaEsimVsEtaEsimHists = new TH2D(*(TH2D*)inputsigmaEsimVsEtaEsim[0]);
    
    for(int iter=1;iter<inputsigmaEsimVsEtaEsim.size();iter++){
      mergedSigmaEsimVsEtaEsimHists->Add((TH2D*)inputsigmaEsimVsEtaEsim[iter]);
    }
    
    canvas->cd(6);
//    mergedSigmaEsimVsEtaEsimHists->Rebin2D(2,2);
    mergedSigmaEsimVsEtaEsimHists->Draw("colz");
    mergedSigmaEsimVsEtaEsimHists->GetXaxis()->SetTitle("|#eta|");
    mergedSigmaEsimVsEtaEsimHists->GetYaxis()->SetTitle("(E_{rec}-E_{sim})/E_{sim}");
//    mergedSigmaEsimVsEtaEsimHists->GetYaxis()->SetRangeUser(-1.0,0.2);
    
    canvas->cd(7);
    TH1D *sigmaEvsEsim = mergedSigmaEsimVsEtaEsimHists->ProjectionY();
    sigmaEvsEsim->Rebin(2);
    sigmaEvsEsim->Draw();
    sigmaEvsEsim->GetXaxis()->SetRangeUser(-1,1);
  }
  
  vector<TObject*> inputSeparation = getHistsWithName("twoSeparation","two clusters separation");
  if(inputSeparation.size() > 0){
    
    TH1D *mergedSeparationHists = new TH1D(*(TH1D*)inputSeparation[0]);
    
    for(int iter=1;iter<inputSeparation.size();iter++){
      mergedSeparationHists->Add((TH1D*)inputSeparation[iter]);
    }
    
    canvas->cd(8);
    mergedSeparationHists->Draw();
    mergedSeparationHists->GetXaxis()->SetTitle("sep");
  }
  
  vector<TObject*> inputNrecNsim = getHistsWithName("NrecNsim","NrecNsim");
  if(inputNrecNsim.size() > 0){
    
    TH2D *mergedNrecNsim = new TH2D(*(TH2D*)inputNrecNsim[0]);
    
    for(int iter=1;iter<inputNrecNsim.size();iter++){
      mergedNrecNsim->Add((TH2D*)inputNrecNsim[iter]);
    }
    
    canvas->cd(9);
    mergedNrecNsim->Draw("colz");
    mergedNrecNsim->GetXaxis()->SetTitle("Nrec");
    mergedNrecNsim->GetYaxis()->SetTitle("Nsim");
    
    
  }
  
  canvas->SaveAs("/Users/Jeremi/Desktop/erecesim.png");
  
}

