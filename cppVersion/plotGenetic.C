#include <dirent.h>

string basePath = "geneticResults/score_fun_v2_pions/";

void PlotInCanvas(TCanvas *canvas, TCanvas *canvasWgt, TCanvas *canvas1D, int pad, TFile *file, string name)
{
  canvas->cd(pad);
  TH2D *hist = (TH2D*)file->Get(name.c_str());
  if(!hist) return;
  hist->Draw("colz");
  canvasWgt->cd(pad);
  TH2D *histWgt = (TH2D*)file->Get((name+"wgt").c_str());
  if(!histWgt) return;
  histWgt->Divide(hist);
  histWgt->Draw("colz");
  
  canvas1D->cd(pad);
  
//  gPad->SetLogy();
  
  TH1D *histFailed = (TH1D*)file->Get((name+"_failed").c_str());
  TH1D *histPassed = (TH1D*)file->Get((name+"_passed").c_str());
  if(!histFailed || !histPassed) return;
  histFailed->SetLineColor(kRed);
//  histFailed->Rebin(4);
//  histPassed->Rebin(4);
  
//  histFailed->Draw();
//  histPassed->Draw("same");
  
  histFailed->DrawNormalized();
  histPassed->DrawNormalized("same");
  
}

int GetCurrentMaxResultsIndex()
{
  DIR *dir = opendir(basePath.c_str());
  
  if(!dir){
    cout<<"could not open directory:"<<basePath<<endl;
    return -1;
  }
  
  struct dirent *ent;
  int maxIndex = -1;
  
  while( (ent = readdir(dir)) ){
    string fileName = ent->d_name;
    
    size_t pos = fileName.find("results_");
    if(pos != std::string::npos){
      string indexString = fileName.substr(fileName.find("results_") + 8);
      int index = stoi(indexString);
      if(index > maxIndex) maxIndex=index;
    }
  }
  closedir(dir);
  
  if(ent){ delete ent; ent=0;}
  
  return maxIndex;
}

void plotGenetic(string path = "")
{
  if(path==""){
    int index = GetCurrentMaxResultsIndex();
    path = basePath+"results_"+to_string(index)+"/geneticHists.root";
    cout<<"Automatically set path to:"<<path<<endl;
  }
  
  TFile *inFile = TFile::Open(path.c_str());
  
  TCanvas *canvas = new TCanvas("genetic","genetic",2880,1800);
  TCanvas *canvasWgt = new TCanvas("genetic wgt","genetic wgt",2880,1800);
  TCanvas *canvas1D = new TCanvas("genetic 1D","genetic 1D",2880,1800);
  
  canvas->Divide(3,4);
  canvasWgt->Divide(3,4);
  canvas1D->Divide(3,4);

  TGraph *graph = (TGraph*)inFile->Get("scores");
  if(graph){
    graph->SetMarkerSize(1.0);
    graph->SetMarkerStyle(20);
    graph->SetMarkerColor(kGreen+2);
    
    canvas->cd(1);
    graph->Draw("AP");
    canvasWgt->cd(1);
    graph->Draw("AP");
    canvas1D->cd(1);
    graph->Draw("AP");
  }
  
  PlotInCanvas(canvas,canvasWgt, canvas1D,2,inFile,"critDistEE");
  PlotInCanvas(canvas,canvasWgt, canvas1D,3,inFile,"critDistFH");
  PlotInCanvas(canvas,canvasWgt, canvas1D,4,inFile,"critDistBH");
  PlotInCanvas(canvas,canvasWgt, canvas1D,5,inFile,"kernel");
  PlotInCanvas(canvas,canvasWgt, canvas1D,6,inFile,"deltaEE");
  PlotInCanvas(canvas,canvasWgt, canvas1D,7,inFile,"deltaFH");
  PlotInCanvas(canvas,canvasWgt, canvas1D,8,inFile,"deltaBH");
  PlotInCanvas(canvas,canvasWgt, canvas1D,9,inFile,"kappa");
  PlotInCanvas(canvas,canvasWgt, canvas1D,10,inFile,"eMin");
  PlotInCanvas(canvas,canvasWgt, canvas1D,11,inFile,"matchingDist");
}
