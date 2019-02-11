// include particles decaying before EE:
//string baseDir = "clusteringResults/noReachedEE_centerWithinRecCluster/";

// play with overlapping circles:
//string baseDir = "clusteringResults/reachedEE_centersOverlap/";
//string baseDir = "clusteringResults/reachedEE_centerWithinRecCluster_inverted/";
//string baseDir = "clusteringResults/reachedEE_centerWithin_80pc_RecCluster_inverted/";

// from C++ version
//string baseDir = "clusteringResultsCXX/twoPions_Pt80_Eta2_DeltaR0p4/";
//string baseDir = "clusteringResultsCXX/twoPions_jniedzie";
//string baseDir = "clusteringResultsCXX/geneticOptimizer";
//string baseDir = "clusteringResultsCXX/geneticOptimizerBest";
//string baseDir = "clusteringResultsCXX/singleGamma";
//string baseDir = "clusteringResultsCXX/pedja_photon_Pt6_dR0p10";
string baseDir = "clusteringResultsCXX/qcd";
//string baseDir = "clusteringResultsCXX/pedja_photon_Pt10_dR0p02";

//string baseDir = "clusteringResultsCXX/geneticOptimizerTwoPionDeltaR0p4/";


//string baseDir = "clusteringResultsCXX/geneticOptimizerTwoPhotons/";
//string baseDir = "clusteringResultsCXX/clemens_twoPhoton_fixed/";

vector<string> fileNames = {
  "ErecVsEsimUnmatched",
  "ErecVsEsimDetIdMatching",
  "resolutionVsEta",
  "NrecVsNsim",
  "resolution",
  "separation",
  "containment",
  "deltaNclusters",
  "fakeVsLayer"
};

void analyzeClustering(string inputPath="")
{
  if(inputPath == "") inputPath = baseDir;
  
  TF1 *fun = new TF1("fun","x",0,100);
  TCanvas *canvas = new TCanvas("canvas","canvas",2880,1800);
  canvas->Divide(3,3);
  
  int iter=1;
  for(auto names : fileNames){
    TFile *file = TFile::Open((inputPath+"/"+names+".root").c_str());
    if(!file){
      iter++;
      continue;
    }
    TH1D *hist = (TH1D*)file->Get(names.c_str());
    
    canvas->cd(iter);
    if(iter>4) hist->SetFillColorAlpha(kBlue,0.3);
    
    hist->Draw("colz");
    
    if(iter < 3) fun->Draw("same");
    iter++;
  }
}

