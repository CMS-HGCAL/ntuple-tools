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
string baseDir = "clusteringResultsCXX/geneticOptimizerBest";
//string baseDir = "clusteringResultsCXX/singleGamma";
//string baseDir = "clusteringResultsCXX/pedja_photon_Pt6_dR0p10";
//string baseDir = "clusteringResultsCXX/pedja_photon_Pt10_dR0p02";

vector<vector<string>> fileNames = {
  {"energyComparisonNoMatchingHist.root","no matching"},
  {"energyComparisonClosestHist.root","closest rec cluster"},
  {"ErecEsimVsEta.root","ErecEsim vs. eta"},
  {"NrecNsim.root","NrecNsim"},
  {"simgaEvsEtaEsim.root","sigma(E)Esim vs. eta"},
  {"containment.root","containment"},
  {"deltaN.root","numberClusters"},
  {"resolution.root","resolution"},
  {"separation.root","separation"}
};

void analyzeClustering(string inputPath="")
{
  if(inputPath == "") inputPath = baseDir;
  
  TF1 *fun = new TF1("fun","x",0,100);
  TCanvas *canvas = new TCanvas("canvas","canvas",2880,1800);
  canvas->Divide(3,3);
  
  int iter=1;
  for(auto names : fileNames){
    TFile *file = TFile::Open((inputPath+"/"+names[0]).c_str());
    TH1D *hist = (TH1D*)file->Get(names[1].c_str());
    
    canvas->cd(iter);
    hist->Draw("colz");
    
    if(iter < 3) fun->Draw("same");
    iter++;
  }
}

