#include <string>
#include <cstdlib>
#include <iostream>
#include <algorithm>

#include "Event.hpp"
#include "RecHitCalibration.hpp"
#include "ImagineAlgo.hpp"

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>

using namespace std;

//----------------------------------------------------------------------------------------
// HGCal Imaging Algo parameters:
bool dependSensor = true;
double deltac[3] = {2., 2., 5.}; // in cartesian coordiantes in cm, per detector
double multiclusterRadii[3] = {2., 5., 5.}; // in cartesian coordiantes in cm, per detector
int minClusters = 3; // request at least minClusters+1 2D clusters

// cut on energy (also passed to ImagingAlgo):
double energyMin = 3; // relative to the noise

// test only within this layers range:
int minLayer=0;
int maxLayer=40;

// range of ntuples to test (will be appended to the inputPath string below):
int minNtuple = 11;
int maxNtuple = 11;

// base input and output paths:
string inputPath = "../../data/_SingleGammaPt100Eta1p6_2p8_PhaseIITDRFall17DR-noPUFEVT_93X_upgrade2023_realistic_v2-v1_GEN-SIM-RECO/NTUP/_SingleGammaPt100Eta1p6_2p8_PhaseIITDRFall17DR-noPUFEVT_93X_upgrade2023_realistic_v2-v1_GEN-SIM-RECO_NTUP_";

string outDir = "../clusteringResultsCXX/pythonCompare";
//----------------------------------------------------------------------------------------


//# get pandas mask of hits in given layer
//def getLayerMask(hits,layer):
//  return hits["layer"]==layer


vector<Hexel*>* GetHexelsInCluster(vector<Hexel*> &hexels,int clusterIndex)
{
  vector<Hexel*> *hexelsInCluster = new vector<Hexel*>(hexels);
  for(auto hexel : hexels){
    if(hexel->clusterIndex == clusterIndex){
      hexelsInCluster->push_back(hexel);
    }
  }
  return hexelsInCluster;
}

// get pandas mask of hits above noice threshold
RecHits* GetHitsInLayer(RecHits *hits,int layer)
{
  RecHits *hitsInLayer = new RecHits();
  RecHit *hit = nullptr;
  
  for(int iHit=0;iHit<hits->N();iHit++){
    hit = hits->GetHit(iHit);
    if(!hit) continue;
    if(hit->layer == layer){
      hitsInLayer->AddHit(hit);
    }
  }
  return hitsInLayer;
}

RecHits* GetHitsAboveNoice(RecHits *hits,double ecut)
{
  double sigmaNoise = 1.;
  double thickIndex = -1;
  
  RecHits *hitsAboveNoice = new RecHits();
  
  RecHitCalibration *recHitCalib = new RecHitCalibration();
  RecHit *hit;
  
  for(int iHit=0;iHit<hits->N();iHit++){
    hit = hits->GetHit(iHit);
    int layer = hit->layer;
    float thickness = hit->thickness;
    
    if(layer <= 40){  // EE + FH
      if     (thickness > 99.  && thickness < 101.) thickIndex = 0;
      else if(thickness > 199. && thickness < 201.) thickIndex = 1;
      else if(thickness > 299. && thickness < 301.) thickIndex = 2;
      else cout<<"ERROR - silicon thickness has a nonsensical value"<<endl;
    }
    // determine noise for each sensor/subdetector using RecHitCalibration library
    sigmaNoise = 0.001 * recHitCalib->sigmaNoiseMeV(layer, thickIndex);  // returns threshold for EE, FH, BH (in case of BH thickIndex does not play a role)
    if(hit->energy >= ecut * sigmaNoise){
      // checks if rechit energy is above the threshold of ecut (times the sigma noise for the sensor, if that option is set)
      hitsAboveNoice->AddHit(hit);
    }
  }
  return hitsAboveNoice;
}

// groups hits into array of clusters
void getHitsPerCluster(vector<RecHits*> &hitsPerCluster, RecHits *hits,SimClusters *clusters)
{
  RecHits *hitsAboveNoice = GetHitsAboveNoice(hits, energyMin);
  vector<unsigned int> *hitsDetIDs = hitsAboveNoice->detid;
  
  for(int iCluster=0;iCluster<clusters->N();iCluster++){
    vector<unsigned int> hitsInClusterDetIDs = clusters->hits->at(iCluster);
    vector<unsigned int> clusterToHitID;
    
    sort(hitsDetIDs->begin(), hitsDetIDs->end());
    sort(hitsInClusterDetIDs.begin(), hitsInClusterDetIDs.end());
    
    set_intersection(hitsDetIDs->begin(), hitsDetIDs->end(),
                     hitsInClusterDetIDs.begin(), hitsInClusterDetIDs.end(),
                     std::back_inserter(clusterToHitID));
    
    RecHits *hitsInThisCluster = new RecHits();
    
    for(unsigned int i : clusterToHitID){
      ptrdiff_t pos = distance(hits->detid->begin(), find(hits->detid->begin(), hits->detid->end(), i));
      hitsInThisCluster->AddHit(hits->GetHit((int)pos));
    }
    
    hitsPerCluster.push_back(hitsInThisCluster);
  }
}

// groups hits associated with hexels into array of clusters
void getRecHitsPerHexel(vector<RecHits*> &hitsClustered, RecHits *hits,vector<Hexel*> hexels){
  RecHits *hitsAboveNoice = GetHitsAboveNoice(hits, energyMin);
  
  vector<int> clusterIndices;
  vector<unsigned int> *hitDetIDs = hitsAboveNoice->detid; //ok
  vector<int> hexelDetIDs;

  for(Hexel *hexel : hexels){
    if(find(clusterIndices.begin(), clusterIndices.end(), hexel->clusterIndex) == clusterIndices.end()){
      clusterIndices.push_back(hexel->clusterIndex);
    }
    hexelDetIDs.push_back(hexel->detid);//ok
  }
  
  vector<int> hexelToHitID;
  vector<int> new_vec_1 = vector<int>(hexelDetIDs);
  std::sort(std::begin(new_vec_1), std::end(new_vec_1));
  
  for (int i = 0; i < hitDetIDs->size(); ++i) {
    if (std::binary_search(std::begin(new_vec_1),
                           std::end(new_vec_1),
                           hitDetIDs->at(i))) {
      hexelToHitID.push_back(i);
    }
  }
  
  for(int clusterIndex=0;
      clusterIndex<= *max_element(clusterIndices.begin(),clusterIndices.end());
      clusterIndex++){
    hitsClustered.push_back(new RecHits());
  }
  
  for(int iHex=0;iHex<hexels.size();iHex++){
    int hitIndex = hexelToHitID[iHex];
    RecHit *hit = hits->GetHit(hitIndex);
    hitsClustered[hexels[iHex]->clusterIndex]->AddHit(hit);
  }
}
// get clustered hexels by re-running the clustering algorithm

void getRecClustersFromImagingAlgo(vector<Hexel*> &hexelsClustered_rerun, RecHits *hits)
{
  ImagingAlgo *algo = new ImagingAlgo(energyMin, deltac, minClusters, dependSensor, 0);
  
  std::vector<std::vector<std::vector<Hexel*>>> clusters2D_rerun;
  algo->makeClusters(clusters2D_rerun, hits,energyMin); // nested list of "hexels", per layer, per 2D cluster
  
  std::vector<BasicCluster*> clusters2DList_rerun;
  algo->getClusters(clusters2DList_rerun, clusters2D_rerun); // flat list of 2D clusters (as basic clusters)
  
  for(BasicCluster *bClust : clusters2DList_rerun){
    for(Hexel *iNode : bClust->thisCluster){
      if(!iNode->isHalo){
        hexelsClustered_rerun.push_back(iNode);
      }
    }
  }
}

bool pointWithinCircle(double px,double py,double x,double y,double r){
  return pow(r,2) >= pow(px-x,2) + pow(py-y,2);
}


double duration(std::chrono::time_point<std::chrono::system_clock> t0,
                std::chrono::time_point<std::chrono::system_clock> t1)
{
  auto elapsed_secs = t1-t0;
  typedef std::chrono::duration<float> float_seconds;
  auto secs = std::chrono::duration_cast<float_seconds>(elapsed_secs);
  return secs.count();
}

std::chrono::time_point<std::chrono::system_clock> now()
{
  return std::chrono::system_clock::now();
}


int main()
{
  gROOT->ProcessLine(".L loader.C+");
  
  std::system(("mkdir -p "+outDir).c_str());
  
  for(int nTupleIter=minNtuple;nTupleIter<=maxNtuple;nTupleIter++){
    cout<<"\nCurrent ntup: "<<nTupleIter<<endl;

    TFile *inFile = TFile::Open(Form("%s%i.root",inputPath.c_str(),nTupleIter));
    TTree *tree = (TTree*)inFile->Get("ana/hgc");
    long long nEvents = tree->GetEntries();
    
    
    cout<<"\n\nLoading ntuple...";
    auto start = now();
    Event *hgCalEvent = new Event(tree);
    auto end = now();
    cout<<" done ("<<duration(start,end)<<" s)"<<endl;
    
    cout<<"n entries:"<<nEvents<<endl;
    
    // start event loop
    for(int iEvent=0;iEvent<nEvents;iEvent++){
      if(iEvent>100) break;
      auto startEvent = now();
      hgCalEvent->GoToEvent(iEvent);
  
      // check if particles reached EE
      bool skipEvent = false;
      for(auto reachedEE : *(hgCalEvent->genParticles->reachedEE)){
        if(reachedEE==0){
//          cout<<"particle didn't reach EE -- skipping the event!!"<<endl;
          skipEvent = true;
          break;
        }
      }
      if(skipEvent) continue;
      
      string eventDir = outDir+"/ntup"+to_string(nTupleIter)+"/event"+to_string(iEvent);
      std::system(("mkdir -p "+eventDir).c_str());
      
      cout<<"\nCurrent event:"<<iEvent<<endl;
      
      RecHits *recHitsRaw = hgCalEvent->recHits;
      SimClusters *simClusters = hgCalEvent->simClusters;
      
      // get simulated hits associated with a cluster
      cout<<"preparing simulated hits and clusters...";
      start = now();
      vector<RecHits*> simHitsPerClusterArray;
      getHitsPerCluster(simHitsPerClusterArray, recHitsRaw, simClusters);
      end = now();
      cout<<" done ("<<duration(start,end)<<" s)"<<endl;
      
      
      // re-run clustering with HGCalAlgo, save to file
      cout<<"running clustering algorithm...";
      start = now();
      std::vector<Hexel*> recClusters;
      getRecClustersFromImagingAlgo(recClusters, recHitsRaw);
      end = now();
      cout<<" done ("<<duration(start,end)<<" s)"<<endl;
      
      // recClusters -> array of hexel objects
      cout<<"looking for hits associated with hexels...";
      start = now();
      vector<RecHits*> recHitsPerClusterArray;
      getRecHitsPerHexel(recHitsPerClusterArray, recHitsRaw, recClusters);
      end = now();
      cout<<" done ("<<duration(start,end)<<" s)\n"<<endl;
      
      // perform final analysis, fill in histograms and save to files
      
      TH2D *energyComparisonHist = new TH2D("energy comparison","energy comparison",100,0,100,100,0,100);
      TH2D *energyComparisonOverlapHist = new TH2D("energy comparison overlap.","energy comparison overlap.",100,0,100,100,0,100);
      
      cout<<"\n\nGenerating final hists...\n\n";
      start = now();
      
      
      for(int layer=minLayer;layer<maxLayer;layer++){
//        cout<<"layer:"<<layer<<endl;
        for(int recClusterIndex=0;recClusterIndex<recHitsPerClusterArray.size();recClusterIndex++){
          
          RecHits *recCluster = recHitsPerClusterArray[recClusterIndex];
          RecHits *recHitsInLayerInCluster = GetHitsInLayer(recCluster,layer);
          
          if(recHitsInLayerInCluster->N()==0) continue;
          
          double recEnergy = recHitsInLayerInCluster->GetTotalEnergy();
          double xMaxRec   = recHitsInLayerInCluster->GetXmax();
          double xMinRec   = recHitsInLayerInCluster->GetXmin();
          double yMaxRec   = recHitsInLayerInCluster->GetYmax();
          double yMinRec   = recHitsInLayerInCluster->GetYmin();
          
          double recClusterX = xMinRec+(xMaxRec-xMinRec)/2.;
          double recClusterY = yMinRec+(yMaxRec-yMinRec)/2.;
          double recClusterR = max((xMaxRec-xMinRec)/2.,(yMaxRec-yMinRec)/2.);
          
          double assocSimEnergy = 0;
          
          for(int simClusterIndex=0;simClusterIndex<simHitsPerClusterArray.size();simClusterIndex++){
            RecHits *simCluster = simHitsPerClusterArray[simClusterIndex];
            RecHits *simHitsInLayerInCluster = GetHitsInLayer(simCluster,layer);
            
            if(simHitsInLayerInCluster->N()==0) continue;
            
            double simEnergy = simHitsInLayerInCluster->GetTotalEnergy();
            double xMaxSim   = simHitsInLayerInCluster->GetXmax();
            double xMinSim   = simHitsInLayerInCluster->GetXmin();
            double yMaxSim   = simHitsInLayerInCluster->GetYmax();
            double yMinSim   = simHitsInLayerInCluster->GetYmin();
            
            double simClusterX = xMinSim+(xMaxSim-xMinSim)/2.;
            double simClusterY = yMinSim+(yMaxSim-yMinSim)/2.;
            double simClusterR = max((xMaxSim-xMinSim)/2.,(yMaxSim-yMinSim)/2.);
            
            if(recEnergy*simEnergy != 0){
              energyComparisonHist->Fill(recEnergy,simEnergy);
            }
            
            if(pointWithinCircle(simClusterX,simClusterY,recClusterX,recClusterY,recClusterR)){
              assocSimEnergy += simEnergy;
            }
            if(recEnergy*assocSimEnergy != 0){
              energyComparisonOverlapHist->Fill(recEnergy,assocSimEnergy);
            }
          }
        }
      }
      energyComparisonHist->SaveAs(Form("%s/energyComparisonHist.root",eventDir.c_str()));
      energyComparisonOverlapHist->SaveAs(Form("%s/energyComparisonOverlapHist.root",eventDir.c_str()));
      end = now();
      cout<<" done ("<<duration(start,end)<<" s)"<<endl;
      
      auto endEvent = now();
      cout<<"Total event processing time: "<<duration(startEvent,endEvent)<<" s"<<endl;
      
    }
    
  }
    return 0;
}









