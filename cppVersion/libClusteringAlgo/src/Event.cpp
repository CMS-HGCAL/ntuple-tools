//
//  HGCalEvent.cpps
//
//  Created by Jeremi Niedziela on 07/06/2018.
//

#include "Event.hpp"

using namespace std;

Event::Event(TTree* _tree, TTree *_testbeamTracksTree) :
tree(_tree),
testbeamTracksTree(_testbeamTracksTree),
isTestBeam(false),
eventNumber(-1),
runNumber(-1),
pdgID(0),
beamEnergy(0.0),
trueBeamEnergy(0.0)
{
  if(tree->GetBranchStatus("trueBeamEnergy")){ // If we are in the testbeam setup
    isTestBeam = true;
    tree->SetBranchAddress("event",         &eventNumber);
    tree->SetBranchAddress("run",           &runNumber);
    tree->SetBranchAddress("pdgID",         &pdgID);
    tree->SetBranchAddress("beamEnergy",    &beamEnergy);
    tree->SetBranchAddress("trueBeamEnergy",&trueBeamEnergy);
    
    recHits = shared_ptr<RecHits>(new RecHits(tree));
    if(testbeamTracksTree){
      testbeamTrack = shared_ptr<TestbeamTrack>(new TestbeamTrack(testbeamTracksTree));
    }
    else{
      cout<<"WARNING -- no test beam track tree was provided! Track info will not be available."<<endl;
    }
  }
  else{ // If we work with the MC data
    genParticles  = shared_ptr<GenParticles>(new GenParticles(tree));
    recHits       = shared_ptr<RecHits>(new RecHits(tree));
    simClusters   = shared_ptr<SimClusters>(new SimClusters(tree));
    clusters2D    = shared_ptr<Clusters2D>(new Clusters2D(tree));
  }
}

Event::~Event()
{
}

void Event::Print()
{
  cout<<"Event (test beam: "<< (isTestBeam ? "yes" : "no") <<")"<<endl;
  cout<<"\tRun number:"<<runNumber<<"\tevent number:"<<eventNumber<<endl;
  cout<<"\tPDG ID:"<<pdgID<<"\tbeam energy:"<<beamEnergy<<"\ttrue beam energy:"<<trueBeamEnergy<<endl;
}

void Event::GoToEvent(int event)
{
  if(!isTestBeam){
    genParticles->Clean();
    simClusters->Clean();
    clusters2D->Clean();
  }
  
  recHits->Clean();
  tree->GetEntry(event);
  if(!testbeamTracksTree){
//    cout<<"Test beam tree does not exist!"<<endl;
  }
  else if(event > testbeamTracksTree->GetEntries()){
    cout<<"Requested an event beyond test beam tracks tree limits!!"<<endl;
  }
  else{
    testbeamTracksTree->GetEntry(event);
  }
}
