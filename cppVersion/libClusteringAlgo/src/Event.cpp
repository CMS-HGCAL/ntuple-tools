//
//  HGCalEvent.cpps
//
//  Created by Jeremi Niedziela on 07/06/2018.
//

#include "Event.hpp"

using namespace std;

Event::Event(TTree* _tree) :
tree(_tree),
isTestBeam(false),
eventNumber(-1),
runNumber(-1),
pdgID(0),
beamEnergy(0.0),
trueBeamEnergy(0.0)
{
  if(tree->GetBranchStatus("trueBeamEnergy")){ // If we are in the testbeam setup
    isTestBeam = true;
    tree->SetBranchAddress("event", &eventNumber);
    tree->SetBranchAddress("run", &runNumber);
    tree->SetBranchAddress("pdgID", &pdgID);
    tree->SetBranchAddress("beamEnergy", &beamEnergy);
    tree->SetBranchAddress("trueBeamEnergy", &trueBeamEnergy);
    
    recHits = shared_ptr<RecHits>(new RecHits(_tree));
  }
  else{ // If we work with the MC data
    genParticles = shared_ptr<GenParticles>(new GenParticles(_tree));
    recHits = shared_ptr<RecHits>(new RecHits(_tree));
    simClusters = shared_ptr<SimClusters>(new SimClusters(_tree));
    clusters2D = shared_ptr<Clusters2D>(new Clusters2D(_tree));
  }
}

Event::~Event()
{
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
}
