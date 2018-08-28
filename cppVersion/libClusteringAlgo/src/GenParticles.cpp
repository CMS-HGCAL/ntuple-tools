//
//  HGGenParticle.cpp
//
//  Created by Jeremi Niedziela on 07/06/2018.
//

#include "GenParticles.hpp"

#include <iostream>

using namespace std;

GenParticles::GenParticles(TTree *_tree) :
reachedEE(nullptr),
pid(nullptr)
{
  _tree->SetBranchAddress("genpart_reachedEE",&reachedEE);
  _tree->SetBranchAddress("genpart_pid",&pid);
  //... more branches can be added if needed
}

GenParticles::~GenParticles()
{
  if(reachedEE){ reachedEE->clear(); delete reachedEE;}
  if(pid){ pid->clear(); delete pid;}
}

void GenParticles::Print(int index)
{
  cout<<"Gen particle "<<index<<":";
  cout<<"\t reached EE:"<<(reachedEE->at(index) ? "yes" : "no");
  cout<<"\t PDG PID:"<<pid->at(index)<<endl;
}

void GenParticles::Clean()
{
  if(reachedEE){reachedEE->clear(); delete reachedEE; reachedEE = nullptr;}
  if(pid){pid->clear(); delete pid; pid = nullptr;}
}
