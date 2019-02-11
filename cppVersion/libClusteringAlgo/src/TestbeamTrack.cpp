//
//  TestbeamTrack.cpp
//
//  Created by Jeremi Niedziela on 30/11/2018.
//

#include "TestbeamTrack.hpp"

#include <iostream>

using namespace std;

TestbeamTrack::TestbeamTrack(TTree *_tree) :
slopeX(0),
slopeY(0),
offsetX(0),
offsetY(0)
{
  for(int iLayer=0;iLayer<nLayers;iLayer++){
    posX.push_back(nullptr);
    posY.push_back(nullptr);
    
    _tree->SetBranchAddress(("impactX_HGCal_layer_"+to_string(iLayer+1)).c_str(),posX[iLayer]);
    _tree->SetBranchAddress(("impactY_HGCal_layer_"+to_string(iLayer+1)).c_str(),posY[iLayer]);
  }
  _tree->SetBranchAddress("m_x",&slopeX);
  _tree->SetBranchAddress("m_y",&slopeY);
  _tree->SetBranchAddress("b_x",&offsetX);
  _tree->SetBranchAddress("b_y",&offsetY);
}

TestbeamTrack::~TestbeamTrack()
{
}

void TestbeamTrack::Print()
{
  cout<<"Test beam track"<<endl;
  cout<<"\t slope x:"<<slopeX<<"\ty:"<<slopeY<<endl;
  cout<<"\t offset x:"<<offsetX<<"\ty:"<<offsetY<<endl;
}
