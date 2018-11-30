//
//  BasicCluster.cpp
//
//  Created by Jeremi Niedziela on 12/06/2018.
//

#include <memory>
#include <iostream>

#include "BasicCluster.hpp"

#include <TMath.h>

using namespace std;
using namespace TMath;

BasicCluster::BasicCluster(double _energy, double _x, double _y, double _z,
                           vector<shared_ptr<Hexel>> _thisCluster,
                           vector<vector<shared_ptr<Hexel>>> _thisClusters):
energy((_energy>=0) ? _energy : 0.),
x(_x),
y(_y),
z(_z),
usedIn3Dcluster(0),
thisCluster(_thisCluster),
thisClusters(_thisClusters)
{
  double theta = acos(z/sqrt(x*x+y*y+z*z));
  eta = -log(tan(theta/2.));
  phi = atan2(y, x);
  radius = 0;
}

BasicCluster::BasicCluster(double _energy, double _x, double _y, double _z, double _eta, double _radius) :
energy(_energy),
x(_x),
y(_y),
z(_z),
eta(_eta),
radius(_radius)
{
  thisCluster = std::vector<std::shared_ptr<Hexel>>(0);
  phi = 0;
}

BasicCluster::BasicCluster(const BasicCluster &b) :
energy(b.energy),
x(b.x),
y(b.y),
z(b.z),
eta(b.eta),
phi(b.phi),
radius(b.radius),
usedIn3Dcluster(b.usedIn3Dcluster)
{
  for(auto hex : b.thisCluster){
    thisCluster.push_back(hex);
  }
}

BasicCluster::~BasicCluster()
{
  thisCluster.clear();
}

void BasicCluster::Print()
{
  cout<<"Basic cluster:"<<endl;
  cout<<"\tE:"<<energy<<"\tposition:("<<x<<","<<y<<","<<z<<")"<<endl;
  cout<<"\tradius:"<<radius<<"\teta:"<<eta<<"\tphi:"<<phi<<endl;
}

bool BasicCluster::operator==(BasicCluster &b)
{
  if(   fabs(energy - b.energy) < 0.000001
     && fabs(x - b.x)           < 0.000001
     && fabs(y - b.y)           < 0.000001
     && fabs(z - b.z)           < 0.000001
     && fabs(eta - b.eta)       < 0.000001
     && fabs(phi - b.phi)       < 0.000001)
    return true;
  
  return false;
}
