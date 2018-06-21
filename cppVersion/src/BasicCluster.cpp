//
//  BasicCluster.cpp
//
//  Created by Jeremi Niedziela on 12/06/2018.
//

#include <memory>

#include "BasicCluster.hpp"

#include <TMath.h>

using namespace std;
using namespace TMath;

BasicCluster::BasicCluster(double _energy, double _x, double _y, double _z,
                           vector<shared_ptr<Hexel>> _thisCluster):
energy((_energy>=0) ? _energy : 0.),
x(_x),
y(_y),
z(_z),
thisCluster(_thisCluster)
{
  double theta = acos(z/sqrt(x*x+y*y+z*z));
  eta = -log(tan(theta/2.));
  phi = atan2(y, x);
}

BasicCluster::~BasicCluster()
{
  thisCluster.clear();
}
