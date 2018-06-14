//
//  Hexel.cpp
//
//  Created by Jeremi Niedziela on 12/06/2018.
//

#include "Hexel.hpp"

Hexel::Hexel(double _eta, double _phi, double _x, double _y, double _z, double _weight,double _thickness,double _time,
             int _detid, int _layer, int _clusterRECOindex,
             bool _isHalf) :
eta(_eta),
phi(_phi),
x(_x),
y(_y),
z(_z),
weight(_weight),
thickness(_thickness),
time(_time),
detid(_detid),
layer(_layer),
clusterRECOIndex(_clusterRECOindex),
isHalf(_isHalf)
{
  rho = 0;
  fraction = 1;
  delta = 0;
  sigmaNoise = 0.;
  isBorder = false;
  isHalo = false;
  nearestHigher = -1;
  clusterIndex = -1;
}

Hexel::~Hexel()
{
  
}
