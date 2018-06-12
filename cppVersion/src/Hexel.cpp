//
//  Hexel.cpp
//
//  Created by Jeremi Niedziela on 12/06/2018.
//

#include "Hexel.hpp"

Hexel::Hexel(RecHit *hit,double _sigmaNoise)
{
  eta = 0;
  phi = 0;
  x = 0;
  y = 0;
  z = 0;
  time = -1;
  isHalfCell = false;
  weight = 0;
  fraction = 1;
  detid = -1;
  rho = 0;
  delta = 0;
  nearestHigher = -1;
  isBorder = false;
  isHalo = false;
  clusterIndex = -1;
  clusterRECOIndex = -1;
  sigmaNoise = 0.;
  thickness = 0.;
  
  if(hit){
    eta = hit->eta;
    phi = hit->phi;
    x = hit->x;
    y = hit->y;
    z = hit->z;
    weight = hit->energy;
    detid = hit->detid;
    layer = hit->layer;
    isHalfCell = hit->isHalf;
    thickness = hit->thickness;
    time = hit->time;
    clusterRECOIndex = hit->cluster2d;
  }
  if(_sigmaNoise>=0) sigmaNoise = _sigmaNoise;
}
