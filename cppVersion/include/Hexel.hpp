//
//  Hexel.hpp
//
//  Created by Jeremi Niedziela on 12/06/2018.
//

#ifndef Hexel_hpp
#define Hexel_hpp

#include "RecHits.hpp"

class Hexel{
public:
  Hexel(RecHit *hit=nullptr,double _sigmaNoise=-1);
  
  bool __gt__(double _rho){return rho > _rho;}
  
  double eta, phi, x, y, z;
  double time, rho;
  double weight, fraction;
  double delta, sigmaNoise, thickness;
  bool isHalfCell, isBorder, isHalo;
  int layer, detid, clusterIndex, clusterRECOIndex, nearestHigher;
  
};

#endif /* Hexel_hpp */
