//
//  Hexel.hpp
//
//  Created by Jeremi Niedziela on 12/06/2018.
//

#ifndef Hexel_hpp
#define Hexel_hpp

class Hexel{
public:
  Hexel(double _eta=0, double _phi=0, double _x=0, double _y=0, double _z=0, double _weight=0,double _thickness=0,double _time=0,
        int _detid=-1, int _layer=-1, int _clusterRECOindex=-1,
        bool _isHalf = false);
  
  // those can be read from RecHit:
  double eta, phi, x, y, z, weight, thickness, time;
  int layer, detid, clusterRECOIndex;
  bool isHalf;
  
  // those must be calculated and set later
  double rho, fraction, delta, sigmaNoise;
  bool isBorder, isHalo;
  int clusterIndex, nearestHigher;
};

#endif /* Hexel_hpp */
