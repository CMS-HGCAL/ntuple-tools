//
//  Hexel.hpp
//
//  Created by Jeremi Niedziela on 12/06/2018.
//

#ifndef Hexel_hpp
#define Hexel_hpp

/// Stores basic information about reconstructed hexel
class Hexel{
public:
  /// Default constructor, takes hexel parameters as an input
  Hexel(double _eta=0, double _phi=0, double _x=0, double _y=0, double _z=0, double _weight=0,double _thickness=0,double _time=0,
        int _detid=-1, int _layer=-1, int _clusterRECOindex=-1,
        bool _isHalf = false);
  ~Hexel();

  // those can be read from RecHit:
  double eta, phi, x, y, z, weight, thickness, time;
  int detid, layer, clusterRECOIndex;
  bool isHalf;

  // those must be calculated and set later
  double rho, fraction, delta, sigmaNoise;
  bool isBorder, isHalo;
  int clusterIndex, nearestHigher;
};

#endif /* Hexel_hpp */
