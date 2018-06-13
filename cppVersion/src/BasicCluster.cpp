//
//  BasicCluster.cpp
//
//  Created by Jeremi Niedziela on 12/06/2018.
//

#include "BasicCluster.hpp"

#include <TMath.h>

using namespace std;
using namespace TMath;

BasicCluster::BasicCluster(double _energy,
                           std::tuple<double,double,double> position,
                           std::vector<Hexel *> _thisCluster,
                           int _algoId,int _caloId):
energy((_energy>=0) ? _energy : 0.),
thisCluster(_thisCluster),
algoId(_algoId),
caloId(_caloId)

{
  eta = 0;
  phi = 0;
  x = 0;
  y = 0;
  z = 0;
  
  if(std::get<0>(position) != 0 && std::get<1>(position) != 0 && std::get<2>(position) != 0){
    x = std::get<0>(position);
    y = std::get<1>(position);
    z = std::get<2>(position);
    
    double theta = acos(z/sqrt(x*x+y*y+z*z));
    eta = -log(tan(theta/2.));
    phi = atan2(y, x);
  }
}
