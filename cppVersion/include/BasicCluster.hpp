//
//  BasicCluster.hpp
//  xHGCalClustering
//
//  Created by Jeremi Niedziela on 12/06/2018.
//  Copyright © 2018 Jeremi Niedziela. All rights reserved.
//

#ifndef BasicCluster_hpp
#define BasicCluster_hpp

#include "Hexel.hpp"

#include <vector>

class BasicCluster{
public:
  BasicCluster(double _energy=-1,
               std::tuple<double,double,double> position=std::make_tuple(0,0,0),
               std::vector<Hexel *> _thisCluster=std::vector<Hexel*>(0),
               int _algoId=-1,int _caloId=-1);
  
  double x,y,z,eta,phi,energy;
  std::vector<Hexel*> thisCluster;
  int algoId, caloId;
  
};


#endif /* BasicCluster_hpp */
