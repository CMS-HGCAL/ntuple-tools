//
//  RecHitCalibration.hpp
//
//  Created by Jeremi Niedziela on 07/06/2018.
//
//  Copy of the python script HGCalImagineAlgo.py:
//
//  Implementation of (stand-alone) functionalities of HGCalImagingAlgo,
//  HGCal3DClustering, and HGCalDepthPreClusterer based on
//  their CMSSW implementations mainly in RecoLocalCalo/HGCalRecAlgos
//

#ifndef HGCalImagineAlgo_hpp
#define HGCalImagineAlgo_hpp

#include "RecHitCalibration.hpp"
#include "RecHits.hpp"
#include "Hexel.hpp"
#include "BasicCluster.hpp"

#include <iostream>
#include <vector>
#include <numeric>
#include <functional>

// definition of the HGCalImagingAlgo class's methods & variables
class ImagingAlgo {
public:
  ImagingAlgo(double _ecut=-1,double _deltac[3]=0, int _minClusters=-1, int _dependSensor=false,int _verbosityLevel=0);
  
  // calculate distance to the nearest hit with higher density (still does not use KDTree)
  void calculateDistanceToHigher(std::vector<Hexel*> &nodes);
 
  // find cluster centers that satisfy delta & maxdensity/kappa criteria, and assign coresponding hexels
  void findAndAssignClusters(std::vector<std::vector<Hexel*>> &current_clusters,
                             std::vector<Hexel*> &nodes,
                             std::vector<double> points_0,std::vector<double> points_1,
                             double maxdensity,int layer);
  
  // make list of Hexels out of rechits
  void populate(std::vector<std::vector<Hexel*>> &points, RecHits *hits,double _ecut=-1);
  
  // make 2D clusters out of rechits (need to introduce class with input params: delta_c, kappa, ecut, ...)
  void makeClusters(std::vector<std::vector<std::vector<Hexel*>>> &clusters, RecHits *hits,double _ecut=-1);
  
    // get basic clusters from the list of 2D clusters
  void getClusters(std::vector<BasicCluster*> &clusters_v,
                   std::vector<std::vector<std::vector<Hexel*>>> &clusters);
  
  std::tuple<double,double,double>  calculatePosition(std::vector<Hexel*> &cluster);
  std::tuple<bool,double>           recHitAboveThreshold(RecHit *hit,double _ecut,bool dependSensor=true);
  
  std::vector<int> sortIndicesRhoInverted(const std::vector<Hexel*> &v);
  std::vector<int> sortIndicesDeltaInverted(const std::vector<Hexel*> &v);
  std::vector<int> query_ball_point(std::vector<double> lpX, std::vector<double> lpY, double x, double y, double r);
  
  // distance squared (in x/y) between the two objects (hexels, clusters)
  double distanceReal2(double x1, double y1, double x2, double y2);
  
  // calculate max local density in a 2D plane of hexels
  double calculateLocalDensity(std::vector<Hexel*> &nd,std::vector<double> lpX, std::vector<double> lpY, int layer);
  
  inline void SetVerbosityLevel(int level){verbosityLevel=level;}
  
private:
  // detector layers to consider
  int lastLayerEE = 28;  // last layer of EE
  int lastLayerFH = 40;  // last layer of FH
  int maxlayer = 52;  // last layer of BH
  
  bool dependSensor;
  double deltac[3];
  double kappa;
  double ecut;
  int minClusters;
  int verbosityLevel;
  
  RecHitCalibration *recHitCalib;
};


#endif
