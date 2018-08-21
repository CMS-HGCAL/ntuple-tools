//
//  RecHitCalibration.hpp
//
//  Created by Jeremi Niedziela on 07/06/2018.
//
//  Copy of the python script HGCalImagineAlgo.py

#ifndef HGCalImagineAlgo_hpp
#define HGCalImagineAlgo_hpp

#include "RecHitCalibration.hpp"
#include "RecHits.hpp"
#include "Hexel.hpp"
#include "BasicCluster.hpp"
#include "ConfigurationManager.hpp"

#include <TF1.h>

#include <iostream>
#include <vector>
#include <numeric>
#include <functional>

///  ImagingAlgo performs 2D clusterization of hits.
///
///  Implementation of (stand-alone) functionalities of HGCalImagingAlgo
///  based on their CMSSW implementations mainly in RecoLocalCalo/HGCalRecAlgos.

class ImagingAlgo {
public:
  /// Default constructor
  ImagingAlgo();
  
  /// Constructor with config path. The global configuration will be ignored and values will be read directly from this config
  ImagingAlgo(std::string _configPath);
  ~ImagingAlgo();
  
  /// Get clustered hexels by re-running the clustering algorithm
  /// \param hexelsClustered Will be filled with non-halo 2D hexels containing info about cluster index and layer
  /// \param hits Rec hits to be clusterized
  void getRecClusters(std::vector<std::shared_ptr<Hexel>> &hexelsClustered, std::shared_ptr<RecHits> &hits);
  
private:
  std::string configPath;
  
  bool dependSensor;     ///< Should algo depend on the sensor type
  double kappa;          ///< ??
  double ecut;           ///< Minimum hit energy
//  int minClusters;       ///< Request at least minClusters+1 2D clusters
  int verbosityLevel;    ///< Current verbosity level
  double criticalDistanceEE;
  double criticalDistanceFH;
  double criticalDistanceBH;
  double deltacEE;
  double deltacFH;
  double deltacBH;
  
  TF1 *energyDensityFunction; ///< Function that will be used to determine energy density for each hit
  
  std::unique_ptr<RecHitCalibration> recHitCalib; ///< Contains calibration of rec hits
  ConfigurationManager *config;                   ///< Manager keeping current configuration
  
  /// Calculates distance to the nearest hit with higher density
  /// \param nodes Nodes will be filled with the information about nearest higher-density hit
  void calculateDistanceToHigher(std::vector<std::unique_ptr<Hexel>> &nodes);
  
  /// Find cluster centers that satisfy delta & maxdensity/kappa criteria, and assign coresponding hexels
  /// \param clusters Will be filled with hexels grouped by clusters
  /// \param nodes Input hexels (in given layer) to be assigned to clusters
  /// \param points_0 X coordinates of hexels all hexels in this layer
  /// \param points_1 Y coordinates of hexels all hexels in this layer
  /// \param maxDensity Maximum local density
  /// \param layer Current layer
  void findAndAssignClusters(std::vector<std::vector<std::unique_ptr<Hexel>>> &clusters,
                             std::vector<std::unique_ptr<Hexel>> &nodes,
                             std::vector<double> points_0,std::vector<double> points_1,
                             double maxDensity,int layer);
  
  /// Make list of Hexels out of rechits
  /// \param points Will be filled with hexels grouped by layer
  /// \param hits Input hits to be grouped by layer
  void populate(std::vector<std::vector<std::unique_ptr<Hexel>>> &points, std::shared_ptr<RecHits> &hits);
  
  /// Make 2D clusters out of recHits (need to introduce class with input params: delta_c, kappa, ecut, ...)
  /// \param clusters Will be filled with hits grouped by layer by cluster
  /// \param hits Input hits to be clusterized
  void makeClusters(std::vector<std::vector<std::vector<std::unique_ptr<Hexel>>>> &clusters, std::shared_ptr<RecHits> &hits);
  
  /// Get flat list of BasicClusters from the list of hexels grouped by layer by cluster
  /// \param clustersFlat Will be filled with BasicClusters
  /// \param clusters Input hexels grouped by layer by cluster to be flatten
  void getBasicClusters(std::vector<std::unique_ptr<BasicCluster>> &clustersFlat,
                   std::vector<std::vector<std::vector<std::unique_ptr<Hexel>>>> &clusters);

  /// Calculates XYZ position of a cluster
  /// \param cluster Input cluster of hexels
  /// \return Returns a tuple containing X,Y,Z coordinates of the cluster
  std::tuple<double,double,double>  calculatePosition(std::vector<std::unique_ptr<Hexel>> &cluster);
  
  /// Sorts hexels by rho
  /// \param v Input vector of hexels to sort (will not be modified)
  /// \return Returns vector of sorted indices
  std::vector<int> sortIndicesRhoInverted(const std::vector<std::unique_ptr<Hexel>> &v);
  
  /// Sorts hexels by delta
  /// \param v Input vector of hexels to sort (will not be modified)
  /// \return Returns vector of sorted indices
  std::vector<int> sortIndicesDeltaInverted(const std::vector<std::unique_ptr<Hexel>> &v);
  
  /// Calculate max local density in a 2D plane of hexels
  /// \param hexels A vector of hexels for which max local density will be calculated
  /// \param lpX X coordinates of all hexels in this layer
  /// \param lpY Y coordinates of all hexels in this layer
  /// \param layer Current layer
  /// \return Returns max local density
  double calculateLocalDensity(std::vector<std::unique_ptr<Hexel>> &hexels,
                               std::vector<double> lpX, std::vector<double> lpY, int layer);
  
  
};


#endif
