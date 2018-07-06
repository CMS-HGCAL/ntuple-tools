//
//  RecHits.hpp
//
//  Created by Jeremi Niedziela on 07/06/2018.
//

#ifndef RecHits_hpp
#define RecHits_hpp

#include "RecHitCalibration.hpp"
#include "SimClusters.hpp"
#include "Hexel.hpp"
#include "ConfigurationManager.hpp"

#include <TTree.h>
#include <vector>

class RecHit;

/// RecHits is a collection of RecHit objects
class RecHits {
public:
  /// Default constructor, returned object will contain no hits
  RecHits();
  /// Constructor with TTree as an input, returned object will contain hits representation of rechit_* branches of TTree
  RecHits(TTree *_tree);
  ~RecHits();
  
  /// Returns number of stored RecHit objects
  inline int N(){return (int)eta->size();}
  
  /// Returns sum of the energies of all hits in the collection
  double GetTotalEnergy();
  
  /// Returns the smallest X among all hits stored in the collection
  double GetXmin();
  
  /// Returns the greatest X among all hits stored in the collection
  double GetXmax();
  
  /// Returns the smallest Y among all hits stored in the collection
  double GetYmin();
  
  /// Returns the greatest Y among all hits stored in the collection
  double GetYmax();
  
  /// Returns the average of min and max pseudorapidity among all hits stored in the collection
  double GetCenterEta();
  
  /// Returns a single RecHit
  /// \param index Index of the desired RecHit
  /// \return Pointer to the RecHit object
  std::unique_ptr<RecHit> GetHit(int index);
  
  /// Adds a single hit to the colleciton
  /// \param hit Pointer to RecHit object to be added to the collection
  void AddHit(std::unique_ptr<RecHit> &hit);
  
  /// Get subset of RecHits that are above the noise threshold
  /// \param energyMin Noise threshold - minimum energy of hit to be included
  /// \return Pointer to a collection of hits above specified noise threshold
  std::unique_ptr<RecHits> GetHitsAboveNoise(double energyMin);
  
  /// Groups hits in clusters
  /// \param hitsPerCluster A vector that will be filled with RecHits collections, one for each cluster
  /// \param clusters Collection of simulated clusters to which hits will be assigned
  /// \param energyMin Noise threshold - minimum energy of hit to be included
  void GetHitsPerSimCluster(std::vector<RecHits*> &hitsPerCluster,
                            std::shared_ptr<SimClusters> clusters, double energyMin);
  
  /// Groups hits associated with hexels into array of clusters
  /// \param hitsClustered Will be filled with RecHits collections, one per each cluster
  /// \param hexels Vector of hexels to which hits are associated
  /// \param energyMin Noise threshold - minimum energy of hit to be included
  void GetRecHitsPerHexel(std::vector<RecHits*> &hitsClustered,
                          std::vector<std::shared_ptr<Hexel>> &hexels, double energyMin);
  
  /// Returns subset of hits that are within given layer
  /// \param layer Layer index in which to look for hits
  /// \return Hits filtered by layer
  std::unique_ptr<RecHits> GetHitsInLayer(int layer);

  /// Returns layer index of i-th hit
  /// \param i Index of the hit
  inline int GetLayerOfHit(int i){return layer->at(i);}
  
  /// Checks if i-th hit is above the noise threshold and calculates sigma noise
  /// \return Returns a tuple: aboveThreshold, sigmaNoise
  std::tuple<bool,double> RecHitAboveThreshold(double iHit);
  
private:
  std::unique_ptr<RecHitCalibration> recHitCalib; ///< Stores current hits calibration
  
  std::vector<float> *eta;        ///< Pseudorapidity values of hit
  std::vector<float> *phi;        ///< Polar angle values of hit
  std::vector<float> *energy;     ///< Energy values of hit (GeV)
  std::vector<float> *x;          ///< X coordinates of hit
  std::vector<float> *y;          ///< Y coordinates of hit
  std::vector<float> *z;          ///< Z coordinates of hit
  std::vector<int> *layer;        ///< Layer indices of hit
  std::vector<unsigned int> *detid; ///< Unique detector IDs of hits stored in the collection
  std::vector<float> *thickness;  ///< Thickness of the sensor to which hit belongs
  std::vector<bool> *isHalf;      ///< Is the sensor to which hit belongs half of the full size
  std::vector<float> *time;       ///< Time signature of the hits
  std::vector<int> *cluster2d;    ///< Index of 2D clusters to which hit belongs
};

/// Class representing a single rec hit
class RecHit {
public:
  /// Default constructor, sets all internal values to zero
  RecHit();
  /// Constructor taking hit parameters as an input
  RecHit(float _eta, float _phi, float _energy, float _x, float _y, float _z, int _layer, unsigned int _detid, float _thickness, bool _isHalf, float _time, int _cluster2d);
  ~RecHit();
  
  /// Returns a hexel with parameters identical to those of this hit
  std::unique_ptr<Hexel> GetHexel();
  
  /// Checks if hit is above the noise threshold and calculates sigma noise
  /// \return Returns a tuple: aboveThreshold, sigmaNoise
  std::tuple<bool,double> RecHitAboveThreshold();
 
 float eta;        ///< Pseudorapidity values of hit
 float phi;        ///< Polar angle values of hit
 float energy;     ///< Energy values of hit (GeV)
 float x;          ///< X coordinates of hit
 float y;          ///< Y coordinates of hit
 float z;          ///< Z coordinates of hit
 int layer;        ///< Layer indices of hit
 unsigned int detid; ///< Unique detector IDs of hits stored in the collection
 float thickness;  ///< Thickness of the sensor to which hit belongs
 bool isHalf;      ///< Is the sensor to which hit belongs half of the full size
 float time;       ///< Time signature of the hits
 int cluster2d;    ///< Index of 2D clusters to which hit belongs
  
private:
  std::unique_ptr<RecHitCalibration> recHitCalibration;
  
};

#endif /* RecHits_hpp */
