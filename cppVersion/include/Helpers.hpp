//
//  Helpers.hpp
//
//  Created by Jeremi Niedziela on 13/06/2018.
//

#ifndef Helpers_h
#define Helpers_h

#include "BasicCluster.hpp"

#include <chrono>
#include <iostream>

// Define constants
const int lastLayerEE = 28;  ///< Last layer of EE
const int lastLayerFH = 40;  ///< Last layer of FH
const int maxlayer = 52;     ///< Last layer of BH


/// Check if point (px,py) is withing circle (x,y,r)
inline bool pointWithinCircle(double px,double py,double x,double y,double r){
  return pow(r,2) >= pow(px-x,2) + pow(py-y,2);
}

/// Check if two circles (x1,x2,r1 and x2,y2,r2) overlap
inline bool circlesOverlap(double x1,double y1,double r1, double x2,double y2,double r2){
  return pow(r1+r2,2) <= pow(x1-x2,2) + pow(y1-y2,2);
}

/// Calculate distance between two points (x1,y1) and (x2,y2)
inline double distanceReal2(double x1, double y1, double x2, double y2)
{
  return pow(x2-x1, 2) + pow(y2-y1, 2);
}

/// Finds a circle among a set of them that is closest to the given point
/// \param Xs X coordinates of circles
/// \param Ys Y coordinates of circles
/// \param Rs Radii of the circles
/// \param x X coordinate of the point
/// \param y Y coordinate of the point
/// \return Returns index of the circle that is the closest to the point
inline int findClosestCircle(std::vector<double> Xs, std::vector<double> Ys, std::vector<double> Rs, double x, double y)
{
  int closestCircleIndex = -1;
  double currentLowestDistance = 99999999;
  
  for(uint i=0;i<Xs.size();i++){
    double distance = sqrt(pow(Xs[i]-x,2)+pow(Ys[i]-y,2)) - Rs[i];
    if(distance < currentLowestDistance){
      currentLowestDistance = distance;
      closestCircleIndex = i;
    }
  }
  return closestCircleIndex;
}

/// Finds all points in a circle
/// \param lpX X coordinates of points to verify
/// \param lpY Y coordinates of points to verify
/// \param x X coordinate of the search region center
/// \param y Y coordinate of the search region center
/// \param r Radius of the search region
/// \return Returns indices of points that are within the search circle (including the edge)
inline std::vector<int> queryBallPoint(std::vector<double> lpX, std::vector<double> lpY, double x, double y, double r)
{
  std::vector<int> foundIndices;

  for(uint i=0;i<lpX.size();i++){
    if(pow(lpX[i]-x,2)+pow(lpY[i]-y,2) <= pow(r,2)){
      foundIndices.push_back(i);
    }
  }
  return foundIndices;
}

/// Calculate duration between two events
/// \param t0 Start time
/// \param t1 End time
/// \return Difference between events t0 and t1 in seconds
inline double duration(std::chrono::time_point<std::chrono::system_clock> t0,
                std::chrono::time_point<std::chrono::system_clock> t1)
{
  auto elapsed_secs = t1-t0;
  typedef std::chrono::duration<float> float_seconds;
  auto secs = std::chrono::duration_cast<float_seconds>(elapsed_secs);
  return secs.count();
}

/// Returns current time
inline std::chrono::time_point<std::chrono::system_clock> now()
{
  return std::chrono::system_clock::now();
}

/// Struct containing one rec cluster and a vector of sim clusters matched to it
struct MatchedClusters {
  
  MatchedClusters(){
    recCluster = nullptr;
    simClusters = new std::vector<BasicCluster*>;
  }
  
  BasicCluster *recCluster;
  std::vector<BasicCluster*> *simClusters;
  
  inline double GetTotalSimEnergy(){
    double energy = 0;
    for(BasicCluster *cluster : *simClusters){
      energy += cluster->GetEnergy();
    }
    return energy;
  }
};

inline BasicCluster* GetBasicClusterFromRecHits(std::unique_ptr<RecHits> &hits)
{
  double recEnergy = hits->GetTotalEnergy();
  double xMax   = hits->GetXmax();
  double xMin   = hits->GetXmin();
  double yMax   = hits->GetYmax();
  double yMin   = hits->GetYmin();
  
  double clusterX = (xMax+xMin)/2.;
  double clusterY = (yMax+yMin)/2.;
  double clusterEta = hits->GetCenterEta();
  double clusterR = std::max((xMax-xMin)/2.,(yMax-yMin)/2.);
  
  BasicCluster *basicCluster = new BasicCluster(recEnergy,clusterX,clusterY,0,clusterEta,clusterR);
  return basicCluster;
}

/// Fills a vector of matched rec and sim clusters finding the nearest rec cluster for given sim cluster
/// \param matched Vector that will be filled with rec and sim clusters
/// \param recHitsPerCluster Vector of rec clusters
/// \param simHitsPerCluster Vector or sim clusters
/// \param layer Layer index
inline void matchClustersClosest(std::vector<MatchedClusters*> &matched, std::vector<RecHits*> &recHitsPerCluster, std::vector<RecHits*> &simHitsPerCluster, int layer)
{
  std::vector<int> alreadyAssociatedClusters;
  std::vector<RecHits*> hitsMatchedToRecClusters;
  std::vector<double> Xs,Ys,Rs, Es;
  
  for(uint recClusterIndex=0;recClusterIndex < recHitsPerCluster.size();recClusterIndex++){
    
    RecHits *recCluster = recHitsPerCluster[recClusterIndex];
    std::unique_ptr<RecHits> recHitsInLayerInCluster = recCluster->GetHitsInLayer(layer);

    if(recHitsInLayerInCluster->N()==0) continue;

    BasicCluster *basicCluster = GetBasicClusterFromRecHits(recHitsInLayerInCluster);
    
    Xs.push_back(basicCluster->GetX());
    Ys.push_back(basicCluster->GetY());
    Rs.push_back(basicCluster->GetRadius());
    Es.push_back(basicCluster->GetEnergy());
    
    MatchedClusters *matchedCluster = new MatchedClusters();
    
    matchedCluster->recCluster = basicCluster;
    matched.push_back(matchedCluster);
  }
  
  for(uint simClusterIndex=0;simClusterIndex<simHitsPerCluster.size();simClusterIndex++){
    
    RecHits *simCluster = simHitsPerCluster[simClusterIndex];
    std::unique_ptr<RecHits> simHitsInLayerInCluster = simCluster->GetHitsInLayer(layer);
    
    if(simHitsInLayerInCluster->N()==0) continue;
    
    BasicCluster *basicCluster = GetBasicClusterFromRecHits(simHitsInLayerInCluster);
    int parentRecCluster = findClosestCircle(Xs, Ys, Rs, basicCluster->GetX(), basicCluster->GetY());
    
    if(parentRecCluster < 0){
      if(ConfigurationManager::Instance()->GetVerbosityLevel() >= 1){
        std::cout<<"No rec cluster found for a sim cluster!!"<<std::endl;
      }
      continue;
    }
    
    
    matched[parentRecCluster]->simClusters->push_back(basicCluster);
  }
}

/// Fills a vector of unmatched rec and sim clusters (simply assigns all sim clusters to all rec each rec cluster)
/// \param matched Vector that will be filled with rec and sim clusters
/// \param recHitsPerCluster Vector of rec clusters
/// \param simHitsPerCluster Vector or sim clusters
/// \param layer Layer index
inline void matchClustersAllToAll(std::vector<MatchedClusters*> &matched, std::vector<RecHits*> &recHitsPerCluster, std::vector<RecHits*> &simHitsPerCluster, int layer)
{
  for(uint recClusterIndex=0;recClusterIndex < recHitsPerCluster.size();recClusterIndex++){
    
    RecHits *recCluster = recHitsPerCluster[recClusterIndex];
    std::unique_ptr<RecHits> recHitsInLayerInCluster = recCluster->GetHitsInLayer(layer);
    
    if(recHitsInLayerInCluster->N()==0) continue;
    
    for(uint simClusterIndex=0;simClusterIndex < simHitsPerCluster.size();simClusterIndex++){
      
      RecHits *simCluster = simHitsPerCluster[simClusterIndex];
      std::unique_ptr<RecHits> simHitsInLayerInCluster = simCluster->GetHitsInLayer(layer);
      
      if(simHitsInLayerInCluster->N()==0) continue;
      
      MatchedClusters *matchedCluster = new MatchedClusters();
      matchedCluster->recCluster = GetBasicClusterFromRecHits(recHitsInLayerInCluster);
    
      BasicCluster *basicCluster = GetBasicClusterFromRecHits(simHitsInLayerInCluster);
      matchedCluster->simClusters->push_back(basicCluster);
      
      matched.push_back(matchedCluster);
    }
  }
}
  
#endif /* Helpers_h */
