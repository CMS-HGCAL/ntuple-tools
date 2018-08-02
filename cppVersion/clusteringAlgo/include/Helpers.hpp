//
//  Helpers.hpp
//
//  Created by Jeremi Niedziela on 13/06/2018.
//

#ifndef Helpers_h
#define Helpers_h

#include <TMath.h>

#include <chrono>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>


// Define constants
const int lastLayerEE = 28;  ///< Last layer of EE
const int lastLayerFH = 40;  ///< Last layer of FH
const int maxlayer = 52;     ///< Last layer of BH

enum EDet {
  kEE,  ///< electromagneric endcap (silicon)
  kFH,  ///< front hadronic endcap (silicon)
  kBH   ///< back hadronic endcap (plastic)
};

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

inline int findMostDetIDsharingCluster(std::vector<std::vector<unsigned int>> &recDetIDs,
                                       std::vector<unsigned int> &detIDsInSimCluster)
{
  int mostSharingIndex = -1;
  
  size_t currentMaxSharing = 0;
  int index=0;
  
  for(std::vector<unsigned int> recDetIDsInCluster : recDetIDs){
    std::vector<unsigned int> common;
    std::set_intersection(detIDsInSimCluster.begin(), detIDsInSimCluster.end(),
                          detIDsInSimCluster.begin(), detIDsInSimCluster.end(), std::back_inserter(common));
    
    if(common.size() > currentMaxSharing){
      currentMaxSharing = common.size();
      mostSharingIndex = index;
    }
    index++;
  }
  return mostSharingIndex;
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
  
#endif /* Helpers_h */
