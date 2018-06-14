//
//  Helpers.hpp
//
//  Created by Jeremi Niedziela on 13/06/2018.
//

#ifndef Helpers_h
#define Helpers_h


// Define constants
const int lastLayerEE = 28;  ///< Last layer of EE
const int lastLayerFH = 40;  ///< Last layer of FH
const int maxlayer = 52;     ///< Last layer of BH


/// Check if point (px,py) is withing circle (x,y,r)
inline bool pointWithinCircle(double px,double py,double x,double y,double r){
  return pow(r,2) >= pow(px-x,2) + pow(py-y,2);
}

/// Calculate distance between two points (x1,y1) and (x2,y2)
inline double distanceReal2(double x1, double y1, double x2, double y2)
{
  return pow(x2-x1, 2) + pow(y2-y1, 2);
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
  
  for(int i=0;i<lpX.size();i++){
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
