//
//  TestbeamTrack.hpp
//
//  Created by Jeremi Niedziela on 30/11/2018.
//

#ifndef TestbeamTrack_hpp
#define TestbeamTrack_hpp

#include <vector>

#include <TTree.h>

/// Collection of test beam tracks
class TestbeamTrack{
public:
  /// Default constructor
  /// \param _tree Input tree from which test beam tracks will be read
  TestbeamTrack(TTree *_tree);
  
  /// Defauld destructor
  ~TestbeamTrack();
  
  /// Prints basic information about track
  void Print();
  
  /// Returns slope in X for given track
  inline double GetSlopeX(){return slopeX;}
  
  /// Returns slope in Y for given track
  inline double GetSlopeY(){return slopeY;}
  
  /// Returns offset in X for given track
  inline double GetOffsetX(){return offsetX;}
  
  /// Returns offset in Y for given track
  inline double GetOffsetY(){return offsetY;}
  
  /// Returns position in X for given track in given layer
  inline double GetX(int iLayer){return posX[iLayer] ? *(posX[iLayer]) : 0;}
  
  /// Returns position in Y for given track in given layer
  inline double GetY(int iLayer){return posY[iLayer] ? *(posY[iLayer]) : 0;}
  
private:
  std::vector<float*> posX; ///< Position X in given layer
  std::vector<float*> posY; ///< Position Y in given layer
  
  double slopeX; ///< Slope in X
  double slopeY; ///< Slope in Y
  
  double offsetX; ///< Offset in X
  double offsetY; ///< Offset in Y
  
  const int nLayers = 28;
};

#endif /* TestbeamTrack_hpp */
