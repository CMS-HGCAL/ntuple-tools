//
//  RecHits.hpp
//
//  Created by Jeremi Niedziela on 07/06/2018.
//

#ifndef RecHits_hpp
#define RecHits_hpp

#include <TTree.h>

#include <vector>

class RecHit;

class RecHits {
public:
  RecHits();
  RecHits(TTree *_tree);
  ~RecHits();
  
  int N(){return (int)eta->size();}

  double GetTotalEnergy();
  double GetXmin();
  double GetXmax();
  double GetYmin();
  double GetYmax();
  
  std::vector<float> *eta;
  std::vector<float> *phi;
  std::vector<float> *energy;
  std::vector<float> *x;
  std::vector<float> *y;
  std::vector<float> *z;
  std::vector<int> *layer;
  std::vector<unsigned int> *detid;
  std::vector<float> *thickness;
  std::vector<bool> *isHalf;
  std::vector<float> *time;
  std::vector<int> *cluster2d;
  
  RecHit* GetHit(int index);
  void AddHit(RecHit *hit);
  
private:
  
};

class RecHit {
public:
  RecHit();
  RecHit(float _eta, float _phi, float _energy, float _x, float _y, float _z, int _layer, unsigned int _detid, float _thickness, bool _isHalf, float _time, int _cluster2d);
  ~RecHit();
  
  float eta;
  float phi;
  float energy;
  float x;
  float y;
  float z;
  int layer;
  unsigned int detid;
  float thickness;
  bool isHalf;
  float time;
  int cluster2d;
};

#endif /* RecHits_hpp */
