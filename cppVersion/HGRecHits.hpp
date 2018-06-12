//
//  HGRecHits.hpp
//
//  Created by Jeremi Niedziela on 07/06/2018.
//

#ifndef HGRecHits_hpp
#define HGRecHits_hpp

#include <TTree.h>

#include <vector>

class HGRecHit;

class HGRecHits {
public:
  HGRecHits();
  HGRecHits(TTree *_tree);
  ~HGRecHits();
  
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
  
  HGRecHit* GetHit(int index);
  void AddHit(HGRecHit *hit);
  
private:
  
};

class HGRecHit {
public:
  HGRecHit();
  HGRecHit(float _eta, float _phi, float _energy, float _x, float _y, float _z, int _layer, unsigned int _detid, float _thickness, bool _isHalf, float _time, int _cluster2d);
  ~HGRecHit();
  
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

#endif /* HGRecHits_hpp */
