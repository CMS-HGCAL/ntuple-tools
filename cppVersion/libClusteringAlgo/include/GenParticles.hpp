//
//  HGGenParticle.hpp
//
//  Created by Jeremi Niedziela on 07/06/2018.
//

#ifndef HGGenParticle_hpp
#define HGGenParticle_hpp

#include <TTree.h>

#include <vector>

/// Collection of generated particles
class GenParticles{
public:
  /// Default constructor
  /// \param _tree Input tree from which generated partiles will be read
  GenParticles(TTree *_tree);
  
  /// Defauld destructor
  ~GenParticles();
  
  /// Prints basic information about gen particle at given index
  void Print(int index);
  
  /// Deletes all objects and leaves collection with size zero
  void Clean();
  
  /// Returns number of stored gen particle objects
  inline int N(){return (int)reachedEE->size();}
  
  /// Returns pointer to vector telling for each gen particle if it reached the EE
  inline std::vector<int>* GetReachedEE(){return reachedEE;}
  
  /// Returns pointer to vector telling for each gen particle if it reached the EE
  inline std::vector<int>* GetPid(){return pid;}
  
private:
  std::vector<int> *reachedEE; ///< Info whether or not particle reached EE
  std::vector<int> *pid; ///< PDG PID
};

#endif /* HGGenParticle_hpp */
