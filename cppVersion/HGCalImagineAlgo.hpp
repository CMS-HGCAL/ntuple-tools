//
//  RecHitCalibration.hpp
//
//  Created by Jeremi Niedziela on 07/06/2018.
//
//  Simple copy of the python script HGCalImagineAlgo.py:
//
//Implementation of (stand-alone) functionalities of HGCalImagingAlgo,
// HGCal3DClustering, and HGCalDepthPreClusterer based on
// their CMSSW implementations mainly in RecoLocalCalo/HGCalRecAlgos
//

// ??
// sys.setrecursionlimit(100000)
// noise thresholds and MIPs

#ifndef HGCalImagineAlgo_hpp
#define HGCalImagineAlgo_hpp

#include "RecHitCalibration.hpp"
#include "HGRecHits.hpp"

#include <iostream>
#include <vector>
#include <numeric>
#include <functional>

// definition of Hexel element

class Hexel{
public:
  Hexel(HGRecHit *hit=nullptr,double _sigmaNoise=-1){
    eta = 0;
    phi = 0;
    x = 0;
    y = 0;
    z = 0;
    time = -1;
    isHalfCell = false;
    weight = 0;
    fraction = 1;
    detid = -1;
    rho = 0;
    delta = 0;
    nearestHigher = -1;
    isBorder = false;
    isHalo = false;
    clusterIndex = -1;
    clusterRECOIndex = -1;
    sigmaNoise = 0.;
    thickness = 0.;
    
    if(hit){
      eta = hit->eta;
      phi = hit->phi;
      x = hit->x;
      y = hit->y;
      z = hit->z;
      weight = hit->energy;
      detid = hit->detid;
      layer = hit->layer;
      isHalfCell = hit->isHalf;
      thickness = hit->thickness;
      time = hit->time;
      clusterRECOIndex = hit->cluster2d;
    }
    if(_sigmaNoise>=0) sigmaNoise = _sigmaNoise;
  }
  
  bool __gt__(double _rho){return rho > _rho;}
  
  double eta, phi, x, y, z;
  double time, rho;
  double weight, fraction;
  double delta, sigmaNoise, thickness;
  bool isHalfCell, isBorder, isHalo;
  int layer, detid, clusterIndex, clusterRECOIndex, nearestHigher;
  
};
  
// definition of basic cluster (based on a set of sub-clusters or set of hexels)


class BasicCluster{
public:
  BasicCluster(double _energy=-1,std::tuple<double,double,double> position=std::make_tuple(0,0,0),std::vector<Hexel *> _thisCluster=std::vector<Hexel*>(0), int _algoId=-1,int _caloId=-1){
    eta = 0;
    phi = 0;
    x = 0;
    y = 0;
    z = 0;
    algoId = _algoId;
    caloId = _caloId;
    energy = (_energy>=0) ? _energy : 0.;
    
    if(std::get<0>(position) != 0 && std::get<1>(position) != 0 && std::get<2>(position) != 0){
      x = std::get<0>(position);
      y = std::get<1>(position);
      z = std::get<2>(position);

      double theta = acos(z/sqrt(x*x+y*y+z*z));
      eta = -log(tan(theta/2.));
      phi = atan2(y, x);
    }
  
    thisCluster = _thisCluster;
  }
  double x,y,z,eta,phi,energy;
  std::vector<Hexel*> thisCluster;
  int algoId, caloId;
  
};

// definition of the HGCalImagingAlgo class's methods & variables
class HGCalImagingAlgo {
private:
    // depth of the KDTree before brute force is applied
  int leafsize = 100000;
    // detector layers to consider
  int lastLayerEE = 28;  // last layer of EE
  int lastLayerFH = 40;  // last layer of FH
  int maxlayer = 52;  // last layer of BH

  bool dependSensor;
  double deltac[3];
  double kappa;
  double ecut;
  bool realSpaceCone;
  double multiclusterRadii[3];
  int minClusters;
  int verbosityLevel;
public:
  HGCalImagingAlgo(double _ecut=-1,double _deltac[3]=0,double _multiclusterRadii[3]=0,int _minClusters=-1, int _dependSensor=false,int _verbosityLevel=0){
        // sensor dependance or not
    dependSensor = _dependSensor;
    deltac[0] = 2.;
    deltac[1] = 2.;
    deltac[2] = 2.;
    
        // (multi)clustering parameters
    if(!dependSensor){  // (no sensor dependence, eta/phi coordinates for multi-clustering)
      // 2D clustering
      kappa = 10.;
      ecut = 0.060;  // in absolute units
      
      // multi-clustering
      realSpaceCone = false;
      multiclusterRadii[0] = 0.015;  // it's in eta/phi coordinates, per detector
      multiclusterRadii[1] = 0.015;
      multiclusterRadii[2] = 0.015;
      minClusters = 3;
    }
    else{  // (with sensor dependence, cartesian coordinates for multi-clustering)
      // 2D clustering
      kappa = 9.;
      ecut = 3;  // relative to the noise
      
      // multi-clustering
      realSpaceCone = true;
      multiclusterRadii[0] = 2.;  // it's in cartesian coordiantes, per detector
      multiclusterRadii[1] = 2.;
      multiclusterRadii[2] = 2.;
      minClusters = 3;
    }
        // adjust params according to inputs, if necessary
    if(_ecut>=0) ecut = _ecut;
    if(_deltac[0]>=0) for(int i=0;i<3;i++){deltac[i] = _deltac[i];};
    if(_minClusters>=0) minClusters = _minClusters;
    if(_multiclusterRadii[0]>=0) for(int i=0;i<3;i++){multiclusterRadii[i] = _multiclusterRadii[i];};

    // others
    verbosityLevel = 0;  // 0 - only basic info (default); 1 - additional info; 2 - detailed info printed
    if(_verbosityLevel>0) verbosityLevel = _verbosityLevel;

    // print out the setup
    if(verbosityLevel >= 1){
      std::cout<<"HGCalImagingAlgo setup: "<<std::endl;
      std::cout<<"   dependSensor: "<<dependSensor<<std::endl;
      std::cout<<"   deltac: "<<deltac<<std::endl;
      std::cout<<"   kappa: "<<kappa<<std::endl;
      std::cout<<"   ecut: "<<ecut<<std::endl;
      std::cout<<"   realSpaceCone: "<<realSpaceCone<<std::endl;
      std::cout<<"   multiclusterRadii: "<<multiclusterRadii<<std::endl;
      std::cout<<"   minClusters: "<<minClusters<<std::endl;
      std::cout<<"   verbosityLevel: "<<verbosityLevel<<std::endl;
    }
  }
  std::vector<int> query_ball_point(std::vector<double> lpX, std::vector<double> lpY, double x, double y, double r)
  {
    std::vector<int> foundIndices;
    
    for(int i=0;i<lpX.size();i++){
      if( pow(lpX[i]-x,2)+pow(lpY[i]-y,2) <= pow(r,2) ){
        foundIndices.push_back(i);
      }
    }
    return foundIndices;
  }
  
  // distance squared (in x/y) between the two objects (hexels, clusters)
  double distanceReal2(double x1, double y1, double x2, double y2)
  {
    return (pow(x2-x1, 2) + pow(y2-y1, 2));
  }
  
  // calculate max local density in a 2D plane of hexels
  double calculateLocalDensity(std::vector<Hexel*> &nd,std::vector<double> lpX, std::vector<double> lpY, int layer)
  {
    double maxdensity = 0;
    double delta_c = 0;
    
    if(layer <= lastLayerEE)       delta_c = deltac[0];
    else if(layer <= lastLayerFH)  delta_c = deltac[1];
    else                           delta_c = deltac[2];
    
    for(Hexel *iNode : nd){
      // search in a circle of radius delta_c or delta_c*sqrt(2) (not identical to search in the box delta_c)
      auto found = query_ball_point(lpX, lpY, iNode->x, iNode->y, delta_c);
      for(int j : found){
        double dist = distanceReal2(iNode->x,iNode->y, nd[j]->x, nd[j]->y);
        if(dist < delta_c * delta_c){
          iNode->rho += nd[j]->weight;
          if(iNode->rho > maxdensity) maxdensity = iNode->rho;
        }
      }
    }
    return maxdensity;
  }
    
  
  bool compareRhoInverted(int a, int b, std::vector<Hexel*> data)
  {
    return data[a]->rho > data[b]->rho;
  }
  
  std::vector<int> sortIndicesRhoInverted(const std::vector<Hexel*> &v) {
    std::vector<int> idx(v.size());
    iota(idx.begin(), idx.end(), 0);
    std::stable_sort(idx.begin(), idx.end(),
              [&v](int i1, int i2) {return v[i1]->rho > v[i2]->rho;});
    return idx;
  }
  
  std::vector<int> sortIndicesDeltaInverted(const std::vector<Hexel*> &v) {
    std::vector<int> idx(v.size());
    iota(idx.begin(), idx.end(), 0);
    std::stable_sort(idx.begin(), idx.end(),
              [&v](int i1, int i2) {return v[i1]->delta > v[i2]->delta;});
    return idx;
  }
  
  // calculate distance to the nearest hit with higher density (still does not use KDTree)
  void calculateDistanceToHigher(std::vector<Hexel*> &nodes)
  {
    // sort vector of Hexels by decreasing local density
    std::vector<int> sortedIndices = sortIndicesRhoInverted(nodes);
    
    // intial values, and check if there are any hits
    double maxdensity = 0.0;
    double nearestHigher = -1;
    if(nodes.size() > 0) maxdensity = nodes[sortedIndices[0]]->rho;
    else return;
        
    //   start by setting delta for the highest density hit to the most distant hit - this is a convention
    double dist2 = 0.;
    for(Hexel *jNode : nodes){
      double tmp = distanceReal2(nodes[sortedIndices[0]]->x,nodes[sortedIndices[0]]->y, jNode->x, jNode->y);
      if(tmp > dist2) dist2 = tmp;
    }
    nodes[sortedIndices[0]]->delta = pow(dist2, 0.5);
    nodes[sortedIndices[0]]->nearestHigher = nearestHigher;
    
    // now we save the largest distance as a starting point
    double max_dist2 = dist2;
    // calculate all remaining distances to the nearest higher density
    for(int oi=1;oi<nodes.size();oi++){  // start from second-highest density
      dist2 = max_dist2;
      // we only need to check up to oi since hits are ordered by decreasing density
      // and all points coming BEFORE oi are guaranteed to have higher rho and the ones AFTER to have lower rho
      for(int oj=0;oj<oi;oj++){
        double tmp = distanceReal2(nodes[sortedIndices[oi]]->x,nodes[sortedIndices[oi]]->y, nodes[sortedIndices[oj]]->x,nodes[sortedIndices[oj]]->y);
        if(tmp <= dist2){  // this "<=" instead of "<" addresses the (rare) case when there are only two hits
          dist2 = tmp;
          nearestHigher = sortedIndices[oj];
        }
      }
      nodes[sortedIndices[oi]]->delta = pow(dist2, 0.5);
      nodes[sortedIndices[oi]]->nearestHigher = nearestHigher;  // this uses the original unsorted hitlist
    }
  }
  
     
  // find cluster centers that satisfy delta & maxdensity/kappa criteria, and assign coresponding hexels
  void findAndAssignClusters(std::vector<std::vector<Hexel*>> &current_clusters,
                             std::vector<Hexel*> &nodes,
                             std::vector<double> points_0,std::vector<double> points_1,
                             double maxdensity,int layer,int _verbosityLevel=-1){
    
    // adjust verbosityLevel if necessary
    if(_verbosityLevel<0) _verbosityLevel = verbosityLevel;
    int clusterIndex = 0;
    
    std::vector<int> rs = sortIndicesRhoInverted(nodes);
    std::vector<int> ds = sortIndicesDeltaInverted(nodes);
    
    double _delta_c;
    if(layer <= lastLayerEE)      _delta_c = deltac[0];
    else if(layer <= lastLayerFH) _delta_c = deltac[1];
    else                          _delta_c = deltac[2];
    
    for(int i=0; i<nodes.size();i++){
      if(nodes[ds[i]]->delta < _delta_c) break;  // no more cluster centers to be looked at
      // skip this as a potential cluster center because it fails the density cut
      if(dependSensor){
        if(nodes[ds[i]]->rho < kappa * nodes[ds[i]]->sigmaNoise){
          continue;  // set equal to kappa times noise threshold
        }
      }
      else{
        if(nodes[ds[i]]->rho < maxdensity / kappa) continue;
      }
      // store cluster index
      nodes[ds[i]]->clusterIndex = clusterIndex;
      
      if(verbosityLevel >= 2){
        std::cout<<"Adding new cluster with index "<<clusterIndex<<std::endl;
        std::cout<<"Cluster center is hit "<<nodes[ds[i]]<<" with density rho: "<<nodes[ds[i]]->rho<<"and delta: "<< nodes[ds[i]]->delta<<std::endl;
      }
      clusterIndex++;
    }
    // at this point clusterIndex is equal to the number of cluster centers - if it is zero we are done
    if(clusterIndex == 0){
      return;
    }
    
    for(int i=0;i<clusterIndex;i++){
      current_clusters.push_back(std::vector<Hexel*>());
    }
    
    // assign to clusters, using the nearestHigher set from previous step (always set except for top density hit that is skipped)...
    for(int oi=1;oi<nodes.size();oi++){
      int ci = nodes[rs[oi]]->clusterIndex;
      if(ci == -1) nodes[rs[oi]]->clusterIndex = nodes[nodes[rs[oi]]->nearestHigher]->clusterIndex;
    }
    
    // assign points closer than dc to other clusters to border region and find critical border density
    double rho_b[clusterIndex];
    for(int i=0;i<clusterIndex;i++){rho_b[i]=0.;}
    
//    lp = spatial.KDTree(list(zip(points_0, points_1)), leafsize=self.leafsize);  // new KDTree
    // now loop on all hits again :( and check: if there are hits from another cluster within d_c -> flag as border hit
    
    for(Hexel *iNode : nodes){
      
      int ci = iNode->clusterIndex;
      bool flag_isolated = true;
      if(ci != -1){// search in a circle of radius delta_c or delta_c*sqrt(2) (not identical to search in the box delta_c)

        auto found = query_ball_point(points_0,points_1, iNode->x, iNode->y, _delta_c);
        for(int j : found){
          // check if the hit is not within d_c of another cluster
          if(nodes[j]->clusterIndex != -1){
            double dist2 = distanceReal2(nodes[j]->x, nodes[j]->y , iNode->x, iNode->y);
            if(dist2 < _delta_c * _delta_c && nodes[j]->clusterIndex != ci){
              // in which case we assign it to the border
              iNode->isBorder = true;
              break;
            }
            // because we are using two different containers, we have to make sure that we don't unflag the
            // hit when it finds *itself* closer than delta_c
            
            if(dist2 < _delta_c * _delta_c && dist2 != 0. && nodes[j]->clusterIndex == ci){
              // this is not an isolated hit
              flag_isolated = false;
            }
          }
        }
        if(flag_isolated) iNode->isBorder = true;  // the hit is more than delta_c from any of its brethren
      }
      // check if this border hit has density larger than the current rho_b and update
      if(iNode->isBorder && rho_b[ci] < iNode->rho)
        rho_b[ci] = iNode->rho;
    }
    // flag points in cluster with density < rho_b as halo points, then fill the cluster vector
    for(Hexel *iNode : nodes){
      int ci = iNode->clusterIndex;
      if(ci != -1 && iNode->rho <= rho_b[ci]){
//        std::cout<<"rho:"<<iNode->rho<<"\trho_b:"<<rho_b[ci]<<std::endl;
        iNode->isHalo = true;  // some issues to be debugged?
      }
      if(ci != -1){
        current_clusters[ci].push_back(iNode);
        if (verbosityLevel >= 2){
          std::cout<<"Pushing hit "<<iNode<<" into cluster with index "<<ci<<std::endl;
          std::cout<<"   rho_b[ci]: "<<rho_b[ci]<<", iNode.rho: "<<iNode->rho<<" iNode.isHalo: "<<iNode->isHalo<<std::endl;
        }
      }
    }
    return;
  }
  // make list of Hexels out of rechits
  void populate(std::vector<std::vector<Hexel*>> &points, HGRecHits *hits,double _ecut=-1)
  {
    std::cout<<"Populate:"<<hits->N()<<std::endl;
    // adjust ecut if necessary
    if(_ecut<0) _ecut = ecut;
    
    // init 2D hexels
//    std::vector<std::vector<Hexel*>> points;
    for(int iLayer=0;iLayer<2*(maxlayer+1);iLayer++){
      points.push_back(std::vector<Hexel*>());
    }
    
    std::cout<<"N hits:"<<hits->N()<<std::endl;
    int skipIter=0;
    int skipIter2=0;
    // loop over all hits and create the Hexel structure, skip energies below ecut
    for(int iHit=0;iHit<hits->N();iHit++){
      HGRecHit *hit = hits->GetHit(iHit);
      
      if (hit->layer > maxlayer){
        skipIter++;
        continue;  // current protection
      }
      // energy treshold dependent on sensor
      auto thresholdResult = recHitAboveThreshold(hit, _ecut, dependSensor);
      if(!std::get<0>(thresholdResult)){
        skipIter2++;
        continue;
      }
      // organise layers accoring to the sgn(z)
      int layerID = hit->layer + (hit->z > 0) * (maxlayer + 1);  // +1 - yes or no?
      points[layerID].push_back(new Hexel(hit,std::get<1>(thresholdResult)));
    }
    std::cout<<"Skipped:"<<skipIter<<std::endl;
    std::cout<<"Skipped2:"<<skipIter2<<std::endl;
    std::cout<<"N points in populate:"<<points.size()<<std::endl;
//    return &points;
  }
  // make 2D clusters out of rechits (need to introduce class with input params: delta_c, kappa, ecut, ...)
  void makeClusters(std::vector<std::vector<std::vector<Hexel*>>> &clusters, HGRecHits *hits,double _ecut=-1){
    std::cout<<"Make clusters:"<<hits->N()<<std::endl;
    // adjust ecut if necessary
    if(_ecut<0) _ecut = ecut;
    
    // init 2D cluster lists
//    std::vector<std::vector<std::vector<Hexel*>>> clusters; = new std::vector<std::vector<std::vector<Hexel*>>>();
    
    // initialise list of per-layer-clusters
    for(int i=0;i<2 * (maxlayer + 1);i++){
      clusters.push_back(std::vector<std::vector<Hexel*>>());
    }
    
    // get the list of Hexels out of raw rechits
    std::vector<std::vector<Hexel*>> points;
    populate(points, hits, _ecut);
    std::cout<<"N points:"<<points.size()<<std::endl;
    // loop over all layers, and for each layer create a list of clusters. layers are organised according to the sgn(z)
    for(int layerID=0;layerID<2*(maxlayer + 1);layerID++){
      if (points[layerID].size() == 0) continue;  // protection
      int layer = layerID - (points[layerID][0]->z > 0) * (maxlayer + 1);  // map back to actual layer
      
      std::vector<double> pointX; // list of hexels'coordinate 0 for current layer
      std::vector<double> pointY; // list of hexels'coordinate 1 for current layer
      
      for(Hexel *hex : points[layerID]){
        pointX.push_back(hex->x);
        pointY.push_back(hex->y);
      }
      
      double maxdensity = calculateLocalDensity(points[layerID],pointX,pointY, layer);
      calculateDistanceToHigher(points[layerID]);// get distances to the nearest higher density
      std::vector<std::vector<Hexel*>> tmp;
      findAndAssignClusters(tmp, points[layerID],pointX,pointY, maxdensity, layer);  // get clusters per layer
      clusters[layerID] = tmp;
      
    }
    // return the clusters list
//    return clusters;
  }
    // get basic clusters from the list of 2D clusters
  void getClusters(std::vector<BasicCluster*> &clusters_v,
                   std::vector<std::vector<std::vector<Hexel*>>> &clusters,
                   int _verbosityLevel=-1)
  {
    std::cout<<"Get clusters n clusters:"<<clusters.size()<<std::endl;
    int notHalo = 0;
    // loop over all layers and all clusters in each layer
    for(std::vector<std::vector<Hexel*>> clist_per_layer : clusters){
      for(std::vector<Hexel*> cluster : clist_per_layer){
        auto position = calculatePosition(cluster);
        
        // skip the clusters where position could not be computed (either all weights are 0, or all hexels are tagged as Halo)
        if(std::get<0>(position)==0 && std::get<1>(position)==0 && std::get<2>(position)==0){
          std::cout<<"skip"<<std::endl;
          continue;
        }
          
        double energy = 0;
        for(Hexel *iNode : cluster){
          if(!iNode->isHalo){
            notHalo++;
            energy += iNode->weight;
          }
        }
        clusters_v.push_back(new BasicCluster(energy, position, cluster));
      }
      
      std::sort(clusters_v.begin( ), clusters_v.end( ), [ ](const BasicCluster *lhs,const BasicCluster *rhs){
        return lhs->energy > rhs->energy;
      });
    }
    
    std::cout<<"Not halo:"<<notHalo<<std::endl;
  }
  /*
    // make multi-clusters starting from the 2D clusters, without KDTree
    def makePreClusters(self, clusters, multiclusterRadii=None, minClusters=None, verbosityLevel=None):
        // adjust multiclusterRadii, minClusters and/or verbosityLevel if necessary
        if multiclusterRadii is None:
            multiclusterRadii = self.multiclusterRadii
        if minClusters is None:
            minClusters = self.minClusters
        if verbosityLevel is None:
            verbosityLevel = self.verbosityLevel
        // get clusters in one list (just following original approach)
        thecls = self.getClusters(clusters)

        // init lists and vars
        thePreClusters = []
        vused = [0.] * len(thecls)
        used = 0
        // indices sorted by decreasing energy
        es = sorted(range(len(thecls)), key=lambda k: thecls[k].energy, reverse=True)
        // loop over all clusters
        index = 0
        for i in range(0, len(thecls)):
            if(vused[i] == 0):
                temp = [thecls[es[i]]]
                if (thecls[es[i]].z > 0):
                    vused[i] = 1
                else:
                    vused[i] = -1
                used += 1
                for j in range(i + 1, len(thecls)):
                    if(vused[j] == 0):
                        distanceCheck = 9999.
                        if(self.realSpaceCone):
                            distanceCheck = distanceReal2(thecls[es[i]], thecls[es[j]])
                        else:
                            distanceCheck = distanceDR2(thecls[es[i]], thecls[es[j]])
                        layer = thecls[es[j]].thisCluster[0].layer
                        multiclusterRadius = 9999.
                        multiclusterRadius = multiclusterRadii[0]
                        if(layer > self.lastLayerEE and layer <= self.lastLayerFH):
                            multiclusterRadius = multiclusterRadii[1]
                        else:
                            multiclusterRadius = multiclusterRadii[2]
                        if(distanceCheck < multiclusterRadius * multiclusterRadius and int(thecls[es[i]].z * vused[i]) > 0):
                            temp.append(thecls[es[j]])
                            vused[j] = vused[i]
                            used += 1
                if(len(temp) > minClusters):
                    position = getMultiClusterPosition(temp)
                    energy = getMultiClusterEnergy(temp)
                    thePreClusters.append(BasicCluster(energy=energy, position=position, thisCluster=temp))
                    if (verbosityLevel >= 1):
                        print("Multi-cluster index: ", index, ", No. of 2D-clusters = ", len(temp), ", Energy  = ",
                              energy, ", Phi = ", position.phi(), ", Eta = ", position.eta(), ", z = ", position.z())
                    index += 1
        return thePreClusters
*/
  
  /*
    // make multi-clusters starting from the 2D clusters, with KDTree
    def make3DClusters(self, clusters, multiclusterRadii=None, minClusters=None, verbosityLevel=None):
        // adjust multiclusterRadii, minClusters and/or verbosityLevel if necessary
        if multiclusterRadii is None:
            multiclusterRadii = self.multiclusterRadii
        if minClusters is None:
            minClusters = self.minClusters
        if verbosityLevel is None:
            verbosityLevel = self.verbosityLevel
        // get clusters in one list (just following original approach)
        thecls = self.getClusters(clusters)

        // init "points" of 2D clusters for KDTree serach and zees of layers (check if it is really needed)
        points = [[] for i in range(0, 2 * (self.maxlayer + 1))]  // initialise list of per-layer-lists of clusters
        zees = [0. for layer in range(0, 2 * (self.maxlayer + 1))]
        for cls in thecls:  // organise layers accoring to the sgn(z)
            layerID = cls.thisCluster[0].layer
            layerID += (cls.z > 0) * (self.maxlayer + 1)  // +1 - yes or no?
            points[layerID].append(cls)
            zees[layerID] = cls.z

        // init lists and vars
        thePreClusters = []
        vused = [0.] * len(thecls)
        used = 0

        // indices sorted by decreasing energy
        es = sorted(range(len(thecls)), key=lambda k: thecls[k].energy, reverse=True)
        // loop over all clusters
        index = 0
        for i in range(0, len(thecls)):
            // if(vused[i]==0):
            if (thecls[es[i]]._usedIn3DClust == 0):
                temp = [thecls[es[i]]]
                if (thecls[es[i]].z > 0):
                    thecls[es[i]]._usedIn3DClust = 1
                else:
                    thecls[es[i]]._usedIn3DClust = -1
                used += 1
                from_ = [thecls[es[i]].x, thecls[es[i]].y, thecls[es[i]].z]
                firstlayer = (thecls[es[i]].z > 0) * (self.maxlayer + 1)
                lastlayer = firstlayer + self.maxlayer + 1
                for j in range(firstlayer, lastlayer):
                    if(zees[j] == 0.):
                        continue
                    to_ = [0., 0., zees[j]]
                    to_[0] = (from_[0] / from_[2]) * to_[2]
                    to_[1] = (from_[1] / from_[2]) * to_[2]
                    layer = j - (zees[j] > 0) * (self.maxlayer + 1)  // maps back from index used for KD trees to actual layer
                    multiclusterRadius = 9999.
                    if(layer <= self.lastLayerEE):
                        multiclusterRadius = multiclusterRadii[0]
                    elif(layer <= self.lastLayerFH):
                        multiclusterRadius = multiclusterRadii[1]
                    elif(layer <= self.maxlayer):
                        multiclusterRadius = multiclusterRadii[2]
                    else:
                        print("ERROR: Nonsense layer value - cannot assign multicluster radius")
                    // KD-tree search in layer j
                    points_0 = [cls.x for cls in points[j]]  // list of cls' coordinate 0 for layer j
                    points_1 = [cls.y for cls in points[j]]  // list of cls' coordinate 1 for layer j
                    hit_kdtree = spatial.KDTree(list(zip(points_0, points_1)), leafsize=self.leafsize)  // create KDTree
                    found = hit_kdtree.query_ball_point([to_[0], to_[1]], multiclusterRadius)
                    for k in found:
                        h_to = Hexel()
                        h_to.x = to_[0]
                        h_to.y = to_[1]  // dummy object
                        if((points[j][k]._usedIn3DClust == 0) and (distanceReal2(points[j][k], h_to) < multiclusterRadius**2)):
                            temp.append(points[j][k])
                            points[j][k]._usedIn3DClust = thecls[es[i]]._usedIn3DClust
                            used += 1
                if(len(temp) > minClusters):
                    position = getMultiClusterPosition(temp)
                    energy = getMultiClusterEnergy(temp)
                    thePreClusters.append(BasicCluster(energy=energy, position=position, thisCluster=temp))
                    if (verbosityLevel >= 1):
                        print ("Multi-cluster index: ", index, ", No. of 2D-clusters = ", len(temp), ", Energy  = ",
                               energy, ", Phi = ", position.phi(), ", Eta = ", position.eta(), ", z = ", position.z())
                    index += 1
        return thePreClusters
*/

  /*
// distance squared (in eta/phi) between the two objects (hexels, clusters)
def distanceDR2(Hex1, Hex2):
    return (pow(Hex2.eta - Hex1.eta, 2) + pow(Hex2.phi - Hex1.phi, 2))
*/
// position of the cluster, based on hexels positions weighted by the energy


  std::tuple<double,double,double> calculatePosition(std::vector<Hexel*> &cluster){
    double total_weight = 0.;
    double x = 0.;
    double y = 0.;
    double z = 0.;
    bool haloOnlyCluster = true;
    
    // check if haloOnlyCluster
    for(Hexel *iNode : cluster){
      if(!iNode->isHalo) haloOnlyCluster = false;
    }
    if(!haloOnlyCluster){
      for(Hexel *iNode : cluster){
        if(!iNode->isHalo){
          total_weight += iNode->weight;
          x += iNode->x * iNode->weight;
          y += iNode->y * iNode->weight;
          z += iNode->z * iNode->weight;
        }
      }
      if(total_weight != 0.) return std::make_tuple(x / total_weight, y / total_weight, z / total_weight);
      else return std::make_tuple(0,0,0);
    }
    
    double maxenergy = -1.0;
    double maxenergy_x=0., maxenergy_y=0., maxenergy_z=0.;
    for(Hexel *iNode : cluster){
      if(iNode->weight > maxenergy){
        maxenergy = iNode->weight;
        maxenergy_x = iNode->x;
        maxenergy_y = iNode->y;
        maxenergy_z = iNode->z;
      }
    }
    return std::make_tuple(maxenergy_x, maxenergy_y, maxenergy_z);
  }

  /*
// get position of the multi-cluster, based on the positions of its 2D clusters weighted by the energy
def getMultiClusterPosition(multi_clu):
    if(len(multi_clu) == 0):
        return ROOT.Math.XYZPoint()
    mcenergy = getMultiClusterEnergy(multi_clu)
    if (mcenergy == 0):
        return ROOT.Math.XYZPoint()

    // compute weighted mean x/y/z position
    acc_x = 0.0
    acc_y = 0.0
    acc_z = 0.0
    totweight = 0.0
    for layer_clu in multi_clu:
        if(layer_clu.energy < 0.01 * mcenergy):
            continue  // cutoff < 1% layer energy contribution
        weight = layer_clu.energy  // weight each corrdinate only by the total energy of the layer cluster
        acc_x += layer_clu.x * weight
        acc_y += layer_clu.y * weight
        acc_z += layer_clu.z * weight
        totweight += weight
    if (totweight != 0):
        acc_x /= totweight
        acc_y /= totweight
        acc_z /= totweight

    return ROOT.Math.XYZPoint(acc_x, acc_y, acc_z)  // return x/y/z in absolute coordinates

// get energy of the multi-cluster, based on its 2D clusters
*/

  /*
def getMultiClusterEnergy(multi_clu):
    acc = 0.
    for layer_clu in multi_clu:
        acc += layer_clu.energy
    return acc
*/
// determine if the rechit energy is above the desired treshold

  
  std::tuple<bool,double> recHitAboveThreshold(HGRecHit *hit,double _ecut,bool dependSensor=true)
  {
    double sigmaNoise = 1.;
    int layer = hit->layer;
    double thickness = hit->thickness;
    double energy = hit->energy;
    
    if(dependSensor){
      int thickIndex = -1;
      
      if(layer <= lastLayerFH){  // EE + FH
        if(thickness > 99. and thickness < 101.)        thickIndex = 0;
        else if(thickness > 199. and thickness < 201.)  thickIndex = 1;
        else if(thickness > 299. and thickness < 301.)  thickIndex = 2;
        else std::cout<<"ERROR - silicon thickness has a nonsensical value"<<std::endl;
        // determine noise for each sensor/subdetector using RecHitCalibration library
      }
      RecHitCalibration *recHitCalib = new RecHitCalibration();
      sigmaNoise = 0.001 * recHitCalib->sigmaNoiseMeV(layer, thickIndex);  // returns threshold for EE, FH, BH (in case of BH thickIndex does not play a role)
    }
    bool aboveThreshold = energy >= _ecut * sigmaNoise;  // this checks if rechit energy is above the threshold of ecut (times the sigma noise for the sensor, if that option is set)
    return std::make_tuple(aboveThreshold,sigmaNoise);
  }
};
  
  /*
def getEnergy(item):
     eturn item.energy
   */

#endif
