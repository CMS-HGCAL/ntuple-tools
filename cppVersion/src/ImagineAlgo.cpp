//
//  HGCalImagineAlgo.cpp
//
//  Created by Jeremi Niedziela on 12/06/2018.
//

#include "ImagineAlgo.hpp"

using namespace std;

ImagingAlgo::ImagingAlgo(double _ecut,double _deltac[3],int _minClusters, int _dependSensor,int _verbosityLevel)
{
  recHitCalib = new RecHitCalibration();
  
  dependSensor = _dependSensor;
  for(int i=0;i<3;i++){deltac[i] = 2.0;};
  minClusters = 3;
  
  if(!dependSensor){  // (no sensor dependence, eta/phi coordinates for multi-clustering)
    kappa = 10.;
    ecut = 0.060;  // in absolute units
  }
  else{  // (with sensor dependence, cartesian coordinates for multi-clustering)
    kappa = 9.;
    ecut = 3;  // relative to the noise
  }
  // adjust params according to inputs, if necessary
  if(_ecut>=0) ecut = _ecut;
  if(_deltac[0]>=0) for(int i=0;i<3;i++){deltac[i] = _deltac[i];};
  if(_minClusters>=0) minClusters = _minClusters;
  
  // others
  verbosityLevel = 0;  // 0 - only basic info (default); 1 - additional info; 2 - detailed info printed
  if(_verbosityLevel>0) verbosityLevel = _verbosityLevel;
  
  // print out the setup
  if(verbosityLevel >= 1){
    cout<<"HGCalImagingAlgo setup: "<<endl;
    cout<<"   dependSensor: "<<dependSensor<<endl;
    cout<<"   deltac: "<<deltac<<endl;
    cout<<"   kappa: "<<kappa<<endl;
    cout<<"   ecut: "<<ecut<<endl;
    cout<<"   minClusters: "<<minClusters<<endl;
    cout<<"   verbosityLevel: "<<verbosityLevel<<endl;
  }
}

vector<int> ImagingAlgo::query_ball_point(vector<double> lpX, vector<double> lpY,
                                               double x, double y, double r)
{
  vector<int> foundIndices;
  
  for(int i=0;i<lpX.size();i++){
    if( pow(lpX[i]-x,2)+pow(lpY[i]-y,2) <= pow(r,2) ){
      foundIndices.push_back(i);
    }
  }
  return foundIndices;
}

double ImagingAlgo::distanceReal2(double x1, double y1, double x2, double y2)
{
  return pow(x2-x1, 2) + pow(y2-y1, 2);
}

double ImagingAlgo::calculateLocalDensity(vector<Hexel*> &nd,
                                               vector<double> lpX, vector<double> lpY,
                                               int layer)
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

vector<int> ImagingAlgo::sortIndicesRhoInverted(const vector<Hexel*> &v)
{
  vector<int> idx(v.size());
  iota(idx.begin(), idx.end(), 0);
  stable_sort(idx.begin(), idx.end(),[&v](int i1, int i2) {return v[i1]->rho > v[i2]->rho;});
  return idx;
}

vector<int> ImagingAlgo::sortIndicesDeltaInverted(const vector<Hexel*> &v)
{
  vector<int> idx(v.size());
  iota(idx.begin(), idx.end(), 0);
  stable_sort(idx.begin(), idx.end(),[&v](int i1, int i2) {return v[i1]->delta > v[i2]->delta;});
  return idx;
}

void ImagingAlgo::calculateDistanceToHigher(vector<Hexel*> &nodes)
{
  // sort vector of Hexels by decreasing local density
  vector<int> sortedIndices = sortIndicesRhoInverted(nodes);
  
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

void ImagingAlgo::findAndAssignClusters(vector<vector<Hexel*>> &current_clusters,
                                             vector<Hexel*> &nodes,
                                             vector<double> points_0,vector<double> points_1,
                                             double maxdensity,int layer)
{
  int clusterIndex = 0;
  
  vector<int> rs = sortIndicesRhoInverted(nodes);
  vector<int> ds = sortIndicesDeltaInverted(nodes);
  
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
      cout<<"Adding new cluster with index "<<clusterIndex<<endl;
      cout<<"Cluster center is hit "<<nodes[ds[i]]<<" with density rho: "<<nodes[ds[i]]->rho<<"and delta: "<< nodes[ds[i]]->delta<<endl;
    }
    clusterIndex++;
  }
  // at this point clusterIndex is equal to the number of cluster centers - if it is zero we are done
  if(clusterIndex == 0){
    return;
  }
  
  for(int i=0;i<clusterIndex;i++){
    current_clusters.push_back(vector<Hexel*>());
  }
  
  // assign to clusters, using the nearestHigher set from previous step (always set except for top density hit that is skipped)...
  for(int oi=1;oi<nodes.size();oi++){
    int ci = nodes[rs[oi]]->clusterIndex;
    if(ci == -1) nodes[rs[oi]]->clusterIndex = nodes[nodes[rs[oi]]->nearestHigher]->clusterIndex;
  }
  
  // assign points closer than dc to other clusters to border region and find critical border density
  double rho_b[clusterIndex];
  for(int i=0;i<clusterIndex;i++){rho_b[i]=0.;}
  
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
      //        cout<<"rho:"<<iNode->rho<<"\trho_b:"<<rho_b[ci]<<endl;
      iNode->isHalo = true;  // some issues to be debugged?
    }
    if(ci != -1){
      current_clusters[ci].push_back(iNode);
      if(verbosityLevel >= 2){
        cout<<"Pushing hit "<<iNode<<" into cluster with index "<<ci<<endl;
        cout<<"   rho_b[ci]: "<<rho_b[ci]<<", iNode.rho: "<<iNode->rho<<" iNode.isHalo: "<<iNode->isHalo<<endl;
      }
    }
  }
  return;
}

void ImagingAlgo::populate(vector<vector<Hexel*>> &points, RecHits *hits,double _ecut)
{
  // adjust ecut if necessary
  if(_ecut<0) _ecut = ecut;
  
  // init 2D hexels
  //    vector<vector<Hexel*>> points;
  for(int iLayer=0;iLayer<2*(maxlayer+1);iLayer++){
    points.push_back(vector<Hexel*>());
  }
  RecHit *hit;
  
  // loop over all hits and create the Hexel structure, skip energies below ecut
  for(int iHit=0;iHit<hits->N();iHit++){
    hit = hits->GetHit(iHit);
    if (hit->layer > maxlayer){
      delete hit;
      continue;  // current protection
    }
    // energy treshold dependent on sensor
    auto thresholdResult = recHitAboveThreshold(hit, _ecut, dependSensor);
    if(!get<0>(thresholdResult)){
      delete hit;
      continue;
    }
    // organise layers accoring to the sgn(z)
    int layerID = hit->layer + (hit->z > 0) * (maxlayer + 1);  // +1 - yes or no?
    double sigmaNoice = get<1>(thresholdResult);
    Hexel *hexel = hit->GetHexel();
//    Hexel *hexel = new Hexel(hit);
//    Hexel *hexel = new Hexel(hit->eta,hit->phi,hit->x,hit->y,hit->z,hit->energy,hit->thickness,hit->time,hit->detid,hit->layer,hit->cluster2d,hit->isHalf);
    
    hexel->sigmaNoise = sigmaNoice;
    points[layerID].push_back(hexel);
    delete hit;
  }
}

void ImagingAlgo::makeClusters(vector<vector<vector<Hexel*>>> &clusters, RecHits *hits,double _ecut)
{
  // adjust ecut if necessary
  if(_ecut<0) _ecut = ecut;
  
  // initialise list of per-layer-clusters
  for(int i=0;i<2 * (maxlayer + 1);i++){
    clusters.push_back(vector<vector<Hexel*>>());
  }
  
  // get the list of Hexels out of raw rechits
  vector<vector<Hexel*>> points;
  populate(points, hits, _ecut);
  
  // loop over all layers, and for each layer create a list of clusters. layers are organised according to the sgn(z)
  for(int layerID=0;layerID<2*(maxlayer + 1);layerID++){
    if (points[layerID].size() == 0) continue;  // protection
    int layer = layerID - (points[layerID][0]->z > 0) * (maxlayer + 1);  // map back to actual layer
    
    vector<double> pointX; // list of hexels'coordinate 0 for current layer
    vector<double> pointY; // list of hexels'coordinate 1 for current layer
    
    for(Hexel *hex : points[layerID]){
      pointX.push_back(hex->x);
      pointY.push_back(hex->y);
    }
    
    double maxdensity = calculateLocalDensity(points[layerID],pointX,pointY, layer);
    calculateDistanceToHigher(points[layerID]);// get distances to the nearest higher density
    vector<vector<Hexel*>> tmp;
    findAndAssignClusters(tmp, points[layerID],pointX,pointY, maxdensity, layer);  // get clusters per layer
    clusters[layerID] = tmp;
  }
}

void ImagingAlgo::getClusters(vector<BasicCluster*> &clusters_v,
                                   vector<vector<vector<Hexel*>>> &clusters)
{
  // loop over all layers and all clusters in each layer
  for(vector<vector<Hexel*>> clist_per_layer : clusters){
    for(vector<Hexel*> cluster : clist_per_layer){
      auto position = calculatePosition(cluster);
      
      // skip the clusters where position could not be computed (either all weights are 0, or all hexels are tagged as Halo)
      if(get<0>(position)==0 && get<1>(position)==0 && get<2>(position)==0) continue;
      
      double energy = 0;
      for(Hexel *iNode : cluster){
        if(!iNode->isHalo) energy += iNode->weight;
      }
      clusters_v.push_back(new BasicCluster(energy, position, cluster));
    }
    
    sort(clusters_v.begin( ), clusters_v.end( ), [ ](const BasicCluster *lhs,const BasicCluster *rhs){
      return lhs->energy > rhs->energy;
    });
  }
}

tuple<double,double,double> ImagingAlgo::calculatePosition(vector<Hexel*> &cluster)
{
  double total_weight=0, x=0, y=0, z=0;
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
    if(total_weight != 0.) return make_tuple(x / total_weight, y / total_weight, z / total_weight);
    else return make_tuple(0,0,0);
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
  return make_tuple(maxenergy_x, maxenergy_y, maxenergy_z);
}

tuple<bool,double> ImagingAlgo::recHitAboveThreshold(RecHit *hit,double _ecut,bool dependSensor)
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
      else cout<<"ERROR - silicon thickness has a nonsensical value"<<endl;
      // determine noise for each sensor/subdetector using RecHitCalibration library
    }
    sigmaNoise = 0.001 * recHitCalib->sigmaNoiseMeV(layer, thickIndex);  // returns threshold for EE, FH, BH (in case of BH thickIndex does not play a role)
  }
  bool aboveThreshold = energy >= _ecut * sigmaNoise;  // this checks if rechit energy is above the threshold of ecut (times the sigma noise for the sensor, if that option is set)
  return make_tuple(aboveThreshold,sigmaNoise);
}
