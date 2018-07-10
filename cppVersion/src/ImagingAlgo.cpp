//
//  HGCalImagineAlgo.cpp
//
//  Created by Jeremi Niedziela on 12/06/2018.
//

#include "ImagingAlgo.hpp"
#include "Helpers.hpp"

#include <algorithm>

using namespace std;

ImagingAlgo::ImagingAlgo()
{
  recHitCalib = unique_ptr<RecHitCalibration>(new RecHitCalibration());
  config = ConfigurationManager::Instance();
  
  dependSensor = config->GetDependSensor();
  kappa = config->GetKappa();
  ecut = config->GetEnergyMin();
  minClusters = config->GetMinClusters();
  verbosityLevel = config->GetVerbosityLevel();

  if(config->GetEnergyDensityFunction() == "step"){
    // param [0] says what's the limit to include or reject hit (critical distance)
    // param [1] says how much should be added to the energy density if hit is accepted
    energyDensityFunction = new TF1("step function", "((x < [0]) ? [1] : 0)", -1000, 1000);
  }
  else if(config->GetEnergyDensityFunction() == "gaus"){
    // param [0] is the distribution width
    // param [1] scales the distribution (should be set to something proportional to the energy of the hit)
    energyDensityFunction = new TF1("gaussian", "[1]/(sqrt(2*TMath::Pi())*[0])*exp(-x*x/(2*[0]*[0]))", -1000, 1000);
  }
  else if(config->GetEnergyDensityFunction() == "exp"){
    // param [0] is the critical distance (further than that we don't include hits)
    // param [1] scales the distribution (should be set to something proportional to the energy of the hit)
    energyDensityFunction = new TF1("exp", "((x < [0]) ? [1]*exp(-x/[0]) : 0)", -1000, 1000);
  }
  else{
    cout<<"ERROR -- unknown energy density function:"<<config->GetEnergyDensityFunction()<<endl;
    exit(0);
  }
  
  // print out the setup
  if(verbosityLevel >= 1){
    cout<<"HGCalImagingAlgo setup: "<<endl;
    cout<<"   dependSensor: "<<dependSensor<<endl;
    cout<<"   deltac_EE: "<<config->GetDeltac(kEE)<<endl;
    cout<<"   deltac_FH: "<<config->GetDeltac(kFH)<<endl;
    cout<<"   deltac_BH: "<<config->GetDeltac(kBH)<<endl;
    cout<<"   kappa: "<<kappa<<endl;
    cout<<"   ecut: "<<ecut<<endl;
    cout<<"   minClusters: "<<minClusters<<endl;
    cout<<"   verbosityLevel: "<<verbosityLevel<<endl;
  }
}

ImagingAlgo::~ImagingAlgo()
{
}

double ImagingAlgo::calculateLocalDensity(vector<unique_ptr<Hexel>> &hexels,
                                          vector<double> lpX, vector<double> lpY,
                                          int layer)
{
  double maxdensity = 0;
  double criticalDistance = 0;

  if(layer <= lastLayerEE)       criticalDistance = config->GetCriticalDistance(kEE);
  else if(layer <= lastLayerFH)  criticalDistance = config->GetCriticalDistance(kFH);
  else                           criticalDistance = config->GetCriticalDistance(kBH);

  for(unique_ptr<Hexel> &iNode : hexels){
    // search in a circle of radius "criticalDistance" or "criticalDistance"*sqrt(2) (not identical to search in the box "criticalDistance")
    auto found = queryBallPoint(lpX, lpY, iNode->x, iNode->y, criticalDistance);
    for(int j : found){
      double dist = sqrt(distanceReal2(iNode->x,iNode->y, hexels[j]->x, hexels[j]->y));
      
      energyDensityFunction->SetParameter(0,criticalDistance);
      energyDensityFunction->SetParameter(1,hexels[j]->weight);
      
      iNode->rho += energyDensityFunction->Eval(dist);
      if(iNode->rho > maxdensity) maxdensity = iNode->rho;
    }
  }
  return maxdensity;
}

vector<int> ImagingAlgo::sortIndicesRhoInverted(const vector<unique_ptr<Hexel>> &v)
{
  vector<int> idx(v.size());
  iota(idx.begin(), idx.end(), 0);
  stable_sort(idx.begin(), idx.end(),[&v](int i1, int i2) {return v[i1]->rho > v[i2]->rho;});
  return idx;
}

vector<int> ImagingAlgo::sortIndicesDeltaInverted(const vector<unique_ptr<Hexel>> &v)
{
  vector<int> idx(v.size());
  iota(idx.begin(), idx.end(), 0);
  stable_sort(idx.begin(), idx.end(),[&v](int i1, int i2) {return v[i1]->delta > v[i2]->delta;});
  return idx;
}

void ImagingAlgo::calculateDistanceToHigher(vector<unique_ptr<Hexel>> &nodes)
{
  // sort vector of Hexels by decreasing local density
  vector<int> sortedIndices = sortIndicesRhoInverted(nodes);

  // intial values, and check if there are any hits
  double maxdensity = 0.0; // this value is never used - do we need it? !!!
  double nearestHigher = -1;
  if(nodes.size() > 0) maxdensity = nodes[sortedIndices[0]]->rho;
  else return;

  //   start by setting delta for the highest density hit to the most distant hit - this is a convention
  double dist2 = 0.;
  for(auto &jNode : nodes){
    double tmp = distanceReal2(nodes[sortedIndices[0]]->x,nodes[sortedIndices[0]]->y, jNode->x, jNode->y);
    if(tmp > dist2) dist2 = tmp;
  }
  nodes[sortedIndices[0]]->delta = pow(dist2, 0.5);
  nodes[sortedIndices[0]]->nearestHigher = nearestHigher;

  // now we save the largest distance as a starting point
  double max_dist2 = dist2;
  // calculate all remaining distances to the nearest higher density
  for(uint oi=1;oi<nodes.size();oi++){  // start from second-highest density
    dist2 = max_dist2;
    // we only need to check up to oi since hits are ordered by decreasing density
    // and all points coming BEFORE oi are guaranteed to have higher rho and the ones AFTER to have lower rho
    for(uint oj=0;oj<oi;oj++){
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

void ImagingAlgo::findAndAssignClusters(vector<vector<unique_ptr<Hexel>>> &clusters,
                                             vector<unique_ptr<Hexel>> &nodes,
                                             vector<double> points_0,vector<double> points_1,
                                             double maxDensity,int layer)
{
  int clusterIndex = 0;

  vector<int> rs = sortIndicesRhoInverted(nodes);
  vector<int> ds = sortIndicesDeltaInverted(nodes);

  double delta_c=0;
  if(layer <= lastLayerEE)      delta_c = config->GetDeltac(kEE);
  else if(layer <= lastLayerFH) delta_c = config->GetDeltac(kFH);
  else                          delta_c = config->GetDeltac(kBH);

  for(uint i=0; i<nodes.size();i++){
    if(nodes[ds[i]]->delta < delta_c) break;  // no more cluster centers to be looked at
    // skip this as a potential cluster center because it fails the density cut
    if(dependSensor){
      if(nodes[ds[i]]->rho < kappa * nodes[ds[i]]->sigmaNoise){
        continue;  // set equal to kappa times noise threshold
      }
    }
    else{
      if(nodes[ds[i]]->rho < maxDensity / kappa) continue;
    }
    // store cluster index
    nodes[ds[i]]->clusterIndex = clusterIndex;

    if(verbosityLevel >= 2){
      cout<<"Adding new cluster with index "<<clusterIndex<<endl;
      cout<<"Cluster center is hit "<<&nodes[ds[i]]<<" with density rho: "<<nodes[ds[i]]->rho<<"and delta: "<< nodes[ds[i]]->delta<<endl;
    }
    clusterIndex++;
  }
  // at this point clusterIndex is equal to the number of cluster centers - if it is zero we are done
  if(clusterIndex == 0) return;

  for(int i=0;i<clusterIndex;i++){
    clusters.push_back(vector<unique_ptr<Hexel>>());
  }

  // assign to clusters, using the nearestHigher set from previous step (always set except for top density hit that is skipped)...
  for(uint oi=1;oi<nodes.size();oi++){
    int ci = nodes[rs[oi]]->clusterIndex;
    if(ci == -1) nodes[rs[oi]]->clusterIndex = nodes[nodes[rs[oi]]->nearestHigher]->clusterIndex;
  }

  // assign points closer than dc to other clusters to border region and find critical border density
  double rho_b[clusterIndex];
  for(int i=0;i<clusterIndex;i++){rho_b[i]=0.;}

  double criticalDistance = 0;
  
  if(layer <= lastLayerEE)       criticalDistance = config->GetCriticalDistance(kEE);
  else if(layer <= lastLayerFH)  criticalDistance = config->GetCriticalDistance(kFH);
  else                           criticalDistance = config->GetCriticalDistance(kBH);
  
  // now loop on all hits again :( and check: if there are hits from another cluster within d_c -> flag as border hit
  for(auto &iNode : nodes){
    int ci = iNode->clusterIndex;
    bool flag_isolated = true;
    if(ci != -1){// search in a circle of radius criticalDistance or criticalDistance*sqrt(2) (not identical to search in the box criticalDistance)

      auto found = queryBallPoint(points_0,points_1, iNode->x, iNode->y, criticalDistance);
      for(int j : found){
        // check if the hit is not within d_c of another cluster
        if(nodes[j]->clusterIndex != -1){
          double dist2 = distanceReal2(nodes[j]->x, nodes[j]->y , iNode->x, iNode->y);
          if(dist2 < criticalDistance * criticalDistance && nodes[j]->clusterIndex != ci){
            // in which case we assign it to the border
            iNode->isBorder = true;
            break;
          }
          // because we are using two different containers, we have to make sure that we don't unflag the
          // hit when it finds *itself* closer than criticalDistance

          if(dist2 < criticalDistance * criticalDistance && dist2 != 0. && nodes[j]->clusterIndex == ci){
            // this is not an isolated hit
            flag_isolated = false;
          }
        }
      }
      if(flag_isolated) iNode->isBorder = true;  // the hit is more than criticalDistance from any of its brethren
    }
    // check if this border hit has density larger than the current rho_b and update
    if(iNode->isBorder && rho_b[ci] < iNode->rho) // here rho_b[ci] is always zero - this seems wrong!!
      rho_b[ci] = iNode->rho;
  }
  // flag points in cluster with density < rho_b as halo points, then fill the cluster vector
  for(auto &iNode : nodes){
    int ci = iNode->clusterIndex;
    if(ci != -1 && iNode->rho <= rho_b[ci]){
      iNode->isHalo = true;  // some issues to be debugged?
    }
    if(ci != -1){
      clusters[ci].push_back(move(iNode));
      if(verbosityLevel >= 2){
        cout<<"Pushing hit "<<&iNode<<" into cluster with index "<<ci<<endl;
        cout<<"   rho_b[ci]: "<<rho_b[ci]<<", iNode.rho: "<<iNode->rho<<" iNode.isHalo: "<<iNode->isHalo<<endl;
      }
    }
  }
  return;
}

void ImagingAlgo::populate(vector<vector<unique_ptr<Hexel>>> &points, shared_ptr<RecHits> &hits)
{
  // init 2D hexels
  for(int iLayer=0;iLayer<2*(maxlayer+1);iLayer++){
    points.push_back(vector<unique_ptr<Hexel>>());
  }
  unique_ptr<RecHit> hit;
  // loop over all hits and create the Hexel structure, skip energies below ecut
  for(int iHit=0;iHit<hits->N();iHit++){
    
    if (hits->GetLayerOfHit(iHit) > maxlayer){
      continue;  // current protection
    }
    
    // energy treshold dependent on sensor
    auto thresholdResult = hits->RecHitAboveThreshold(iHit);
    if(!get<0>(thresholdResult)){
      continue;
    }
    // organise layers accoring to the sgn(z)
    hit = hits->GetHit(iHit);
    int layerID = hit->layer + (hit->z > 0) * (maxlayer + 1);  // +1 - yes or no?
    double sigmaNoise = get<1>(thresholdResult);
    unique_ptr<Hexel> hexel = hit->GetHexel();
    hexel->sigmaNoise = sigmaNoise;
    points[layerID].push_back(move(hexel));
  }
}

void ImagingAlgo::makeClusters(vector<vector<vector<unique_ptr<Hexel>>>> &clusters, shared_ptr<RecHits> &hits)
{
  // initialise list of per-layer-clusters
  for(int i=0;i<2 * (maxlayer + 1);i++){
    clusters.push_back(vector<vector<unique_ptr<Hexel>>>());
  }

  // get the list of Hexels out of raw rechits
  vector<vector<unique_ptr<Hexel>>> points;
  populate(points, hits);

  // loop over all layers, and for each layer create a list of clusters. layers are organised according to the sgn(z)
  for(int layerID=0;layerID<2*(maxlayer + 1);layerID++){
    if (points[layerID].size() == 0) continue;  // protection
    int layer = layerID - (points[layerID][0]->z > 0) * (maxlayer + 1);  // map back to actual layer

    vector<double> pointX; // list of hexels'coordinate 0 for current layer
    vector<double> pointY; // list of hexels'coordinate 1 for current layer

    for(auto &hex : points[layerID]){
      pointX.push_back(hex->x);
      pointY.push_back(hex->y);
    }

    double maxdensity = calculateLocalDensity(points[layerID],pointX,pointY, layer);
    calculateDistanceToHigher(points[layerID]);// get distances to the nearest higher density
    vector<vector<unique_ptr<Hexel>>> tmp;
    findAndAssignClusters(tmp, points[layerID],pointX,pointY, maxdensity, layer);  // get clusters per layer
    clusters[layerID] = move(tmp);
  }
}

void ImagingAlgo::getBasicClusters(vector<unique_ptr<BasicCluster>> &clustersFlat,
                              vector<vector<vector<unique_ptr<Hexel>>>> &clusters)
{
  // loop over all layers and all clusters in each layer
  for(vector<vector<unique_ptr<Hexel>>> &clustersInLayer : clusters){
    for(vector<unique_ptr<Hexel>> &cluster : clustersInLayer){
      auto position = calculatePosition(cluster);

      // skip the clusters where position could not be computed (either all weights are 0, or all hexels are tagged as Halo)
      if(get<0>(position)==0 && get<1>(position)==0 && get<2>(position)==0) continue;

      double energy = 0;
      vector<shared_ptr<Hexel>> sharedCluster;
      for(auto &iNode : cluster){
        if(!iNode->isHalo) energy += iNode->weight;
        sharedCluster.push_back(shared_ptr<Hexel>(new Hexel(*iNode)));
      }
      clustersFlat.push_back(unique_ptr<BasicCluster>(new BasicCluster(energy,get<0>(position),get<1>(position),get<2>(position), sharedCluster)));
    }

    sort(clustersFlat.begin( ), clustersFlat.end( ), [ ](const unique_ptr<BasicCluster> &lhs,const  unique_ptr<BasicCluster> &rhs){
      return lhs->GetEnergy() > rhs->GetEnergy();
    });
  }
}

void ImagingAlgo::getRecClusters(vector<shared_ptr<Hexel>> &hexelsClustered, shared_ptr<RecHits> &hits)
{
  // get 3D array of hexels (per layer, per 2D cluster)
  vector<vector<vector<std::unique_ptr<Hexel>>>> clusters2D;
  makeClusters(clusters2D, hits);
  
  // get flat list of 2D clusters (as basic clusters)
  std::vector<unique_ptr<BasicCluster>> clusters2Dflat;
  getBasicClusters(clusters2Dflat, clusters2D);
  
  // keep only non-halo hexels
  for(auto &basicCluster : clusters2Dflat){
    for(auto hexel : basicCluster->GetHexelsInThisCluster()){
      if(!hexel->isHalo){
        hexelsClustered.push_back(hexel);
      }
    }
  }
}

tuple<double,double,double> ImagingAlgo::calculatePosition(vector<unique_ptr<Hexel>> &cluster)
{
  double total_weight=0, x=0, y=0, z=0;
  bool haloOnlyCluster = true;

  // check if haloOnlyCluster
  for(auto &iNode : cluster){
    if(!iNode->isHalo) haloOnlyCluster = false;
  }
  if(!haloOnlyCluster){
    for(auto &iNode : cluster){
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
  for(auto &iNode : cluster){
    if(iNode->weight > maxenergy){
      maxenergy = iNode->weight;
      maxenergy_x = iNode->x;
      maxenergy_y = iNode->y;
      maxenergy_z = iNode->z;
    }
  }
  return make_tuple(maxenergy_x, maxenergy_y, maxenergy_z);
}
