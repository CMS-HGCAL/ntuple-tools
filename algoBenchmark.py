from __future__ import print_function
import ROOT
import os
import time
import numpy as np
import pandas as pd
from HGCalImagingAlgo import *
from NtupleDataFormat import HGCalNtuple
from RecHitCalibration import RecHitCalibration

#----------------------------------------------------------------------------------------
# HGCal Imaging Algo parameters:
dependSensor = True
deltac = [2., 2., 5.] # in cartesian coordiantes in cm, per detector
multiclusterRadii = [2., 5., 5.] # in cartesian coordiantes in cm, per detector
minClusters = 3 # request at least minClusters+1 2D clusters

# cut on energy (also passed to HGCalImagingAlgo):
energyMin = 3 # relative to the noise

# other cuts
clusterAcceptScale = 1.0 # for E_sim vs. E_rec, accept only sim clusters with clusterAcceptScale*clusterRadius

# test only within this layers range:
minLayer=0
maxLayer=40

# range of ntuples to test (will be appended to the inputPath string below):
minNtuple = 11
maxNtuple = 11

# base input and output paths:
inputPath = "../data/_SingleGammaPt100Eta1p6_2p8_PhaseIITDRFall17DR-noPUFEVT_93X_upgrade2023_realistic_v2-v1_GEN-SIM-RECO/NTUP/_SingleGammaPt100Eta1p6_2p8_PhaseIITDRFall17DR-noPUFEVT_93X_upgrade2023_realistic_v2-v1_GEN-SIM-RECO_NTUP_"

outDir = "clusteringResults"
#----------------------------------------------------------------------------------------


# get pandas mask of hits in given layer
def getLayerMask(hits,layer):
  return hits["layer"]==layer

# get pandas mask of hits above noice threshold
def getMaskAboveNoice(hits,ecut):
  sigmaNoise = 1.
  thickIndex = -1
  mask=[]
  RecHitCalib = RecHitCalibration()
  for hit in hits.itertuples():
    layer=hit.layer
    thickness=hit.thickness
    if(layer <= 40):  # EE + FH
      if  (thickness > 99.  and thickness < 101.): thickIndex = 0
      elif(thickness > 199. and thickness < 201.): thickIndex = 1
      elif(thickness > 299. and thickness < 301.): thickIndex = 2
      else: print("ERROR - silicon thickness has a nonsensical value")
      # determine noise for each sensor/subdetector using RecHitCalibration library
    sigmaNoise = 0.001 * RecHitCalib.sigmaNoiseMeV(layer, thickIndex)  # returns threshold for EE, FH, BH (in case of BH thickIndex does not play a role)
    mask.append(hit.energy >= ecut * sigmaNoise)  # checks if rechit energy is above the threshold of ecut (times the sigma noise for the sensor, if that option is set)
  return pd.Series(mask)

# groups hits into array of clusters
def getHitsPerCluster(hits, clusters):
  
  noiceMask = getMaskAboveNoice(hits,energyMin)
  hitsAboveNoice = hits[noiceMask]

  hitsDetIDs = hitsAboveNoice["detid"]
  hitsPerCluster = []
  
  for cluster in clusters.itertuples():
    clusterToHitID = np.nonzero(np.in1d(hitsDetIDs, cluster.hits))
    hitsInThisCluster = hitsAboveNoice.iloc[clusterToHitID[0].tolist()]
    hitsPerCluster.append(hitsInThisCluster)

  return hitsPerCluster

# groups hits associated with hexels into array of clusters
def getRecHitsPerHexel(hits, hexels):
  noiceMask = getMaskAboveNoice(hits,energyMin)
  hitsAboveNoice = hits[noiceMask]

  clusterIndices = []
  hexelDetIDs = []

  for hexel in hexels:
    if hexel.clusterIndex not in clusterIndices:
      clusterIndices.append(hexel.clusterIndex)
    hexelDetIDs.append(hexel.detid)

  hitDetIDs = hitsAboveNoice["detid"]
  hexelToHitID = np.nonzero(np.in1d(hitDetIDs, hexelDetIDs))

  hitsClustered = []

  for clusterIndex in range(0,np.max(clusterIndices)+1):
    hitsClustered.append(pd.DataFrame(data=None, columns=hitsAboveNoice.columns, index=hitsAboveNoice.index))

  for hexelIndex, hexel in enumerate(hexels):
    hitIndex = hexelToHitID[0][hexelIndex]
    hitsClustered[hexel.clusterIndex].loc[hitsAboveNoice.index[hitIndex]] = hitsAboveNoice.iloc[hitIndex]

  for clusterIndex in range(0,np.max(clusterIndices)+1):
    hitsClustered[clusterIndex].dropna(subset=["eta"], inplace=True)

  return hitsClustered

# get clustered hexels by re-running the clustering algorithm
def getRecClustersFromImagingAlgo(recHitsRaw):
  HGCalAlgo = HGCalImagingAlgo(energyMin, deltac, multiclusterRadii, minClusters, dependSensor, verbosityLevel = 0)
  clusters2D_rerun = HGCalAlgo.makeClusters(recHitsRaw,energyMin,True) # nested list of "hexels", per layer, per 2D cluster
  clusters2DList_rerun = HGCalAlgo.getClusters(clusters2D_rerun, verbosityLevel = 0) # flat list of 2D clusters (as basic clusters)
  hexelsClustered_rerun = [iNode for bClust in clusters2DList_rerun for iNode in bClust.thisCluster if not iNode.isHalo]  # flat list of clustered "hexeles", without the "halo" hexels
  return hexelsClustered_rerun

# check is two circles overlap
def circlesOverlap(x1,y1,r1,x2,y2,r2,scale=1.0):
  return (scale*(r1+r2))**2 >= (x1-x2)**2 + (y1-y2)**2

# check if point is within a circle
def pointWithinCircle(px,py,x,y,r,scale=1.0):
  return (scale*r)**2 >= (px-x)**2 + (py-y)**2


def main():
  if not os.path.exists(outDir): os.makedirs(outDir)
  
  for ntupleNumber in range(minNtuple,maxNtuple+1):
    print("\nCurrent ntup: ", ntupleNumber)

    ntuple = HGCalNtuple(inputPath+"{}.root".format(ntupleNumber));

    # start event loop
    for event in ntuple:
      startEvent = time.time()
      eventID = event.entry()
      startEvent = time.time()
      
      print("\nCurrent event: ", eventID)

      # check if particles reached EE
      genParticles = event.genParticles()
      skipEvent = False
      for particle in genParticles:
        if not particle.reachedEE():
#          print("particle didn't reach EE -- skipping the event!!")
          skipEvent = True
          break
      if skipEvent: continue
    
      eventDir = outDir+"/ntup{}/event{}".format(ntupleNumber,eventID)
      if not os.path.exists(eventDir): os.makedirs(eventDir)
    
      # get raw rec hits
      print("\n\npreparing raw recHits...",end='')
      start = time.time()
      recHitsRaw = event.getDataFrame("rechit")
      end = time.time()
      print(" done (",end-start," s)")


      # get simulated hits associated with a cluster
      print("preparing simulated hits and clusters...",end='')
      start = time.time()
      simClusters = event.getDataFrame("simcluster")
      simHitsPerClusterArray = getHitsPerCluster(recHitsRaw, simClusters)
      end = time.time()
      print(" done (",end-start," s)")


      # re-run clustering with HGCalAlgo, save to file
      print("running clustering algorithm...",end='')
      start = time.time()
      recClusters = getRecClustersFromImagingAlgo(recHitsRaw)
      end = time.time()
      print(" done (",end-start," s)")
      
      
      # recClusters -> array of hexel objects
      print("looking for hits associated with hexels...",end='')
      start = time.time()
      recHitsPerClusterArray = getRecHitsPerHexel(recHitsRaw, recClusters)
      end = time.time()
      print(" done (",end-start," s)")


      # perform final analysis, fill in histograms and save to files
      print("\nGenerating final hists...")
      start = time.time()
      energyComparisonHist = ROOT.TH2D("energy comparison","energy comparison",100,0,100,100,0,100)
      energyComparisonOverlapHist = ROOT.TH2D("energy comparison overlap.","energy comparison overlap.",100,0,100,100,0,100)

      for layer in range(minLayer,maxLayer):
#        print("layer:",layer)
        for recClusterIndex, recCluster in enumerate(recHitsPerClusterArray):
#          print("rec cluster:",recCluster)

          recHitsInLayerInCluster = recCluster[getLayerMask(recCluster,layer)]
      
          recEnergy = recHitsInLayerInCluster["energy"].sum()
          xMaxRec   = recHitsInLayerInCluster["x"].max()
          xMinRec   = recHitsInLayerInCluster["x"].min()
          yMaxRec   = recHitsInLayerInCluster["y"].max()
          yMinRec   = recHitsInLayerInCluster["y"].min()
          
          recClusterX = xMinRec+(xMaxRec-xMinRec)/2.
          recClusterY = yMinRec+(yMaxRec-yMinRec)/2.
          recClusterR = max((xMaxRec-xMinRec)/2.,(yMaxRec-yMinRec)/2.)
        
          assocSimEnergy = 0

          for simClusterIndex, simCluster in enumerate(simHitsPerClusterArray):
#            print("sim cluster:",simCluster)

            simHitsInLayerInCluster = simCluster[getLayerMask(simCluster,layer)]
          
            simEnergy = simHitsInLayerInCluster["energy"].sum()
            xMaxSim   = simHitsInLayerInCluster["x"].max()
            xMinSim   = simHitsInLayerInCluster["x"].min()
            yMaxSim   = simHitsInLayerInCluster["y"].max()
            yMinSim   = simHitsInLayerInCluster["y"].min()
          
            simClusterX = xMinSim+(xMaxSim-xMinSim)/2.
            simClusterY = yMinSim+(yMaxSim-yMinSim)/2.
            simClusterR = max((xMaxSim-xMinSim)/2.,(yMaxSim-yMinSim)/2.)
          
            if recEnergy*simEnergy != 0:  
              energyComparisonHist.Fill(recEnergy,simEnergy)
#              if circlesOverlap(recClusterX,recClusterY,recClusterR,simClusterX,simClusterY,simClusterR):
#                energyComparisonOverlapHist.Fill(recEnergy,simEnergy)

            if pointWithinCircle(simClusterX,simClusterY,recClusterX,recClusterY,recClusterR,clusterAcceptScale):
#            if circlesOverlap(recClusterX,recClusterY,recClusterR,simClusterX,simClusterY,simClusterR,clusterAcceptScale):
              assocSimEnergy += simEnergy
            
          if recEnergy*assocSimEnergy != 0:
            energyComparisonOverlapHist.Fill(recEnergy,assocSimEnergy)

      energyComparisonHist.SaveAs("{}/energyComparisonHist.root".format(eventDir))
      energyComparisonOverlapHist.SaveAs("{}/energyComparisonOverlapHist.root".format(eventDir))
      end = time.time()
      print(" done (",end-start," s)")

      endEvent = time.time()
      print("Total event processing time: ",endEvent-startEvent," s")

if __name__ == '__main__':
  main()









