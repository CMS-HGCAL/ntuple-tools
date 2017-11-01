# investigate shower development based on RecHits and SimClusters
# import ROOT
# import os
import optparse
# from array import array
# from HGCalImagingAlgo import recHitAboveThreshold
from NtupleDataFormat import HGCalNtuple
# from GeoUtils import GeoUtil
# import math
import hgcalHelpers
# import hgcalHistHelpers
# import numpy as np
import pandas as pd
from itertools import repeat
maxlayer = 52
energyWeights = list(repeat(1.02, 28)) + list(repeat(0.86, 12)) + list(repeat(1.12, 12))

def getConeRadius(frontRadius, backRadius, z, maxval=9999.):
    depthTerm = backRadius * (abs(z)-320.7)/(407.8-320.7)
    val = frontRadius + depthTerm
    if val > maxval:
        return maxval
    return val


def getMegaClusters(genParticles, multiClusters, layerClusters, recHits, gun_type, GEN_engpt, pidSelected, energyRadius=6, frontRadius=3, backRadius=8):
    """
    get the actual mega clusters.
    frontRadius: cone at front of EE
    backRadius: cone at back of FH (frontRadius to be added to it)
    returns a dataframe containing 4-vectors
    """

    # use genParticles with generated pT/energy that reach EE before converting
    selectedGen = genParticles[(abs(genParticles.pid) == pidSelected) & (genParticles.reachedEE > 0)]
    if gun_type == "pt":
        selectedGen = selectedGen[(selectedGen.pt >= GEN_engpt*.999)]
    else:
        selectedGen = selectedGen[(selectedGen.energy >= GEN_engpt*.999)]
    # print selectedGen

    # for the mega cluster axis, take highest energy multicluster within dR = 0.1
    bestMultiClusterIndices = hgcalHelpers.getHighestEnergyObjectIndex(selectedGen[['eta', 'phi']], multiClusters[['eta', 'phi']], multiClusters['energy'], 0.1)
    # print bestMultiClusterIndices

    megaClusters = []

    for idx, genPart in selectedGen.iterrows():
        matchedMultiCluster = multiClusters.iloc[[bestMultiClusterIndices[idx]]]
        energySum = 0
        pTSum = 0
        # now find layer clusters within the multicluster cone
        # maybe it's good to do this per layer to save some computing time
        for layer in range(1, maxlayer+1):
            # match only in same detector side
            selectedLayerClusters = layerClusters[(layerClusters.layer == layer) & (layerClusters.eta*genPart.eta > 0)]
            if selectedLayerClusters.shape[0] == 0:
                # continue of no layer clusters selected
                continue
            # take first layer cluster z value
            layer_z = selectedLayerClusters.head(1).z.item()
            # get multi cluster x and y coordinates
            multiClusPosDF = hgcalHelpers.convertToXY(matchedMultiCluster.eta, matchedMultiCluster.phi, layer_z)
            # calculate radius based on current layer's z position
            coneRadius = getConeRadius(frontRadius, backRadius, layer_z)
            # mind that we need only the first index since there is only one multiCluster
            layerClusterIndices = hgcalHelpers.getIndicesWithinRadius(multiClusPosDF[['x', 'y']], selectedLayerClusters[['x', 'y']], coneRadius)
            # now we need to recalculate the layer cluster energies using associated RecHits
            for layerClusterIndex in layerClusterIndices[0]:
                associatedRecHits = recHits.iloc[selectedLayerClusters.iloc[layerClusterIndex].rechits]
                # find maximum energy RecHit
                maxEnergyRecHitIndex = associatedRecHits['energy'].argmax()
                # considering only associated RecHits within a radius of energyRadius (6 cm)
                matchedRecHitIndices = hgcalHelpers.getIndicesWithinRadius(associatedRecHits.loc[[maxEnergyRecHitIndex]][['x', 'y']], associatedRecHits[['x', 'y']], energyRadius)[maxEnergyRecHitIndex]
                # sum up energies and pT
                selectedRecHits = associatedRecHits.iloc[matchedRecHitIndices]
                # correct energy by subdetector weights
                energySum += selectedRecHits[["energy"]].sum()[0]*energyWeights[layer-1]
                pTSum += selectedRecHits[["pt"]].sum()[0]*energyWeights[layer-1]

        # use as coordinates eta and phi of matched multi cluster
        megaCluster = [pTSum, matchedMultiCluster['eta'].item(), matchedMultiCluster['phi'].item(), energySum]
        megaClusters.append(megaCluster)

    megaClustersDF = pd.DataFrame(megaClusters, columns=['pt', 'eta', 'phi', 'energy'])

    return megaClustersDF


def getCollections(event):
    """
    get the collections to be fed to the mega clustering.
    need genParticles, multiClusters, layerClusters, recHits
    """
    genParticles = event.getDataFrame(prefix="genpart")
    multiClusters = event.getDataFrame(prefix="multiclus")
    layerClusters = event.getDataFrame(prefix="cluster2d")
    recHits = event.getDataFrame(prefix="rechit")
    return genParticles, multiClusters, layerClusters, recHits


def main():
    global opt, args

    usage = ('usage: %prog [options]\n' + '%prog -h for help')
    parser = optparse.OptionParser(usage)

    # input options
    # parser.add_option('', '--files', dest='fileString', type='string',  default='root://eoscms.cern.ch//eos/cms/store/cmst3/group/hgcal/CMG_studies/Production/_SinglePiPt50Eta1p6_2p8_PhaseIITDRFall17DR-noPUFEVT_93X_upgrade2023_realistic_v2-v1_GEN-SIM-RECO/NTUP/_SinglePiPt50Eta1p6_2p8_PhaseIITDRFall17DR-noPUFEVT_93X_upgrade2023_realistic_v2-v1_GEN-SIM-RECO_NTUP_1.root', help='comma-separated file list')
    parser.add_option('', '--files', dest='fileString', type='string',  default='/afs/cern.ch/work/e/escott/public/HGCStudies/Ntuples/partGun_Pion_Pt25_93X_PU140.root', help='comma-separated file list')
    parser.add_option('', '--gunType', dest='gunType', type='string',  default='pt', help='pt or e')
    parser.add_option('', '--pid', dest='pid', type='int',  default=211, help='pdgId int')
    parser.add_option('', '--genValue', dest='genValue', type='float',  default=25, help='generated pT or energy')

    # store options and arguments as global variables
    global opt, args
    (opt, args) = parser.parse_args()

    print "files:", opt.fileString
    print "gunType:", opt.gunType
    print "pid:", opt.pid
    print "GEN_engpt:", opt.genValue

    # set sample/tree - for photons
    gun_type = opt.gunType
    pidSelected = opt.pid
    GEN_engpt = opt.genValue

    fileList = opt.fileString.split(",")

    for fileName in fileList:
        ntuple = HGCalNtuple(opt.fileString)

        for event in ntuple:
            if (event.entry() > 11):
                break
            # get collections
            genParticles, multiClusters, layerClusters, recHits = getCollections(event)
            megaClusters = getMegaClusters(genParticles, multiClusters, layerClusters, recHits, gun_type, GEN_engpt, pidSelected)
            print megaClusters


if __name__ == '__main__':
    main()
