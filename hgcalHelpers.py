import numpy as np
from HGCalImagingAlgo import recHitAboveThreshold
import math
from scipy.spatial import cKDTree


def getRecHitDetIds(rechits):
    recHitsList = []
    for rHit in rechits:
        recHitsList.append(rHit.detid())
    # print "RecHits -"*10
    # print recHitsList
    recHits = np.array(recHitsList)
    return recHits


def getHitList(simClus, recHitDetIds):
    sClusHitsList = []
    for DetId in simClus.hits():
        sClusHitsList.append(DetId)
    sClusHits = np.array(sClusHitsList)
    # thanks to http://stackoverflow.com/questions/11483863/python-intersection-indices-numpy-array
    recHitIndices = np.nonzero(np.in1d(recHitDetIds, sClusHits))
    return recHitIndices

# get list of rechist associated to sim-cluster hits


def getRecHitsSimAssoc(rechits_raw, simcluster, dependSensor=True, ecut=3, verbosityLevel=0):
    # get sim-cluster associations
    nSimClus = 0
    simClusHitAssoc = []
    recHitDetIds = getRecHitDetIds(rechits_raw)
    for simClusIndex, simClus in enumerate(simcluster):
        simClusHitAssoc.append(getHitList(simClus, recHitDetIds))
        nSimClus += 1
    # get list of rechist associated to simhits
    rHitsSimAssoc = [[] for k in range(0, nSimClus)]
    for simClusIndex, simCl in enumerate(simcluster):
        if (verbosityLevel >= 1):
            print "Sim-cluster index: ", simClusIndex, ", pT: ", simCl.pt(), ", E: ", simCl.energy(), ", phi: ", simCl.phi(), ", eta: ", simCl.eta()
        # loop over sim clusters and then rechits
        rHitsSimAssocTemp = []
        for hitIndexArray in simClusHitAssoc[simClusIndex]:
            for hitIndex in hitIndexArray:
                thisHit = rechits_raw[hitIndex]
                if(not recHitAboveThreshold(thisHit, ecut, dependSensor)[1]):
                    continue
                # independent of sim cluster, after cleaning
                rHitsSimAssocTemp.append(thisHit)
        rHitsSimAssoc[simClusIndex] = rHitsSimAssocTemp
    return rHitsSimAssoc

# get list of rechist associated to sim-cluster which is within dR<0.2 of the genParticle


def deltaRSquared(obj1, obj2):
    # return (obj1.eta() - obj2.eta())**2 + (obj1.phi() - obj2.phi())**2
    return (obj1.eta - obj2.eta())**2 + (obj1.phi - obj2.phi())**2


def deltaR(obj1, obj2):
    return math.sqrt(deltaRSquared(obj1, obj2))


def getRecHitsSimAssocPUP(rechits_raw, simcluster, genparticles, pidSelected, GEN_engpt):
    # get sim-cluster associations
    # nSimClus = 0
    # simClusHitAssoc = []
    recHitDetIds = getRecHitDetIds(rechits_raw)
    rHitsSimAssoc = [[] for k in range(0, len(simcluster))]
    for simClusIndex, simClus in enumerate(simcluster):
        for genPartIndex, genPart in enumerate(genparticles):
            if (genPart.pid() == pidSelected and math.sqrt(math.fabs(genPart.eta() - simClus.eta())**2 + math.fabs(genPart.phi() - simClus.phi())**2) < 0.1 and simClus.pt() > 0.7 * GEN_engpt):
                rHitsSimAssocTemp = []
                for hitIndexArray in getHitList(simClus, recHitDetIds):
                    for hitIndex in hitIndexArray:
                        thisHit = rechits_raw[hitIndex]
                        # if(not recHitAboveThreshold(thisHit, ecut, dependSensor)[1]): continue
                        rHitsSimAssocTemp.append(thisHit)
                rHitsSimAssoc[simClusIndex] = rHitsSimAssocTemp
    return rHitsSimAssoc


def getClosestObjectIndices(ref_etaphi, obj_etaphi, deltaR=0.1):
    '''Match object with smallest DeltaR'''
    kdtree = cKDTree(obj_etaphi)
    matched_indices = {}

    for index, row in ref_etaphi.iterrows():
        # matched = kdtree.query_ball_point([row.eta, row.phi], deltaR)
        closest = kdtree.query([row.eta, row.phi], 1)
        # Handle the -pi pi transition, do the matching again
        # matched_sym = kdtree.query_ball_point([row.eta, row.phi-np.sign(row.phi)*2.*math.pi], deltaR)
        closest_sym = kdtree.query([row.eta, row.phi-np.sign(row.phi)*2.*math.pi], 1)
        # keep only unique indices
        matched = [closest, closest_sym]
        # Choose the match with minimum deltaR
        best_match = min(matched, key=lambda k: k[0])
        # print best_match
        if (best_match[0] < deltaR):
            matched_indices[index] = best_match[1]
    return matched_indices
