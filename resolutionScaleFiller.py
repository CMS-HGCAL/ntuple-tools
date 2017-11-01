# investigate shower development based on RecHits and SimClusters
import ROOT
import os
import optparse
# from array import array
# from HGCalImagingAlgo import recHitAboveThreshold
from NtupleDataFormat import HGCalNtuple
# from GeoUtils import GeoUtil
import math
import hgcalHelpers
import hgcalHistHelpers
import numpy as np
import timeit

# filtering parameters
dependSensor = True
ecut = 3  # relative to the noise
# verbosity etc.
verbosityLevel = 0  # 0 - only basic info (default); 1 - additional info; 2 - detailed info printed, histograms produced

# basic settings
# names and pid mapping
pidmap = {11: "electron", 13: "muon", 22: "photon", 211: "pion"}

etaBins = {"eta1p479to1p6": (1.479, 1.6), "eta1p6to1p8": (1.6, 1.8), "eta1p8to2p0": (1.8, 2.0), "eta2p0to2p2": (2.0, 2.2), "eta2p2to2p4": (
    2.2, 2.4), "eta2p4to2p6": (2.4, 2.6), "eta2p6to2p8": (2.6, 2.8), "eta2p8to3p0": (2.8, 3.0), "eta1p479to3p0": (1.479, 3.0), "eta1p6to2p8": (1.6, 2.8)}
# phiBins = {"phi0to0p5pi":(0.*math.pi, 0.5*math.pi), "phi0p5to1p0pi":(0.5*math.pi, 1.0*math.pi), "phim1p0pitom0p5pi":(-1.0*math.pi, -0.5*math.pi), "phim0p5pito0":(-0.5*math.pi, 0.*math.pi),"phim1p0pito1p0pi":(-1.0*math.pi, 1.0*math.pi) }

# these are to run only inclusive bins
# etaBins = {"eta1p479to3p0":(1.479, 3.0)}
phiBins = {"phim1p0pito1p0pi": (-1.0 * math.pi, 1.0 * math.pi)}

deltaRMaxRef = 0.1
relativeFractionRef = 0.05


class ResolutionScaleObject:
    """Class to simplify typical values to calculate for scale and resolution"""

    __slots__ = ['refTLV', 'objTLV']

    def __init__(self, ref, obj):
        self.refTLV = ROOT.TLorentzVector()
        self.refTLV.SetPtEtaPhiE(ref.pt, ref.eta, ref.phi, ref.energy)
        self.objTLV = ROOT.TLorentzVector()
        self.objTLV.SetPtEtaPhiE(obj.pt, obj.eta, obj.phi, obj.energy)

    def deltaR(self):
        return self.refTLV.DeltaR(self.objTLV)

    def deltaEta(self):
        return self.refTLV.Eta() - self.objTLV.Eta()

    def deltaPhi(self):
        return self.refTLV.DeltaPhi(self.objTLV)

    def scale(self):
        return

    def dEoverE(self):
        return (self.objTLV.E() - self.refTLV.E()) / self.refTLV.E()

    def dPtoverPt(self):
        return (self.objTLV.Pt() - self.refTLV.Pt()) / self.refTLV.Pt()


def eventLoop(ntuple, refName, objName, gun_type, pidOfInterest, GEN_engpt, histDict):
    """
    Loop over ntuple,
    for the collection of interest, match with genPart to select relevant objects,
    then pass selected objects to calculate scale and resolution
    """
    # common strings
    GEN_pTEng = "{0}={1:.1f} GeV".format(gun_type, GEN_engpt)
    GEN_partId = pidmap[pidOfInterest]

    # define some global lists and dictionaries
    # obj_Eng_EngRelDiff = {pid: [] for pid in s_all_pids}
    resolutionScaleObjects = []

    # initialisation of GeoUtils
    # gu = GeoUtil()

    # loop over the events
    print "Total events to process (PID:", GEN_partId, ",", GEN_pTEng, "):", ntuple.nevents()
    # for event in ntuple:
    for event in ntuple:
        if (event.entry() > 10):
            break
        if (verbosityLevel >= 0):
            if (event.entry() % 1 == 0):
                print "Event: ", event.entry()
        # get collections
        # genParticles = event.genParticles()
        genParticles = event.getDataFrame(prefix=refName)
        # pfClusters = event.pfClusters()
        pfClusters = event.getDataFrame(prefix=objName)

        # let's look at pfClusters compared to genParticles
        collectionOfInterest = pfClusters
        referenceCollection = genParticles
        # print "collections:", len(collectionOfInterest), len(referenceCollection)
        # filter reference Collection for faster matching
        if (gun_type == "e"):
            referenceCollection = filterReferenceCollection(referenceCollection, pidOfInterest, refMinE=GEN_engpt*.9)
        else:
            referenceCollection = filterReferenceCollection(referenceCollection, pidOfInterest, refMinPt=GEN_engpt*.9)

        pairs = getReferencePairs(referenceCollection, collectionOfInterest)
        resolutionScaleObjects += getResolutionScaleObjects(pairs)
    fillComparisonHistograms(resolutionScaleObjects, GEN_engpt, histDict)


def filterReferenceCollection(referenceCollection, pidOfInterest, refMinPt=0, refMinE=0):
    """cut on reference pdgId and minPt or minE
    probably need a new function for other collections"""

    referencePIDselected = referenceCollection[(abs(referenceCollection.pid) == pidOfInterest)]
    if refMinPt > 0:
        # referencePairs = [x for x in referencePairs if x[0].pt() > refMinPt]
        referencePIDselected = referencePIDselected[referencePIDselected.pt > refMinPt]
    if refMinE > 0:
        # referencePairs = [x for x in referencePairs if x[0].energy() > refMinE]
        referencePIDselected = referencePIDselected[referencePIDselected.energy > refMinE]
    return referencePIDselected


def getReferencePairs(referenceCollection, collectionOfInterest):
    """
    - match collectionOfInterest with closest reference (including DeltaR cut)
    - return list of pairs
    """

    # start_time = timeit.default_timer()
    referencePair = []
    # for reference in referenceCollection.itertuples():
    # # for reference in referenceCollection:
    #     # first check for relevant pdgId
    #     # if abs(reference.pid()) != pidOfInterest:
    #     if abs(reference.pid) != pidOfInterest:
    #         continue
    #     # sort indices of collectionOfInterest according to the distance of
    #     # clusters to the genParticles (index 0 is the nearest collectionOfInterest object)
    #     nearestObjectIndex = min(range(len(collectionOfInterest)),
    #                              key=lambda k:
    #                              hgcalHelpers.deltaRSquared(reference, collectionOfInterest[k]))
    #     nearestObject = collectionOfInterest[nearestObjectIndex]
    #
    #     # require nearest object to be reasonably matched
    #     # if (hgcalHelpers.deltaR(reference, nearestObject) < deltaRMaxRef and nearestObject.energy() > reference.energy() * relativeFractionRef):
    #     if (hgcalHelpers.deltaR(reference, nearestObject) < deltaRMaxRef and nearestObject.energy() > reference.energy * relativeFractionRef):
    #         # store the pair
    #         referencePair.append((reference, nearestObject))
    matched_indices = hgcalHelpers.getClosestObjectIndices(referenceCollection[['eta', 'phi']], collectionOfInterest[['eta', 'phi']], deltaR=deltaRMaxRef)
    for idx1, idx2 in matched_indices.iteritems():
        if collectionOfInterest.iloc[idx2].energy > referenceCollection.iloc[idx1].energy * relativeFractionRef:
            referencePair.append((referenceCollection.iloc[idx1], collectionOfInterest.iloc[idx2]))
    # elapsed = timeit.default_timer() - start_time
    # print "Time:", elapsed
    return referencePair


def getResolutionScaleObjects(referencePairs):
    """
    referencePairs: (reference, collection of interest)
    """
    objectsForHists = []
    for ref, obj in referencePairs:
        resolutionScaleObject = ResolutionScaleObject(ref, obj)
        if verbosityLevel > 0:
            relE = resolutionScaleObject.dEoverE()
            deltaR = resolutionScaleObject.deltaR()
            deltaEta = resolutionScaleObject.deltaEta()
            deltaPhi = resolutionScaleObject.deltaPhi()
            print "relative energy:", relE, "reference energy:", ref.energy, "reference pT:", ref.pt
            print "deltaR:", deltaR, "deltaEta:", deltaEta, "deltaPhi:", deltaPhi
        objectsForHists.append(resolutionScaleObject)
    return objectsForHists


def fillComparisonHistograms(resolutionScaleObjects, GEN_engpt, histDict):
    """fill lists from resolutionScaleObjects, then histograms"""

    for etaBinName in etaBins:
        GEN_eta = "[{0:.3f} - {1:.1f}]".format(etaBins[etaBinName][0], etaBins[etaBinName][1])
        if verbosityLevel > 0:
            print "Extracting info for eta range ", GEN_eta
        histDict[etaBinName] = {}
        for phiBinName in phiBins:
            GEN_phi = "[{0:.2f} - {1:.2f}]".format(phiBins[phiBinName][0], phiBins[phiBinName][1])
            # print "Extracting info for phi range ", GEN_phi
            histDict[etaBinName][phiBinName] = {}
            # print some info, fill dE/E values for current eta/phi bin
            # print "Mean dE/E (%)", "\t", "\t", "eta", "\t\t", "phi"
            # get the 1D lists
            valueLists = np.array([(x.refTLV.E(), x.refTLV.Pt(), x.objTLV.E(), x.objTLV.Pt(), x.dEoverE(), x.dPtoverPt())
                          for x in resolutionScaleObjects
                          if ((math.fabs(x.refTLV.Eta()) >= etaBins[etaBinName][0] and math.fabs(x.refTLV.Eta()) < etaBins[etaBinName][1])
                          and (x.refTLV.Phi() >= phiBins[phiBinName][0] and x.refTLV.Phi() < phiBins[phiBinName][1]))], dtype=float)
            if len(valueLists) > 0 and verbosityLevel > 0:
                print valueLists[:, 1]
            # fill the hists
            if len(valueLists) > 0:
                rangeGeV = GEN_engpt * 1.6
                nbins = int(rangeGeV)
                if (rangeGeV < 50.):
                    nbins = int(10 * rangeGeV)
                binsBoundariesX_eng = [nbins, 0, rangeGeV]
                # reference energy
                histDict[etaBinName][phiBinName] = hgcalHistHelpers.histValue1D(valueLists[:, 0], histDict[etaBinName][phiBinName], tag="ref_Energy_eta" + etaBinName + "_phi" + phiBinName, title="Energy of reference object, #eta=" + GEN_eta + ", #phi=" + GEN_phi + ")",   axunit="E [GeV]", binsBoundariesX=binsBoundariesX_eng, ayunit="N(clusters)")
                # reference Pt
                histDict[etaBinName][phiBinName] = hgcalHistHelpers.histValue1D(valueLists[:, 1], histDict[etaBinName][phiBinName], tag="ref_Pt_eta" + etaBinName + "_phi" + phiBinName, title="Pt of reference object, #eta=" + GEN_eta + ", #phi=" + GEN_phi + ")",   axunit="p_{T} [GeV]", binsBoundariesX=binsBoundariesX_eng, ayunit="N(clusters)")
                # object of interest energy
                histDict[etaBinName][phiBinName] = hgcalHistHelpers.histValue1D(valueLists[:, 2], histDict[etaBinName][phiBinName], tag="obj_Energy_eta" + etaBinName + "_phi" + phiBinName, title="Energy of object of interest, #eta=" + GEN_eta + ", #phi=" + GEN_phi + ")",   axunit="E [GeV]", binsBoundariesX=binsBoundariesX_eng, ayunit="N(clusters)")
                # object of interest Pt
                histDict[etaBinName][phiBinName] = hgcalHistHelpers.histValue1D(valueLists[:, 3], histDict[etaBinName][phiBinName], tag="obj_Pt_eta" + etaBinName + "_phi" + phiBinName, title="Pt of object of interest, #eta=" + GEN_eta + ", #phi=" + GEN_phi + ")",   axunit="p_{T} [GeV]", binsBoundariesX=binsBoundariesX_eng, ayunit="N(clusters)")
                # resolution
                binsBoundariesX_engDiff = [[800, -100, 60], [650, -80, 50]]["1p" in etaBinName]
                histDict[etaBinName][phiBinName] = hgcalHistHelpers.histValue1D(valueLists[:, 4]*100, histDict[etaBinName][phiBinName], tag="obj_dEoverE_eta" + etaBinName + "_phi" + phiBinName, title="dEoverE, #eta=" + GEN_eta + ", #phi=" + GEN_phi + ")", axunit="#Delta E_{clust}/E_{clust}[%]", binsBoundariesX=binsBoundariesX_engDiff, ayunit="N(clusters)")
                histDict[etaBinName][phiBinName] = hgcalHistHelpers.histValue1D(valueLists[:, 5]*100, histDict[etaBinName][phiBinName], tag="obj_dPtoverPt_eta" + etaBinName + "_phi" + phiBinName, title="dPtoverPt, #eta=" + GEN_eta + ", #phi=" + GEN_phi + ")", axunit="#Delta p_{T, clust}/p_{T, clust}[%]", binsBoundariesX=binsBoundariesX_engDiff, ayunit="N(clusters)")


def main():

    global opt, args

    usage = ('usage: %prog [options]\n' + '%prog -h for help')
    parser = optparse.OptionParser(usage)

    # input options
    parser.add_option('', '--files', dest='fileString', type='string',  default='root://eoscms.cern.ch//eos/cms/store/cmst3/group/hgcal/CMG_studies/Production/_RelValSingleGammaPt25Eta1p7_2p7_CMSSW_9_3_2-PU25ns_93X_upgrade2023_realistic_v2_2023D17PU200-v1_GEN-SIM-RECO/NTUP/_RelValSingleGammaPt25Eta1p7_2p7_CMSSW_9_3_2-PU25ns_93X_upgrade2023_realistic_v2_2023D17PU200-v1_GEN-SIM-RECO_NTUP_1.root', help='comma-separated file list')
    parser.add_option('', '--gunType', dest='gunType', type='string',  default='pt', help='pt or e')
    parser.add_option('', '--pid', dest='pid', type='int',  default=22, help='pdgId int')
    parser.add_option('', '--genValue', dest='genValue', type='int',  default=25, help='generated pT or energy')
    parser.add_option('', '--tag', dest='tag', type='string',  default='noPU', help='some tag, best used for PU and other info')
    parser.add_option('', '--ref', dest='refName', type='string',  default='genpart', help='reference collection')
    parser.add_option('', '--obj', dest='objName', type='string',  default='pfcluster', help='object of interest collection')

    # store options and arguments as global variables
    global opt, args
    (opt, args) = parser.parse_args()

    print "files:", opt.fileString
    print "gunType:", opt.gunType
    print "pid:", opt.pid
    print "GEN_engpt:", opt.genValue
    print "refName:", opt.refName
    print "objName:", opt.objName

    # set sample/tree - for photons
    gun_type = opt.gunType
    pidSelected = opt.pid
    GEN_engpt = opt.genValue
    tag = opt.tag
    refName = opt.refName
    objName = opt.objName

    histDict = {}
    fileList = opt.fileString.split(",")

    start_time = timeit.default_timer()

    for fileName in fileList:
        ntuple = HGCalNtuple(opt.fileString)
        eventLoop(ntuple, refName, objName, gun_type, pidSelected, GEN_engpt, histDict)

    f = ROOT.TFile("{}_{}_{}GeV_{}_{}_{}.root".format(gun_type, pidSelected, GEN_engpt, refName, objName, tag), "recreate")
    for etaBinName in etaBins:
        for phiBinName in phiBins:
            if "ref_Energy_eta" + etaBinName + "_phi" + phiBinName in histDict[etaBinName][phiBinName]:
                histDict[etaBinName][phiBinName]["ref_Energy_eta" + etaBinName + "_phi" + phiBinName].Write()
                histDict[etaBinName][phiBinName]["ref_Pt_eta" + etaBinName + "_phi" + phiBinName].Write()
                histDict[etaBinName][phiBinName]["obj_Energy_eta" + etaBinName + "_phi" + phiBinName].Write()
                histDict[etaBinName][phiBinName]["obj_Pt_eta" + etaBinName + "_phi" + phiBinName].Write()
                histDict[etaBinName][phiBinName]["obj_dEoverE_eta" + etaBinName + "_phi" + phiBinName].Write()
                histDict[etaBinName][phiBinName]["obj_dPtoverPt_eta" + etaBinName + "_phi" + phiBinName].Write()

    f.Write()
    f.Close()
    elapsed = timeit.default_timer() - start_time
    print "Time:", elapsed


if __name__ == '__main__':
    main()
