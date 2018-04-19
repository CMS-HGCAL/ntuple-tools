# Translate the HGCAL N-tuples to the format used for standalone calibration (SAC) by Pedro
# import ROOT
# import os
import optparse
from NtupleDataFormat import HGCalNtuple
import hgcalHelpers
import numpy as np
import pandas as pd
from itertools import repeat
from collections import defaultdict
from ROOT import TVector2,TMath

class MatchedGenHits:
        """A wrapper with the same format as the trees used for calibration in the standalone setup"""
        def __init__(self,genId,en,pt,eta,phi,si_sumen=[],sci_sumen=[]):
                self.genId=genId
                self.genEta=eta
                self.genEt=pt
                self.genEn=en
                self.genPhi=phi
                self.si_sumen=[x for x in si_sumen]
                self.sci_sumen=[x for x in sci_sumen]

class SACevent:
        """Parses the reco-ntuples, associates the recHits to the genParticles, stores info as used in the standalone calibration setup"""
	def __init__(self, event,nlayers):
		self.hgcEvent = event
                self.nlayers=nlayers
                self.matchedHitSums = self.getMatchedRecHitSums()

	def isNoise(self, iRecHit,recHits,genParticles):
		iMatch=self.isMatched(recHits.eta[iRecHit],recHits.phi[iRecHit],genParticles,xrange(0,len(genParticles)))
                return True if iMatch<0 else False

	def isMatched(self, eta,phi, genParticles, gpIdx=[],minDR = 0.5):
		belong = -1
                for iGen in gpIdx:
                        geta,gphi=genParticles.eta[iGen],genParticles.phi[iGen]
                        deta=geta-eta
                        dphi=TVector2.Phi_mpi_pi(gphi-phi)
                        dR=TMath.Sqrt(deta**2+dphi**2)
                        if dR<minDR:
                                belong=iGen
                                break
                return belong

	def getMatchedRecHitSums(self):
                matchedHitSums=[]
		genParticles = self.hgcEvent.getDataFrame(prefix="genpart")
		recHits      = self.hgcEvent.getDataFrame(prefix="rechit")

                #filter gen particles for those reaching endcap
                #and not suffering nuclear interactions in the tracker
                gpIdx=[]
		for iGen in range(0, len(genParticles)):
                        if genParticles.reachedEE[iGen]==0: continue
                        if genParticles.gen[iGen]<0 : continue                        
                        gpIdx.append(iGen)
                if len(gpIdx)==0 : return matchedHitSums
                        
                #associate the rec hits to each selected genParticle by deltaR
                si_sumen=defaultdict(lambda: [0.] * self.nlayers)
                sci_sumen=defaultdict(lambda: [0.] * self.nlayers)
                for iRecHit in range (0,len(recHits)):
                        iMatch = self.isMatched(recHits.eta[iRecHit],recHits.phi[iRecHit],genParticles,gpIdx)
                        if iMatch<0 : continue
                        isSi = True if int(recHits.thickness[iRecHit]) in [100,200,300] else False
                        layerIdx=recHits.layer[iRecHit]-1
                        if layerIdx>=self.nlayers: continue 
                        if isSi:
                                si_sumen[iMatch][layerIdx]+=recHits.energy[iRecHit]
                        else:
                                sci_sumen[iMatch][layerIdx]+=recHits.energy[iRecHit]

                #finalize list of matched hit sums
		for iGen in gpIdx:
                        genId=genParticles.pid[iGen]
                        en=genParticles.energy[iGen]
                        pt=genParticles.pt[iGen]
                        eta=genParticles.eta[iGen]
                        phi=genParticles.phi[iGen]
                        matchedHitSums.append(  MatchedGenHits(genId,en,pt,eta,phi,si_sumen[iGen],sci_sumen[iGen]) )
                return matchedHitSums

	def Print(self):
                for h in self.matchedHitSums:
			print "gen id = %d, energy = %f, Et = %f, eta = %f, phi = %f" %(h.genId, h.genEn, h.genEt, h.genEta, h.genPhi)
			for l in range(0,len(h.si_sumen)):
				print "Layer number %d: Si energy sum = %f, Sci energy sum = %f" %(l,h.si_sumen[l],h.sci_sumen[l])



def main():
    global opt, args

    usage = ('usage: %prog [options]\n' + '%prog -h for help')
    parser = optparse.OptionParser(usage)

    # input options
    parser.add_option('', '--files', dest='fileString', type='string',  default='/eos/cms/store/cmst3/group/hgcal/CMG_studies/Production/FlatRandomPtGunProducer_SinglePiPt2Eta1p6_2p8_Fall17DR-NoPUFEVT_clange_20180129/NTUP/partGun_PDGid211_x60_Pt2.0To2.0_NTUP_6.root', help='comma-separated file list')
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
            SACEvt = SACevent(event,60)
            SACEvt.Print()


if __name__ == '__main__':
    main()
