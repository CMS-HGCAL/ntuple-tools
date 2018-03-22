# Translate the HGCAL N-tuples to the format used for standalone calibration (SAC) by Pedro
# import ROOT
# import os
import optparse
from NtupleDataFormat import HGCalNtuple
import hgcalHelpers
import numpy as np
import pandas as pd
from itertools import repeat
from ROOT import TLorentzVector

class SACevent:
	def __init__(self, event):
		self.hgcEvent = event
	def getGenInfo(self):
		if hasattr(self, "genParticles"):
			return self.genId, self.genEn, self.genEt, self.genEta, self.genPhi
		self.genParticles = self.hgcEvent.getDataFrame(prefix="genpart")
		self.genId = self.genParticles.pid
		self.genEta = self.genParticles.eta
		self.genEt = self.genParticles.pt
		self.genEn = self.genParticles.energy
		self.genPhi = self.genParticles.phi
		return self.genId, self.genEn, self.genEt, self.genEta, self.genPhi

	def isNoise(self, iRecHit):
		myrechit = TLorentzVector()
		IsNoise = True
		myrechit.SetPtEtaPhiE(self.recHits.pt[iRecHit],self.recHits.eta[iRecHit],self.recHits.phi[iRecHit], self.recHits.energy[iRecHit])
		for iGen in range(0, len(self.genParticles)):
			if isMatched(iRecHit, iGen, 0.5):
				IsNoise = False
				break
		return IsNoise

	def isMatched(self, iRecHit, iGen, DR = 0.5):
		myrechit = TLorentzVector()
		belong_ = False
		myrechit.SetPtEtaPhiE(self.recHits.pt[iRecHit],self.recHits.eta[iRecHit],self.recHits.phi[iRecHit], self.recHits.energy[iRecHit])
		mygen = TLorentzVector()
		mygen.SetPtEtaPhiE(self.genParticles.pt[iGen],self.genParticles.eta[iGen],self.genParticles.phi[iGen], self.genParticles.energy[iGen])
		if myrechit.DeltaR(mygen) < DR:
			belong_ = True
		return belong_

	def getRecHitInfo(self):
		if not hasattr(self,"genParticles"):
			self.getGenInfo()
		if hasattr(self, "recHits"):
			return self.array_si_sim_sumen, self.array_sci_sim_sumen, self.array_layers
		self.recHits = self.hgcEvent.getDataFrame(prefix="rechit")
		self.array_layers = []
		self.array_si_sim_sumen = []
		self.array_sci_sim_sumen = []
		for iGen in range(0, len(self.genParticles)):
			layers = []
			si_sim_sumen = []
			sci_sim_sumen = []
			for iRecHit in range (0,len(self.recHits)):
				if not self.isMatched(iRecHit,iGen):
					continue
				isSi = (int(self.recHits.thickness[iRecHit]) == 100 or int(self.recHits.thickness[iRecHit]) == 200 or int(self.recHits.thickness[iRecHit]) == 300)
				if not self.recHits.layer[iRecHit] in layers:
					layers.append(self.recHits.layer[iRecHit])
					si_sim_sumen.append(0)
					sci_sim_sumen.append(0)
				index = layers.index(self.recHits.layer[iRecHit])
				if isSi:
					si_sim_sumen[index]+=self.recHits.energy[iRecHit]
				else:
					sci_sim_sumen[index]+=self.recHits.energy[iRecHit]
			self.array_si_sim_sumen.append(si_sim_sumen)
			self.array_sci_sim_sumen.append(sci_sim_sumen)	
			self.array_layers.append(layers)
		return self.array_si_sim_sumen, self.array_sci_sim_sumen,self.array_layers

	def Print(self):
		self.getGenInfo()
		self.getRecHitInfo()
		for gen in range(0,len(self.genId)):
			print "gen id = %d, energy = %f, Et = %f, eta = %f, phi = %f" %(self.genId[gen], self.genEn[gen], self.genEt[gen], self.genEta[gen], self.genPhi[gen])
			for l in range(0,len(self.array_layers[gen])):
				print "Layer number %d: Si energy sum = %f, Sci energy sum = %f" %(self.array_layers[gen][l], self.array_si_sim_sumen[gen][l], self.array_sci_sim_sumen[gen][l])



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
            SACEvt = SACevent(event)
            SACEvt.Print()

            break

if __name__ == '__main__':
    main()
