# investigate shower development based on RecHits and SimClusters
import ROOT
import os
import optparse
from array import array
import math
import hgcalHelpers
import hgcalHistHelpers
import timeit

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


def runCalibrationScaleResolution(pidSelected, GEN_engpt, ntuple, relFractionE, gun_type, histDict):

    # common strings
    GEN_pTEng = "{}={0:.1f} GeV".format(gun_type, GEN_engpt)
    GEN_partId = pidmap[pidSelected]

    # init output stuff
    outDir = "simClusterHitsContained_ScaleResolution_" + GEN_partId + "s_" + "{0:.1f}GeV".format(GEN_engpt) + "_testSpatial"
    if not os.path.exists(outDir):
        os.makedirs(outDir)

    # define some global lists and dictionaries
    obj_Eng_EngRelDiff = {pid: [] for pid in s_all_pids}
    # define some global variables
    counterDivError = 0
    counterClustError = 0

    # initialisation of GeoUtils
    gu = GeoUtil()

    # loop over the events
    print "Total events to process (for multiClusters,", GEN_partId, ", ", GEN_pTEng, "): ", ntuple.nevents()
    hitsSelected = []
    # for event in ntuple:
    for event in ntuple:
        if (event.entry() > 3):
            break
        if (verbosityLevel >= 0):
            if (event.entry() % 1 == 0):
                print "Event: ", event.entry()
        # get collections
        genParticles = event.genParticles()
        simClusters = event.simClusters()
        # pfClusters = event.pfClusters()
        # multiClusters = event.multiClusters()
        recHits = event.recHits()
        # check for cases where number of Sim clusters is incorrect
        # if (len(simClusters)==0 or len(simClusters)>2): # skip if no sim clusters or skip if more then 2 sim clusters
        # if (verbosityLevel>=1): print "WARNING (Event:",iEvt,"): Sim cluster length is improper."
        # counterClustError+=1
        # continue

        # fill simhits based on rechits mapping and filtering ###
        # get rechit-simclus assoc.
        # energy sums
        # rHitsSimAssoc = getRecHitsSimAssoc(recHits, simClusters)
        rHitsSimAssoc = hgcalHelpers.getRecHitsSimAssocPUP(recHits, simClusters, genParticles, pidSelected, GEN_engpt)
        simClusters_energyAboveThreshold = []
        simClusters_energyAboveThresholdContained = []
        simClusters_energyAll = []

        # pt sums
        simClusters_ptAboveThreshold = []
        simClusters_ptAboveThresholdContained = []
        simClusters_ptAll = []
        # store the rechit info
        hitsSelected.extend([(10 * rechit.x(), 10 * rechit.y(), rechit.energy(), recHitAboveThreshold(rechit, ecut, dependSensor)[1]) for k in range(len(simClusters)) for rechit in rHitsSimAssoc[k]])
        # store in list the energy with filtering of simi hits
        for k in range(len(simClusters)):
            # values for energies
            sum_simHitsEnergyAll = 0
            sum_simHitsEnergyAbove = 0
            sum_simHitsEnergyContained = 0

        # values for pt
            sum_simHitsPtAll = 0
            sum_simHitsPtAbove = 0
            sum_simHitsPtContained = 0
            for rechitind, rechit in enumerate(rHitsSimAssoc[k]):
                # energy for rechits above noise
                sum_simHitsEnergyAll += rechit.energy() * simClusters[k].fractions()[rechitind]
                sum_simHitsPtAll += rechit.pt() * simClusters[k].fractions()[rechitind]
                # energy for rechits above noise
                if recHitAboveThreshold(rechit, ecut, dependSensor)[1]:
                    sum_simHitsEnergyAbove += rechit.energy() * simClusters[k].fractions()[rechitind]
                    sum_simHitsPtAbove += rechit.pt() * simClusters[k].fractions()[rechitind]
                # energy of rechits above noise and after masking
                if (recHitAboveThreshold(rechit, ecut, dependSensor)[1] and (rechit.layer() > 35 or (rechit.layer() <= 35 and gu.planes[rechit.layer() - 1].contains(10 * rechit.x(), 10 * rechit.y())))):
                    sum_simHitsEnergyContained += rechit.energy() * simClusters[k].fractions()[rechitind]
                    sum_simHitsPtContained += rechit.pt() * simClusters[k].fractions()[rechitind]
            simClusters_energyAll.append(sum_simHitsEnergyAll)
            simClusters_energyAboveThreshold.append(sum_simHitsEnergyAbove)
            simClusters_energyAboveThresholdContained.append(sum_simHitsEnergyContained)

            simClusters_ptAll.append(sum_simHitsPtAll)
            simClusters_ptAboveThreshold.append(sum_simHitsPtAbove)
            simClusters_ptAboveThresholdContained.append(sum_simHitsPtContained)

        # loop over the gen particles in current event
        collectionOfInterest = simClusters  # or collectionOfInterest = pfClusters # or collectionOfInterest = multiClusters or ...
        if (verbosityLevel >= 1):
            print "PDG ID", "\t", "Sim E(GeV)", "\t", "GEN E(GeV)", "\t", "Eng.Loss (%)", "\t", "DR(GEN, Sim)"
        for genParticle in genParticles:
            # sort indecies of PF clusters according to the distance of clusters to the GEN particles (index 0 is the nearest PF cluster)
            pfs = sorted(range(len(collectionOfInterest)), key=lambda k: math.sqrt(math.fabs(genParticle.eta() -
                                                                                             collectionOfInterest[k].eta())**2 + math.fabs(genParticle.phi() - collectionOfInterest[k].phi())**2), reverse=False)
            # check pathological cases
            if (len(collectionOfInterest) != len(simClusters)):
                if (verbosityLevel >= 1):
                    print "ERROR (Event:", event.entry(), "): multiClusters != simClusters."
                counterClustError += 1
                continue
            if (simClusters[pfs[0]].energy() == 0.):
                if (verbosityLevel >= 1):
                    print "ERROR (Event:", event.entry(), "): simClusters[pfs[0]].energy() == 0."
                counterDivError += 1
                continue
            # prepare some observables and store it for every eta/phi
            #    for simClusters in geometry: energy_to_calibrate = simClusters_energyAboveThresholdContained[pfs[0]]
            #    for simClusters noise: energy_to_calibrate = simClusters_energyAboveThreshold[pfs[0]]
            #    for multiClusters: energy_to_calibrate = multiClusters[pfs[0]].energy()
            #    for pfClusters: energy_to_calibrate = pfClusters[pfs[0]].correctedEnergy()
            energy_to_calibrate = simClusters_energyAboveThresholdContained[pfs[0]]
            energy_to_compare = simClusters_energyAboveThreshold[pfs[0]]  # or pfClusters[pfs[0]].correctedEnergy() or pfClusters[pfs[0]].energy() or ...
            simulated_energy = simClusters[pfs[0]].energy()

            pt_to_calibrate = simClusters_ptAboveThresholdContained[pfs[0]]
            pt_to_compare = simClusters_ptAboveThreshold[pfs[0]]  # or pfClusters[pfs[0]].correctedEnergy() or pfClusters[pfs[0]].energy() or ...
            simulated_pt = simClusters[pfs[0]].pt()

            DR_gen_sim = math.sqrt((genParticle.eta() - collectionOfInterest[pfs[0]].eta())**2 + (genParticle.phi() - collectionOfInterest[pfs[0]].phi())**2)
            if (gun_type == "e"):
                relDiff_sim = (energy_to_calibrate - simulated_energy) / simulated_energy  # relative diff. in energy between Sim clusters with/without filtering
                # print some info
                if (verbosityLevel >= 2):
                    print genParticle.pid(), "\t", energy_to_calibrate, "\t", genParticle.energy(), "\t", relDiff_sim, "\t", DR_gen_sim
                # store info in lists for cases where gen particles and Sim cluster are matched, and where reco. E for Sim cluster is close to the GEN-one
                # and simulated_energy >= (1-relFractionE)*GEN_engpt and simulated_energy <= (1+relFractionE)*GEN_engpt):
                if (math.fabs(DR_gen_sim) < 0.1 and energy_to_calibrate > simulated_energy / 20.):
                    if (abs(genParticle.pid()) in s_all_pids):
                        if (verbosityLevel == 1):
                            print genParticle.pid(), "\t", energy_to_calibrate, "\t", genParticle.energy(), "\t", relDiff_sim, "\t", DR_gen_sim
                        obj_Eng_EngRelDiff[abs(genParticle.pid())].append((simClusters[pfs[0]].phi(), simClusters[pfs[0]].eta(),
                                                                           simulated_energy, energy_to_compare, energy_to_calibrate, 100 * relDiff_sim))
            elif (gun_type == "pt"):
                relDiff_sim = (pt_to_calibrate - simulated_pt) / simulated_pt  # relative diff. in energy between Sim clusters with/without filtering
                # print some info
                if (verbosityLevel >= 2):
                    print genParticle.pid(), "\t", pt_to_calibrate, "\t", genParticle.pt(), "\t", relDiff_sim, "\t", DR_gen_sim
                # store info in lists for cases where gen particles and Sim cluster are matched, and where reco. E for Sim cluster is close to the GEN-one
                if (math.fabs(DR_gen_sim) < 0.1 and pt_to_calibrate > simulated_pt / 20.):  # and simulated_energy >= (1-relFractionE)*GEN_engpt and simulated_energy <= (1+relFractionE)*GEN_engpt):
                    if (abs(genParticle.pid()) in s_all_pids):
                        if (verbosityLevel == 1):
                            print genParticle.pid(), "\t", pt_to_calibrate, "\t", genParticle.pt(), "\t", relDiff_sim, "\t", DR_gen_sim
                        obj_Eng_EngRelDiff[abs(genParticle.pid())].append((simClusters[pfs[0]].phi(), simClusters[pfs[0]].eta(), simulated_pt, pt_to_compare, pt_to_calibrate, 100 * relDiff_sim))

    # prepare lists for individual eta/phi bins
    dEOverE_pidSelected = []
    dEOverE_pidSelected_values = {}
    for etaBinName in etaBins:
        GEN_eta = "[{0:.3f} - {1:.1f}]".format(etaBins[etaBinName][0], etaBins[etaBinName][1])
        print "Extracting info for eta range ", GEN_eta
        histDict[etaBinName] = {}
        dEOverE_pidSelected_values[etaBinName] = {}
        for phiBinName in phiBins:
            GEN_phi = "[{0:.2f} - {1:.2f}]".format(phiBins[phiBinName][0], phiBins[phiBinName][1])
            print "Extracting info for phi range ", GEN_phi
            histDict[etaBinName][phiBinName] = {}
            dEOverE_pidSelected_values[etaBinName][phiBinName] = {}
            # histograms per pid, print some info, fill dE/E values for current eta/phi bin
            print "PDG ID", "\t", "Mean dE/E (%)", "\t", "\t", "eta", "\t\t", "phi"
            for pid in s_all_pids:
                if (pid == 211 or pid == 22 or pid == 11 or pid == 13):
                    # get the 1D lists
                    # for energy of sim cluster
                    sim_Energy = [eng for (phi, eta, eng, engToCompare, engToCalibrate, relDiff) in obj_Eng_EngRelDiff[pid] if ((math.fabs(eta) >= etaBins[etaBinName]
                                                                                                                                 [0] and math.fabs(eta) < etaBins[etaBinName][1]) and (phi >= phiBins[phiBinName][0] and phi < phiBins[phiBinName][1]))]
                    # for energy of pf clusters before corrections
                    obj_EnergyToCompare = [engToCompare for (phi, eta, eng, engToCompare, engToCalibrate, relDiff) in obj_Eng_EngRelDiff[pid] if (
                        (math.fabs(eta) >= etaBins[etaBinName][0] and math.fabs(eta) < etaBins[etaBinName][1]) and (phi >= phiBins[phiBinName][0] and phi < phiBins[phiBinName][1]))]
                    # for energy of pf clusters after corrections
                    obj_EnergyToCalibrate = [engToCalibrate for (phi, eta, eng, engToCompare, engToCalibrate, relDiff) in obj_Eng_EngRelDiff[pid] if (
                        (math.fabs(eta) >= etaBins[etaBinName][0] and math.fabs(eta) < etaBins[etaBinName][1]) and (phi >= phiBins[phiBinName][0] and phi < phiBins[phiBinName][1]))]
                    # for energy loss
                    obj_EnergyDiffFilter = [relDiff for (phi, eta, eng, engToCompare, engToCalibrate, relDiff) in obj_Eng_EngRelDiff[pid] if (
                        (math.fabs(eta) >= etaBins[etaBinName][0] and math.fabs(eta) < etaBins[etaBinName][1]) and (phi >= phiBins[phiBinName][0] and phi < phiBins[phiBinName][1]))]
                    # fill the hists
                    rangeGeV = GEN_engpt * 1.6
                    nbins = int(rangeGeV)
                    if (rangeGeV < 50.):
                        nbins = int(10 * rangeGeV)
                    binsBoundariesX_eng = [nbins, 0, rangeGeV]
                    # for energy of sim cluster
                    histDict[etaBinName][phiBinName] = hgcalHistHelpers.histValue1D(sim_Energy, histDict[etaBinName][phiBinName], tag="sim_Energy_eta" + etaBinName + "_phi" + phiBinName + "_pid" + str(pid), title="Energy of the Sim cluster (matched to pid=" + str(
                        pid) + ", GEN-level " + GEN_partId + " " + GEN_pTEng + ", #eta=" + GEN_eta + ", #phi=" + GEN_phi + ")",   axunit="E[GeV]", binsBoundariesX=binsBoundariesX_eng, ayunit="N(clusters)")
                    # for energy of pf clusters before corrections
                    histDict[etaBinName][phiBinName] = hgcalHistHelpers.histValue1D(obj_EnergyToCompare, histDict[etaBinName][phiBinName], tag="obj_EnergyToCompare_eta" + etaBinName + "_phi" + phiBinName + "_pid" + str(
                        pid), title="Energy of the PF cluster before correction (cluster matched to pid=" + str(pid) + ", GEN-level " + GEN_partId + " " + GEN_pTEng + ", #eta=" + GEN_eta + ", #phi=" + GEN_phi + ")",   axunit="E[GeV]", binsBoundariesX=binsBoundariesX_eng, ayunit="N(clusters)")
                    # for energy of pf clusters after corrections
                    histDict[etaBinName][phiBinName] = hgcalHistHelpers.histValue1D(obj_EnergyToCalibrate, histDict[etaBinName][phiBinName], tag="obj_EnergyToCalibrate_eta" + etaBinName + "_phi" + phiBinName + "_pid" + str(
                        pid), title="Energy of the PF cluster after correction (cluster matched to pid=" + str(pid) + ", GEN-level " + GEN_partId + " " + GEN_pTEng + ", #eta=" + GEN_eta + ", #phi=" + GEN_phi + ")",   axunit="E[GeV]", binsBoundariesX=binsBoundariesX_eng, ayunit="N(clusters)")
                    # for energy loss from filter
                    binsBoundariesX_engDiff = [[800, -100, 60], [650, -80, 50]]["1p" in etaBinName]
                    histDict[etaBinName][phiBinName] = hgcalHistHelpers.histValue1D(obj_EnergyDiffFilter, histDict[etaBinName][phiBinName], tag="obj_EnergyDiffFilter_eta" + etaBinName + "_phi" + phiBinName + "_pid" + str(pid), title="Energy loss from PF hits filtering (cluster matched to pid=" + str(
                        pid) + ", GEN-level " + GEN_partId + " " + GEN_pTEng + ", #eta=" + GEN_eta + ", #phi=" + GEN_phi + ")",   axunit="#DeltaE_{clust}/E_{clust}[%]", binsBoundariesX=binsBoundariesX_engDiff, ayunit="N(clusters)")
                    # extract mean and std for eng diff
                    h_tempEdiff = histDict[etaBinName][phiBinName]["obj_EnergyDiffFilter_eta" + etaBinName + "_phi" + phiBinName + "_pid" + str(pid)]
                    hEntries = h_tempEdiff.GetEntries()
                    gMean = h_tempEdiff.GetMean()
                    gMeanError = h_tempEdiff.GetMeanError()
                    gStd = h_tempEdiff.GetRMS()
                    if (hEntries > 100):  # extract mean/error from fit if enough statistics
                        (h_tempEdiff, gMean, gStd) = hgcalHistHelpers.fitGauss(h_tempEdiff)
                        gMeanError = gStd / (hEntries**0.5)
                    effSigma = hgcalHistHelpers.getEffSigma(h_tempEdiff)
                    # print table for rel.diff.
                    print pid, "\t", "{0:.2f}".format(gMean) + " +/- " + "{0:.2f}".format(gMeanError), "\t", GEN_eta, "\t", GEN_phi, "\t", "(pf clusters: " + str(int(hEntries)) + ")", "\t", "effSigma: ", effSigma, "\t", "gStd: ", gStd
                    # for pidSelected only: overlap hists together, fill list with phi, eta and dE/E values
                    if (pid == pidSelected):
                        # overlap three hists together
                        h_tempE = histDict[etaBinName][phiBinName]["sim_Energy_eta" + etaBinName + "_phi" + phiBinName + "_pid" + str(pid)]
                        h_tempEToCalibrate = histDict[etaBinName][phiBinName]["obj_EnergyToCalibrate_eta" + etaBinName + "_phi" + phiBinName + "_pid" + str(pid)]
                        h_tempEToCompare = histDict[etaBinName][phiBinName]["obj_EnergyToCompare_eta" + etaBinName + "_phi" + phiBinName + "_pid" + str(pid)]

                        # h_tempE.Write()
                        # h_tempEToCalibrate.Write()
                        # h_tempEToCompare.Write()

                        gMeanCalibEnergy = (gMean / 100. + 1.) * GEN_engpt
                        gMeanCalibEnergyError = (gMeanError / 100.) * GEN_engpt
                        # prepare histograms and their properties. this should depend on the input option which energies shpould be extracted/calibrated/compared
                        # this is for PF clusters
#                        histsAndProps = {h_tempE:{"leg":"Sim cluster ("+GEN_pTEng+")","color":ROOT.kRed}, h_tempEToCalibrate:{"leg":"Corrected PF cluster","color":ROOT.kBlue}, h_tempEToCompare:{"leg":"Non-corrected PF cluster","color":ROOT.kGreen-6}}
                        # this is for PF multiclusters
#                        histsAndProps = {h_tempE:{"leg":"Sim cluster ("+GEN_pTEng+")","color":ROOT.kRed}, h_tempEToCalibrate:{"leg":"Multicluster (uncorrected)","color":ROOT.kBlue}, h_tempEToCompare:{"leg":"PF cluster (corrected)","color":ROOT.kGreen-6}}
                        # this is for SIM clusters (3 sigma, 5 sigma, realistic, etc.)
#                        histsAndProps = {h_tempE:{"leg":"Sim cluster ("+GEN_pTEng+")","color":ROOT.kRed}, h_tempEToCalibrate:{"leg":"SIM cluster (detids matched to rechits, 5#sigma above noise)","color":ROOT.kBlue}, h_tempEToCompare:{"leg":"PF cluster (corrected)","color":ROOT.kGreen-6}}
                        # this is for different geometries
                        histsAndProps = {h_tempE: {"leg": "Sim cluster (" + GEN_pTEng + ")", "color": ROOT.kRed}, h_tempEToCalibrate: {"leg": "SIM cluster (3#sigma noise filtering, realistic Si geometry)",
                                                                                                                                       "color": ROOT.kBlue}, h_tempEToCompare: {"leg": "SIM cluster (3#sigma above noise filtering, idealistic/CMSSW Si geometry)", "color": ROOT.kGreen - 6}}
                        # plot these histograms on top of each other (for each eta/phi bin)
                        hgcalHistHelpers.histsPrintSaveSameCanvas(histsAndProps, outDir, tag="obj_EnergyHistsOverlapped_eta" + etaBinName + "_phi" + phiBinName + "_pid" + str(pid), latexComment=["E_{fit} = " + "{0:.2f}".format(
                            gMeanCalibEnergy) + " #pm " + "{0:.2f}".format(gMeanCalibEnergyError) + " GeV", "#sigma_{fit} = " + "{0:.1f}".format(gStd) + " %, #sigma_{eff} = " + "{0:.1f}".format(effSigma) + " %"], funcsAndProps=None)
                        # fill list with phi, eta and dE/E values
                        # skip filling the lists for overll hists in case of repeating inclusive bins
                        if not (len(etaBins.keys()) > 1 and etaBinName == "eta1p479to3p0") and not (len(phiBins.keys()) > 1 and phiBinName == "phim1p0pito1p0pi"):
                            dEOverE_pidSelected.append((sum(phiBins[phiBinName]) / 2., sum(etaBins[etaBinName]) / 2., gMean))
                        dEOverE_pidSelected_values[etaBinName][phiBinName]['mean'] = gMean
                        dEOverE_pidSelected_values[etaBinName][phiBinName]['meanError'] = gMeanError
                        dEOverE_pidSelected_values[etaBinName][phiBinName]['effSigma'] = effSigma
            # plos histograms for all pids in current eta/phi bin (if needed)
            # histPrintSaveAll(histDict[etaBinName][phiBinName], outDir, "_eta"+etaBinName+"_phi"+phiBinName)

    # global 2D histogram
    histDictGlob = {}
    # get edge boudaries for eta and phi
    bns_phi = list(set([item for sublist in [[phiBins[phiBin][0], phiBins[phiBin][1]] for phiBin in phiBins] for item in sublist]))
    bns_phi.sort(key=lambda x: x)
    bns_eta = list(set([item for sublist in [[etaBins[etaBin][0], etaBins[etaBin][1]] for etaBin in etaBins] for item in sublist]))
    bns_eta.sort(key=lambda x: x)
    # plot histograms for 3D list (phi, eta, relDiff)
    # simClust_EnergyDiffFilter_FullpidSelected = [(phi, eta, relDiff) for (phi, eta, eng, engToCalibrate, engAbove, relDiff) in simClust_Eng_EngRelDiff[pidSelected]]
    # histDictGlob = histValues2D(simClust_EnergyDiffFilter_FullpidSelected,
    # histDictGlob, tag = "dEOverE_pidSelected_weigthed2D_full", title =
    # "Relative energy loss #DeltaE_{clust}/E_{clust} from hits filtering
    # (%)", axunit = "#phi", binsBoundariesX = [len(bns_phi)-1, bns_phi],
    # ayunit = "#eta", binsBoundariesY = [len(bns_eta)-1, bns_eta], weighted2D
    # = True)
    histDictGlob = hgcalHistHelpers.histValues2D(dEOverE_pidSelected, histDictGlob, tag="dEOverE_pidSelected_weigthed2D", title="Relative energy loss #DeltaE_{clust}/E_{clust} from hits filtering (%)", axunit="#phi", binsBoundariesX=[
        len(bns_phi) - 1, bns_phi], ayunit="#eta", binsBoundariesY=[len(bns_eta) - 1, bns_eta], weighted2D=True)
    recHits_Filtered = [(x, y, int(filt)) for (x, y, eng, filt) in hitsSelected]
    histDictGlob = hgcalHistHelpers.histValues2D(recHits_Filtered, histDictGlob, tag="recHits_filteredByNoise", title="Spatial distribution of hits filtered by noise",
                                                 axunit="x[mm]", binsBoundariesX=[3200, -1600, 1600], ayunit="y[mm]", binsBoundariesY=[3200, -1600, 1600], weighted2D=True)
    hgcalHistHelpers.histPrintSaveAll(histDictGlob, outDir, "_glob")

    # print 2D eta-phi table
    print "Effect on the cluster energy (relative loss in %):"
    print "eta\phi", "\t", "\t".join(["[{0:.2f} - {1:.2f}]".format(phiBins[phiBinName][0], phiBins[phiBinName][1]) for phiBinName in phiBins])
    for etaBinName in etaBins:
        print "[{0:.3f} - {1:.1f}]".format(etaBins[etaBinName][0], etaBins[etaBinName][1]), "\t", "\t".join(["{0:.2f}".format(dEOverE_pidSelected_values[etaBinName][phiBinName]['mean']) + " +/- " + "{0:.2f}".format(dEOverE_pidSelected_values[etaBinName][phiBinName]['meanError']) for phiBinName in phiBins])

    # print 2D eta-phi table
    print "Effect on the cluster energy (sigma effective in %):"
    print "eta\phi", "\t", "\t".join(["[{0:.2f} - {1:.2f}]".format(phiBins[phiBinName][0], phiBins[phiBinName][1]) for phiBinName in phiBins])
    for etaBinName in etaBins:
        print "[{0:.3f} - {1:.1f}]".format(etaBins[etaBinName][0], etaBins[etaBinName][1]), "\t", "\t".join(["{0:.2f}".format(dEOverE_pidSelected_values[etaBinName][phiBinName]['effSigma']) for phiBinName in phiBins])

    # print 2D eta-phi table
    print "Effect on the cluster energy (value of calibration const.):"
    print "eta\phi", "\t", "\t".join(["[{0:.2f} - {1:.2f}]".format(phiBins[phiBinName][0], phiBins[phiBinName][1]) for phiBinName in phiBins])
    for etaBinName in etaBins:
        print "[{0:.3f} - {1:.1f}]".format(etaBins[etaBinName][0], etaBins[etaBinName][1]), "\t", "\t".join(["{0:.2f}".format(1 / (1 + dEOverE_pidSelected_values[etaBinName][phiBinName]['mean'] / 100.)) + " -/+ " + "{0:.2f}".format((1 / ((1 + dEOverE_pidSelected_values[etaBinName][phiBinName]['mean'] / 100.)**2)) * dEOverE_pidSelected_values[etaBinName][phiBinName]['meanError']) for phiBinName in phiBins])


# plot graphs with energy resolution for various scenarios (currently needs to be invoked manually, one all the energy points have been processed)
def graphTest():
    graphsAndProps = {}
    xEng = array('f', [5., 20., 50., 100., 300.])

# photons, ideal vs. realistic geometry (matched/filtred simHits)
    # eta [1.479 - 3.0] pfClusters
    yRes = array('f', [9.80, 5.40, 3.60, 2.80, 2.00])
    graphsAndProps[ROOT.TGraph(len(xEng), xEng, yRes)] = {"leg": "PF clusters (realistic SIM clusters, ideal geometry)", "color": ROOT.kBlue, "MarkerStyle": 21, "LineStyle": 4}
    # eta [1.479 - 3.0] matched/filtred simHits, 3 sigma, ideal geometry
    yRes = array('f', [9.90, 5.40, 3.50, 2.60, 1.50])
    graphsAndProps[ROOT.TGraph(len(xEng), xEng, yRes)] = {"leg": "SIM clusters (noise filtered @3#sigma, ideal geometry)", "color": ROOT.kGreen - 6, "MarkerStyle": 21, "LineStyle": 1}
    # eta [1.479 - 3.0] matched/filtred simHits, 3 sigma, realistic geometry
    yRes = array('f', [11.50, 6.30, 4.20, 3.10, 1.90])
    graphsAndProps[ROOT.TGraph(len(xEng), xEng, yRes)] = {"leg": "SIM clusters (noise filtered @3#sigma, realistic geometry)", "color": ROOT.kRed - 6, "MarkerStyle": 26, "LineStyle": 1}
    hgcalHistHelpers.drawGraphsTest(graphsAndProps, title="Energy resolution (photons, #eta #in [1.479 - 3.0], PU=0)", tag="test_photonsResolutionComparisonGeo")


def main():

    global opt, args

    usage = ('usage: %prog [options]\n' + '%prog -h for help')
    parser = optparse.OptionParser(usage)

    # input options
    parser.add_option('', '--gunType', dest='gunType', type='string',  default='pt', help='pt or e')
    parser.add_option('', '--pid', dest='pid', type='int',  default=22, help='pdgId int')
    parser.add_option('', '--genValue', dest='genValue', type='float',  default=25, help='generated pT or energy')
    parser.add_option('', '--tag', dest='tag', type='string',  default='noPU', help='some tag, best used for PU and other info')
    parser.add_option('', '--ref', dest='refName', type='string',  default='genpart', help='reference collection')
    parser.add_option('', '--obj', dest='objName', type='string',  default='pfcluster', help='object of interest collection')

    # store options and arguments as global variables
    global opt, args
    (opt, args) = parser.parse_args()

    print "gunType:", opt.gunType
    print "pid:", opt.pid
    print "GEN_engpt:", opt.genValue
    print "tag:", opt.tag
    print "refName:", opt.refName
    print "objName:", opt.objName

    # set sample/tree - for photons
    gun_type = opt.gunType
    pidSelected = opt.pid
    GEN_engpt = opt.genValue
    tag = opt.tag
    refName = opt.refName
    objName = opt.objName

    start_time = timeit.default_timer()

    histList = ["ref_Energy", "ref_Pt", "obj_Energy", "obj_Pt"]
    resolutionHistList = ["obj_dEoverE", "obj_dPtoverPt"]

    f = ROOT.TFile.Open("{}_{}_{}GeV_{}_{}_{}.root".format(gun_type, pidSelected, GEN_engpt, refName, objName, tag), "read")
    for etaBinName in etaBins:
        for phiBinName in phiBins:
            for histPrefix in histList:
                histo = f.Get(histPrefix + "_eta" + etaBinName + "_phi" + phiBinName)
                # histsAndProps = {histo: {"leg": "Sim cluster (" + GEN_pTEng + ")", "color": ROOT.kRed}, h_tempEToCalibrate: {"leg": "SIM cluster (3#sigma noise filtering, realistic Si geometry)", "color": ROOT.kBlue}, h_tempEToCompare: {"leg": "SIM cluster (3#sigma above noise filtering, idealistic/CMSSW Si geometry)", "color": ROOT.kGreen - 6}}
                # # plot these histograms on top of each other (for each eta/phi bin)
                # hgcalHistHelpers.histsPrintSaveSameCanvas(histsAndProps, outDir, tag="obj_EnergyHistsOverlapped_eta" + etaBinName + "_phi" + phiBinName + "_pid" + str(pid), latexComment=["E_{fit} = " + "{0:.2f}".format(
                #     gMeanCalibEnergy) + " #pm " + "{0:.2f}".format(gMeanCalibEnergyError) + " GeV", "#sigma_{fit} = " + "{0:.1f}".format(gStd) + " %, #sigma_{eff} = " + "{0:.1f}".format(effSigma) + " %"], funcsAndProps=None)
            for histPrefix in resolutionHistList:
                histo = f.Get(histPrefix + "_eta" + etaBinName + "_phi" + phiBinName)
                # do the fits etc.

    f.Close()
    elapsed = timeit.default_timer() - start_time
    print "Time:", elapsed


if __name__ == '__main__':
    main()
