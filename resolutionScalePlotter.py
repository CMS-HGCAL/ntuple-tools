# investigate shower development based on RecHits and SimClusters
import ROOT
import os
import sys
import optparse
from array import array
import math
# import hgcalHelpers
import hgcalHistHelpers
import timeit

# verbosity etc.
verbosityLevel = 0  # 0 - only basic info (default); 1 - additional info; 2 - detailed info printed, histograms produced

# basic settings
# names and pid mapping
pidmap = {11: "electron", 13: "muon", 22: "photon", 211: "pion"}

etaBins = {"eta1p6to1p8": (1.6, 1.8), "eta1p8to2p0": (1.8, 2.0), "eta2p0to2p2": (2.0, 2.2), "eta2p2to2p4": (2.2, 2.4), "eta2p4to2p6": (2.4, 2.6), "eta2p6to2p8": (2.6, 2.8), "eta1p6to2p8": (1.6, 2.8)}
# phiBins = {"phi0to0p5pi":(0.*math.pi, 0.5*math.pi), "phi0p5to1p0pi":(0.5*math.pi, 1.0*math.pi), "phim1p0pitom0p5pi":(-1.0*math.pi, -0.5*math.pi), "phim0p5pito0":(-0.5*math.pi, 0.*math.pi),"phim1p0pito1p0pi":(-1.0*math.pi, 1.0*math.pi) }

# these are to run only inclusive bins
# etaBins = {"eta1p6to2p8": (1.6, 2.8)}
phiBins = {"phim1p0pito1p0pi": (-1.0 * math.pi, 1.0 * math.pi)}


# format entries for the 2D eta-phi tables
def tableEntriesFormating(type='mean', resScale_valuesPerBin={}):
    if (type == 'mean'):
        entryValue = "{0:.2f}".format(resScale_valuesPerBin['mean']) + " +/- " + "{0:.2f}".format(resScale_valuesPerBin['meanError'])
    if (type == 'effSigma'):
        entryValue = "{0:.2f}".format(resScale_valuesPerBin['effSigma'])
    if (type == 'calib'):
        entryValue = "{0:.2f}".format(1 / (1 + resScale_valuesPerBin['mean'] / 100.)) + " -/+ " + "{0:.2f}".format((1 /
                                                                                                                    ((1 + resScale_valuesPerBin['mean'] / 100.)**2)) * resScale_valuesPerBin['meanError'])
    return entryValue


# print 2D eta-phi tables
def printEtaPhiTable(resScale_values, type='mean'):
    print {'mean': "Effect on the cluster energy (relative loss in %):", 'effSigma': "Effect on the cluster energy (sigma effective in %):", 'calib': "Effect on the cluster energy (value of calibration const.):"}[type]
    print "eta\phi", "\t", "\t".join(["[{0:.2f} - {1:.2f}]".format(phiBins[phiBinName][0], phiBins[phiBinName][1]) for phiBinName in phiBins])
    for etaBinName in etaBins:
        print "[{0:.3f} - {1:.1f}]".format(etaBins[etaBinName][0], etaBins[etaBinName][1]), "\t", "\t".join([tableEntriesFormating(type, resScale_values[etaBinName][phiBinName]) for phiBinName in phiBins])


# extract resolution/scale info for all eta and phi bins
def extractResolutionScale(resolutionFileAndInfoMap, GEN_engpt):
    if (verbosityLevel >= 1):
        print "Mean dE/E (%)", "\t", "\t", "eta", "\t\t", "phi"
    resScale_values = {}  # to keep values for printing the final 2D eta=phi tables
    for etaBinName in etaBins:
        if (verbosityLevel >= 1):
            GEN_eta = "[{0:.3f} - {1:.1f}]".format(etaBins[etaBinName][0], etaBins[etaBinName][1])
            print "Exrtacting info for eta range ", GEN_eta
        resScale_values[etaBinName] = {}
        for phiBinName in phiBins:
            if (verbosityLevel >= 1):
                GEN_phi = "[{0:.2f} - {1:.2f}]".format(phiBins[phiBinName][0], phiBins[phiBinName][1])
                print "Exrtacting info for phi range ", GEN_phi
            resScale_values[etaBinName][phiBinName] = {}

            # extract renormalisation factor (to noPU values)
            hist_renorm = resolutionFileAndInfoMap['renormNoPU']['file'].Get(resolutionFileAndInfoMap['renormNoPU']['hist_prefix'] + "_eta" + etaBinName + "_phi" + phiBinName)  # get histogram for renormalisation
            (hist_renorm, hEntriesTmp, gMeanTmp, gMeanErrorTmp, gStdTmp, effSigmaTmp) = hgcalHistHelpers.getHistMeanStd(hist_renorm)  # get histo stat. properties
            renorm = GEN_engpt / gMeanTmp

            # resolution/scale histograms and info (e, pt, etc.)
            histo = resolutionFileAndInfoMap['file'].Get(resolutionFileAndInfoMap['hist_prefix'] + "_eta" + etaBinName + "_phi" + phiBinName)  # get resolution/scale histogram
            (histo, hEntries, gMean, gMeanError, gStd, effSigma) = hgcalHistHelpers.getHistMeanStd(histo)  # get histo stat. properties

            # save values for printing the 2D tables
            resScale_values[etaBinName][phiBinName]['mean'] = gMean
            resScale_values[etaBinName][phiBinName]['meanError'] = gMeanError
            resScale_values[etaBinName][phiBinName]['effSigma'] = effSigma * renorm
            resScale_values[etaBinName][phiBinName]['calib'] = resScale_values[etaBinName][phiBinName]['mean'] / 100. + 1.
            # print basic info on scale/resolution
            if (verbosityLevel >= 1):
                "{0:.2f}".format(gMean) + " +/- " + "{0:.2f}".format(gMeanError), "\t", GEN_eta, "\t", GEN_phi, "\t", "(pf clusters: " + \
                    str(int(hEntries)) + ")", "\t", "effSigma * renorm: ", effSigma * renorm, "\t", "gStd: ", gStd
    return resScale_values


# prepare/plot/save histograms for comparison (e, pt, etc.)
def plotComparisons(histsFilesAndInfoMap, resScale_values, pidSelected, GEN_engpt, outDir):
    for etaBinName in etaBins:
        if (verbosityLevel >= 1):
            GEN_eta = "[{0:.3f} - {1:.1f}]".format(etaBins[etaBinName][0], etaBins[etaBinName][1])
            print "Exrtacting info for eta range ", GEN_eta
        for phiBinName in phiBins:
            if (verbosityLevel >= 1):
                GEN_phi = "[{0:.2f} - {1:.2f}]".format(phiBins[phiBinName][0], phiBins[phiBinName][1])
                print "Extracting info for phi range ", GEN_phi
            # get histograms to plot overalpped
            histsAndProps = {histsFilesAndInfoMap[obj]['file'].Get(histsFilesAndInfoMap[obj]['hist_prefix'] + "_eta" + etaBinName + "_phi" + phiBinName): histsFilesAndInfoMap[obj] for obj in histsFilesAndInfoMap.keys() if len(histsFilesAndInfoMap[obj]['file'].GetListOfKeys()) > 0}
            # prepare basic info and plot these histograms on top of each other
            gMeanCalibEnergy = (resScale_values[etaBinName][phiBinName]['mean'] / 100. + 1.) * GEN_engpt
            gMeanCalibEnergyError = (resScale_values[etaBinName][phiBinName]['meanError'] / 100.) * GEN_engpt
            plotComments = ["p_{T,fit} = " + "{0:.2f}".format(gMeanCalibEnergy) + " #pm " + "{0:.2f}".format(gMeanCalibEnergyError) + " GeV",
                            "#frac{#sigma_{eff}}{p_{T,Mean}} = " + "{0:.1f}".format(resScale_values[etaBinName][phiBinName]['effSigma']) + " %"]
            plotFileTag = "obj_histsOverlayed_" + etaBinName + "_" + phiBinName + "_pid" + str(pidSelected) + "_engpt" + str(int(GEN_engpt))
            # plot and save comparison hists
            hgcalHistHelpers.histsPrintSaveSameCanvas(histsAndProps, outDir, tag=plotFileTag, latexComment=plotComments)


# prepare summary resolution graphs (needs one graph per scenario, multiple graphs could be on one plot)
def setupSummaryGraphs(pidSelected, resScale_values, scenarios, type='effSigma', vsDep = 'pt'):
    # common stuff
    graphsAndProps = {}
    grOptions = {}
    if (vsDep not in ['pt','eta']):
        print "Error: Required depedance (" + vsDep + ") not implemented. No graph added."
        return (graphsAndProps, grOptions)
    grOptions['title'] = {'effSigma':"Energy resolution (",'calib':"Energy scale ("}[type] + pidmap[pidSelected] + "s)"
    grOptions['Xaxis'] = {'pt':{'effSigma':"p_{T}[GeV]",'calib':"p_{T}[GeV]"}[type], 'eta':{'effSigma':"|#eta|",'calib':"|#eta|"}[type]}[vsDep]
    grOptions['Yaxis'] = {'effSigma':"#sigma_{eff}(E)/E [%]",'calib':"E_{mean}/E_{true}"}[type]
    # buils list of graphs to be compared/plotted together, one graph for one scenario (for each eta bin)
    etaBinCenters = [0.5*(etaBins[etaBinName][0] + etaBins[etaBinName][1]) for etaBinName in etaBins]
    for scenario in scenarios:
        grIDKeys = {'pt':etaBins.keys(), 'eta':resScale_values[scenario].keys()}[vsDep]
        xBinCenters = {'pt':resScale_values[scenario].keys(), 'eta':etaBinCenters}[vsDep]
        for indx, grID in enumerate(grIDKeys):
            xEng = array('f', xBinCenters)
            if (vsDep == 'pt'): # for pt depedance graphs
                sGraphTag = "|#eta| #in [{0:.3f} - {1:.1f}]".format(etaBins[grID][0], etaBins[grID][1])
                yRes = array('f', [resScale_values[scenario][engpt][grID]['phim1p0pito1p0pi'][type] for engpt in xBinCenters])
            if (vsDep == 'eta'): # for eta depedance graphs
                sGraphTag = "p_{T} = " + "{0:.1f} GeV".format(grID)
                yRes = array('f', [resScale_values[scenario][grID][etaBinName]['phim1p0pito1p0pi'][type] for etaBinName in etaBins])
            if (scenario == "PF_PU200"):  # for scenario PF_PU200
                graphsAndProps[ROOT.TGraph(len(xEng), xEng, yRes)] = {"leg": "PF cluster, pileup 200, "+sGraphTag, "color": ROOT.kBlue, "MarkerStyle": 21+indx, "LineStyle": 4+indx}
            elif (scenario == "PF_noPU"):  # for scenario PF_noPU
                graphsAndProps[ROOT.TGraph(len(xEng), xEng, yRes)] = {"leg": "PF cluster, no pileup, "+sGraphTag, "color": ROOT.kGreen - 6, "MarkerStyle": 21+indx, "LineStyle": 4+indx}
            elif (scenario == "Mega_PU200"):  # for scenario PF_PU200
                graphsAndProps[ROOT.TGraph(len(xEng), xEng, yRes)] = {"leg": "Megacluster, pileup 200, "+sGraphTag, "color": ROOT.kRed - 6, "MarkerStyle": 25+indx, "LineStyle": 1+indx}
            elif (scenario == "Mega_noPU"):  # for scenario PF_noPU
                graphsAndProps[ROOT.TGraph(len(xEng), xEng, yRes)] = {"leg": "Megacluster, no pileup, "+sGraphTag, "color": ROOT.kBlack, "MarkerStyle": 20+indx, "LineStyle": 1+indx}
            else:  # here implement graphs for additional scenarios...
                print "Warning: Required scenario (" + scenario + ") not implemented. No graph added."
    return (graphsAndProps, grOptions)


# setup scenario for resolution and scale extraction and comparions
def setupResScaleScenario(inputdir, gun_type, pidSelected, GEN_engpt, refName, scenario="PF_noPU"):
    if (scenario == "PF_PU200"):  # scenario: PU, resolutoin from PF, comparison "PF vs. PF corrected.
        # list of files and corresponding info
        filePF = ROOT.TFile.Open(inputdir + "/{}_{}_{}GeV_{}_{}_{}.root".format(gun_type, pidSelected, int(GEN_engpt), refName, "pfcluster", "PU200"), "read")  # info based on PF energy
        filePFuncalib = ROOT.TFile.Open(inputdir + "/{}_{}_{}GeV_{}_{}_{}.root".format(gun_type, pidSelected, int(GEN_engpt), refName, "pfcluster_uncalib", "PU200"), "read")  # info based on PF energy
        fileRenormNoPU = ROOT.TFile.Open(inputdir + "/{}_{}_{}GeV_{}_{}_{}.root".format(gun_type, pidSelected, int(GEN_engpt), refName, "pfcluster", "noPU"), "read")  # info based on megacluster energy, noPU
        # map of histograms and files
        histsFilesAndInfoMap = {"obj_Pt": {'file': filePF, 'hist_prefix': "obj_Pt", 'leg': "PF cluster (calibrated)", 'color': ROOT.kBlue},
                                "ref_Pt": {'file': filePF, 'hist_prefix': "ref_Pt", 'leg': "GEN particle (" + gun_type + "={0:.1f} GeV".format(GEN_engpt) + ")", 'color': ROOT.kRed},
                                "cmp_Pt": {'file': filePFuncalib, 'hist_prefix': "obj_Pt", 'leg': "PF (non-calibrated) cluster", 'color': ROOT.kGreen - 6}}
        resolutionFileAndInfoMap = {'file': filePF, 'hist_prefix': "obj_dEoverE", 'leg': "PF (calibrated) cluster", 'color': ROOT.kBlue,
                                    'renormNoPU': {'file': fileRenormNoPU,  'hist_prefix': "obj_Pt"}}
    elif (scenario == "PF_noPU"):  # scenario: noPU, resolutoin from PF, comparison "PF vs. PF corrected
        # list of files and corresponding info
        filePF = ROOT.TFile.Open(inputdir + "/{}_{}_{}GeV_{}_{}_{}.root".format(gun_type, pidSelected, int(GEN_engpt), refName, "pfcluster", "noPU"), "read")  # info based on PF energy
        filePFuncalib = ROOT.TFile.Open(inputdir + "/{}_{}_{}GeV_{}_{}_{}.root".format(gun_type, pidSelected, int(GEN_engpt), refName, "pfcluster_uncalib", "PU200"), "read")  # info based on PF energy
        fileRenormNoPU = ROOT.TFile.Open(inputdir + "/{}_{}_{}GeV_{}_{}_{}.root".format(gun_type, pidSelected, int(GEN_engpt), refName, "pfcluster", "noPU"), "read")  # info based on megacluster energy, noPU
        # map of histograms and files
        histsFilesAndInfoMap = {"obj_Pt": {'file': filePF, 'hist_prefix': "obj_Pt", 'leg': "PF cluster (calibrated)", 'color': ROOT.kBlue},
                                "ref_Pt": {'file': filePF, 'hist_prefix': "ref_Pt", 'leg': "GEN particle (" + gun_type + "={0:.1f} GeV".format(GEN_engpt) + ")", 'color': ROOT.kRed},
                                "cmp_Pt": {'file': filePFuncalib, 'hist_prefix': "obj_Pt", 'leg': "PF (non-calibrated) cluster", 'color': ROOT.kGreen - 6}}
        resolutionFileAndInfoMap = {'file': filePF, 'hist_prefix': "obj_dEoverE", 'leg': "PF (calibrated) cluster", 'color': ROOT.kBlue,
                                    'renormNoPU': {'file': fileRenormNoPU,  'hist_prefix': "obj_Pt"}}
    elif (scenario == "Mega_PU200"):  # scenario: PU, resolutoin from Mega cluster, comparison "Mega vs. PF corrected.
        # list of files and corresponding info
        filePF = ROOT.TFile.Open(inputdir + "/{}_{}_{}GeV_{}_{}_{}.root".format(gun_type, pidSelected, int(GEN_engpt), refName, "pfcluster", "PU200"), "read")  # info based on PF energy
        fileMega = ROOT.TFile.Open(inputdir + "/{}_{}_{}GeV_{}_{}_{}.root".format(gun_type, pidSelected, int(GEN_engpt), refName, "megacluster", "PU200"), "read")  # info based on megacluster energy
        fileRenormNoPU = ROOT.TFile.Open(inputdir + "/{}_{}_{}GeV_{}_{}_{}.root".format(gun_type, pidSelected, int(GEN_engpt), refName, "megacluster", "noPU"), "read")  # info based on megacluster energy, noPU
        # map of histograms and files
        histsFilesAndInfoMap = {"obj_Pt": {'file': fileMega, 'hist_prefix': "obj_Pt", 'leg': "Megacluster (non-calibrated)", 'color': ROOT.kBlue},
                                "ref_Pt": {'file': filePF,   'hist_prefix': "ref_Pt", 'leg': "GEN particle (" + gun_type + "={0:.1f} GeV".format(GEN_engpt) + ")", 'color': ROOT.kRed},
                                "cmp_Pt": {'file': filePF,   'hist_prefix': "obj_Pt", 'leg': "PF (calibrated) cluster", 'color': ROOT.kGreen - 6}}
        resolutionFileAndInfoMap = {'file': fileMega, 'hist_prefix': "obj_dEoverE", 'leg': "PF (calibrated) cluster", 'color': ROOT.kBlue,
                                    'renormNoPU': {'file': fileRenormNoPU,  'hist_prefix': "obj_Pt"}}
    elif (scenario == "Mega_noPU"):  # scenario: noPU, resolutoin from Mega cluster, comparison "Mega vs. PF corrected
        # list of files and corresponding info
        filePF = ROOT.TFile.Open(inputdir + "/{}_{}_{}GeV_{}_{}_{}.root".format(gun_type, pidSelected, int(GEN_engpt), refName, "pfcluster", "noPU"), "read")  # info based on PF energy
        fileMega = ROOT.TFile.Open(inputdir + "/{}_{}_{}GeV_{}_{}_{}.root".format(gun_type, pidSelected, int(GEN_engpt), refName, "megacluster", "noPU"), "read")  # info based on megacluster energy
        fileRenormNoPU = ROOT.TFile.Open(inputdir + "/{}_{}_{}GeV_{}_{}_{}.root".format(gun_type, pidSelected, int(GEN_engpt), refName, "megacluster", "noPU"), "read")  # info based on megacluster energy, noPU
        # map of histograms and files
        histsFilesAndInfoMap = {"obj_Pt": {'file': fileMega, 'hist_prefix': "obj_Pt", 'leg': "Megacluster (non-calibrated)", 'color': ROOT.kBlue},
                                "ref_Pt": {'file': filePF,   'hist_prefix': "ref_Pt", 'leg': "GEN particle (" + gun_type + "={0:.1f} GeV".format(GEN_engpt) + ")", 'color': ROOT.kRed},
                                "cmp_Pt": {'file': filePF,   'hist_prefix': "obj_Pt", 'leg': "PF (calibrated) cluster", 'color': ROOT.kGreen - 6}}
        resolutionFileAndInfoMap = {'file': fileMega, 'hist_prefix': "obj_dEoverE", 'leg': "Mega (non-calibrated) cluster", 'color': ROOT.kBlue,
                                    'renormNoPU': {'file': fileRenormNoPU,  'hist_prefix': "obj_Pt"}}
    elif (scenario == "Mega_noPU_PU200"):  # scenario: noPU, resolutoin from Mega cluster, comparison "Mega vs. PF corrected
        # list of files and corresponding info
        fileMega_PU200 = ROOT.TFile.Open(inputdir + "/{}_{}_{}GeV_{}_{}_{}.root".format(gun_type, pidSelected, int(GEN_engpt), refName, "megacluster", "PU200"), "read")  # info based on megacluster energy
        fileMega_noPU = ROOT.TFile.Open(inputdir + "/{}_{}_{}GeV_{}_{}_{}.root".format(gun_type, pidSelected, int(GEN_engpt), refName, "megacluster", "noPU"), "read")  # info based on megacluster energy
        fileRenormNoPU = ROOT.TFile.Open(inputdir + "/{}_{}_{}GeV_{}_{}_{}.root".format(gun_type, pidSelected, int(GEN_engpt), refName, "megacluster", "noPU"), "read")  # info based on megacluster energy, noPU
        # map of histograms and files
        histsFilesAndInfoMap = {"obj_Pt": {'file': fileMega_PU200, 'hist_prefix': "obj_Pt", 'leg': "Megacluster, pile-up 200", 'color': ROOT.kBlue},
                                "cmp_Pt": {'file': fileMega_noPU,  'hist_prefix': "obj_Pt", 'leg': "Megacluster, no pile-up", 'color': ROOT.kGreen - 6}}
        resolutionFileAndInfoMap = {'file': fileMega_PU200, 'hist_prefix': "obj_dEoverE", 'leg': "Megacluster, pile-up 200", 'color': ROOT.kBlue,
                                    'renormNoPU': {'file': fileRenormNoPU,  'hist_prefix': "obj_Pt"}}
    else:
        print "Error: Required scenario (" + scenario + ") not implemented..."
        sys.exit()
    return histsFilesAndInfoMap, resolutionFileAndInfoMap


def main():

    global opt, args

    usage = ('usage: %prog [options]\n' + '%prog -h for help')
    parser = optparse.OptionParser(usage)

    # input options
    parser.add_option('', '--inputDir', dest='inputDir', type='string',  default='./', help='path to the input directory')
    parser.add_option('', '--gunType', dest='gunType', type='string',  default='pt', help='pt or e')
    parser.add_option('', '--pids', dest='pid', type='string',  default='22', help='pdgId string (comma-separated list)')
    parser.add_option('', '--genValues', dest='genValue', type='string',  default='25', help='generated pT or energy (comma-separated list)')
    parser.add_option('', '--ref', dest='refName', type='string',  default='genpart', help='reference collection')
    parser.add_option('', '--scen', dest='scenarios', type='string',  default='PF_noPU', help='scenario for res/scale/comparisons (PF_noPU, PF_PU200, Mega_noPU, Mega_PU200, Mega_noPU_PU200)')
    parser.add_option('', '--tag', dest='tag', type='string',  default='test', help='some tag, to be attached to the results of processing given scenarios')

    # store options and arguments as global variables
    global opt, args
    (opt, args) = parser.parse_args()

    # prepare list of E/Pt points, PIDs, scenarios
    pids = [int(p.strip(" ")) for p in opt.pid.split(",")]
    GEN_engpts = [float(engpt.strip(" ")) for engpt in opt.genValue.split(",")]
    scenarios = [scen.strip(" ") for scen in opt.scenarios.split(",")]

    print "inputdir: ", opt.inputDir
    print "gunType: ", opt.gunType
    print "pids: ", pids
    print "GEN_engpts: ", GEN_engpts
    print "refName: ", opt.refName
    print "Scenarios for resolution/scale/comparisons: ", scenarios
    print "output tag: ", opt.tag

    # set some variables based on input options
    inputDir = opt.inputDir
    gun_type = opt.gunType
    tag = opt.tag
    refName = opt.refName

    # time stamp - start
    start_time = timeit.default_timer()

    # run over all pid and engpt points, extract resolution/scale info, plot comparison plots
    resScale_values = {}
    for pidSelected in pids:
        print "Processing info for PID " + str(pidSelected)
        resScale_values[pidSelected] = {}
        for scenario in scenarios:
            print "Processing info for scenario " + scenario + "..."
            resScale_values[pidSelected][scenario] = {}
            # prepare some basics
            outDir = "ScaleResolutionPlots_pid" + str(pidSelected) + "_" + scenario + "_" + tag
            if not os.path.exists(outDir):
                os.makedirs(outDir)
            for GEN_engpt in GEN_engpts:
                print "Processing info for E/Pt point " + "{0:.1f}GeV".format(GEN_engpt)
                # setup scenario for resolution/scale/comparisons
                (histsFilesAndInfoMap, resolutionFileAndInfoMap) = setupResScaleScenario(inputDir, gun_type, pidSelected, GEN_engpt, refName, scenario)
                if len(resolutionFileAndInfoMap['file'].GetListOfKeys()) > 0:
                    # extract resolution/scale info for all eta and phi bins
                    resScale_values[pidSelected][scenario][GEN_engpt] = extractResolutionScale(resolutionFileAndInfoMap, GEN_engpt)
                    # prepare/plot/save histograms for comparison (e, pt, etc.)
                    plotComparisons(histsFilesAndInfoMap, resScale_values[pidSelected][scenario][GEN_engpt], pidSelected, GEN_engpt, outDir)
                    # print the 2D eta-phi tables
                    printEtaPhiTable(resScale_values[pidSelected][scenario][GEN_engpt], type='mean')
                    printEtaPhiTable(resScale_values[pidSelected][scenario][GEN_engpt], type='effSigma')
                    printEtaPhiTable(resScale_values[pidSelected][scenario][GEN_engpt], type='calib')
                else:
                    print "Warning: Empty file. Skipping GEN_engpt point ", GEN_engpt, "..."
                # close all the files
                resolutionFileAndInfoMap['file'].Close()
                for obj in histsFilesAndInfoMap.keys():
                    histsFilesAndInfoMap[obj]['file'].Close()

    # produce summary resolution/calibration plots for set of scenarios
    for pidSelected in pids:
        print "Preparing summary plots for PID " + str(pidSelected)
        outDir = "ScaleResolutionPlots_pid" + str(pidSelected) + "_" + tag
        if not os.path.exists(outDir):
            os.makedirs(outDir)
        for type in ['calib', 'effSigma']:
            for vsDep in ['pt','eta']:
                graphsAndProps, grOptions = setupSummaryGraphs(pidSelected, resScale_values[pidSelected], scenarios, type, vsDep)
                if type == 'effSigma' and vsDep == 'pt':
                    for gr in graphsAndProps:
                        (gr, stohasticTerm, constantTerm, noiseTerm) = hgcalHistHelpers.fitResolution(gr, graphsAndProps[gr]['color'], graphsAndProps[gr]['LineStyle'])
                        # graphsAndProps[gr]['latexComment'] = "#frac{#sigma(p_{T})}{p_{T}} = #frac{" +"{0:.1f}".format(stohasticTerm) + "%}{#sqrt{p_{T}}} #oplus " + "{0:.1f}".format(constantTerm) + "%"
                hgcalHistHelpers.drawGraphs(graphsAndProps, grOptions, outDir, tag=type+"_vs_"+vsDep+"_"+ pidmap[pidSelected] + "_scenarios_" + "_".join(scenarios) + "_" + tag)

    # time stamp - end
    elapsed = timeit.default_timer() - start_time
    print "Time:", elapsed


if __name__ == '__main__':
    main()
