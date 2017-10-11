# investigate shower development based on RecHits and SimClusters
import ROOT
import os
import numpy as np
from array import array
from HGCalImagingAlgo import *
from NtupleDataFormat import HGCalNtuple
from GeoUtils import *

# filtering parameters
dependSensor = True
ecut = 3 # relative to the noise
# verbosity etc.
verbosityLevel = 0 # 0 - only basic info (default); 1 - additional info; 2 - detailed info printed, histograms produced

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
def getRecHitsSimAssoc(rechits_raw, simcluster):
    # get sim-cluster associations
    nSimClus = 0
    simClusHitAssoc = []
    recHitDetIds = getRecHitDetIds(rechits_raw)
    for simClusIndex, simClus in enumerate(simcluster):
        simClusHitAssoc.append(getHitList(simClus, recHitDetIds))
        nSimClus += 1
    # get list of rechist associated to simhits
    rHitsSimAssoc = [[] for k in range(0,nSimClus)]
    for simClusIndex, simCl in enumerate(simcluster):
        if (verbosityLevel>=1): print "Sim-cluster index: ",simClusIndex, ", pT: ",simCl.pt(), ", E: ",simCl.energy(), ", phi: ",simCl.phi(), ", eta: ",simCl.eta()
        # loop over sim clusters and then rechits
        rHitsSimAssocTemp = []
        for hitIndexArray in simClusHitAssoc[simClusIndex]:
            for hitIndex in hitIndexArray:
                thisHit = rechits_raw[hitIndex]
                if(not recHitAboveTreshold(thisHit, ecut, dependSensor)[1]): continue
                # independent of sim cluster, after cleaning
                rHitsSimAssocTemp.append(thisHit)
        rHitsSimAssoc[simClusIndex]= rHitsSimAssocTemp
    return rHitsSimAssoc

# 1D histograming of given list of values
def histValue1D(fValues, histDict, tag = "hist1D_", title = "hist 1D", axunit = "a.u.", binsBoundariesX = [10, -1, 1], ayunit = "a.u."):
    # sanity check for hists
    if (histDict == None): return
    # sanity check for boundaries
    if (len(binsBoundariesX) != 3 and len(binsBoundariesX) != 2): return
    # define event-level hists
    elif len(binsBoundariesX) == 3: # bondaries in format [nbins, low, high]
        histDict[tag]  = ROOT.TH1F(tag, title+";"+axunit+";"+ayunit, binsBoundariesX[0], binsBoundariesX[1], binsBoundariesX[2])
    elif len(binsBoundariesX) == 2: # bondaries in format [nbins, list_boundaries]
        histDict[tag]  = ROOT.TH1F(tag, title+";"+axunit+";"+ayunit, binsBoundariesX[0], array('f', binsBoundariesX[1]))
    # set some properties
    histDict[tag].GetYaxis().SetTitleOffset(histDict[tag].GetYaxis().GetTitleOffset()*3.0)
    # loop over all values
    if (verbosityLevel>=3): print "tag: ", tag, ", fValues: ", fValues
    for value in fValues:
        histDict[tag].Fill(value)
    return histDict

# 2D histograming of given list of values
def histValues2D(fValues, histDict, tag = "hist2D_", title = "hist 2D", axunit = "a.u.", binsBoundariesX = [10, -1, 1], ayunit = "a.u.", binsBoundariesY = [10, -1, 1], weighted2D = False):
    # sanity check for hists
    if (histDict == None): return
    # sanity check for boundaries
    if (len(binsBoundariesX) != len(binsBoundariesY)): return
    if (len(binsBoundariesX) != 3 and len(binsBoundariesX) != 2): return
    # define event-level hists
    elif len(binsBoundariesX) == 3: # bondaries in format [nbins, low, high]
        histDict[tag]  = ROOT.TH2F(tag, title+";"+axunit+";"+ayunit, binsBoundariesX[0], binsBoundariesX[1], binsBoundariesX[2], binsBoundariesY[0], binsBoundariesY[1], binsBoundariesY[2])
    elif len(binsBoundariesY) == 2: # bondaries in format [nbins, list_boundaries]
        histDict[tag]  = ROOT.TH2F(tag, title+";"+axunit+";"+ayunit, binsBoundariesX[0], array('f', binsBoundariesX[1]), binsBoundariesY[0], array('f', binsBoundariesY[1]))
    # set some properties
    histDict[tag].GetXaxis().SetTitleOffset(histDict[tag].GetXaxis().GetTitleOffset()*1.0)
    histDict[tag].GetYaxis().SetTitleOffset(histDict[tag].GetYaxis().GetTitleOffset()*3.0)
    # loop over all values
    if (verbosityLevel>=3): print "tag: ", tag, ", fValues: ", fValues
    if (not weighted2D):
        for (valueX, valueY) in fValues:
            histDict[tag].Fill(valueX, valueY)
    else:
        for (valueX, valueY, valueZ) in fValues:
            histDict[tag].Fill(valueX, valueY, valueZ)
    return histDict

# implement simHits
class simHit:
    def __init__(self, energy = None, pt = None, layer = None, eta = None, phi = None, thickness = None):
        self._eta = 0
        self._phi = 0
        self._energy = 0
        self._pt = 0
        self._thickness = 100.
        self._layer = int(1)
        if energy is not None:
            self._energy = energy
        if pt is not None:
            self._pt = pt
        if layer is not None:
            self._layer = int(layer)
        if eta is not None:
            self._eta = eta
        if phi is not None:
            self._phi = phi
        if thickness is not None:
            self._thickness = thickness
        lv = ROOT.TLorentzVector()
        lv.SetPtEtaPhiE(self._pt,self._eta,self._phi,self._energy)
        self._x = lv.X()
        self._y = lv.Y()
        self._z = lv.Z()
    def thickness(self):
        return self._thickness
    def layer(self):
        return self._layer
    def energy(self):
        return self._energy
    def pt(self):
        return self._pt
    def eta(self):
        return self._eta
    def phi(self):
        return self._phi
    def x(self):
        return self._x
    def y(self):
        return self._y
    def z(self):
        return self._z

# print/save list of histograms with their properties on one canvas
def histsPrintSaveSameCanvas(histsAndProps, outDir, tag = "hists1D_", latexComment = "", funcsAndProps = None):
    # supress info messages
    ROOT.gErrorIgnoreLevel = ROOT.kInfo+1
    # set default style values
    ROOT.gStyle.SetPalette(ROOT.kBird)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPadTopMargin(0.08)
    ROOT.gStyle.SetPadBottomMargin(0.12)
    ROOT.gStyle.SetPadLeftMargin(0.12)
    ROOT.gStyle.SetPadRightMargin(0.05)
    # create canvas
    canvas = ROOT.TCanvas(outDir+tag, outDir+tag, 500, 500)
    # prepare the legend
    leg = ROOT.TLegend(0.15,0.75,0.62,0.9)
    #leg.SetHeader("Energy of the clusters before/after filtering")
    leg.SetBorderSize(0)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    #prepare latex comment
    ltx = ROOT.TLatex()
    ltx.SetNDC(ROOT.kTRUE)
    ltx.SetTextFont(42);
    ltx.SetTextSize(0.03)
    # set image extensions
    imgTypes = ["pdf","png"]
    if (verbosityLevel>=3):
        print "histsAndProps: ", histsAndProps
        print "funcsAndProps: ", funcsAndProps
    # loop over all histograms to get max
    y_maxs = [1.]
    for hist in histsAndProps:
        # do not print/save empty histograms
        if (type(hist) == ROOT.TH1F) or (type(hist) == ROOT.TH2F) or (type(hist) == ROOT.TH3F):
            if hist.GetEntries() == 0:
                continue
        curr_max = hist.GetMaximum()
        if (curr_max < hist.GetEntries()/2.):
            y_maxs.append(curr_max)
    #print "y_maxs: ", y_maxs
    # loop over all histograms
    first = True
    for hist in histsAndProps:
        # do not print/save empty histograms
        if (type(hist) == ROOT.TH1F) or (type(hist) == ROOT.TH2F) or (type(hist) == ROOT.TH3F):
            if hist.GetEntries() == 0:
                continue
        # print and save
        if type(hist) == ROOT.TH1F:
            hist.SetLineColor(histsAndProps[hist]["color"])
            leg.AddEntry(hist, histsAndProps[hist]["leg"], "L")
            hist.GetXaxis().SetTitleOffset(hist.GetXaxis().GetTitleOffset()*1.2)
            hist.GetYaxis().SetTitleOffset(hist.GetYaxis().GetTitleOffset()*3.0)
            if (first):
                hist.GetYaxis().SetRangeUser(0,max(y_maxs)*1.4)
                hist.Draw("hist0 goff")
                first = False
            else:
                hist.Draw("hist0 same goff")
    # check if any function should be drawn
    if (funcsAndProps!=None):
        for func in funcsAndProps:
            func.SetLineColor(funcsAndProps[func]["color"])
            leg.AddEntry(func, funcsAndProps[func]["leg"], "L")
            func.Draw("same goff")
    # draw the rest
    leg.Draw("same")
    ltx.DrawLatex(0.170,0.71,latexComment[0])
    ltx.DrawLatex(0.170,0.66,latexComment[1])
    for imgType in imgTypes:
        canvas.SaveAs("{}/{}.{}".format(outDir, tag, imgType))
    return canvas

def drawGraphsTest(graphsAndProps, title = "Resolution" , tag = "graphTest_"):
    # supress info messages
    ROOT.gErrorIgnoreLevel = ROOT.kInfo+1
    # set default style values
    ROOT.gStyle.SetPalette(ROOT.kBird)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPadTopMargin(0.08)
    ROOT.gStyle.SetPadBottomMargin(0.12)
    ROOT.gStyle.SetPadLeftMargin(0.10)
    ROOT.gStyle.SetPadRightMargin(0.05)
    # create canvas
    canvas = ROOT.TCanvas(tag, tag, 800, 600)
    # prepare the legend
    leg = ROOT.TLegend(0.55,0.45,0.9,0.9)
    leg.SetHeader(title)
    leg.SetBorderSize(0)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    #prepare latex comment
    ltx = ROOT.TLatex()
    ltx.SetNDC(ROOT.kTRUE)
    ltx.SetTextFont(42);
    ltx.SetTextSize(0.03)
    # set image extensions
    imgTypes = ["pdf","png","root"]
    if (verbosityLevel>=3):
        print "graphsAndProps: ", graphsAndProps
    # loop over all graphs to get max
    y_maxs = [gr.GetYaxis().GetXmax() for gr in graphsAndProps]
    y_mins = [gr.GetYaxis().GetXmin() for gr in graphsAndProps]
    print "y_mins: ", y_mins
    print "y_maxs: ", y_maxs
    # loop over all histograms
    first = True
    for gr in graphsAndProps:
        gr.GetXaxis().SetTitle( 'E[GeV]' )
        gr.GetYaxis().SetTitle( '#sigma_{eff}[%]' )
        colour = graphsAndProps[gr]["color"]
        gr.SetLineColor( colour )
        gr.SetLineWidth( 1 )
        gr.SetLineStyle( graphsAndProps[gr]["LineStyle"] )
        gr.SetMarkerColor( colour )
        gr.SetMarkerStyle( graphsAndProps[gr]["MarkerStyle"] )
        gr.SetMarkerSize( 0.7 )
        gr.SetFillColor( 0 )
        gr.SetFillStyle( 0 )
        leg.AddEntry(gr, graphsAndProps[gr]["leg"] )
        if (first):
            gr.SetMaximum(max(y_maxs)*1.5)
            gr.SetMinimum(0)
            gr.Draw("AP goff")
            first = False
        else:
            gr.Draw("P same goff")
    # legend
    leg.Draw("same")
    # save
    for imgType in imgTypes:
        canvas.SaveAs("{}.{}".format(tag, imgType))
    return canvas

# print/save all histograms
def histPrintSaveAll(histDict, outDir, tag = "_test"):
    # supress info messages
    ROOT.gErrorIgnoreLevel = ROOT.kInfo+1
    # set image extensions
    imgType = "pdf"
    if (verbosityLevel>=3): print "histDict.items(): ", histDict.items()
    for key, item in histDict.items():
        # do not save empty histograms
        if (type(item) == ROOT.TH1F) or (type(item) == ROOT.TH2F) or (type(item) == ROOT.TH3F):
            if item.GetEntries() == 0:
                continue
        # set default style values
        ROOT.gStyle.SetPalette(ROOT.kBird)
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetPadTopMargin(0.10)
        ROOT.gStyle.SetPadBottomMargin(0.12)
        ROOT.gStyle.SetPadLeftMargin(0.15)
        ROOT.gStyle.SetPadRightMargin(0.05)
        # print and save
        if type(item) == ROOT.TH1F:
            ROOT.gStyle.SetPadRightMargin(0.05)
            canvas = ROOT.TCanvas(outDir+tag+key, outDir+tag+key, 500, 500)
            item.Draw("hist0 goff")
            canvas.SaveAs("{}/{}.{}".format(outDir, key, imgType))
        if type(item) == ROOT.TH2F:
            ROOT.gStyle.SetPadRightMargin(0.15)
            canvas = ROOT.TCanvas(outDir+tag+key, outDir+tag+key, 500, 500)
            item.Draw("colz goff")
            canvas.SaveAs("{}/{}.{}".format(outDir, key, imgType))
        elif type(item) == ROOT.TH3F:
            ROOT.gStyle.SetPadRightMargin(0.05)
            canvas = ROOT.TCanvas(outDir+tag+key, outDir+tag+key, 500, 500)
            item.Draw("box goff")
            canvas.SaveAs("{}/{}.{}".format(outDir, key, imgType))
        else:
            continue
    #del canvas
    return

def fitGauss(hist, paramRangeFactor = 1.8):
    if (hist.GetEntries()==0): return (hist, 0, 0)
    hist.GetXaxis().SetTitleOffset(hist.GetXaxis().GetTitleOffset()*1.2)
    hist.GetYaxis().SetTitleOffset(hist.GetYaxis().GetTitleOffset()*3.0)
    # define the range of the fit from the hist mean and RMS
    meanLimitDn = hist.GetMean() - paramRangeFactor*hist.GetRMS()
    meanLimitUp = hist.GetMean() + paramRangeFactor*hist.GetRMS()
    sigmaLimitDn = hist.GetRMS() / paramRangeFactor
    sigmaLimitUp = hist.GetRMS() * paramRangeFactor
    # define the fitting gausian and range of its parameters
    fGauss = ROOT.TF1("f","[0]*TMath::Gaus(x,[1],[2])",meanLimitDn,meanLimitUp)
    fGauss.SetParLimits(1,meanLimitDn,meanLimitUp)
    fGauss.SetParLimits(2,sigmaLimitDn,sigmaLimitUp)
    # perform fit and extract params
    hist.Fit(fGauss,"Q","",meanLimitDn,meanLimitUp)
    gaussMean = fGauss.GetParameter(1)
    gaussStd = fGauss.GetParameter(2)
    #print "hist ent.: ", hist.GetEntries()
    #print "hist mean: ", hist.GetMean(), ", hist std: ", hist.GetRMS()
    #print "gaussMean: ", gaussMean,      ", gaussStd: ", gaussStd
    return (hist, gaussMean, gaussStd)

# taken from Hgg framework (by Ed)
def getEffSigma( theHist, wmin=-100, wmax=100, epsilon=0.01 ):
    # initialise
    weight = 0.
    points = []
    thesum = theHist.Integral()
    # return -1 in case of empty histogram
    if (thesum==0): return -1.
    # compute the cumulative distr. points
    for i in range(theHist.GetNbinsX()):
        weight += theHist.GetBinContent(i)
        if weight/thesum > epsilon:
            points.append( [theHist.GetBinCenter(i),weight/thesum] )
#    print "thesum: ", thesum
#    print "points: ", points
    # initialise
    low = wmin
    high = wmax
    width = wmax-wmin
    # find minimal 0.683 interval
    for i in range(len(points)):
        for j in range(i,len(points)):
            wy = points[j][1] - points[i][1]
            if abs(wy-0.683) < epsilon:
                wx = points[j][0] - points[i][0]
                if wx < width:
                    low = points[i][0]
                    high = points[j][0]
                    width=wx
                    #print "width: ", width, "low: ", low, ", high: ", high, ", wx: ", wx, ", wy: ", wy,
    return 0.5*(high-low)

def main():
    
    relFractionE=0.001
    # set sample/tree - for pions
    pidSelected = 211
#    GEN_eng = 5.
#    ntuple = HGCalNtuple("/eos/cms/store/cmst3/group/hgcal/CMG_studies/Production/FlatRandomEGunProducer_pdgid211_E5_cmssw93X_withPRs_20170817/NTUP/partGun_PDGid211_x100_E5.0To5.0_NTUP.root")
#    GEN_eng = 20.
#    ntuple = HGCalNtuple("/eos/cms/store/cmst3/group/hgcal/CMG_studies/Production/FlatRandomEGunProducer_pdgid211_E20_cmssw93X_withPRs_20170817/NTUP/partGun_PDGid211_x100_E20.0To20.0_NTUP.root")
#    GEN_eng = 50.
#    ntuple = HGCalNtuple("/eos/cms/store/cmst3/group/hgcal/CMG_studies/Production/FlatRandomEGunProducer_pdgid211_E50_cmssw93X_withPRs_20170809/NTUP/partGun_PDGid211_x100_E50.0To50.0_NTUP.root")
#    GEN_eng = 100.
#    ntuple = HGCalNtuple("/eos/cms/store/cmst3/group/hgcal/CMG_studies/Production/FlatRandomEGunProducer_pdgid211_E100_cmssw93X_withPRs_20170817/NTUP/partGun_PDGid211_x100_E100.0To100.0_NTUP.root")
#    GEN_eng = 300.
#    ntuple = HGCalNtuple("/eos/cms/store/cmst3/group/hgcal/CMG_studies/Production/FlatRandomEGunProducer_pdgid211_E300_cmssw93X_withPRs_20170817/NTUP/partGun_PDGid211_x100_E300.0To300.0_NTUP.root")

    # set sample/tree - for photons
    pidSelected = 22
    GEN_eng = 5.
    ntuple = HGCalNtuple("/eos/cms/store/cmst3/group/hgcal/CMG_studies/Production/FlatRandomEGunProducer_pdgid22_E5_cmssw93X_withPRs_20170809/NTUP/partGun_PDGid22_x100_E5.0To5.0_NTUP.root")
#    GEN_eng = 20.
#    ntuple = HGCalNtuple("/eos/cms/store/cmst3/group/hgcal/CMG_studies/Production/FlatRandomEGunProducer_pdgid22_E20_cmssw93X_withPRs_20170817/NTUP/partGun_PDGid22_x100_E20.0To20.0_NTUP.root")
#    GEN_eng = 50.
#    ntuple = HGCalNtuple("/eos/cms/store/cmst3/group/hgcal/CMG_studies/Production/FlatRandomEGunProducer_pdgid22_E50_cmssw93X_withPRs_20170817/NTUP/partGun_PDGid22_x100_E50.0To50.0_NTUP.root")
#    GEN_eng = 100.
#    ntuple = HGCalNtuple("/eos/cms/store/cmst3/group/hgcal/CMG_studies/Production/FlatRandomEGunProducer_pdgid22_E100_cmssw93X_withPRs_20170817/NTUP/partGun_PDGid22_x100_E100.0To100.0_NTUP.root")
#    GEN_eng = 300.
#    ntuple = HGCalNtuple("/eos/cms/store/cmst3/group/hgcal/CMG_studies/Production/FlatRandomEGunProducer_pdgid22_E300_cmssw93X_withPRs_20170817/NTUP/partGun_PDGid22_x100_E300.0To300.0_NTUP.root")

    runCalibrationScaleResolution(pidSelected, GEN_eng, ntuple, relFractionE)

def runCalibrationScaleResolution(pidSelected, GEN_eng, ntuple, relFractionE):
    # names and pid mapping
    pidmap = {11:"electron", 13:"mion", 22:"photon", 211:"pion"}

    # common strings
    GEN_pTEng = "E={0:.1f} GeV".format(GEN_eng)
    GEN_partId = pidmap[pidSelected]
    
    # init output stuff
    outDir = "simClusterHitsContained_ScaleResolution_"+GEN_partId+"s_"+"{0:.1f}GeV".format(GEN_eng)+"_testSpatial"
    if not os.path.exists(outDir): os.makedirs(outDir)
    histDict = {}
    
    # define eta and phi bins
    etaBins = {"eta1p479to1p6":(1.479, 1.6), "eta1p6to1p8":(1.6, 1.8), "eta1p8to2p0":(1.8, 2.0), "eta2p0to2p2":(2.0, 2.2), "eta2p2to2p4":(2.2, 2.4), "eta2p4to2p6":(2.4, 2.6), "eta2p6to2p8":(2.6, 2.8), "eta2p8to3p0":(2.8, 3.0), "eta1p479to3p0":(1.479, 3.0), "eta1p6to2p8":(1.6, 2.8)}
    #phiBins = {"phi0to0p5pi":(0.*math.pi, 0.5*math.pi), "phi0p5to1p0pi":(0.5*math.pi, 1.0*math.pi), "phim1p0pitom0p5pi":(-1.0*math.pi, -0.5*math.pi), "phim0p5pito0":(-0.5*math.pi, 0.*math.pi),"phim1p0pito1p0pi":(-1.0*math.pi, 1.0*math.pi) }
    # these are to run only inclusive bins
    #etaBins = {"eta1p479to3p0":(1.479, 3.0)}
    phiBins = {"phim1p0pito1p0pi":(-1.0*math.pi, 1.0*math.pi)}
    
    ## define some global lists and dictionaries
    s_all_pids = [211, 22] # pids to consider
    obj_Eng_EngRelDiff = {pid:[] for pid in s_all_pids}
    # define some global variables
    counterDivError = 0
    counterClustError = 0
    
    # initialisation of GeoUtils
    gu = GeoUtil()

    # loop over the events
    print "Total events to process (for multiClusters, ",GEN_partId,", ",GEN_pTEng,"): ", ntuple.nevents()
    hitsSelected = []
    for event in ntuple:
        #if (event.entry()>1000): break
        if (verbosityLevel>=0):
            if (event.entry()%1000 == 0): print "Event: ", event.entry()
        # get collections
        genParticles = event.genParticles()
        simClusters  = event.simClusters()
        pfClusters  = event.pfClusters()
        multiClusters = event.multiClusters()
        recHits = event.recHits()
        # check for cases where number of Sim clusters is incorrect
        if (len(simClusters)==0 or len(simClusters)>2): # skip if no sim clusters or skip if more then 2 sim clusters
            if (verbosityLevel>=1): print "WARNING (Event:",iEvt,"): Sim cluster length is improper."
            counterClustError+=1
            continue
        
        ### fill simhits based on rechits mapping and filtering ###
        # get rechit-simclus assoc.
        rHitsSimAssoc = getRecHitsSimAssoc(recHits, simClusters)
        simClusters_energyAboveThreshold = []
        simClusters_energyAboveThresholdContained = []
        simClusters_energyAll = []
        # store the rechit info
        hitsSelected.extend([(10*rechit.x(), 10*rechit.y(), rechit.energy(), recHitAboveTreshold(rechit, ecut, dependSensor)[1]) for k in range(len(simClusters)) for rechit in rHitsSimAssoc[k]])
        # store in list the energy with filtering of simi hits
        for k in range(len(simClusters)):
            sum_simHitsEnergyAll = 0
            sum_simHitsEnergyAbove = 0
            sum_simHitsEnergyContained = 0
            for rechit in rHitsSimAssoc[k]:
                # energy for rechits above noise
                sum_simHitsEnergyAll += rechit.energy()
                # energy for rechits above noise
                if recHitAboveTreshold(rechit, ecut, dependSensor)[1]:
                    sum_simHitsEnergyAbove += rechit.energy()
                # energy of rechits above noise and after masking
                if (recHitAboveTreshold(rechit, ecut, dependSensor)[1] and (rechit.layer()>35 or (rechit.layer()<=35 and gu.planes[rechit.layer()-1].contains(10*rechit.x(), 10*rechit.y())))):
                    sum_simHitsEnergyContained += rechit.energy()
            simClusters_energyAll.append(sum_simHitsEnergyAll)
            simClusters_energyAboveThreshold.append(sum_simHitsEnergyAbove)
            simClusters_energyAboveThresholdContained.append(sum_simHitsEnergyContained)

        # loop over the gen particles in current event
        objectsToCalibrate = simClusters # or objectsToCalibrate = pfClusters # or objectsToCalibrate = multiClusters or ...
        if (verbosityLevel>=1): print "PDG ID", "\t", "Sim E(GeV)", "\t", "GEN E(GeV)", "\t", "Eng.Loss (%)", "\t", "DR(GEN, Sim)"
        for genParticle in genParticles:
            # sort indecies of PF clusters according to the distance of clusters to the GEN particles (index 0 is the nearest PF cluster)
            pfs = sorted(range(len(objectsToCalibrate)), key=lambda k: math.sqrt( math.fabs(genParticle.eta()-objectsToCalibrate[k].eta())**2 + math.fabs(genParticle.phi()-objectsToCalibrate[k].phi())**2 ), reverse=False)
            # check patological cases
            if (len(objectsToCalibrate) != len(simClusters)):
                if (verbosityLevel>=1): print "ERROR (Event:",event.entry(),"): multiClusters != simClusters."
                counterClustError+=1
                continue
            if (simClusters[pfs[0]].energy() == 0.):
                if (verbosityLevel>=1): print "ERROR (Event:",event.entry(),"): simClusters[pfs[0]].energy() == 0."
                counterDivError+=1
                continue
            # prepare some observables and store it for every eta/phi
            #    for simClusters in geometry: energy_to_calibrate = simClusters_energyAboveThresholdContained[pfs[0]]
            #    for simClusters noise: energy_to_calibrate = simClusters_energyAboveThreshold[pfs[0]]
            #    for multiClusters: energy_to_calibrate = multiClusters[pfs[0]].energy()
            #    for pfClusters: energy_to_calibrate = pfClusters[pfs[0]].correctedEnergy()
            energy_to_calibrate = simClusters_energyAboveThresholdContained[pfs[0]]
            energy_to_compare = simClusters_energyAboveThreshold[pfs[0]] # or pfClusters[pfs[0]].correctedEnergy() or pfClusters[pfs[0]].energy() or ...
            simulated_energy = simClusters[pfs[0]].energy()
            relDiff_sim = (energy_to_calibrate - simulated_energy)/simulated_energy # relative diff. in energy between Sim clusters with/without filtering
            DR_gen_sim = math.sqrt( (genParticle.eta() - objectsToCalibrate[pfs[0]].eta())**2 + (genParticle.phi() - objectsToCalibrate[pfs[0]].phi())**2 )
            # print some info
            if (verbosityLevel>=2): print genParticle.pid(), "\t", energy_to_calibrate, "\t", genParticle.energy(), "\t", relDiff_sim, "\t", DR_gen_sim
            # store info in lists for cases where gen particles and Sim cluster are matched, and where reco. E for Sim cluster is close to the GEN-one
            if (math.fabs(DR_gen_sim) < 0.1 and energy_to_calibrate > simulated_energy/20. and simulated_energy >= (1-relFractionE)*GEN_eng and simulated_energy <= (1+relFractionE)*GEN_eng):
                if (abs(genParticle.pid()) in s_all_pids):
                    if (verbosityLevel==1): print genParticle.pid(), "\t", energy_to_calibrate, "\t", genParticle.energy(), "\t", relDiff_sim, "\t", DR_gen_sim
                    obj_Eng_EngRelDiff[abs(genParticle.pid())].append((simClusters[pfs[0]].phi(), simClusters[pfs[0]].eta(), simulated_energy, energy_to_compare, energy_to_calibrate, 100*relDiff_sim ))

    # prepare lists for individual eta/phi bins
    dEOverE_pidSelected = []
    dEOverE_pidSelected_values = {}
    for etaBinName in etaBins:
        GEN_eta = "[{0:.3f} - {1:.1f}]".format(etaBins[etaBinName][0],etaBins[etaBinName][1])
        print "Exrtacting info for eta range ",GEN_eta
        histDict[etaBinName] = {}
        dEOverE_pidSelected_values[etaBinName] = {}
        for phiBinName in phiBins:
            GEN_phi = "[{0:.2f} - {1:.2f}]".format(phiBins[phiBinName][0],phiBins[phiBinName][1])
            print "Exrtacting info for phi range ",GEN_phi
            histDict[etaBinName][phiBinName] = {}
            dEOverE_pidSelected_values[etaBinName][phiBinName] = {}
            # histograms per pid, print some info, fill dE/E values for current eta/phi bin
            print "PDG ID", "\t", "Mean dE/E (%)", "\t", "\t", "eta", "\t\t", "phi"
            for pid in s_all_pids:
                if (pid==211 or pid==22 or pid==11 or pid==13):
                    ### get the 1D lists
                    # for energy of sim cluster
                    sim_Energy = [eng for (phi, eta, eng, engToCompare, engToCalibrate, relDiff) in obj_Eng_EngRelDiff[pid] if ((math.fabs(eta) >= etaBins[etaBinName][0] and math.fabs(eta) < etaBins[etaBinName][1]) and (phi >= phiBins[phiBinName][0] and phi < phiBins[phiBinName][1]))]
                    # for energy of pf clusters before corrections
                    obj_EnergyToCompare = [engToCompare for (phi, eta, eng, engToCompare, engToCalibrate, relDiff) in obj_Eng_EngRelDiff[pid] if ((math.fabs(eta) >= etaBins[etaBinName][0] and math.fabs(eta) < etaBins[etaBinName][1]) and (phi >= phiBins[phiBinName][0] and phi < phiBins[phiBinName][1]))]
                    # for energy of pf clusters after corrections
                    obj_EnergyToCalibrate = [engToCalibrate for (phi, eta, eng, engToCompare, engToCalibrate, relDiff) in obj_Eng_EngRelDiff[pid] if ((math.fabs(eta) >= etaBins[etaBinName][0] and math.fabs(eta) < etaBins[etaBinName][1]) and (phi >= phiBins[phiBinName][0] and phi < phiBins[phiBinName][1]))]
                    # for energy loss
                    obj_EnergyDiffFilter = [relDiff for (phi, eta, eng, engToCompare, engToCalibrate, relDiff) in obj_Eng_EngRelDiff[pid] if ((math.fabs(eta) >= etaBins[etaBinName][0] and math.fabs(eta) < etaBins[etaBinName][1]) and (phi >= phiBins[phiBinName][0] and phi < phiBins[phiBinName][1]))]
                    ### fill the hists
                    rangeGeV = GEN_eng*1.6
                    nbins = int(rangeGeV)
                    if (rangeGeV<50.): nbins = int(10*rangeGeV)
                    binsBoundariesX_eng = [nbins, 0, rangeGeV]
                    # for energy of sim cluster
                    histDict[etaBinName][phiBinName] = histValue1D(sim_Energy, histDict[etaBinName][phiBinName], tag = "sim_Energy_eta"+etaBinName+"_phi"+phiBinName+"_pid"+str(pid), title = "Energy of the Sim cluster (matched to pid="+str(pid)+", GEN-level "+GEN_partId+" "+GEN_pTEng+", #eta="+GEN_eta+", #phi="+GEN_phi+")",   axunit = "E[GeV]", binsBoundariesX = binsBoundariesX_eng, ayunit = "N(clusters)")
                    # for energy of pf clusters before corrections
                    histDict[etaBinName][phiBinName] = histValue1D(obj_EnergyToCompare, histDict[etaBinName][phiBinName], tag = "obj_EnergyToCompare_eta"+etaBinName+"_phi"+phiBinName+"_pid"+str(pid), title = "Energy of the PF cluster before correction (cluster matched to pid="+str(pid)+", GEN-level "+GEN_partId+" "+GEN_pTEng+", #eta="+GEN_eta+", #phi="+GEN_phi+")",   axunit = "E[GeV]", binsBoundariesX = binsBoundariesX_eng, ayunit = "N(clusters)")
                    # for energy of pf clusters after corrections
                    histDict[etaBinName][phiBinName] = histValue1D(obj_EnergyToCalibrate, histDict[etaBinName][phiBinName], tag = "obj_EnergyToCalibrate_eta"+etaBinName+"_phi"+phiBinName+"_pid"+str(pid), title = "Energy of the PF cluster after correction (cluster matched to pid="+str(pid)+", GEN-level "+GEN_partId+" "+GEN_pTEng+", #eta="+GEN_eta+", #phi="+GEN_phi+")",   axunit = "E[GeV]", binsBoundariesX = binsBoundariesX_eng, ayunit = "N(clusters)")
                    # for energy loss from filter
                    binsBoundariesX_engDiff = [[800, -100, 60],[650, -80, 50]]["1p" in etaBinName]
                    histDict[etaBinName][phiBinName] = histValue1D(obj_EnergyDiffFilter, histDict[etaBinName][phiBinName], tag = "obj_EnergyDiffFilter_eta"+etaBinName+"_phi"+phiBinName+"_pid"+str(pid), title = "Energy loss from PF hits filtering (cluster matched to pid="+str(pid)+", GEN-level "+GEN_partId+" "+GEN_pTEng+", #eta="+GEN_eta+", #phi="+GEN_phi+")",   axunit = "#DeltaE_{clust}/E_{clust}[%]", binsBoundariesX = binsBoundariesX_engDiff, ayunit = "N(clusters)")
                    # extract mean and std for eng diff
                    h_tempEdiff = histDict[etaBinName][phiBinName]["obj_EnergyDiffFilter_eta"+etaBinName+"_phi"+phiBinName+"_pid"+str(pid)]
                    hEntries = h_tempEdiff.GetEntries()
                    gMean = h_tempEdiff.GetMean()
                    gMeanError = h_tempEdiff.GetMeanError()
                    gStd = h_tempEdiff.GetRMS()
                    if (hEntries > 100): # extract mean/error from fit if enough statistics
                        (h_tempEdiff, gMean, gStd) = fitGauss(h_tempEdiff)
                        gMeanError = gStd/(hEntries**0.5)
                    effSigma = getEffSigma(h_tempEdiff)
                    #print table for rel.diff.
                    print pid, "\t", "{0:.2f}".format(gMean) + " +/- " + "{0:.2f}".format(gMeanError), "\t", GEN_eta, "\t", GEN_phi, "\t", "(pf clusters: "+str(int(hEntries))+")", "\t", "effSigma: ", effSigma, "\t", "gStd: ", gStd
                    # for pidSelected only: overlap hists together, fill list with phi, eta and dE/E values
                    if (pid == pidSelected):
                        # overlap three hists together
                        h_tempE = histDict[etaBinName][phiBinName]["sim_Energy_eta"+etaBinName+"_phi"+phiBinName+"_pid"+str(pid)]
                        h_tempEToCalibrate = histDict[etaBinName][phiBinName]["obj_EnergyToCalibrate_eta"+etaBinName+"_phi"+phiBinName+"_pid"+str(pid)]
                        h_tempEToCompare = histDict[etaBinName][phiBinName]["obj_EnergyToCompare_eta"+etaBinName+"_phi"+phiBinName+"_pid"+str(pid)]
                        gMeanCalibEnergy = (gMean/100. + 1.) * GEN_eng
                        gMeanCalibEnergyError = (gMeanError/100.) * GEN_eng
                        # prepare histograms and their properties. this should depend on the input option which energies shpould be extracted/calibrated/compared
                        # this is for PF clusters
#                        histsAndProps = {h_tempE:{"leg":"Sim cluster ("+GEN_pTEng+")","color":ROOT.kRed}, h_tempEToCalibrate:{"leg":"Corrected PF cluster","color":ROOT.kBlue}, h_tempEToCompare:{"leg":"Non-corrected PF cluster","color":ROOT.kGreen-6}}
                        # this is for PF multiclusters
#                        histsAndProps = {h_tempE:{"leg":"Sim cluster ("+GEN_pTEng+")","color":ROOT.kRed}, h_tempEToCalibrate:{"leg":"Multicluster (uncorrected)","color":ROOT.kBlue}, h_tempEToCompare:{"leg":"PF cluster (corrected)","color":ROOT.kGreen-6}}
                        # this is for SIM clusters (3 sigma, 5 sigma, realistic, etc.)
#                        histsAndProps = {h_tempE:{"leg":"Sim cluster ("+GEN_pTEng+")","color":ROOT.kRed}, h_tempEToCalibrate:{"leg":"SIM cluster (detids matched to rechits, 5#sigma above noise)","color":ROOT.kBlue}, h_tempEToCompare:{"leg":"PF cluster (corrected)","color":ROOT.kGreen-6}}
                        # this is for different geometries
                        histsAndProps = {h_tempE:{"leg":"Sim cluster ("+GEN_pTEng+")","color":ROOT.kRed}, h_tempEToCalibrate:{"leg":"SIM cluster (3#sigma noise filtering, realistic Si geometry)","color":ROOT.kBlue}, h_tempEToCompare:{"leg":"SIM cluster (3#sigma above noise filtering, idealistic/CMSSW Si geometry)","color":ROOT.kGreen-6}}
                        # plot these histograms on top of each other (for each eta/phi bin)
                        histsPrintSaveSameCanvas(histsAndProps, outDir, tag = "obj_EnergyHistsOverlapped_eta"+etaBinName+"_phi"+phiBinName+"_pid"+str(pid), latexComment = ["E_{fit} = "+"{0:.2f}".format(gMeanCalibEnergy) + " #pm " + "{0:.2f}".format(gMeanCalibEnergyError)+" GeV","#sigma_{fit} = "+ "{0:.1f}".format(gStd)+" %, #sigma_{eff} = "+ "{0:.1f}".format(effSigma)+" %"], funcsAndProps = None)
                        #fill list with phi, eta and dE/E values
                        #skip filling the lists for overll hists in case of repeating inclusive bins
                        if not (len(etaBins.keys())>1 and etaBinName=="eta1p479to3p0") and not (len(phiBins.keys())>1 and phiBinName=="phim1p0pito1p0pi"):
                            dEOverE_pidSelected.append((sum(phiBins[phiBinName])/2.,sum(etaBins[etaBinName])/2.,gMean))
                        dEOverE_pidSelected_values[etaBinName][phiBinName]['mean'] = gMean
                        dEOverE_pidSelected_values[etaBinName][phiBinName]['meanError'] = gMeanError
                        dEOverE_pidSelected_values[etaBinName][phiBinName]['effSigma'] = effSigma
            # plos histograms for all pids in current eta/phi bin (if needed)
            # histPrintSaveAll(histDict[etaBinName][phiBinName], outDir, "_eta"+etaBinName+"_phi"+phiBinName)

    # global 2D histogram
    histDictGlob = {}
    # get edge boudaries for eta and phi
    bns_phi = list(set([item for sublist in [[phiBins[phiBin][0],phiBins[phiBin][1]] for phiBin in phiBins] for item in sublist]))
    bns_phi.sort(key=lambda x: x)
    bns_eta = list(set([item for sublist in [[etaBins[etaBin][0],etaBins[etaBin][1]] for etaBin in etaBins] for item in sublist]))
    bns_eta.sort(key=lambda x: x)
    # plot histograms for 3D list (phi, eta, relDiff)
    #simClust_EnergyDiffFilter_FullpidSelected = [(phi, eta, relDiff) for (phi, eta, eng, engToCalibrate, engAbove, relDiff) in simClust_Eng_EngRelDiff[pidSelected]]
    #histDictGlob = histValues2D(simClust_EnergyDiffFilter_FullpidSelected, histDictGlob, tag = "dEOverE_pidSelected_weigthed2D_full", title = "Relative energy loss #DeltaE_{clust}/E_{clust} from hits filtering (%)", axunit = "#phi", binsBoundariesX = [len(bns_phi)-1, bns_phi], ayunit = "#eta", binsBoundariesY = [len(bns_eta)-1, bns_eta], weighted2D = True)
    histDictGlob = histValues2D(dEOverE_pidSelected, histDictGlob, tag = "dEOverE_pidSelected_weigthed2D", title = "Relative energy loss #DeltaE_{clust}/E_{clust} from hits filtering (%)", axunit = "#phi", binsBoundariesX = [len(bns_phi)-1, bns_phi], ayunit = "#eta", binsBoundariesY = [len(bns_eta)-1, bns_eta], weighted2D = True)
    recHits_Filtered = [(x, y, int(filt)) for (x, y, eng, filt) in hitsSelected]
    histDictGlob = histValues2D(recHits_Filtered, histDictGlob, tag = "recHits_filteredByNoise", title = "Spatial distribution of hits filtered by noise", axunit = "x[mm]", binsBoundariesX = [3200, -1600, 1600], ayunit = "y[mm]", binsBoundariesY = [3200, -1600, 1600], weighted2D = True)
    histPrintSaveAll(histDictGlob, outDir, "_glob")
    
    # print 2D eta-phi table
    print "Effect on the cluster energy (relative loss in %):"
    print "eta\phi","\t","\t".join(["[{0:.2f} - {1:.2f}]".format(phiBins[phiBinName][0],phiBins[phiBinName][1]) for phiBinName in phiBins])
    for etaBinName in etaBins:
        print "[{0:.3f} - {1:.1f}]".format(etaBins[etaBinName][0],etaBins[etaBinName][1]), "\t", "\t".join(["{0:.2f}".format(dEOverE_pidSelected_values[etaBinName][phiBinName]['mean']) + " +/- " + "{0:.2f}".format(dEOverE_pidSelected_values[etaBinName][phiBinName]['meanError']) for phiBinName in phiBins])

    # print 2D eta-phi table
    print "Effect on the cluster energy (sigma effective in %):"
    print "eta\phi","\t","\t".join(["[{0:.2f} - {1:.2f}]".format(phiBins[phiBinName][0],phiBins[phiBinName][1]) for phiBinName in phiBins])
    for etaBinName in etaBins:
        print "[{0:.3f} - {1:.1f}]".format(etaBins[etaBinName][0],etaBins[etaBinName][1]), "\t", "\t".join(["{0:.2f}".format(dEOverE_pidSelected_values[etaBinName][phiBinName]['effSigma']) for phiBinName in phiBins])

    # print 2D eta-phi table
    print "Effect on the cluster energy (value of calibration const.):"
    print "eta\phi","\t","\t".join(["[{0:.2f} - {1:.2f}]".format(phiBins[phiBinName][0],phiBins[phiBinName][1]) for phiBinName in phiBins])
    for etaBinName in etaBins:
        print "[{0:.3f} - {1:.1f}]".format(etaBins[etaBinName][0],etaBins[etaBinName][1]), "\t", "\t".join(["{0:.2f}".format(1/(1+dEOverE_pidSelected_values[etaBinName][phiBinName]['mean']/100.)) + " -/+ " + "{0:.2f}".format((1/((1+dEOverE_pidSelected_values[etaBinName][phiBinName]['mean']/100.)**2)) * dEOverE_pidSelected_values[etaBinName][phiBinName]['meanError']) for phiBinName in phiBins])

if __name__ == '__main__':
    main()

# plot graphs with energy resolution for various scenarios (currently needs to be invoked manually, one all the energy points have been processed)
def graphTest():
    graphsAndProps = {}
    xEng = array('f',[5., 20., 50., 100., 300.])

# photons, ideal vs. realistic geometry (matched/filtred simHits)
    # eta [1.479 - 3.0] pfClusters
    yRes = array('f',[9.80, 5.40, 3.60, 2.80, 2.00])
    graphsAndProps[ROOT.TGraph( len(xEng), xEng, yRes )] = {"leg":"PF clusters (realistic SIM clusters, ideal geometry)", "color":ROOT.kBlue, "MarkerStyle":21, "LineStyle":4}
    # eta [1.479 - 3.0] matched/filtred simHits, 3 sigma, ideal geometry
    yRes = array('f',[9.90, 5.40, 3.50, 2.60, 1.50])
    graphsAndProps[ROOT.TGraph( len(xEng), xEng, yRes )] = {"leg":"SIM clusters (noise filtered @3#sigma, ideal geometry)", "color":ROOT.kGreen-6, "MarkerStyle":21, "LineStyle":1}
    # eta [1.479 - 3.0] matched/filtred simHits, 3 sigma, realistic geometry
    yRes = array('f',[11.50, 6.30, 4.20, 3.10, 1.90])
    graphsAndProps[ROOT.TGraph( len(xEng), xEng, yRes )] = {"leg":"SIM clusters (noise filtered @3#sigma, realistic geometry)", "color":ROOT.kRed-6, "MarkerStyle":26, "LineStyle":1}
    drawGraphsTest(graphsAndProps, title = "Energy resolution (photons, #eta #in [1.479 - 3.0], PU=0)", tag="test_photonsResolutionComparisonGeo")

