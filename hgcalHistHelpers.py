import ROOT
from array import array


def histValue1D(fValues, histDict, tag="hist1D_", title="hist 1D", axunit="a.u.", binsBoundariesX=[10, -1, 1], ayunit="a.u.", verbosityLevel=0):
    """1D histograming of given list of values."""
    # sanity check for hists
    if histDict is None:
        return
    # sanity check for boundaries
    if (len(binsBoundariesX) != 3 and len(binsBoundariesX) != 2):
        return
    # define event-level hists
    elif len(binsBoundariesX) == 3:  # bondaries in format [nbins, low, high]
        histDict[tag] = ROOT.TH1F(tag, title + ";" + axunit + ";" + ayunit, binsBoundariesX[0], binsBoundariesX[1], binsBoundariesX[2])
    elif len(binsBoundariesX) == 2:  # bondaries in format [nbins, list_boundaries]
        histDict[tag] = ROOT.TH1F(tag, title + ";" + axunit + ";" + ayunit, binsBoundariesX[0], array('f', binsBoundariesX[1]))
    # set some properties
    histDict[tag].GetYaxis().SetTitleOffset(histDict[tag].GetYaxis().GetTitleOffset() * 3.0)
    # loop over all values
    if (verbosityLevel >= 3):
        print "tag: ", tag, ", fValues: ", fValues
    for value in fValues:
        histDict[tag].Fill(value)
    return histDict


def histValues2D(fValues, histDict, tag="hist2D_", title="hist 2D", axunit="a.u.", binsBoundariesX=[10, -1, 1], ayunit="a.u.", binsBoundariesY=[10, -1, 1], weighted2D=False, verbosityLevel=0):
    """2D histograming of given list of values"""
    # sanity check for hists
    if histDict is None:
        return
    # sanity check for boundaries
    if (len(binsBoundariesX) != len(binsBoundariesY)):
        return
    if (len(binsBoundariesX) != 3 and len(binsBoundariesX) != 2):
        return
    # define event-level hists
    elif len(binsBoundariesX) == 3:  # bondaries in format [nbins, low, high]
        histDict[tag] = ROOT.TH2F(tag, title + ";" + axunit + ";" + ayunit, binsBoundariesX[0], binsBoundariesX[1], binsBoundariesX[2], binsBoundariesY[0], binsBoundariesY[1], binsBoundariesY[2])
    elif len(binsBoundariesY) == 2:  # bondaries in format [nbins, list_boundaries]
        histDict[tag] = ROOT.TH2F(tag, title + ";" + axunit + ";" + ayunit, binsBoundariesX[0], array('f', binsBoundariesX[1]), binsBoundariesY[0], array('f', binsBoundariesY[1]))
    # set some properties
    histDict[tag].GetXaxis().SetTitleOffset(histDict[tag].GetXaxis().GetTitleOffset() * 1.0)
    histDict[tag].GetYaxis().SetTitleOffset(histDict[tag].GetYaxis().GetTitleOffset() * 3.0)
    # loop over all values
    if (verbosityLevel >= 3):
        print "tag: ", tag, ", fValues: ", fValues
    if (not weighted2D):
        for (valueX, valueY) in fValues:
            histDict[tag].Fill(valueX, valueY)
    else:
        for (valueX, valueY, valueZ) in fValues:
            histDict[tag].Fill(valueX, valueY, valueZ)
    return histDict


def histsPrintSaveSameCanvas(histsAndProps, outDir, tag="hists1D_", latexComment="", funcsAndProps=None, verbosityLevel=0):
    """print/save list of histograms with their properties on one canvas"""
    # supress info messages
    ROOT.gErrorIgnoreLevel = ROOT.kInfo + 1
    # set default style values
    ROOT.gStyle.SetPalette(ROOT.kBird)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPadTopMargin(0.08)
    ROOT.gStyle.SetPadBottomMargin(0.12)
    ROOT.gStyle.SetPadLeftMargin(0.12)
    ROOT.gStyle.SetPadRightMargin(0.05)
    # create canvas
    canvas = ROOT.TCanvas(outDir + tag, outDir + tag, 500, 500)
    # prepare the legend
    leg = ROOT.TLegend(0.15, 0.75, 0.62, 0.9)
    # leg.SetHeader("Energy of the clusters before/after filtering")
    leg.SetBorderSize(0)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    # prepare latex comment
    ltx = ROOT.TLatex()
    ltx.SetNDC(ROOT.kTRUE)
    ltx.SetTextFont(42)
    ltx.SetTextSize(0.03)
    # set image extensions
    imgTypes = ["pdf", "png"]
    if (verbosityLevel >= 3):
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
        if (curr_max < hist.GetEntries() / 2.):
            y_maxs.append(curr_max)
    # print "y_maxs: ", y_maxs
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
            hist.GetXaxis().SetTitleOffset(hist.GetXaxis().GetTitleOffset() * 1.2)
            hist.GetYaxis().SetTitleOffset(hist.GetYaxis().GetTitleOffset() * 3.0)
            if (first):
                hist.GetYaxis().SetRangeUser(0, max(y_maxs) * 1.4)
                hist.Draw("hist0 goff")
                first = False
            else:
                hist.Draw("hist0 same goff")
    # check if any function should be drawn
    if funcsAndProps is not None:
        for func in funcsAndProps:
            func.SetLineColor(funcsAndProps[func]["color"])
            leg.AddEntry(func, funcsAndProps[func]["leg"], "L")
            func.Draw("same goff")
    # draw the rest
    leg.Draw("same")
    ltx.DrawLatex(0.170, 0.71, latexComment[0])
    ltx.DrawLatex(0.170, 0.66, latexComment[1])
    for imgType in imgTypes:
        canvas.SaveAs("{}/{}.{}".format(outDir, tag, imgType))
    return canvas


def drawGraphsTest(graphsAndProps, title="Resolution", tag="graphTest_", verbosityLevel=0):
    # supress info messages
    ROOT.gErrorIgnoreLevel = ROOT.kInfo + 1
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
    leg = ROOT.TLegend(0.55, 0.45, 0.9, 0.9)
    leg.SetHeader(title)
    leg.SetBorderSize(0)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    # prepare latex comment
    ltx = ROOT.TLatex()
    ltx.SetNDC(ROOT.kTRUE)
    ltx.SetTextFont(42)
    ltx.SetTextSize(0.03)
    # set image extensions
    imgTypes = ["pdf", "png", "root"]
    if (verbosityLevel >= 3):
        print "graphsAndProps: ", graphsAndProps
    # loop over all graphs to get max
    y_maxs = [gr.GetYaxis().GetXmax() for gr in graphsAndProps]
    y_mins = [gr.GetYaxis().GetXmin() for gr in graphsAndProps]
    print "y_mins: ", y_mins
    print "y_maxs: ", y_maxs
    # loop over all histograms
    first = True
    for gr in graphsAndProps:
        gr.GetXaxis().SetTitle('E[GeV]')
        gr.GetYaxis().SetTitle('#sigma_{eff}[%]')
        colour = graphsAndProps[gr]["color"]
        gr.SetLineColor(colour)
        gr.SetLineWidth(1)
        gr.SetLineStyle(graphsAndProps[gr]["LineStyle"])
        gr.SetMarkerColor(colour)
        gr.SetMarkerStyle(graphsAndProps[gr]["MarkerStyle"])
        gr.SetMarkerSize(0.7)
        gr.SetFillColor(0)
        gr.SetFillStyle(0)
        leg.AddEntry(gr, graphsAndProps[gr]["leg"])
        if (first):
            gr.SetMaximum(max(y_maxs) * 1.5)
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


def histPrintSaveAll(histDict, outDir, tag="_test", verbosityLevel=0):
    # supress info messages
    ROOT.gErrorIgnoreLevel = ROOT.kInfo + 1
    # set image extensions
    imgType = "pdf"
    if (verbosityLevel >= 3):
        print "histDict.items(): ", histDict.items()
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
            canvas = ROOT.TCanvas(outDir + tag + key, outDir + tag + key, 500, 500)
            item.Draw("hist0 goff")
            canvas.SaveAs("{}/{}.{}".format(outDir, key, imgType))
        if type(item) == ROOT.TH2F:
            ROOT.gStyle.SetPadRightMargin(0.15)
            canvas = ROOT.TCanvas(outDir + tag + key, outDir + tag + key, 500, 500)
            item.Draw("colz goff")
            canvas.SaveAs("{}/{}.{}".format(outDir, key, imgType))
        elif type(item) == ROOT.TH3F:
            ROOT.gStyle.SetPadRightMargin(0.05)
            canvas = ROOT.TCanvas(outDir + tag + key, outDir + tag + key, 500, 500)
            item.Draw("box goff")
            canvas.SaveAs("{}/{}.{}".format(outDir, key, imgType))
        else:
            continue
    # del canvas
    return


def fitGauss(hist, paramRangeFactor=1.8):
    if (hist.GetEntries() == 0):
        return (hist, 0, 0)
    hist.GetXaxis().SetTitleOffset(hist.GetXaxis().GetTitleOffset() * 1.2)
    hist.GetYaxis().SetTitleOffset(hist.GetYaxis().GetTitleOffset() * 3.0)
    # define the range of the fit from the hist mean and RMS
    meanLimitDn = hist.GetMean() - paramRangeFactor * hist.GetRMS()
    meanLimitUp = hist.GetMean() + paramRangeFactor * hist.GetRMS()
    sigmaLimitDn = hist.GetRMS() / paramRangeFactor
    sigmaLimitUp = hist.GetRMS() * paramRangeFactor
    # define the fitting gausian and range of its parameters
    fGauss = ROOT.TF1("f", "[0]*TMath::Gaus(x,[1],[2])", meanLimitDn, meanLimitUp)
    fGauss.SetParLimits(1, meanLimitDn, meanLimitUp)
    fGauss.SetParLimits(2, sigmaLimitDn, sigmaLimitUp)
    # perform fit and extract params
    hist.Fit(fGauss, "Q", "", meanLimitDn, meanLimitUp)
    gaussMean = fGauss.GetParameter(1)
    gaussStd = fGauss.GetParameter(2)
    # print "hist ent.: ", hist.GetEntries()
    # print "hist mean: ", hist.GetMean(), ", hist std: ", hist.GetRMS()
    # print "gaussMean: ", gaussMean,      ", gaussStd: ", gaussStd
    return (hist, gaussMean, gaussStd)


def getEffSigma(theHist, wmin=-100, wmax=100, epsilon=0.01):
    """taken from Hgg framework (by Ed)"""
    # initialise
    weight = 0.
    points = []
    thesum = theHist.Integral()
    # return -1 in case of empty histogram
    if (thesum == 0):
        return -1.
    # compute the cumulative distr. points
    for i in range(theHist.GetNbinsX()):
        weight += theHist.GetBinContent(i)
        if weight / thesum > epsilon:
            points.append([theHist.GetBinCenter(i), weight / thesum])
#    print "thesum: ", thesum
#    print "points: ", points
    # initialise
    low = wmin
    high = wmax
    width = wmax - wmin
    # find minimal 0.683 interval
    for i in range(len(points)):
        for j in range(i, len(points)):
            wy = points[j][1] - points[i][1]
            if abs(wy - 0.683) < epsilon:
                wx = points[j][0] - points[i][0]
                if wx < width:
                    low = points[i][0]
                    high = points[j][0]
                    width = wx
                    # print "width: ", width, "low: ", low, ", high: ", high, ", wx: ", wx, ", wy: ", wy,
    return 0.5 * (high - low)
