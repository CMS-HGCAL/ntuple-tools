import ROOT
import os
import sys

eventsPerFile = 250


def splitFile(inDir, oldFileName):

    # Get old file, old tree and set top branch address
    dirname = "ana"
    prefix = "root://eoscms.cern.ch/"
    newFileName = oldFileName.strip(".root")
    oldFileFullName = "{}/{}/{}".format(prefix, inDir, oldFileName)

    print "Splitting", oldFileFullName

    oldfile = ROOT.TFile(oldFileFullName)
    oldtree = oldfile.Get("{}/hgc".format(dirname))
    nentries = oldtree.GetEntries()
    if nentries > eventsPerFile:
        oldtree.SetBranchStatus("*", 1)
        # Create a new file + a clone of old tree in new file
        index = 0
        for i in range(0, nentries, eventsPerFile):
            print "Event", i, "of", nentries
            newfile = ROOT.TFile("{}/{}/{}_{}.root".format(prefix, inDir, newFileName, index), "recreate")
            newfile.mkdir(dirname)
            newfile.cd(dirname)
            newtree = oldtree.CopyTree("", "", eventsPerFile, i)
            newtree.AutoSave()
            newfile.Write()
            newfile.Close()
            index += 1

        print "Deleting {}/{}".format(inDir, oldFileName)
        os.remove("{}/{}".format(inDir, oldFileName))
    else:
        print "No need to split", oldFileFullName


def main():

    inDir = sys.argv[1]
    # /eos/cms/store/cmst3/group/hgcal/CMG_studies/Production/_SinglePiPt10Eta1p6_2p8_PhaseIITDRFall17DR-noPUFEVT_93X_upgrade2023_realistic_v2-v1_GEN-SIM-RECO/NTUP
    for rootFile in os.listdir(inDir):
        splitFile(inDir, rootFile)


if __name__ == '__main__':
    main()
