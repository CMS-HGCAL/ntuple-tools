"""batch submission wrapping script"""
import os
# from SampleHelper import SampleManager
import logging
import commands
import time
import optparse


def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in xrange(0, len(l), n):
        yield l[i:i + n]


def processCmd(cmd, quiet=True):
    if not quiet:
        print cmd
    status, output = commands.getstatusoutput(cmd)
    if (status != 0 and not quiet):
        print 'Error in processing command:\n   ['+cmd+']'
        print 'Output:\n   ['+output+'] \n'
    return output


def main():
    """Main function and settings."""

    global opt, args

    usage = ('usage: %prog [options]\n' + '%prog -h for help')
    parser = optparse.OptionParser(usage)

    # input options
    parser.add_option('', '--inputdir', dest='inputdir', type='string',  default='root://eoscms.cern.ch//eos/cms/store/cmst3/group/hgcal/CMG_studies/Production/_RelValSingleGammaPt25Eta1p7_2p7_CMSSW_9_3_2-PU25ns_93X_upgrade2023_realistic_v2_2023D17PU200-v1_GEN-SIM-RECO/NTUP/', help='path to the input directory')
    parser.add_option('', '--outdir', dest='outDir', type='string',  default=os.getcwd(), help='output directory')
    parser.add_option('', '--gunType', dest='gunType', type='string',  default='pt', help='pt or e')
    parser.add_option('', '--pid', dest='pid', type='int',  default=22, help='pdgId int')
    parser.add_option('', '--genValue', dest='genValue', type='float',  default=25, help='generated pT or energy')
    parser.add_option('', '--tag', dest='tag', type='string',  default='noPU', help='some tag, best used for PU and other info')
    parser.add_option('', '--ref', dest='refName', type='string',  default='genpart', help='reference collection')
    parser.add_option('', '--obj', dest='objName', type='string',  default='pfcluster', help='object of interest collection')
    parser.add_option('', '--queue', dest='queue', type='string',  default='8nh', help='queue to run on')

    # store options and arguments as global variables
    global opt, args
    (opt, args) = parser.parse_args()

    queue = opt.queue
    # save working dir
    outDir = opt.outDir
    currentDir = os.getcwd()
    CMSSW_BASE = os.getenv('CMSSW_BASE')
    CMSSW_VERSION = os.getenv('CMSSW_VERSION')
    SCRAM_ARCH = os.getenv('SCRAM_ARCH')
    # geometryFile = currentDir + "/v33-withBH.txt"

    print "inputdir:", opt.inputdir
    print "outdir:", opt.outDir
    print "gunType:", opt.gunType
    print "pid:", opt.pid
    print "genValue:", opt.genValue
    print "tag:", opt.tag
    print "ref:", opt.refName
    print "obj:", opt.objName
    eosdir = opt.inputdir[21:]
    onlyfiles = [f for f in os.listdir(eosdir) if os.path.isfile(os.path.join(eosdir, f))]
    # nEvents = -1
    defaultFilesPerJob = 1
    batchScript = currentDir + '/' + 'resolutionScaleFiller_batchScript.sh'

    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(message)s')
    ch.setFormatter(formatter)

    sampleName = opt.gunType + str(opt.genValue)+'_pid_'+str(opt.pid) + '_ref_'+opt.refName+'_obj_'+opt.objName
    if not os.path.exists('./' + sampleName):
        os.makedirs('./' + sampleName)
    if not os.path.exists('./' + sampleName + '/std'):
        os.makedirs('./' + sampleName + '/std')

    for fileIndex, filename in enumerate(onlyfiles):
        filename = opt.inputdir + filename
        outfilename = opt.gunType + str(opt.genValue)+'_pid_'+str(opt.pid) + '_ref_'+opt.refName+'_obj_'+opt.objName+'_'+str(fileIndex)+'.root'
        subSampleName = opt.gunType + str(opt.genValue)+'_pid_'+str(opt.pid) + '_ref_'+opt.refName+'_obj_'+opt.objName+'_'+str(fileIndex)

        stdOut = "./{sampleName}/std/{subSampleName}.out".format(sampleName=sampleName, subSampleName=subSampleName)
        print stdOut
        stdErr = "./{sampleName}/std/{subSampleName}.err".format(sampleName=sampleName, subSampleName=subSampleName)
        print stdErr
        print filename
        jobName = str(fileIndex)+'_'+str(time.time())
        tag = "{}_{}".format(opt.tag, fileIndex)
        cmd = 'bsub -o {stdOut} -e {stdErr} -q {queue} -J {jobName} {batchScript} {p1} {p2} {p3} {p4} \"{p5}\" {p6} {p7} {p8} {p9} {p10} {p11} {p12} {p13} {p14}'.format(stdOut=stdOut, stdErr=stdErr, queue=queue, jobName=jobName, batchScript=batchScript, p1=currentDir, p2=CMSSW_BASE, p3=CMSSW_VERSION, p4=SCRAM_ARCH, p5=filename, p6=outDir, p7=sampleName, p8=outfilename, p9=opt.gunType, p10=opt.pid, p11=opt.genValue, p12=tag, p13=opt.refName, p14=opt.objName)
        print cmd
        processCmd(cmd)
        time.sleep(1)


if __name__ == '__main__':
    main()
