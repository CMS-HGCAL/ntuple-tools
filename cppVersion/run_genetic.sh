#!/bin/bash

echo "Starting Genetic optimizer"

export LD_LIBRARY_PATH=/cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.16.00/x86_64-centos7-gcc48-opt/lib:/cvmfs/sft.cern.ch/lcg/contrib/gcc/4.8/x86_64-centos7-gcc48-opt/lib64:/cvmfs/cms.cern.ch/slc6_amd64_gcc630/external/llvm/5.0.0-omkpbe2/lib64:/cvmfs/cms.cern.ch/slc6_amd64_gcc630/external/gcc/6.3.0/lib64:/cvmfs/cms.cern.ch/slc6_amd64_gcc630/external/gcc/6.3.0/lib:/cvmfs/sft.cern.ch/lcg/releases/pcre/8.38-a0fda/x86_64-centos7-gcc62-opt/lib/

#/afs/cern.ch/work/j/jniedzie/private/CMSSW_10_2_0_pre1/biglib/slc6_amd64_gcc630:
#/afs/cern.ch/work/j/jniedzie/private/CMSSW_10_2_0_pre1/lib/slc6_amd64_gcc630:
#/afs/cern.ch/work/j/jniedzie/private/CMSSW_10_2_0_pre1/external/slc6_amd64_gcc630/lib:
#/cvmfs/cms.cern.ch/slc6_amd64_gcc630/cms/cmssw/CMSSW_10_2_0_pre1/biglib/slc6_amd64_gcc630:
#/cvmfs/cms.cern.ch/slc6_amd64_gcc630/cms/cmssw/CMSSW_10_2_0_pre1/lib/slc6_amd64_gcc630:
#/cvmfs/cms.cern.ch/slc6_amd64_gcc630/cms/cmssw/CMSSW_10_2_0_pre1/external/slc6_amd64_gcc630/lib:


source /cvmfs/sft.cern.ch/lcg/contrib/gcc/4.8/x86_64-centos7-gcc48-opt/setup.sh
source /cvmfs/sft.cern.ch/lcg/releases/pcre/8.38-a0fda/x86_64-centos7-gcc62-opt/pcre-env.sh
source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.16.00/x86_64-centos7-gcc48-opt/bin/thisroot.sh

echo "Which root:`which root`"
echo "Lib path:`echo $LD_LIBRARY_PATH`"

cd /afs/cern.ch/work/j/jniedzie/private/hgcal/ntuple-tools/cppVersion/

echo "Im in `pwd`"
./geneticOptimizer
