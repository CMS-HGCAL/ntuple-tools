#!/bin/bash
## External vars
currentDir=${1}
CMSSWDIR=${2}
CMSSWVER=${3}
CMSSWARCH=${4}
fileList=${5}
outDir=${6}
sampleName=${7}
gunType=${8}
pid=${9}
genValue=${10}
tag=${11}
refName=${12}
objName=${13}

##Create Work Area
cd $TMPDIR || exit
export SCRAM_ARCH=${CMSSWARCH}
source /afs/cern.ch/cms/cmsset_default.sh
eval `scramv1 project CMSSW ${CMSSWVER}`
cd ${CMSSWVER}/ || exit
# this is dangerous!
rm -rf ./*
cp -r -d ${CMSSWDIR}/* ./
cd src || exit
eval `scramv1 runtime -sh`
edmPluginRefresh -p ../lib/$SCRAM_ARCH

## Execute job and retrieve the outputs
echo "Job running on `hostname` at `date`"
cd $TMPDIR || exit
cp -r $currentDir analysis
cd analysis || exit


#cd RecoNtuples/HGCalAnalysis/test/HGCalAnalysis
pwd
echo python resolutionScaleFiller.py --files $fileList --gunType $gunType --pid $pid --genValue $genValue --tag $tag --ref $refName --obj $objName
python resolutionScaleFiller.py --files $fileList --gunType $gunType --pid $pid --genValue $genValue --tag $tag --ref $refName --obj $objName
ls -l

# copy to outDir
mkdir -p $outDir/${sampleName%_*}
cp $outfilename $outDir/${sampleName%_*}/
