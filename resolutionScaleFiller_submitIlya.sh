#!/bin/bash

# gorbunov_PiPt25_PUP200

SAMPLEFILE=$1

OUTDIR="/eos/user/c/clange/HGCal/ScaleResolution/_new/Ilya/"
PREFIX="root://eoscms.cern.ch//eos/cms/store/cmst3/group/hgcal/CMG_studies/Production/"
GUNTYPE="pt"
POSTFIX="NTUP"
REFNAME="genpart"
OBJNAMES="pfcluster megacluster"
QUEUE="8nh"

for SAMPLE in `cat $SAMPLEFILE`; do
  SAMPLEDIR="${PREFIX}${SAMPLE}/${POSTFIX}/"
  PTVAL=`echo ${SAMPLE} | gawk 'match($0, /.*gorbunov_PiPt(.*)_PUP.*/, arr) { print arr[1] }'`
  PARTICLE=`echo ${SAMPLE} | gawk 'match($0, /.*gorbunov_(.*)Pt.*/, arr) { print arr[1] }'`
  TAG="PU200"
  PID=0
  if [ "$PARTICLE" == "Pi" ]; then
    PID=211
  elif [ "$PARTICLE" == "Photon" ]; then
    PID=22
  fi;
  for OBJNAME in $OBJNAMES; do
    echo python resolutionScaleFiller_batchWrapper.py --inputdir ${SAMPLEDIR} --outdir ${OUTDIR} --gunType $GUNTYPE --pid $PID --genValue $PTVAL --tag $TAG --ref $REFNAME --obj $OBJNAME --queue $QUEUE
    python resolutionScaleFiller_batchWrapper.py --inputdir ${SAMPLEDIR} --outdir ${OUTDIR} --gunType $GUNTYPE --pid $PID --genValue $PTVAL --tag $TAG --ref $REFNAME --obj $OBJNAME --queue $QUEUE
  done
done
