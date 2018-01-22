#!/bin/bash

SAMPLEFILE=$1

OUTDIR="/eos/user/c/clange/HGCal/ScaleResolution/_new/"
PREFIX="root://eoscms.cern.ch//eos/cms/store/cmst3/group/hgcal/CMG_studies/Production/"
GUNTYPE="pt"
POSTFIX="NTUP"
REFNAME="genpart"
OBJNAMES="pfcluster megacluster"
QUEUE="8nh"

for SAMPLE in `cat $SAMPLEFILE`; do
  SAMPLEDIR="${PREFIX}${SAMPLE}/${POSTFIX}/"
  PARTICLE=`echo ${SAMPLE} | gawk 'match($0, /.*Single(.*)Pt.*/, arr) { print arr[1] }'`
  PTVAL=`echo ${SAMPLE} | gawk -v part=".*Single${PARTICLE}Pt(.*)Eta.*" 'match($0, part, arr) { print arr[1] }'`
  TAG=`echo ${SAMPLE} | gawk 'match($0, /.*Fall17DR-(.*)FEVT.*/, arr) { print arr[1] }'`
  PID=0
  if [ "$PARTICLE" == "Pi" ]; then
    PID=211
  elif [ "$PARTICLE" == "Gamma" ]; then
    PID=22
  fi;
  for OBJNAME in $OBJNAMES; do
    echo python resolutionScaleFiller_batchWrapper.py --inputdir ${SAMPLEDIR} --outdir ${OUTDIR} --gunType $GUNTYPE --pid $PID --genValue $PTVAL --tag $TAG --ref $REFNAME --obj $OBJNAME --queue $QUEUE
    python resolutionScaleFiller_batchWrapper.py --inputdir ${SAMPLEDIR} --outdir ${OUTDIR} --gunType $GUNTYPE --pid $PID --genValue $PTVAL --tag $TAG --ref $REFNAME --obj $OBJNAME --queue $QUEUE
  done
done
