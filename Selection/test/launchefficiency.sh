#!/bin/sh

cwd=`pwd`

#list=`find /gpfs/grid/srm/cms/store/data/Run2011A/SingleElectron/AOD/May10ReReco-v1/ -type f`
#list=`find /gpfs/grid/srm/cms/store/data/Run2011A/SingleElectron/RAW-RECO/WElectron-May10ReReco-v1/ -type f`
list=`find /gpfs/grid/srm/cms/store/data/Run2011A/DoubleElectron/RAW-RECO/ZElectron-May10ReReco-v1 -type f`
#list=`find /gpfs/grid/srm/cms/store/data/Run2011A/DoubleElectron/RAW-RECO/ZElectron-05Aug2011-v1 -type f`
#list='/gpfs/grid/srm/cms/store/data/Run2011A/DoubleElectron/RAW-RECO/ZElectron-May10ReReco-v1/0000/EC04DB0E-D97D-E011-B60B-00261894396F.root'
#list=`find /gpfs/grid/srm/cms/store/data/Run2011A/DoubleElectron/RAW-RECO/ZElectron-PromptSkim-v4 -type f`
#list=`find /gpfs/grid/srm/cms/store/data/Run2011A/DoubleElectron/RAW-RECO/ZElectron-PromptSkim-v6 -type f`

i=0

for f in $list; do

  i=`expr $i + 1`

  f=`echo $f | sed -e 's;/gpfs/grid/srm/cms;;'`

  echo $i $f

#  echo $workdir

  echo "process.source.fileNames = ['"$f"']" >> $cwd/zefficiency.py
  sleep 5
  #queue=normal
  #queue=ts_cms
  queue=short
  bsub -q $queue -J $f -e job.log -o job.log cmsRun $cwd/zefficiency.py

  cd $cwd

done

