```
cmsrel CMSSW_12_1_1 
cd CMSSW_12_1_1/src/; 
cmsenv 
git cms-init 
git remote add ppalit-cmssw git@github.com:pritampalit/cmssw.git 
git remote -v 
git pull ppalit-cmssw dataMCRun3 
git cms-addpkg DQM/TrackingMonitorSource 
scram b

```
