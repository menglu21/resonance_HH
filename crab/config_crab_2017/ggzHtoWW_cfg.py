from WMCore.Configuration import Configuration

config = Configuration()

config.section_("General")
config.General.requestName = 'ggZHtoWW'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.allowUndistributedCMSSW = True
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'PSet.py'
config.JobType.scriptExe = '2017_script/crab_script.sh'
# hadd nano will not be needed once nano tools are in cmssw
config.JobType.inputFiles = ['crab_script.py', '../scripts/haddnano.py','keep_and_drop.txt']
config.JobType.sendPythonFolder = True

config.section_("Data")
config.Data.inputDataset = '/GluGluZH_HToWWTo2L2Nu_M125_13TeV_powheg_pythia8_TuneCP5_PSweight/sRunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
#config.Data.outLFNDirBase = '/store/group/phys_top/ExtraYukawa/TTC_version9/'
config.Data.unitsPerJob = 1
config.Data.totalUnits = -1
config.Data.publication = False
config.Data.outputDatasetTag = 'ggZHtoWW'

config.section_("Site")
config.Site.storageSite = "T3_CH_CERNBOX"
