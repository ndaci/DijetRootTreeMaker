from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.workArea = 'WORKINGAREA'
config.General.requestName = 'WORKINGDIR'
config.section_('JobType')
config.JobType.psetName = 'CMSSWCFG'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['OUTFILENAME']
config.section_('Data')
config.Data.inputDataset = 'INPUTDATASET'
config.Data.unitsPerJob = FILESPERJOB #without '' since it must be an int
config.Data.splitting = 'LumiBased'
config.Data.publication = False
config.Data.lumiMask = '/afs/cern.ch/work/m/mdiamant/private/CMSSW_9_2_3_patch2/CMSSW_9_2_3_patch2_JetID/src/CMSDIJET/DijetRootTreeMaker/prod/submitJobsWithCrab3/json_7july_DCSONLY.txt'
#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/Cert_254833_13TeV_PromptReco_Collisions15_JSON.txt'
#config.Data.lumiMask = 'JSON_Certified/Cert_271036-278808_13TeV_PromptReco_Collisions16_JSON_NoL1T_19Aug16.txt'
#config.Data.runRange = '254833-254833'#'208306-238354' # '193093-194075'                             
config.Data.outLFNDirBase = '/store/group/phys_jetmet/magda/bigtrees_Data2017_Run2017B-PromptReco-v1_newtriggers/'
#config.Data.outLFNDirBase = '/store/user/santanas/rootTrees/Run2015B_JetHT_10June2015DCSJson_JECV5_5a70fc3/'
config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'
#config.Site.storageSite = 'T2_IT_Rome'
