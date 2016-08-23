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
#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/Cert_246908-248005_13TeV_PromptReco_Collisions15_ZeroTesla_JSON_CaloOnly.txt'
#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/DCSOnly/json_DCSONLY_Run2015B.txt'
config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/Cert_246908-251883_13TeV_PromptReco_Collisions15_JSON.txt'
#config.Data.lumiMask = 'json_new.txt' ## if you downloaded the file in the working directory
config.Data.runRange = '251586-251883' # '193093-194075'                             
#config.Data.outLFNDirBase = '/store/group/phys_exotica/dijet/Dijet13TeV/juska/Spring15/'
config.Data.outLFNDirBase = '/store/user/santanas/rootTrees/Run2015B_JetHT_25July2015-CertJson-251586-251883_JEC-Summer15_50nsV2_501dfb2/'
config.section_('User')
config.section_('Site')
#config.Site.storageSite = 'T2_CH_CERN'
config.Site.storageSite = 'T2_IT_Rome'
