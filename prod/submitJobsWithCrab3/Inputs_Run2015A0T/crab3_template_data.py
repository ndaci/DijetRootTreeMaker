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
config.Data.unitsPerJob = 20
config.Data.splitting = 'LumiBased'
config.Data.publication = False
config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/Cert_246908-248005_13TeV_PromptReco_Collisions15_ZeroTesla_JSON_CaloOnly.txt'
#config.Data.outLFNDirBase = '/store/group/phys_exotica/dijet/Dijet13TeV/juska/Spring15/'
config.Data.outLFNDirBase = '/store/user/santanas/rootTrees/Run2015A0T_Jet_246908-248005_V1_Tagf0be56e/'
config.section_('User')
config.section_('Site')
#config.Site.storageSite = 'T2_CH_CERN'
config.Site.storageSite = 'T2_IT_Rome'
