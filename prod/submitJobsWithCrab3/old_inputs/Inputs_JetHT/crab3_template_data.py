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
config.Data.ignoreLocality = True 
#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/Cert_246908-248005_13TeV_PromptReco_Collisions15_ZeroTesla_JSON_CaloOnly.txt'
#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/DCSOnly/json_DCSONLY_Run2015B.txt'
#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/Cert_246908-255031_13TeV_PromptReco_Collisions15_50ns_JSON.txt'
#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/Cert_246908-255031_13TeV_PromptReco_Collisions15_25ns_JSON_v2.txt'
config.Data.lumiMask = 'inputJson_final.json' ## if you downloaded the file in the working directory
#config.Data.runRange = '251244-251585' # '193093-194075'                             
#config.Data.runRange = '251586-251883' # '193093-194075'                             
##
#config.Data.outLFNDirBase = '/store/group/phys_exotica/dijet/Dijet13TeV/santanas/rootTrees/Run2015B_SingleMuonPromptReco_08Sep2015-CertJson50ns-XXXXXX-YYYYYY_JEC-Summer15_50nsV4_0d7d556/'
config.Data.outLFNDirBase = '/store/group/phys_exotica/dijet/Dijet13TeV/santanas/rootTrees/Run2015D_JetHT_13Nov2015-DCSJson_bc77ee9/'
#config.Data.outLFNDirBase = '/store/user/santanas/rootTrees/Run2015B_SingleElectron_07Aug2015-CertJson-251586-251883_JEC-Summer15_50nsV2_501dfb2/'
#config.Data.outLFNDirBase = '/store/user/santanas/rootTrees/Run2015Dv4_JetHT_12Nov2015-DCSJson_2e3c9db/'
##
config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'
#config.Site.storageSite = 'T2_IT_Rome'
