import FWCore.ParameterSet.Config as cms 

# User options
nEvents  = 5
nReport  = 1
runJetTB = False
saveEDM  = False

process = cms.Process('XJETS')

process.load('PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')

## ----------------- Global Tag ------------------
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
#process.GlobalTag.globaltag = '80X_dataRun2_Prompt_v8'
#process.GlobalTag.globaltag = THISGLOBALTAG
#process.GlobalTag.globaltag = '80X_dataRun2_2016LegacyRepro_v4'
process.GlobalTag.globaltag = '80X_dataRun2_2016SeptRepro_v7' # recommended for Run2016 03Feb2017

#--------------------- Report and output ---------------------------
# Note: in grid runs this parameter is not used.
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(nEvents))

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = nReport


process.TFileService=cms.Service("TFileService",
                                 fileName=cms.string('tree_data.root'),
                                 #fileName=cms.string(THISROOTFILE),
                                 closeFileFast = cms.untracked.bool(True)
                                 )

## --- suppress long output ---> wantSummary = cms.untracked.bool(False) 

process.options = cms.untracked.PSet(
        allowUnscheduled = cms.untracked.bool(True),
        wantSummary = cms.untracked.bool(False),
)

############## output  edm format ###############
process.out = cms.OutputModule(
    'PoolOutputModule',
    fileName = cms.untracked.string('edm_data.root'),
    outputCommands = cms.untracked.vstring([
            #'keep *'
            'drop *',
            'keep *_*AK8*_*_*',
            #'keep *_*ak8PFJetsCHS*_*_*'
            #'keep *_slimmedJets_*_*',
            #'keep *_slimmedJetsAK8_*_*', 
            ])                                                                                           
    )

if saveEDM:
    process.output = cms.EndPath(process.out)

# Added 'vertexRef().isNonnull() &&' check for 80X data compatibility. Juska
process.chs = cms.EDFilter('CandPtrSelector', src = cms.InputTag('packedPFCandidates'), cut = cms.string('vertexRef().isNonnull() && fromPV'))

from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
process.slimmedGenJetsAK8 = ak4GenJets.clone(src = 'packedGenParticles', rParam = 0.8)


#-------------------------------------------------------
# Gen Particles Pruner
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")


process.prunedGenParticlesDijet = cms.EDProducer('GenParticlePruner',
    src = cms.InputTag("prunedGenParticles"),
    select = cms.vstring(
    "drop  *  ", # by default
    "keep ( status = 3 || (status>=21 && status<=29) )", # keep hard process particles
    )
)



#------------- Recluster Gen Jets to access the constituents -------
#already in toolbox, just add keep statements

process.out.outputCommands.append("keep *_slimmedGenJets_*_*")
process.out.outputCommands.append("keep *_slimmedGenJetsAK8_*_*")

##-------------------- Define the source  ----------------------------

#/JetHT/Run2016B-07Aug17_ver1-v1/MINIAOD
#/JetHT/Run2016E-07Aug17-v1/MINIAOD
#/JetHT/Run2016H-07Aug17-v1/MINIAOD
#/JetHT/Run2016D-07Aug17-v1/MINIAOD
#/JetHT/Run2016F-07Aug17-v1/MINIAOD
#/JetHT/Run2016B-07Aug17_ver2-v1/MINIAOD
#/JetHT/Run2016C-07Aug17-v1/MINIAOD
#/JetHT/Run2016G-07Aug17-v1/MINIAOD

# Note: for grid running it does not matter what's here, as input data is
# handled separately there.

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
        #"/store/data/Run2016H/JetHT/MINIAOD/07Aug17-v1/110000/F6440A09-927C-E711-97F2-0242AC110002.root"
        '/store/data/Run2016G/JetHT/MINIAOD/03Feb2017-v1/100000/006E7AF2-AEEC-E611-A88D-7845C4FC3B00.root'
        )
    )

##--- Use JetToolbox ---##
from JMEAnalysis.JetToolbox.jetToolbox_cff import jetToolbox
if runJetTB:
    jetToolbox( process, 'ak8', 'jetSequence', 'out', PUMethod='CHS', miniAOD=True, 
                #Cut='pt > 200 && abs(eta) < 2.5', # Tight
                addPruning = True, addSoftDrop = True, 
                addNsub    = True, 
                #addPrunedSubjets=True, addSoftDropSubjets=True,  
                #addNsubSubjets  =True 
                )
# Note: addNsubSubjets needs Subjets collection; in case both Pruned and SoftDrop are there, will pick by default SoftDrop
##########

##-------------------- User analyzer  --------------------------------

# Residue from AOD and RECO running
calo_collection=''
cluster_collection=''
pfcalo_collection=''
   

process.dijets     = cms.EDAnalyzer('DijetTreeProducer',

  # There's no avoiding this in Consumes era
  isData          = cms.bool(True),
  useJetTB         = cms.bool(False),
  jetsAK8_TB       = cms.InputTag('selectedPatJetsAK8PFCHS'),

  ## JETS/MET ########################################
  jetsAK4             = cms.InputTag('slimmedJets'), 
  jetsAK8             = cms.InputTag('slimmedJetsAK8'),     
  rho              = cms.InputTag('fixedGridRhoFastjetAll'),
  met              = cms.InputTag('slimmedMETs'),
  vtx              = cms.InputTag('offlineSlimmedPrimaryVertices'),
  ptMinAK4         = cms.double(10),
  ptMinAK8         = cms.double(10),
  
  ## MC ########################################
  pu               = cms.untracked.InputTag('slimmedAddPileupInfo'), # Updated from untracked to 80X by Juska
  ptHat            = cms.untracked.InputTag('generator'), # Why do these need to be 'untracked' anyway?
  genParticles     = cms.InputTag('prunedGenParticlesDijet'),
  genJetsAK4             = cms.InputTag('slimmedGenJets'), 
  genJetsAK8             = cms.InputTag('slimmedGenJetsAK8'),  

  ## trigger ###################################
  triggerAlias = cms.vstring('PFHT900','PFHT800','PFHT650','PFHT600','PFHT475','PFHT400','PFHT350','PFHT300','PFHT250','PFHT200',
                             'PFHT650MJJ950','PFHT650MJJ900',
                             'PFJET500','PFJET450','PFJET200',
                             'HT2000','HT2500','Mu45Eta2p1',
                             'AK8DiPFJet280200TrimMass30Btag','AK8PFHT600TriMass50Btag','AK8PFHT700TriMass50','AK8PFJet360TrimMass50',
                             'AK8PFJet450', 'AK8PFJet500','CaloJet500NoJetID','DiPFJetAve300HFJEC','DiPFJetAve500',
                             'PFHT400SixJet30Btag','PFHT450SixJet40Btag','PFHT750FourJetPt50','QuadPFJetVBF'
                             ),

  triggerSelection = cms.vstring(
     'HLT_PFHT900_v*',
     'HLT_PFHT800_v*',
     'HLT_PFHT650_v*',
     'HLT_PFHT600_v*',
     'HLT_PFHT475_v*',
     'HLT_PFHT400_v*',
     'HLT_PFHT350_v*',
     'HLT_PFHT300_v*',
     'HLT_PFHT250_v*',
     'HLT_PFHT200_v*',
     'HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_v*',
     'HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v*',
     'HLT_PFJet500_v*',
     'HLT_PFJet450_v*',
     'HLT_PFJet200_v*',
     'HLT_HT2000_v*',
     'HLT_HT2500_v*',
     'HLT_Mu45_eta2p1_v*',
     'HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV0p45_v*',
     'HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV0p45_v*',
     'HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_v*',
     'HLT_AK8PFJet360_TrimMass30_v*',
     'HLT_AK8PFJet450_v*',
     'HLT_AK8PFJet500_v*',
     'HLT_CaloJet500_NoJetID_v*',
     'HLT_DiPFJetAve300_HFJEC_v*',
     'HLT_DiPFJetAve500_v*',
     'HLT_PFHT400_SixJet30_BTagCSV0p55_2PFBTagCSV0p72_v*',
     'HLT_PFHT450_SixJet40_PFBTagCSV0p72_v*',
     'HLT_PFHT750_4JetPt50_v*',
     'HLT_QuadPFJet_VBF_v*',
     ###
  ),

  triggerConfiguration = cms.PSet(
    hltResults            = cms.InputTag('TriggerResults','','HLT'),
    l1tResults            = cms.InputTag(''),
    daqPartitions         = cms.uint32(1),
    l1tIgnoreMaskAndPrescale = cms.bool(False),
    #l1tIgnoreMask         = cms.bool(False),
   # l1techIgnorePrescales = cms.bool(False),
    throw                 = cms.bool(False)
  ),

  ## Noise Filters ###################################
  noiseFilterSelection_HBHENoiseFilter = cms.string('Flag_HBHENoiseFilter'),
  noiseFilterSelection_CSCTightHaloFilter = cms.string('Flag_CSCTightHaloFilter'),
  noiseFilterSelection_hcalLaserEventFilter = cms.string('Flag_hcalLaserEventFilter'),
  noiseFilterSelection_EcalDeadCellTriggerPrimitiveFilter = cms.string('Flag_EcalDeadCellTriggerPrimitiveFilter'),
  noiseFilterSelection_goodVertices = cms.string('Flag_goodVertices'),
  #noiseFilterSelection_trackingFailureFilter = cms.string('Flag_trackingFailureFilter'),
  noiseFilterSelection_eeBadScFilter = cms.string('Flag_eeBadScFilter'),
  noiseFilterSelection_ecalLaserCorrFilter = cms.string('Flag_ecalLaserCorrFilter'),
  noiseFilterSelection_trkPOGFilters = cms.string('Flag_trkPOGFilters'),
  # and the sub-filters
  noiseFilterSelection_trkPOG_manystripclus53X = cms.string('Flag_trkPOG_manystripclus53X'),
  noiseFilterSelection_trkPOG_toomanystripclus53X = cms.string('Flag_trkPOG_toomanystripclus53X'),
  noiseFilterSelection_trkPOG_logErrorTooManyClusters = cms.string('Flag_trkPOG_logErrorTooManyClusters'),

  noiseFilterConfiguration = cms.PSet(
    hltResults            = cms.InputTag('TriggerResults','','RECO'), #for prompt reco
    l1tResults            = cms.InputTag(''),
    daqPartitions         = cms.uint32(1),
    l1tIgnoreMaskAndPrescale = cms.bool(False),
    #l1tIgnoreMask         = cms.bool(False),
    #l1techIgnorePrescales = cms.bool(False),
    throw                 = cms.bool(False)
  ),


  ## JECs ################
  redoJECs  = cms.bool(True),

  ## Version Summer15_25nsV3 ( https://hypernews.cern.ch/HyperNews/CMS/get/JetMET/ )
  # Note that it hardly matters what is put in here, as these should be overriden in analysis step anyway. Juska.
  # That's also why these JEC's are greatly dated.
  L1corrAK4_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer15_25nsV3_DATA/Summer15_25nsV3_DATA_L1FastJet_AK4PFchs.txt'),
  L2corrAK4_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer15_25nsV3_DATA/Summer15_25nsV3_DATA_L2Relative_AK4PFchs.txt'),
  L3corrAK4_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer15_25nsV3_DATA/Summer15_25nsV3_DATA_L3Absolute_AK4PFchs.txt'),
  ResCorrAK4_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer15_25nsV3_DATA/Summer15_25nsV3_DATA_L2L3Residual_AK4PFchs.txt'),
  L1corrAK8_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer15_25nsV3_DATA/Summer15_25nsV3_DATA_L1FastJet_AK8PFchs.txt'),
  L2corrAK8_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer15_25nsV3_DATA/Summer15_25nsV3_DATA_L2Relative_AK8PFchs.txt'),
  L3corrAK8_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer15_25nsV3_DATA/Summer15_25nsV3_DATA_L3Absolute_AK8PFchs.txt'),
  ResCorrAK8_DATA = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer15_25nsV3_DATA/Summer15_25nsV3_DATA_L2L3Residual_AK4PFchs.txt'),
  L1corrAK4_MC = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer15_25nsV3_MC/Summer15_25nsV3_MC_L1FastJet_AK4PFchs.txt'),
  L2corrAK4_MC = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer15_25nsV3_MC/Summer15_25nsV3_MC_L2Relative_AK4PFchs.txt'),
  L3corrAK4_MC = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer15_25nsV3_MC/Summer15_25nsV3_MC_L3Absolute_AK4PFchs.txt'),
  L1corrAK8_MC = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer15_25nsV3_MC/Summer15_25nsV3_MC_L1FastJet_AK8PFchs.txt'),
  L2corrAK8_MC = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer15_25nsV3_MC/Summer15_25nsV3_MC_L2Relative_AK8PFchs.txt'),
  L3corrAK8_MC = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/Summer15_25nsV3_MC/Summer15_25nsV3_MC_L3Absolute_AK8PFchs.txt')
)

if runJetTB:                                    
    process.dijets.useJetTB   = cms.bool(True)
    process.dijets.jetsAK8_TB = cms.InputTag('selectedPatJetsAK8PFCHS')

# ------------------ path --------------------------


process.p = cms.Path()
process.p +=                      process.chs
process.p +=                      process.dijets
