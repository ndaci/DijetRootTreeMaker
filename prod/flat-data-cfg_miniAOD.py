import FWCore.ParameterSet.Config as cms 

process = cms.Process('jetToolbox')

process.load('PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')

############################################################################
###### Noise Filters -- load here and apply in process path or before ######
############################################################################
#https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFilters

#process.load('RecoMET.METFilters.CSCTightHaloFilter_cfi') #DOES NOT WORK
# Error message:
#Principal::getByToken: Found zero products matching all criteria
#Looking for type: reco::BeamHaloSummary
#Looking for module label: BeamHaloSummary

#process.load('RecoMET.METFilters.eeBadScFilter_cfi') #DOES NOT WORK
# Error message:
#   [2] Calling event method for module EEBadScFilter/'eeBadScFilter'
#Exception Message:
#Principal::getByToken: Found zero products matching all criteria
#Looking for type: edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >
#Looking for module label: reducedEcalRecHitsEE

#process.load('RecoMET.METFilters.trackingFailureFilter_cfi') #DOES NOT WORK
#process.goodVertices = cms.EDFilter(
#  "VertexSelector",
#  filter = cms.bool(False),
#  src = cms.InputTag("offlinePrimaryVertices"),
#  cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.rho < 2")
#)
# Error message:
#   [2] Calling event method for module VertexSelector/'goodVertices'
#Exception Message:
#Principal::getByToken: Found zero products matching all criteria
#Looking for type: std::vector<reco::Vertex>
#Looking for module label: offlinePrimaryVertice

process.load('RecoMET.METFilters.hcalLaserEventFilter_cfi')
process.load('RecoMET.METFilters.EcalDeadCellTriggerPrimitiveFilter_cfi')
process.load('CommonTools/RecoAlgos/HBHENoiseFilter_cfi')

###################################### Run on AOD instead of MiniAOD? ########
runOnAOD=True
###################################### Run on RECO instead of MiniAOD? ########
runOnRECO=True

## ----------------- Global Tag ------------------
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

#process.GlobalTag.globaltag = 'START53_V27::All'
#process.GlobalTag.globaltag = 'START70_V6::All'
#process.GlobalTag.globaltag = 'POSTLS170_V5::All'
#process.GlobalTag.globaltag = 'POSTLS170_V7::All'
#process.GlobalTag.globaltag = 'PLS170_V7AN1::All'
#process.GlobalTag.globaltag = 'PHYS14_25_V1::All'
#process.GlobalTag.globaltag = 'MCRUN2_74_V9A::All'
process.GlobalTag.globaltag = 'GR_P_V56::All'


#--------------------- Report and output ---------------------------

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 5000


process.TFileService=cms.Service("TFileService",
                                 #fileName=cms.string('dijetTree_signal_M1000.root'),
                                 #fileName=cms.string('dijetTree_signal_M8000.root'),
                                 #fileName=cms.string('dijetTree_QstarToJJ_M_3000_PHYS14.root'),
                                 #fileName=cms.string('dijetTree_data_RunQCD.root'),
                                 fileName=cms.string('dijetTree_data_Run245155.root'),
                                 #fileName=cms.string('dijetTree_data_Run246908.root'),
                                 #fileName=cms.string('dijetTree_data_Run247081.root'),
                                 closeFileFast = cms.untracked.bool(True)
                                 )

## --- suppress long output ---> wantSummary = cms.untracked.bool(False) 

process.options = cms.untracked.PSet(
        allowUnscheduled = cms.untracked.bool(True),
        wantSummary = cms.untracked.bool(False),
)

############## output  edm format ###############
process.out = cms.OutputModule('PoolOutputModule',                                                                                                                  
                               fileName = cms.untracked.string('jettoolbox.root'),                                                                              
                               outputCommands = cms.untracked.vstring([
                                                                       # 'keep *_ak4PFJetsCHS_*_*',                                                                    
                                                                       # 'keep *_patJetsAK4PFCHS_*_*',                                                                  
                                                                       # 'keep *_ca8PFJetsCHS_*_*',                                                                     
                                                                       # 'keep *_patJetsCA8PFCHS_*_*',                                                                  
                                                                       # 'keep *_ak8PFJetsCHS_*_*',                                                                     
                                                                       # 'keep *_patJetsAK8PFCHS_*_*',                                                                  
                                                                      'keep *_slimmedJets_*_*',                                                                  
                                                                      'keep *_slimmedJetsAK8_*_*',                                                                  
                                                                       ])                                                                                           
                               )

#### NOT RUNNING OUTPUT MODULE ######                                                                                                                              
# process.endpath = cms.EndPath(process.out)    

### RUN MINIAOD SEQUENCE
if runOnAOD:
  from FWCore.ParameterSet.Utilities import convertToUnscheduled
  process=convertToUnscheduled(process)
  process.load('Configuration.StandardSequences.PAT_cff')
  from PhysicsTools.PatAlgos.slimming.miniAOD_tools import miniAOD_customizeAllData 
  process = miniAOD_customizeAllData(process)

if runOnRECO:
### RUN PFCLUSTERJETS
  process.load("RecoParticleFlow.PFClusterProducer.particleFlowCluster_cff")
  process.load("RecoLocalCalo.HcalRecAlgos.hcalRecAlgoESProd_cfi")
  process.load("RecoLocalCalo.EcalRecAlgos.EcalSeverityLevelESProducer_cfi")
  process.pfClusterRefsForJetsHCAL = cms.EDProducer("PFClusterRefCandidateProducer",
    src          = cms.InputTag('particleFlowClusterHCAL'),
    particleType = cms.string('pi+')
  )
  process.pfClusterRefsForJetsECAL = cms.EDProducer("PFClusterRefCandidateProducer",
    src          = cms.InputTag('particleFlowClusterECAL'),
    particleType = cms.string('pi+')
  )
  process.pfClusterRefsForJetsHF = cms.EDProducer("PFClusterRefCandidateProducer",
    src          = cms.InputTag('particleFlowClusterHF'),
    particleType = cms.string('pi+')
  )
  process.pfClusterRefsForJetsHO = cms.EDProducer("PFClusterRefCandidateProducer",
    src          = cms.InputTag('particleFlowClusterHO'),
    particleType = cms.string('pi+')
  )
  process.pfClusterRefsForJets = cms.EDProducer("PFClusterRefCandidateMerger",
    src = cms.VInputTag("pfClusterRefsForJetsHCAL", "pfClusterRefsForJetsECAL", "pfClusterRefsForJetsHF", "pfClusterRefsForJetsHO")
  )
  process.load("RecoJets.JetProducers.ak4PFClusterJets_cfi")
  process.pfClusterRefsForJets_step = cms.Sequence(
   process.particleFlowRecHitECAL*
   process.particleFlowRecHitHBHE*
   process.particleFlowRecHitHF*
   process.particleFlowRecHitHO*
   process.particleFlowClusterECALUncorrected*
   process.particleFlowClusterECAL*
   process.particleFlowClusterHBHE*
   process.particleFlowClusterHCAL*
   process.particleFlowClusterHF*
   process.particleFlowClusterHO*
   process.pfClusterRefsForJetsHCAL*
   process.pfClusterRefsForJetsECAL*
   process.pfClusterRefsForJetsHF*
   process.pfClusterRefsForJetsHO*
   process.pfClusterRefsForJets*
   process.ak4PFClusterJets
  )
### RUN PFCALOJETS
  ############ need the following setup when running in <CMSSW_7_5_X:
  # git cms-addpkg RecoParticleFlow/PFProducer
  # git cherry-pick af5c1ba33e88b3be627c262eb93d678f9f70e729
  process.hltParticleFlowBlock = cms.EDProducer("PFBlockProducer",
    debug = cms.untracked.bool(False),
    verbose = cms.untracked.bool(False),
    elementImporters = cms.VPSet(
        cms.PSet(
            source = cms.InputTag("particleFlowClusterECAL"),
            #source = cms.InputTag("particleFlowClusterECALUncorrected"), #we use uncorrected
            importerName = cms.string('GenericClusterImporter')
        ),
        cms.PSet(
            source = cms.InputTag("particleFlowClusterHCAL"),
            importerName = cms.string('GenericClusterImporter')
        ),
        cms.PSet(
            source = cms.InputTag("particleFlowClusterHO"),
            importerName = cms.string('GenericClusterImporter')
        ),
        cms.PSet(
            source = cms.InputTag("particleFlowClusterHF"),
            importerName = cms.string('GenericClusterImporter')
        )
    ),
    linkDefinitions = cms.VPSet(
        cms.PSet(
            linkType = cms.string('ECAL:HCAL'),
            useKDTree = cms.bool(False),
            #linkerName = cms.string('ECALAndHCALLinker')
            linkerName = cms.string('ECALAndHCALCaloJetLinker') #new ECal and HCal Linker for PFCaloJets
        ),
        cms.PSet(
            linkType = cms.string('HCAL:HO'),
            useKDTree = cms.bool(False),
            linkerName = cms.string('HCALAndHOLinker')
        ),
        cms.PSet(
            linkType = cms.string('HFEM:HFHAD'),
            useKDTree = cms.bool(False),
            linkerName = cms.string('HFEMAndHFHADLinker')
        ),
        cms.PSet(
            linkType = cms.string('ECAL:ECAL'),
            useKDTree = cms.bool(False),
            linkerName = cms.string('ECALAndECALLinker')
        )
     )
  )
  from RecoParticleFlow.PFProducer.particleFlow_cfi import particleFlowTmp
  process.hltParticleFlow = particleFlowTmp.clone(
    GedPhotonValueMap = cms.InputTag(""),
    useEGammaFilters = cms.bool(False),
    useEGammaElectrons = cms.bool(False), 
    useEGammaSupercluster = cms.bool(False),
    rejectTracks_Step45 = cms.bool(False),
    usePFNuclearInteractions = cms.bool(False),  
    blocks = cms.InputTag("hltParticleFlowBlock"), 
    egammaElectrons = cms.InputTag(""),
    useVerticesForNeutral = cms.bool(False),
    PFEGammaCandidates = cms.InputTag(""),
    useProtectionsForJetMET = cms.bool(False),
    usePFConversions = cms.bool(False),
    rejectTracks_Bad = cms.bool(False),
    muons = cms.InputTag(""),
    postMuonCleaning = cms.bool(False),
    usePFSCEleCalib = cms.bool(False)
    )
  from RecoJets.JetProducers.PFJetParameters_cfi import *
  process.PFCaloJetParameters = PFJetParameters.clone(
    src = cms.InputTag('hltParticleFlow')
   )
  from RecoJets.JetProducers.AnomalousCellParameters_cfi import *
  process.ak4PFCaloJets = cms.EDProducer(
    "FastjetJetProducer",
    process.PFCaloJetParameters,
    AnomalousCellParameters,
    jetAlgorithm = cms.string("AntiKt"),
    rParam       = cms.double(0.4)
    )
  process.pfClusterRefsForJets_step += process.hltParticleFlowBlock
  process.pfClusterRefsForJets_step += process.hltParticleFlow
  process.pfClusterRefsForJets_step += process.ak4PFCaloJets


# ----------------------- Jet Tool Box  -----------------
# ----- giulia test: do not recluster ak4 and ca8 jets to save time --------


# ##Load the toolBoxMiniHelper
# #from RecoJets.JetProducers.jettoolboxMiniHelper_cff import * 
# ##for some reason doesn't work  ?__?
# ## just copy&paste the cff here

process.chs = cms.EDFilter('CandPtrSelector', src = cms.InputTag('packedPFCandidates'), cut = cms.string('fromPV'))

from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
process.slimmedGenJetsAK8 = ak4GenJets.clone(src = 'packedGenParticles', rParam = 0.8)

# process.ak4PFJets.src = 'packedPFCandidates'
# process.ak4PFJets.doAreaFastjet = True

# process.ak4PFJetsCHS = process.ak4PFJetsCHS.clone(src = 'chs', doAreaFastjet = True) #boh ?
# process.ak8PFJetsCHS = process.ak4PFJetsCHS.clone(src = 'chs', doAreaFastjet = True, rParam = 0.8)
# process.ca8PFJetsCHS = process.ca4PFJets.clone(src = 'chs', doAreaFastjet = True, rParam = 0.8)

# process.ak4GenJets.src = 'packedGenParticles'
# process.ak8GenJets = process.ak4GenJets.clone(src = 'packedGenParticles', rParam = 0.8)
# process.ca8GenJets = process.ca4GenJets.clone(src = 'packedGenParticles', rParam = 0.8)

# process.fixedGridRhoFastjetAll.pfCandidatesTag = 'packedPFCandidates'

# process.ak4PFJetsCHSPruned = ak5PFJetsCHSPruned.clone(
#     src='chs'
#     rParam = 0.4, 
#     jetPtMin = 15.0
#     )
# process.ak4PFJetsCHSFiltered = ak5PFJetsCHSFiltered.clone(
#     src='chs'
#     rParam = 0.4, 
#     jetPtMin = 15.0 
#     )                                                                             
# process.ak4PFJetsCHSTrimmed = ak5PFJetsCHSTrimmed.clone(
#     src='chs'
#     rParam = 0.4,
#     jetPtMin = 15.0  
#     )

# process.ak8PFJetsCHSPruned.src = 'chs'
# process.ak8PFJetsCHSTrimmed.src = 'chs'
# process.ak8PFJetsCHSFiltered.src = 'chs'

# process.ca8PFJetsCHSPruned.src = 'chs'
# process.ca8PFJetsCHSTrimmed.src = 'chs'
# process.ca8PFJetsCHSFiltered.src = 'chs'

# process.cmsTopTagPFJetsCHS.src = 'chs'

# from PhysicsTools.PatAlgos.tools.jetTools import addJetCollection
# from PhysicsTools.PatAlgos.tools.jetTools import switchJetCollection

# addJetCollection(
#     process,
#     labelName = 'AK4PFCHS',
#     jetSource = cms.InputTag('ak4PFJetsCHS'),
#     algo = 'ak4',
#     rParam = 0.4,
#     jetCorrections = ('AK5PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
#     trackSource = cms.InputTag('unpackedTracksAndVertices'),
#     pvSource = cms.InputTag('unpackedTracksAndVertices'),
#     btagDiscriminators = ['combinedSecondaryVertexBJetTags'],
#     ) 

# addJetCollection(
#     process,
#     labelName = 'CA8PFCHS',
#     jetSource = cms.InputTag('ca8PFJetsCHS'),
#     algo = 'ca8',
#     rParam = 0.8,
#     jetCorrections = ('AK7PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
#     trackSource = cms.InputTag('unpackedTracksAndVertices'),
#     pvSource = cms.InputTag('unpackedTracksAndVertices'),
#     )                                                                                                                                                                   

# addJetCollection(
#     process,
#     labelName = 'AK8PFCHS',
#     jetSource = cms.InputTag('ak8PFJetsCHS'),
#     algo = 'ak8',
#     rParam = 0.8,
#     jetCorrections = ('AK7PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
#     trackSource = cms.InputTag('unpackedTracksAndVertices'),
#     pvSource = cms.InputTag('unpackedTracksAndVertices'),
#     ) 

# """
# switchJetCollection(
#     process,
#     jetSource = cms.InputTag('ak4PFJets'),
#     algo = 'ak4',
#     rParam = 0.4,
#     jetCorrections = ('AK5PF', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'Type-1'),
#     # btagDiscriminators = ['jetBProbabilityBJetTags',
#     #                       'jetProbabilityBJetTags',
#     #                       'trackCountingHighPurBJetTags',
#     #                       'trackCountingHighEffBJetTags',
#     #                       'simpleSecondaryVertexHighEffBJetTags',
#     #                       'simpleSecondaryVertexHighPurBJetTags',
#     #                       'combinedSecondaryVertexBJetTags'
#     #                       ],
#     trackSource = cms.InputTag('unpackedTracksAndVertices'),
#     pvSource = cms.InputTag('unpackedTracksAndVertices'),
#     )
# """

# process.patJetsAK4PFCHS.addJetCharge   = False
# process.patJetsAK4PFCHS.addBTagInfo    = True
# process.patJetsAK4PFCHS.getJetMCFlavour = False
# process.patJetsAK4PFCHS.addAssociatedTracks = False
# process.patJetPartonMatchAK4PFCHS.matched='prunedGenParticles'
# process.patJetCorrFactorsAK4PFCHS.primaryVertices = 'offlineSlimmedPrimaryVertices'

# process.patJetsCA8PFCHS.addJetCharge   = False
# process.patJetsCA8PFCHS.addBTagInfo    = False   #For some reason this has to be False
# process.patJetsCA8PFCHS.getJetMCFlavour = False
# process.patJetsCA8PFCHS.addAssociatedTracks = False
# process.patJetPartonMatchCA8PFCHS.matched='prunedGenParticles'
# process.patJetCorrFactorsCA8PFCHS.primaryVertices = 'offlineSlimmedPrimaryVertices'

# process.patJetsAK8PFCHS.addJetCharge   = False
# process.patJetsAK8PFCHS.addBTagInfo    = False    #For some reason this has to be False
# process.patJetsAK8PFCHS.getJetMCFlavour = False
# process.patJetsAK8PFCHS.addAssociatedTracks = False
# process.patJetPartonMatchAK8PFCHS.matched='prunedGenParticles'
# process.patJetCorrFactorsAK8PFCHS.primaryVertices = 'offlineSlimmedPrimaryVertices'

# process.patJets.addJetCharge   = False
# process.patJets.addBTagInfo    = True
# process.patJets.getJetMCFlavour = False
# process.patJets.addAssociatedTracks = False
# process.patJetPartonMatch.matched = 'prunedGenParticles'
# process.patJetCorrFactors.primaryVertices = 'offlineSlimmedPrimaryVertices'


# process.load('RecoBTag.Configuration.RecoBTag_cff')
# process.load('RecoJets.Configuration.RecoJetAssociations_cff')
# process.load('PhysicsTools.PatAlgos.slimming.unpackedTracksAndVertices_cfi')

# process.ak4JetTracksAssociatorAtVertexPF.jets = cms.InputTag('ak4PFJetsCHS')
# process.ak4JetTracksAssociatorAtVertexPF.jets = cms.InputTag('ak4PFJetsCHS')
# process.ak4JetTracksAssociatorAtVertexPF.tracks = cms.InputTag('unpackedTracksAndVertices')
# process.ak8JetTracksAssociatorAtVertexPF=process.ak4JetTracksAssociatorAtVertexPF.clone(jets = cms.InputTag('ak8PFJetsCHS'),
#                                                                                         coneSize = 0.8)
# process.ca8JetTracksAssociatorAtVertexPF=process.ak4JetTracksAssociatorAtVertexPF.clone(jets = cms.InputTag('ca8PFJetsCHS'),
#                                                                                         coneSize = 0.8)

# process.impactParameterTagInfos.primaryVertex = cms.InputTag('unpackedTracksAndVertices')
# process.inclusiveSecondaryVertexFinderTagInfos.extSVCollection = cms.InputTag('unpackedTracksAndVertices','secondary','')
# process.combinedSecondaryVertex.trackMultiplicityMin = 1 #silly sv, uses un filtered tracks.. i.e. any pt

# process.load('FWCore.MessageLogger.MessageLogger_cfi')
# #process.MessageLogger.cerr.FwkReport.reportEvery = 10
# process.MessageLogger.suppressWarning = cms.untracked.vstring('ecalLaserCorrFilter','manystripclus53X','toomanystripclus53X')
# #process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
# #process.options.allowUnscheduled = cms.untracked.bool(True)


# #Load the toolbox

# process.load('CMSDIJET.DijetRootTreeMaker.jettoolbox_cff')


#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#Test with Edmund Code... not working
#Load the GenEvent and GenParticles cff files
#process.load('CMSDIJET.DijetRootTreeMaker.RootTupleMakerV2_GenEventInfo_cfi')
#process.load('CMSDIJET.DijetRootTreeMaker.RootTupleMakerV2_GenParticles_cfi')

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

#process.out.outputCommands.append("keep *_ak4GenJets_*_*")
#process.out.outputCommands.append("keep *_ak8GenJets_*_*")
#process.out.outputCommands.append("keep *_ca8GenJets_*_*")

process.out.outputCommands.append("keep *_slimmedGenJets_*_*")
process.out.outputCommands.append("keep *_slimmedGenJetsAK8_*_*")

##-------------------- Define the source  ----------------------------



process.source = cms.Source("PoolSource",
    #fileNames = cms.untracked.vstring('file:miniAOD_RSGravToJJ_kMpl01_M-8000.root')
    #fileNames = cms.untracked.vstring('file:/cmshome/santanas/CMS/data/Spring14miniaod__RSGravToJJ_kMpl01_M-1000_Tune4C_13TeV-pythia8__MINIAODSIM__PU20bx25_POSTLS170_V5-v1__00000__6AACD832-3707-E411-A167-001E672489D5.root')
    #fileNames = cms.untracked.vstring('file:/cmshome/santanas/CMS/data/Spring14drAODSIM__RSGravToJJ_kMpl01_M-1000_Tune4C_13TeV-pythia8__AODSIM__PU20bx25_POSTLS170_V5-v1__00000__0622C950-58E4-E311-A595-0025904B130A.root')
    #fileNames = cms.untracked.vstring('file:2CEB70D6-D918-E411-B814-003048F30422.root')    
    #fileNames = cms.untracked.vstring('file:QstarToJJ_M_4000_TuneCUETP8M1_13TeV_pythia8__MINIAODSIM__Asympt50ns_MCRUN2_74_V9A-v1__70000__AA35D1E7-FEFE-E411-B1C5-0025905B858A.root')    
    #fileNames = cms.untracked.vstring('/store/data/Run2015A/Jet/AOD/PromptReco-v1/000/247/081/00000/804F6C9F-DB0C-E511-B0B6-02163E0143D9.root')
    #fileNames = cms.untracked.vstring('/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/081/00000/EA9C98A0-E20C-E511-B59A-02163E011B58.root')
fileNames = cms.untracked.vstring(
'/store/data/Commissioning2015/MinimumBias/RECO/PromptReco-v1/000/245/155/00000/0A4A9D19-CC00-E511-B426-02163E01454B.root',
'/store/data/Commissioning2015/MinimumBias/RECO/PromptReco-v1/000/245/155/00000/18183AC8-CC00-E511-8313-02163E014338.root',
'/store/data/Commissioning2015/MinimumBias/RECO/PromptReco-v1/000/245/155/00000/34F0B9F6-D400-E511-A707-02163E01422C.root',
'/store/data/Commissioning2015/MinimumBias/RECO/PromptReco-v1/000/245/155/00000/44361EC3-CB00-E511-AEDD-02163E014100.root',
'/store/data/Commissioning2015/MinimumBias/RECO/PromptReco-v1/000/245/155/00000/74BF29D0-CC00-E511-873C-02163E014218.root',
'/store/data/Commissioning2015/MinimumBias/RECO/PromptReco-v1/000/245/155/00000/7C4EDB4E-CE00-E511-B1A9-02163E0135A8.root',
'/store/data/Commissioning2015/MinimumBias/RECO/PromptReco-v1/000/245/155/00000/84DFF492-D700-E511-9625-02163E0146A9.root',
'/store/data/Commissioning2015/MinimumBias/RECO/PromptReco-v1/000/245/155/00000/908B7240-CC00-E511-9590-02163E014682.root',
'/store/data/Commissioning2015/MinimumBias/RECO/PromptReco-v1/000/245/155/00000/B8BCD533-CC00-E511-AAFA-02163E011DDC.root',
'/store/data/Commissioning2015/MinimumBias/RECO/PromptReco-v1/000/245/155/00000/B8C6E2E8-CD00-E511-A8CB-02163E01355F.root',
'/store/data/Commissioning2015/MinimumBias/RECO/PromptReco-v1/000/245/155/00000/BAE71A7C-CB00-E511-B250-02163E0142DC.root',
'/store/data/Commissioning2015/MinimumBias/RECO/PromptReco-v1/000/245/155/00000/C66D433E-CC00-E511-A025-02163E013976.root',
'/store/data/Commissioning2015/MinimumBias/RECO/PromptReco-v1/000/245/155/00000/D4B470AC-CC00-E511-9548-02163E014719.root',
'/store/data/Commissioning2015/MinimumBias/RECO/PromptReco-v1/000/245/155/00000/F6608E71-D000-E511-93A0-02163E011D60.root',
)
)
fileNames = cms.untracked.vstring(
'/store/relval/CMSSW_7_4_4/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/MCRUN2_740TV1_0Tv2-v1/00000/004451DD-E907-E511-8379-0025905AA9CC.root',
'/store/relval/CMSSW_7_4_4/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/MCRUN2_740TV1_0Tv2-v1/00000/1224482F-DC07-E511-88FD-0025905A605E.root',
'/store/relval/CMSSW_7_4_4/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/MCRUN2_740TV1_0Tv2-v1/00000/2ADC4CF1-DB07-E511-B9CD-003048FFD760.root',
'/store/relval/CMSSW_7_4_4/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/MCRUN2_740TV1_0Tv2-v1/00000/346D28AC-DC07-E511-B413-0025905A497A.root',
'/store/relval/CMSSW_7_4_4/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/MCRUN2_740TV1_0Tv2-v1/00000/5AEBFB42-E707-E511-AFDC-0025905A6056.root',
'/store/relval/CMSSW_7_4_4/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/MCRUN2_740TV1_0Tv2-v1/00000/5CE63EA4-DC07-E511-8F07-0025905A6110.root',
'/store/relval/CMSSW_7_4_4/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/MCRUN2_740TV1_0Tv2-v1/00000/6A9EE3C3-ED07-E511-B4C4-0025905A6118.root',
'/store/relval/CMSSW_7_4_4/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/MCRUN2_740TV1_0Tv2-v1/00000/76FBE3E3-DC07-E511-B158-0025905A609A.root',
'/store/relval/CMSSW_7_4_4/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/MCRUN2_740TV1_0Tv2-v1/00000/DCD19734-DD07-E511-A156-0025905A60CE.root',
'/store/relval/CMSSW_7_4_4/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/MCRUN2_740TV1_0Tv2-v1/00000/F6156231-EA07-E511-A1D0-00259059642A.root',
'/store/relval/CMSSW_7_4_4/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/MCRUN2_740TV1_0Tv2-v1/00000/FCFD84E8-DB07-E511-A6C1-003048FFD79C.root',
)
fileNames = cms.untracked.vstring(
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/04C73F4E-320B-E511-9820-02163E01374A.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/04C9C349-2E0B-E511-B2D7-02163E01184E.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/08BE08AF-490B-E511-8357-02163E01410A.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/0C71723A-390B-E511-9EFD-02163E01365C.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/16B05D7A-380B-E511-A925-02163E01446F.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/1816AE9A-1E0B-E511-B799-02163E01410A.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/18402E98-3A0B-E511-A447-02163E013837.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/1AA3B0D1-140B-E511-B89F-02163E011D41.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/1C988818-3B0B-E511-B180-02163E01374A.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/1E7689B5-350B-E511-9302-02163E0141D2.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/1E7F85B7-2E0B-E511-B580-02163E014146.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/1EBE0646-300B-E511-B016-02163E01365F.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/205EFAEA-360B-E511-888F-02163E013837.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/221CEB34-320B-E511-AF60-02163E0136A3.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/22B6EE97-360B-E511-BA44-02163E0121C5.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/24936509-A10B-E511-A095-02163E01432F.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/24CB9CB3-320B-E511-9691-02163E0141D2.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/284EC9D9-2B0B-E511-A74E-02163E01186C.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/2A5E6674-4B0B-E511-AD22-02163E0145D2.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/2C1F0FE6-540B-E511-B661-02163E0143EB.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/2C6E5651-350B-E511-8C76-02163E0138D3.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/2C8AC5B4-3C0B-E511-B783-02163E01429D.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/2E5DBA3B-320B-E511-ADBD-02163E011D25.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/30833ECA-2E0B-E511-912A-02163E012BE0.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/329066C0-2E0B-E511-85E2-02163E013613.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/32BFD97C-2F0B-E511-894F-02163E014110.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/32DA01BB-2B0B-E511-8425-02163E0142E6.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/344D6112-2D0B-E511-9C0E-02163E014349.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/3646A6A8-380B-E511-BCF9-02163E013614.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/3800665D-210B-E511-B8A7-02163E013826.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/3C47EB36-590B-E511-8699-02163E0141B9.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/3C6948BA-330B-E511-8190-02163E01288E.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/3C97665B-370B-E511-9297-02163E011ACE.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/3E233BDA-2D0B-E511-B36B-02163E014613.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/3E5B07C3-450B-E511-9D74-02163E0142BC.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/40712BDC-180B-E511-AC6D-02163E01288E.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/422DBA7B-400B-E511-8EB9-02163E0136B3.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/427E56CF-180B-E511-BB45-02163E0145F4.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/42AD048E-3B0B-E511-A10D-02163E0121DA.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/42DC1231-1A0B-E511-9180-02163E01257B.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/46DD66DC-380B-E511-8691-02163E0136C3.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/481E76B7-2E0B-E511-BE82-02163E0133DC.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/54FC7C61-340B-E511-967F-02163E0142D7.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/5ADB5ED2-2E0B-E511-A88C-02163E013729.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/5C4EA8C8-2E0B-E511-B99F-02163E01365F.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/5EEC6712-300B-E511-99D7-02163E012551.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/6046FE45-320B-E511-B260-02163E012AA9.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/66DF374D-1B0B-E511-B1A5-02163E013860.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/6851A60C-4C0B-E511-8B8E-02163E01383E.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/68992D28-480B-E511-995B-02163E012AA9.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/6AFC8165-2E0B-E511-B97D-02163E011D25.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/721F2A67-2C0B-E511-90C6-02163E012808.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/76460718-320B-E511-B2EE-02163E014338.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/7E18B4A2-2F0B-E511-AC3B-02163E014613.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/7E3EE5D8-420B-E511-9F9E-02163E01180A.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/84A0ADEC-2B0B-E511-84DD-02163E011DC2.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/84CB2C1D-330B-E511-A4FB-02163E01184E.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/86949AAE-4F0B-E511-A69D-02163E0136A7.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/887E40EC-340B-E511-9F9F-02163E013570.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/8C2D92DD-340B-E511-AA8D-02163E013999.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/8CF41181-450B-E511-A781-02163E012AA9.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/8E20D2E3-300B-E511-BBD9-02163E014113.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/94AD05E4-370B-E511-AC46-02163E011DA4.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/94ADF5F2-2E0B-E511-B5FC-02163E0136A3.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/98173C0C-2F0B-E511-BFB0-02163E013604.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/987356C3-380B-E511-8627-02163E01374A.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/989275B5-3A0B-E511-97C7-02163E011A9B.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/9E366DF8-3C0B-E511-A2E3-02163E01475D.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/A0FBC451-330B-E511-B97C-02163E0143B6.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/A2970551-320B-E511-8BD1-02163E011D69.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/A41098EB-340B-E511-917C-02163E014503.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/AC370FAF-3E0B-E511-9957-02163E014565.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/AC4AAA2D-300B-E511-81C3-02163E012462.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/ACB2A196-350B-E511-9CF0-02163E012032.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/B0DB22DE-300B-E511-B743-02163E0136A3.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/B21D56FC-4A0B-E511-832E-02163E012AA9.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/B224107E-440B-E511-92BA-02163E0138EB.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/B890FB49-430B-E511-AE45-02163E012925.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/BCDF57CD-420B-E511-9D59-02163E0135CA.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/C6E518D7-350B-E511-B7BE-02163E01180A.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/C6F94C87-380B-E511-855A-02163E014166.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/C875DB32-300B-E511-B3E0-02163E0141CC.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/CAC98058-320B-E511-B55F-02163E011B55.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/CC633A46-350B-E511-8F9E-02163E0141BD.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/CE644C35-4A0B-E511-80EF-02163E012AA9.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/D2A5C31E-680B-E511-BC35-02163E011D69.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/D407BDCF-180B-E511-8EF6-02163E0141CC.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/D40E8A4A-2E0B-E511-A3AE-02163E011926.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/D4A8F746-340B-E511-8B81-02163E0119DC.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/D62DB0AF-4A0B-E511-B99D-02163E01220A.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/D6FE50D5-300B-E511-B787-02163E0136A3.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/DA9EE188-2F0B-E511-AA54-02163E01365C.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/DAE7EA2B-B70B-E511-BBD3-02163E0143EC.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/DC6C4A3F-630B-E511-81C5-02163E01220A.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/DE876ACA-180B-E511-9813-02163E01475C.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/E0937738-300B-E511-8260-02163E014565.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/E0CE74BD-310B-E511-9655-02163E013729.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/E2CC1316-300B-E511-94AE-02163E014374.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/E40307F1-320B-E511-8758-02163E013604.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/E4CF63D4-2B0B-E511-8D9E-02163E01475D.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/E6ACA1A9-2F0B-E511-87AE-02163E0121DA.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/E6E30CB1-3F0B-E511-B9D3-02163E011ACE.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/E87F8442-340B-E511-A612-02163E01383E.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/ECA05354-330B-E511-A02A-02163E0142E6.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/EE165AE1-630B-E511-B134-02163E011D69.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/F04EE6B0-320B-E511-A3A6-02163E014613.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/F0B756D1-4C0B-E511-A18C-02163E012298.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/F0E92325-450B-E511-B9FE-02163E014565.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/F2D229D9-310B-E511-917C-02163E0134D9.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/F41E61AB-570B-E511-A5C0-02163E014793.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/FAD55223-300C-E511-8E06-02163E0138C6.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/FC1660E1-3B0B-E511-98A1-02163E013518.root',
'/store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/908/00000/FC6A78FB-3C0B-E511-B34F-02163E0143FB.root',
    )


# #Keep statements for valueMaps (link Reco::Jets to associated quantities)
# #You don't have to keep them to deswizzle below

# #pocess.out.outputCommands += ['keep *_pileupJetIdEvaluator_*_*',
# #                       'keep *_QGTagger_*_*']

# process.out.outputCommands += ['keep *_NjettinessCA8_*_*',
#                                'keep *_NjettinessAK8_*_*',
                               
# #                               'keep *_QJetsAdderCA8_*_*',
#                                'keep *_ca8PFJetsCHSPrunedLinks_*_*',
#                                'keep *_ca8PFJetsCHSTrimmedLinks_*_*',
#                                'keep *_ca8PFJetsCHSFilteredLinks_*_*',
#                                'keep *_cmsTopTagPFJetsCHSLinksCA8_*_*',
# #                               'keep *_QJetsAdderAK8_*_*',
#                                'keep *_ak8PFJetsCHSPrunedLinks_*_*',
#                                'keep *_ak8PFJetsCHSTrimmedLinks_*_*',
#                                'keep *_ak8PFJetsCHSFilteredLinks_*_*',
#                                'keep *_cmsTopTagPFJetsCHSLinksAK8_*_*',
#                                ]


# #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# #Deswizzle valueMaps and attach to PAT::Jets as userFloats

# #process.patJetsAK4PFCHS.userData.userFloats.src += ['pileupJetIdEvaluator:fullDiscriminant','QGTagger:qgLikelihood']
# #process.patJetsAK4PFCHS.userData.userInts.src   += ['pileupJetIdEvaluator:cutbasedId','pileupJetIdEvaluator:fullId']

# #run these modules before pat in the sequence

# process.patJetsCA8PFCHS.userData.userFloats.src += ['NjettinessCA8:tau1',
#                                                     'NjettinessCA8:tau2',
#                                                     'NjettinessCA8:tau3',
# #                                                    'QJetsAdderCA8:QjetsVolatility',
#                                                     'ca8PFJetsCHSPrunedLinks',
#                                                     'ca8PFJetsCHSTrimmedLinks',
#                                                     'ca8PFJetsCHSFilteredLinks',
#                                                     'cmsTopTagPFJetsCHSLinksCA8'
#                                                     ]



# process.patJetsAK8PFCHS.userData.userFloats.src += ['NjettinessAK8:tau1',
#                                                     'NjettinessAK8:tau2',
#                                                     'NjettinessAK8:tau3',
#  #                                                   'QJetsAdderAK8:QjetsVolatility',
#                                                     'ak8PFJetsCHSPrunedLinks',
#                                                     'ak8PFJetsCHSTrimmedLinks',
#                                                     'ak8PFJetsCHSFilteredLinks',
#                                                     'cmsTopTagPFJetsCHSLinksAK8'
#                                                     ]





##-------------------- User analyzer  --------------------------------


process.dijets     = cms.EDAnalyzer('DijetTreeProducer',
  ## JETS/MET ########################################
  # jetsAK4             = cms.InputTag('patJetsAK4PFCHS'), 
  # jetsAK8         = cms.InputTag('patJetsAK8PFCHS'),     
  # jetsCA8         = cms.InputTag('patJetsCA8PFCHS'),
  jetsAK4             = cms.InputTag('slimmedJets'), 
  jetsAK4Calo         = cms.InputTag('ak4CaloJets'), 
  jetsAK4PFCluster    = cms.InputTag('ak4PFClusterJets'), 
  jetsAK4PFCalo    = cms.InputTag('ak4PFCaloJets'), 
  jetsAK8             = cms.InputTag('slimmedJetsAK8'),     
  rho              = cms.InputTag('fixedGridRhoFastjetAll'),
  met              = cms.InputTag('slimmedMETs'),
  vtx              = cms.InputTag('offlineSlimmedPrimaryVertices'),
  ptMinAK4         = cms.double(10),
  ptMinAK8         = cms.double(10),
  #ptMinCA8         = cms.double(10),
  #mjjMin           = cms.double(700),
  #dEtaMax          = cms.double(1.3),

  #forse non serve perche' gia' aggiunti al pat::jet -giulia-
  #Cutbasedid        = cms.InputTag("pileupJetIdEvaluator:cutbasedId"),
  #fullDiscriminant  = cms.InputTag("pileupJetIdEvaluator:fullDiscriminant")
  #fullId            = cms.InputTag("pileupJetIdEvaluator:fullId"),
  ## MC ########################################
  pu               = cms.untracked.InputTag('addPileupInfo'),
  ptHat            = cms.untracked.InputTag('generator'),
  genParticles     = cms.InputTag('prunedGenParticlesDijet'),
  # genJetsAK4             = cms.InputTag('ak4GenJets'), 
  # genJetsAK8         = cms.InputTag('ak8GenJets'),     
  # genJetsCA8         = cms.InputTag('ca8GenJets'),
  genJetsAK4             = cms.InputTag('slimmedGenJets'), 
  genJetsAK8             = cms.InputTag('slimmedGenJetsAK8'),     


  ## trigger ###################################
  #triggerAlias     = cms.vstring('Fat','PFHT650','PFNoPUHT650','HT750','HT550'),
  triggerAlias     = cms.vstring('PFHT900'),
  triggerSelection = cms.vstring(
     #'HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v*',
     #'HLT_PFHT650_v*', #giulia : commented because not found in new entuples
     ### giulia
     #'HLT_HT650_v*',
     ### end giulia
     #'HLT_PFNoPUHT650_v*',
     #'HLT_HT750_v*',  
     #'HLT_HT550_v*'
     'HLT_PFHT900_v*'
  ),
  triggerConfiguration = cms.PSet(
    hltResults            = cms.InputTag('TriggerResults','','HLT'),
    l1tResults            = cms.InputTag(''),
    daqPartitions         = cms.uint32(1),
    l1tIgnoreMask         = cms.bool(False),
    l1techIgnorePrescales = cms.bool(False),
    throw                 = cms.bool(False)
  ),
  ## JECs ######################################
  redoJECs  = cms.bool(True),
  L1corrAK4 = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/PHYS14_25_V2_L1FastJet_AK4PFchs.txt'),
  L2corrAK4 = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/PHYS14_25_V2_L2Relative_AK4PFchs.txt'),
  L3corrAK4 = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/PHYS14_25_V2_L3Absolute_AK4PFchs.txt'),
  L1corrAK8 = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/PHYS14_25_V2_L1FastJet_AK8PFchs.txt'),
  L2corrAK8 = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/PHYS14_25_V2_L2Relative_AK8PFchs.txt'),
  L3corrAK8 = cms.FileInPath('CMSDIJET/DijetRootTreeMaker/data/PHYS14_25_V2_L3Absolute_AK8PFchs.txt')
)


# ------------------ path --------------------------
process.filter_step = cms.Path(
              process.EcalDeadCellTriggerPrimitiveFilter*
              process.hcalLaserEventFilter)

process.p = cms.Path()

if runOnRECO:
   process.p += process.pfClusterRefsForJets_step
	      
# Noise filters added first in path as recommended in Twiki

                     #process.CSCTightHaloFilter* does not work
                     #process.eeBadScFilter* does not work
                     #process.goodVertices*process.trackingFailureFilter* does not work
process.p +=                      process.HBHENoiseFilter
                     
                     
                     #process.prunedGenParticlesDijet* #GENPAR REMOVED
process.p +=                      process.chs

                     #process.slimmedGenJetsAK8 * #GENPAR REMOVED
                     
                     #process.ak4PFJetsCHS *
                     #process.ak4GenJets *
                     #process.patJetsAK4PFCHS *
                     ##process.patJetPartonMatchAK4PFCHS *
                     ##process.patJetCorrFactorsAK4PFCHS* 
                     #process.ak8PFJetsCHS *
                     #process.ak8GenJets 
                     #process.NjettinessAK8 *
                     ##process.QJetsAdderAK8 * #questo rallenta - da capire
                     # process.ak8PFJetsCHSPrunedLinks *
                     # process.ak8PFJetsCHSFilteredLinks *
                     # process.ak8PFJetsCHSTrimmedLinks *
                     # process.patJetsAK8PFCHS *
                     # process.ca8PFJetsCHS *
                     # process.ca8GenJets * 
                     # process.NjettinessCA8 *
                     # #process.QJetsAdderCA8 * #questo rallenta - da capire
                     # process.ca8PFJetsCHSPrunedLinks *
                     # process.ca8PFJetsCHSFilteredLinks *
                     # process.ca8PFJetsCHSTrimmedLinks *
                     # # process.patJetsCA8PFCHS *    

                     # #process.pileupJetId        ##recipe not working for now
                     # #process.puJetIdSequence    ##recipe not working for now
                     # #process.pileupJetIdCalculator *  ##recipe not working for now
                     # #process.pileupJetIdEvaluator     ##recipe not working for now
                     # #process.QGTagger *

process.p +=                      process.dijets
