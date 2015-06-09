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

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 5000


process.TFileService=cms.Service("TFileService",
                                 #fileName=cms.string('dijetTree_signal_M1000.root'),
                                 #fileName=cms.string('dijetTree_signal_M8000.root'),
                                 #fileName=cms.string('dijetTree_QstarToJJ_M_3000_PHYS14.root'),
                                 #fileName=cms.string('dijetTree_data_RunQCD.root'),
                                 #fileName=cms.string('dijetTree_data_Run245155.root'),
                                 fileName=cms.string('dijetTree_data_Run246908.root'),
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
fileNames = cms.untracked.vstring()
)
#process.source.fileNames.extend(
[
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/246/960/00000/4CED9FE9-DB0B-E511-A069-02163E012432.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/246/963/00000/104D5CB1-F20B-E511-9F6A-02163E0135C8.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/246/963/00000/1862B0AD-0C0C-E511-B29F-02163E0143B6.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/246/963/00000/308A9778-010C-E511-98BC-02163E012925.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/246/963/00000/40668AC8-040C-E511-9307-02163E0142F0.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/246/963/00000/407857E2-F70B-E511-BB79-02163E013604.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/246/963/00000/48EFAD97-080C-E511-B4BB-02163E011D11.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/246/963/00000/4A48510F-F70B-E511-AA2A-02163E01437C.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/246/963/00000/6289FF2E-FD0B-E511-A78A-02163E012BE0.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/246/963/00000/702F6994-FC0B-E511-9EFF-02163E01440B.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/246/963/00000/7283B16A-FB0B-E511-8D2F-02163E01474E.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/246/963/00000/7EB02C16-FB0B-E511-8218-02163E011C14.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/246/963/00000/864FD55A-E90B-E511-9470-02163E0138BC.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/246/963/00000/8E6D2472-EE0B-E511-8527-02163E0126FB.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/246/963/00000/923A4E48-F50B-E511-885F-02163E01447E.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/246/963/00000/9A01C119-FB0B-E511-A621-02163E01438B.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/246/963/00000/AAF148AD-F80B-E511-B6FF-02163E0134D5.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/246/963/00000/B62C208B-100C-E511-97D5-02163E014698.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/246/963/00000/B8090ADF-EE0B-E511-8528-02163E014196.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/246/963/00000/BA0F8C4D-F10B-E511-8977-02163E0140E8.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/246/963/00000/D053DCB5-F20B-E511-A1DB-02163E0118B8.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/246/963/00000/D25CDDB5-FF0B-E511-858E-02163E0140FA.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/246/963/00000/DE0A13C4-150C-E511-8A6F-02163E0141CE.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/246/963/00000/DE9101ED-EF0B-E511-A9E1-02163E011B6D.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/246/963/00000/E20D2F7F-E70B-E511-A7C2-02163E013399.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/246/963/00000/EE922BAE-110C-E511-98E3-02163E011BDB.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/246/963/00000/FA739B1C-0A0C-E511-9B6F-02163E013768.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/246/963/00000/FCD6BD64-E50B-E511-B5CB-02163E0138C4.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/246/965/00000/1CCEC492-E90B-E511-B57E-02163E0134DE.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/047/00000/5AD344AF-680C-E511-8193-02163E011B55.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/049/00000/F2D7A36A-6B0C-E511-B58C-02163E014565.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/052/00000/B08A4739-740C-E511-BF1C-02163E01392B.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/054/00000/F65CF239-790C-E511-AAA2-02163E011D4A.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/056/00000/44370C3B-7C0C-E511-985D-02163E013976.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/057/00000/4E394495-7E0C-E511-AB67-02163E01442E.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/059/00000/1E0BE554-810C-E511-BB8E-02163E01452F.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/060/00000/CAA465BF-820C-E511-B4CB-02163E012927.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/063/00000/E2E6296F-870C-E511-99C6-02163E012376.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/068/00000/0C271720-B60C-E511-A7EF-02163E0134D5.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/068/00000/1CBE3236-A10C-E511-9950-02163E0142E6.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/068/00000/1CD8C7FA-A00C-E511-A5FC-02163E01475E.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/068/00000/26039859-9C0C-E511-9FE4-02163E0136A3.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/068/00000/30241F14-A30C-E511-8147-02163E01287B.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/068/00000/5C1FB049-A10C-E511-9C4A-02163E0136A7.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/068/00000/6462F0EE-A70C-E511-9D15-02163E01456E.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/068/00000/70DDFE56-9B0C-E511-A071-02163E013671.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/068/00000/72C3221E-9B0C-E511-8676-02163E014212.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/068/00000/76B71834-A10C-E511-9DD1-02163E014565.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/068/00000/86F4731E-9B0C-E511-9A50-02163E014172.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/068/00000/8A5063CE-AA0C-E511-B2D1-02163E011C14.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/068/00000/8E2AF409-D30C-E511-8EAC-02163E013686.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/068/00000/9009EB87-A00C-E511-8AAD-02163E011E07.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/068/00000/94440BF9-9D0C-E511-9E79-02163E01453D.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/068/00000/94C2B6F5-9A0C-E511-8C13-02163E011C2F.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/068/00000/9E16D2C7-9A0C-E511-A8F6-02163E01257B.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/068/00000/A21344B9-D10C-E511-8F09-02163E0143FB.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/068/00000/C435D6B7-A70C-E511-8121-02163E011A9D.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/068/00000/CAAAD242-AA0C-E511-8C92-02163E0121DA.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/068/00000/CE2A43C5-A60C-E511-A4EC-02163E01374A.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/068/00000/E421D2FE-9E0C-E511-94B7-02163E0144F5.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/068/00000/E60B122B-A60C-E511-9FE3-02163E0139A3.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/068/00000/F83C5654-9B0C-E511-9DF5-02163E014220.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/068/00000/FAEEAAAE-AB0C-E511-AE15-02163E0138C6.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/068/00000/FC60BA2D-A90C-E511-B727-02163E0121DA.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/070/00000/026B36BF-AE0C-E511-896C-02163E01369F.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/070/00000/0CD2B7DC-AD0C-E511-9A4C-02163E0134D5.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/070/00000/0E3ECC12-AF0C-E511-9723-02163E011E07.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/070/00000/14700840-AE0C-E511-AF98-02163E01184E.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/070/00000/18C8DC1A-B20C-E511-8EB6-02163E011E08.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/070/00000/1A0794AF-AC0C-E511-8816-02163E012283.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/070/00000/1A51148F-AD0C-E511-86C2-02163E013942.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/070/00000/1CACD869-AC0C-E511-94E0-02163E014565.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/070/00000/229374A6-AF0C-E511-9B6A-02163E011883.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/070/00000/3C9C5206-AE0C-E511-9AC0-02163E012069.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/070/00000/40C7B9D0-B10C-E511-B5E6-02163E01420B.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/070/00000/44C44443-AE0C-E511-8B1A-02163E01184D.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/070/00000/4650F248-B80C-E511-A265-02163E01453C.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/070/00000/5ADBACBD-AD0C-E511-BEF8-02163E01447E.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/070/00000/789AE63B-AF0C-E511-9468-02163E011CFB.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/070/00000/7A0B4A95-E80C-E511-AD4E-02163E0144E5.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/070/00000/8057EA0F-AE0C-E511-AF8F-02163E0145E7.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/070/00000/868DBCBB-AD0C-E511-9F0A-02163E01184E.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/070/00000/8A9893EA-CA0C-E511-9B04-02163E0142D7.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/070/00000/8E95E690-AD0C-E511-AC0A-02163E013952.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/070/00000/8EA65BC1-AB0C-E511-A92F-02163E0134D5.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/070/00000/96CA7FEF-AF0C-E511-8DB7-02163E01432C.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/070/00000/96CD75FD-B20C-E511-A292-02163E0127ED.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/070/00000/A054B7B1-C20C-E511-BD0E-02163E0119D0.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/070/00000/AE7BC17A-AB0C-E511-895D-02163E01383E.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/070/00000/B45B8CAD-AE0C-E511-BB77-02163E011BD8.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/070/00000/B46062E1-B80C-E511-9AC2-02163E011968.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/070/00000/B66B3FEF-AC0C-E511-BD1C-02163E01180A.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/070/00000/B6A72B02-B30C-E511-9430-02163E0143FC.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/070/00000/B827F72E-AC0C-E511-BE09-02163E01369F.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/070/00000/C800B33A-BA0C-E511-B0D1-02163E013985.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/070/00000/CCB8FFD5-060D-E511-BAC5-02163E0145E3.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/070/00000/CE03B529-AC0C-E511-859D-02163E01447E.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/070/00000/D054A4BC-AE0C-E511-A3AD-02163E011D4A.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/070/00000/D2C59C4C-AE0C-E511-9230-02163E01424D.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/070/00000/D2CDFA2D-D00C-E511-8D86-02163E01190E.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/070/00000/D89F0796-BF0C-E511-AE72-02163E014682.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/070/00000/DA73BFAF-AE0C-E511-A69E-02163E01257B.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/070/00000/DE2C5F66-AF0C-E511-A4E9-02163E0138F3.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/070/00000/DE625040-AE0C-E511-A81B-02163E01446F.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/070/00000/E0610614-C40C-E511-9660-02163E013493.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/070/00000/E48E482F-BD0C-E511-A0EE-02163E01444B.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/070/00000/E6C036E7-AB0C-E511-B8F1-02163E011B55.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/070/00000/F063735B-AD0C-E511-B8BF-02163E013675.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/070/00000/F815B817-C50C-E511-9B16-02163E0136B5.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/070/00000/F86AD4F4-AD0C-E511-A312-02163E0139A3.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/070/00000/FC7D9717-B10C-E511-9241-02163E011A91.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/00223189-DA0C-E511-A02B-02163E011C0B.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/0083D6EB-EE0C-E511-8639-02163E013837.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/02591B88-C20C-E511-91A3-02163E0136A3.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/046815A6-D50C-E511-9512-02163E011966.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/04CB83DD-C20C-E511-9B2C-02163E014293.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/069D6F54-CD0C-E511-84A9-02163E014227.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/0A74E858-ED0C-E511-98F0-02163E01458C.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/0AF6AD83-C10C-E511-BF3F-02163E013445.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/0C0E1F1A-BE0C-E511-808F-02163E011A3E.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/0C1A6580-BF0C-E511-9EB0-02163E011C10.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/0E3E2DE0-EE0C-E511-B26A-02163E0140DB.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/10C54686-EC0C-E511-B06E-02163E012980.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/12486302-E80C-E511-BB64-02163E01442E.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/16BAFFA4-D30C-E511-B9F2-02163E0138F6.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/18DC8C8C-C20C-E511-AAD1-02163E011C5E.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/1A595E57-CA0C-E511-9769-02163E011915.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/20125F28-C60C-E511-809C-02163E0136EB.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/206D787E-F10C-E511-987D-02163E011D25.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/209BDE50-E00C-E511-B7E5-02163E011E08.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/20B8E00D-D70C-E511-896C-02163E012AE5.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/20BAC87B-C50C-E511-8625-02163E0143B6.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/225C8564-EA0C-E511-806B-02163E011865.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/22C63360-C70C-E511-8748-02163E0126A4.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/2641DDE5-CD0C-E511-85CB-02163E011A58.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/2A58E771-C70C-E511-A80C-02163E0142EA.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/2E1D4526-D80C-E511-B280-02163E012553.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/2EEAAA99-D20C-E511-A0D4-02163E01453D.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/3050F5F3-CF0C-E511-B476-02163E013642.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/306F46A6-D50C-E511-8DF0-02163E01295D.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/30B80D82-C10C-E511-A504-02163E011D9B.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/32040BCC-C20C-E511-AB41-02163E011A33.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/32737A8B-C50C-E511-8463-02163E0135CA.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/368F5EB7-C20C-E511-B51A-02163E011C9E.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/36E6C04B-D40C-E511-95FF-02163E01474E.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/38D15E7D-CC0C-E511-B774-02163E0141CE.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/3AE5B140-CD0C-E511-957A-02163E01410A.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/3C5EEA00-D20C-E511-ADD7-02163E011D7C.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/3CE9098B-C50C-E511-AFB4-02163E011C04.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/3E3F4814-C40C-E511-917E-02163E0141F9.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/425F05D7-CF0C-E511-B4B6-02163E014601.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/44145386-C10C-E511-B3D7-02163E01429A.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/44BFF3A3-C10C-E511-BDE9-02163E0137C2.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/481D4BA5-F70C-E511-86FE-02163E011D69.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/48A32086-C50C-E511-B138-02163E0141EB.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/48F46651-E00C-E511-92A5-02163E0146A9.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/4C61DF78-C50C-E511-9262-02163E01410A.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/4C71E6B2-E20C-E511-B1C3-02163E01357B.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/4E9BFDD9-FA0C-E511-9E0D-02163E0136B6.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/50D6F19D-C40C-E511-9D6C-02163E013946.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/5234C635-C60C-E511-B091-02163E014484.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/52692F07-C00C-E511-8253-02163E014409.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/52E1CB24-C50C-E511-AEFE-02163E01476E.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/52FFC6C8-CB0C-E511-9D5E-02163E014133.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/5477E8CD-C10C-E511-899D-02163E011A3E.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/5493245F-F60C-E511-8B6E-02163E0145B1.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/56007D56-C90C-E511-8329-02163E0121D4.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/560BE6EE-DA0C-E511-8DAC-02163E0120CB.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/56DD8E19-D10C-E511-99CE-02163E014120.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/5E7A9A61-EC0C-E511-A704-02163E0138AC.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/609E9775-BF0C-E511-909B-02163E011BE3.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/648B70A7-C30C-E511-B686-02163E011D69.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/66898348-D80C-E511-B3FF-02163E011A48.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/6A1534B4-C20C-E511-A7AB-02163E014613.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/6A770605-E70D-E511-A7E1-02163E012069.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/6AED1599-0A0D-E511-96A6-02163E012462.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/6C2587D8-CF0C-E511-9D79-02163E013684.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/6E6C1AAE-C80C-E511-935B-02163E011BE1.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/6EEBBF80-E50C-E511-9F06-02163E014718.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/708386D3-C30C-E511-9BF2-02163E01420B.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/7098C87B-BF0C-E511-B23C-02163E01475F.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/70F1B8A5-BF0C-E511-91BC-02163E0142EA.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/70FA7852-EA0C-E511-BEBA-02163E013463.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/74D7C480-F30C-E511-A675-02163E013729.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/766E3E39-DD0C-E511-A147-02163E011A15.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/820DDDB0-CB0C-E511-9CCA-02163E014204.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/8273DC91-C10C-E511-96E5-02163E01191E.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/843445D5-D40C-E511-9D54-02163E01240C.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/86BC0ED3-EC0C-E511-A57D-02163E0146B8.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/86F235FF-D10C-E511-8D30-02163E01295D.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/88FC8FE3-D70C-E511-A05C-02163E011FE0.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/8AA8FAF6-D80C-E511-86EF-02163E0138C4.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/8C44442E-000D-E511-A0A8-02163E011E08.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/8E886304-C00C-E511-9B45-02163E0146A5.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/904F1FFD-D80C-E511-870F-02163E01389A.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/905F1983-CB0C-E511-B884-02163E0141CC.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/964D0778-C20C-E511-8761-02163E013765.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/9AE3D22E-C40C-E511-B6C2-02163E013491.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/9E0B7388-C10C-E511-9F65-02163E01440F.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/A0ABC608-D30C-E511-B4B5-02163E013940.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/A46B33E4-EE0C-E511-816F-02163E011C8D.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/A6404F6E-E50C-E511-8723-02163E014185.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/A658DF69-F10C-E511-A5D9-02163E0128A1.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/A6B9F65C-BD0C-E511-BA4C-02163E011C04.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/A8A42CAE-C10C-E511-9254-02163E011B61.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/AC3BCB3A-CE0C-E511-BFCC-02163E011C5E.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/AE27AB72-BE0C-E511-8A54-02163E014119.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/B2016B08-CC0C-E511-ADA2-02163E011C21.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/B24AD6CC-C60C-E511-8F97-02163E013502.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/B2F1C1A3-DA0C-E511-8215-02163E01451C.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/BC3C3439-C40C-E511-95DB-02163E011BFD.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/BC6D0387-C10C-E511-BBE0-02163E01460E.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/C2D9DAD9-CA0C-E511-B634-02163E0118C8.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/C45A1047-C30C-E511-8D20-02163E011C14.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/C4E80894-C20C-E511-AA07-02163E0141D2.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/C4EE5663-C70C-E511-9423-02163E01448D.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/C616FE28-EA0C-E511-BA38-02163E014185.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/CA83262E-D80C-E511-8896-02163E014216.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/D2D94335-DD0C-E511-B9C4-02163E0138F2.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/D49B376A-BE0C-E511-9351-02163E014642.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/D8AC53B6-BE0C-E511-B855-02163E0136A7.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/D8C7A343-C60C-E511-AE16-02163E012017.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/D8E99CFF-DA0C-E511-9A58-02163E014484.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/DC3EFB6E-BE0C-E511-91D2-02163E0136E4.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/DC710316-C80C-E511-B142-02163E0142F4.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/DE93840A-C50C-E511-B27C-02163E01377B.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/E0A7CB9C-C10C-E511-9F0E-02163E012B06.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/E0BA7F8C-C40C-E511-AF55-02163E0146A9.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/E2F0B46C-E50C-E511-9197-02163E011CFE.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/E6984BED-C20C-E511-B973-02163E011A33.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/EE1CD15D-C80C-E511-9848-02163E011B3D.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/EEDADDEA-CA0C-E511-BBCD-02163E0143C5.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/F04AA792-C40C-E511-8BE8-02163E013502.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/F07F0541-C90C-E511-9C02-02163E01457D.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/F41BFE09-C20C-E511-A252-02163E014185.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/F424E841-BF0C-E511-92FB-02163E012462.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/F45FD7AA-C20C-E511-9F33-02163E0136B5.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/F633E643-EA0C-E511-9AC7-02163E0128A1.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/F6C26BF5-CE0C-E511-8E04-02163E0133E6.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/F8A7908E-C50C-E511-B49F-02163E011AF5.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/FA320973-C20C-E511-B43B-02163E0144FC.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/FA9D3EF9-D20C-E511-86F9-02163E014658.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/FAE4442B-030D-E511-9D10-02163E01443F.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/FCC5DD27-EA0C-E511-BC4B-02163E014122.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/073/00000/FE654145-D60C-E511-B946-02163E0142BF.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/0282AEB2-CC0C-E511-94EE-02163E0133E6.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/04A1DD5A-D90C-E511-AD7D-02163E011E39.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/04C73F43-CD0C-E511-9092-02163E014402.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/0E279DA2-EC0C-E511-A862-02163E0142F4.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/0E7A5602-FE0C-E511-BBD6-02163E0143EB.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/1424C08A-CF0C-E511-A301-02163E01432C.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/1851215E-DD0C-E511-8444-02163E011AF5.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/1AE5EEAC-DA0C-E511-91CA-02163E0118B8.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/1C4FEF60-D30C-E511-AD0F-02163E013502.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/1CB4A353-DB0C-E511-86C9-02163E0142D7.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/20ED3542-CD0C-E511-BA0A-02163E0135C8.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/2410C26C-E50C-E511-9B5B-02163E011A7E.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/24A534AD-CF0C-E511-9140-02163E011DE4.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/26754B0D-D20C-E511-A651-02163E013729.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/282C1B2B-EA0C-E511-8CD6-02163E011AAA.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/28682C46-DC0C-E511-9A32-02163E0145B7.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/2A77CD47-CE0C-E511-9989-02163E01190E.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/3059B982-CE0C-E511-8C38-02163E01369F.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/30A7C73C-D10C-E511-A9C0-02163E0139A3.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/324DF5AF-E20C-E511-BF85-02163E0140DD.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/365FF49F-E20C-E511-B0D8-02163E013946.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/3858B3C3-F50C-E511-8E50-02163E0140DB.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/3CDD86A7-EC0C-E511-935A-02163E01181B.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/40635AAC-D50C-E511-9766-02163E014216.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/44FE5F42-CD0C-E511-8F22-02163E013879.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/4AA1949F-E70C-E511-93B9-02163E0145CC.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/4E84C7A5-D10C-E511-B0AD-02163E011DE4.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/52450BD7-D30C-E511-8CA4-02163E011AF5.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/54927A41-F10C-E511-99C5-02163E0145DB.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/54BAE837-F10C-E511-84A9-02163E01457A.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/66F200D0-DD0C-E511-BB50-02163E014150.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/683868BE-CC0C-E511-A5A6-02163E01184E.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/68837E3E-CD0C-E511-8CB1-02163E013775.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/701DAB45-F70C-E511-BE3C-02163E0121DA.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/706F960F-D20C-E511-8E0C-02163E011AF5.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/72D7285A-E00C-E511-98E6-02163E01442F.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/76752A75-E50C-E511-89E2-02163E0133AE.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/80DAFDF8-E70C-E511-8C37-02163E013556.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/86D179C8-E50C-E511-96FB-02163E012283.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/92400425-CF0C-E511-8AF2-02163E01410A.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/9428C91C-D30C-E511-B74B-02163E014271.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/96D13AA8-D50C-E511-A5D9-02163E01440B.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/987901B3-040D-E511-95F1-02163E01460C.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/9AABB805-D60C-E511-B849-02163E014601.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/9C276880-CE0C-E511-98AB-02163E01184E.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/9EA462BE-FD0C-E511-9611-02163E0143CF.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/9EC8CF57-F10C-E511-9010-02163E013949.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/A021B80E-E80C-E511-90BE-02163E013454.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/AC4198DF-D30C-E511-839A-02163E0141F7.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/AE7D0B17-DC0C-E511-B41B-02163E0140FA.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/B0C0DCCA-D00C-E511-B5AF-02163E0146E0.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/B2F473AC-D10C-E511-BA61-02163E011E07.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/B4947B64-D40C-E511-9AFB-02163E012124.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/B888FC7A-E50C-E511-A566-02163E012BAA.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/BC615A06-E80C-E511-8A2A-02163E012A4D.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/BC6F938A-DF0D-E511-93B0-02163E012032.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/BCA17886-E60C-E511-94C7-02163E0145E3.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/BEA2099C-EC0C-E511-B086-02163E0139A9.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/C2BB6E74-DD0C-E511-931E-02163E0119DC.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/C63F9D35-E80C-E511-AD8B-02163E011D69.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/C8B3CAC2-D00C-E511-AFA6-02163E0135CA.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/C8D51B80-CD0C-E511-969F-02163E011FD2.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/CA150925-EF0C-E511-AC25-02163E013446.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/CACDF6A8-010D-E511-B7C0-02163E0138F3.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/CC0E155F-E00C-E511-9BF7-02163E0126EB.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/CCCC4AC8-CE0C-E511-9B11-02163E01297C.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/D44C9D06-D40C-E511-B5E1-02163E01257B.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/D4985B77-210D-E511-8511-02163E014384.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/D8E1A759-CE0C-E511-ABFB-02163E014601.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/DE5FC04D-D90C-E511-896F-02163E011BDB.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/DEE4638C-E50C-E511-A251-02163E014659.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/E04A22D5-D50C-E511-8351-02163E0138F6.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/EAD3E96E-C80C-E511-9AA2-02163E013653.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/ECED5C61-DD0C-E511-B463-02163E01369F.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/ECFAA193-CC0C-E511-BCB9-02163E0143B6.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/EE1AA9FF-D20C-E511-9853-02163E01396C.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/EE7938BC-D70C-E511-A232-02163E01478F.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/F203778F-D30C-E511-AA70-02163E014353.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/F66B0975-E50C-E511-BC78-02163E0138D9.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/F86C5510-CE0C-E511-979A-02163E011D9B.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/FAA1751A-DB0C-E511-96FD-02163E012BE1.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/078/00000/FEC2E074-CD0C-E511-A782-02163E01192B.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/079/00000/127F75C1-DA0C-E511-B02B-02163E0119CE.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/079/00000/140E7FF7-D20C-E511-9085-02163E013395.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/079/00000/1AC4EE66-F10C-E511-B73F-02163E0138CA.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/079/00000/34177A10-2C0D-E511-9883-02163E0145CC.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/079/00000/34A34599-DD0C-E511-9667-02163E013949.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/079/00000/3E5E7DAC-D70C-E511-ADB5-02163E0146D0.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/079/00000/4A4DA1F5-D80C-E511-A99E-02163E013899.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/079/00000/4CDF92CC-E50C-E511-8020-02163E013940.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/079/00000/648504D3-D70C-E511-93B9-02163E01410A.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/079/00000/6AD99220-ED0C-E511-91AC-02163E0146AC.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/079/00000/78EFAD1C-DA0C-E511-BF21-02163E013446.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/079/00000/7E226067-D60C-E511-B936-02163E01364D.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/079/00000/7E70089F-DA0C-E511-A6C1-02163E0144CF.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/079/00000/842EC8F4-E20C-E511-A8F1-02163E012298.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/079/00000/AC90FF83-DD0C-E511-A5E6-02163E011C7F.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/079/00000/C6237C37-D80C-E511-B808-02163E011F76.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/079/00000/D28E9704-D30C-E511-A017-02163E014476.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/079/00000/DAA305DF-FE0C-E511-A582-02163E01443F.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/079/00000/E27CBACB-E50C-E511-A0D8-02163E0133E6.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/079/00000/F2D0076D-D40C-E511-A0B0-02163E011E07.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/079/00000/FE669A03-E60C-E511-AF22-02163E013556.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/081/00000/08832611-E00C-E511-AB60-02163E011D9B.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/081/00000/108B2D7A-EC0C-E511-96DD-02163E01433E.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/081/00000/2085B09A-DA0C-E511-BAB6-02163E01344D.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/081/00000/32E86D4E-EA0C-E511-B486-02163E013907.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/081/00000/347EB558-DD0C-E511-A100-02163E01475F.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/081/00000/384ED732-EA0C-E511-99D4-02163E012884.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/081/00000/38FBD36F-E50C-E511-95C4-02163E011BCF.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/081/00000/427BBD49-D90C-E511-AEA9-02163E01450B.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/081/00000/523F422B-E00C-E511-98F2-02163E011A7E.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/081/00000/54F7D661-DD0C-E511-813C-02163E014111.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/081/00000/5859B96B-DB0C-E511-8634-02163E012283.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/081/00000/660D9370-E50C-E511-B8EB-02163E0136F3.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/081/00000/6C6CC3EF-EE0C-E511-A56A-02163E01460E.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/081/00000/72A17CB9-DB0C-E511-830A-02163E014384.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/081/00000/788EEC8C-F50C-E511-9A78-02163E012A4D.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/081/00000/925DE853-EA0C-E511-8704-02163E014607.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/081/00000/92D66C71-D90C-E511-A727-02163E011865.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/081/00000/9626AC7A-EC0C-E511-8437-02163E014147.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/081/00000/9CC6AC3A-EA0C-E511-A824-02163E012BAA.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/081/00000/9CFAC46E-DD0C-E511-944D-02163E01471D.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/081/00000/A4D34C9A-DA0C-E511-8474-02163E0134CC.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/081/00000/AE0E77E4-EE0C-E511-8FF2-02163E0133C2.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/081/00000/B0ECD0C8-DA0C-E511-893B-02163E0145B1.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/081/00000/B8C0A042-FB0C-E511-AADA-02163E01438B.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/081/00000/C488506E-DD0C-E511-947C-02163E014718.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/081/00000/C489D399-EC0C-E511-B7FE-02163E0138F8.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/081/00000/DCDDB426-E00C-E511-8BC9-02163E012A49.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/081/00000/DE0CD87A-EC0C-E511-BC8B-02163E013985.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/081/00000/EA47EBCC-DB0C-E511-9DD5-02163E014245.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/081/00000/EA9C98A0-E20C-E511-B59A-02163E011B58.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/252/00000/F2D7701A-E80D-E511-98B7-02163E014569.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/253/00000/C4CD92A7-E90D-E511-A110-02163E014699.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/262/00000/326B75FE-EF0D-E511-BC34-02163E01457A.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/266/00000/B288F8AD-F70D-E511-AE47-02163E01410A.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/267/00000/060021B4-0C0E-E511-8C75-02163E01476A.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/267/00000/0A313404-100E-E511-9F2E-02163E011D35.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/267/00000/109E9E0F-4F0E-E511-91F1-02163E014295.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/267/00000/120A5355-0A0E-E511-8009-02163E011C14.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/267/00000/18A19E41-120E-E511-A2AD-02163E0142BB.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/267/00000/202DD9B5-0B0E-E511-BF0E-02163E013942.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/267/00000/2C783F45-140E-E511-A38A-02163E014565.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/267/00000/325F1E37-1A0E-E511-BCE9-02163E013835.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/267/00000/36D01395-0B0E-E511-84F3-02163E011B61.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/267/00000/507E52F6-0C0E-E511-89F0-02163E0134D2.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/267/00000/5C6167EF-0A0E-E511-91A9-02163E0134D9.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/267/00000/6029C84C-1C0E-E511-9D12-02163E014166.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/267/00000/6A5562C8-250E-E511-BB5E-02163E01429D.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/267/00000/7C5A2DDF-090E-E511-A35E-02163E011FF2.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/267/00000/9825245B-0A0E-E511-8837-02163E0136A7.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/267/00000/9888085B-360E-E511-A72C-02163E0138F3.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/267/00000/9C1FDFFB-0A0E-E511-B5A6-02163E013502.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/267/00000/B45BC015-0C0E-E511-A932-02163E011AE2.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/267/00000/B8FA8744-1D0E-E511-8B02-02163E014166.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/267/00000/BA7E8B12-0C0E-E511-90D7-02163E011FBE.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/267/00000/C8488552-1D0E-E511-B946-02163E01379B.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/267/00000/CCC4DAA7-0B0E-E511-A0AA-02163E014742.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/267/00000/D83726E5-080E-E511-A059-02163E0134D9.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/267/00000/D89C201A-0B0E-E511-A931-02163E01456E.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/267/00000/DAF0F0E0-090E-E511-92F8-02163E013502.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/267/00000/E458FA32-0A0E-E511-ADF5-02163E014652.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/267/00000/E67460A8-0B0E-E511-9214-02163E0142E6.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/267/00000/E828F5E5-0A0E-E511-B642-02163E011953.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/267/00000/EEF01DCD-0E0E-E511-B75C-02163E014698.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/267/00000/F6524D2A-090E-E511-B044-02163E014338.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/267/00000/FC63835B-120E-E511-999F-02163E01228C.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/267/00000/FE8E0EAA-0B0E-E511-9D91-02163E0138F3.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/373/00000/8062961F-DB0E-E511-9B17-02163E014150.root',
'/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/247/376/00000/DCB61895-DD0E-E511-BD4F-02163E0138B3.root',
]
#process.source.fileNames.extend(
[
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
]
#process.source.fileNames.extend(
[
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
]
process.source.fileNames.extend(
[
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
])


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
