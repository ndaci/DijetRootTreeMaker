#ifndef DijetTreeProducer_h
#define DijetTreeProducer_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionData.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionEvaluator.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionParser.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "TTree.h"
#include "TH1F.h"
// For JECs
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"

// For 80X updates (for consumes system)
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "DataFormats/PatCandidates/interface/MET.h"


class DijetTreeProducer : public edm::EDAnalyzer 
{
 public:
  typedef reco::Particle::LorentzVector LorentzVector;
  explicit DijetTreeProducer(edm::ParameterSet const& cfg);
  virtual void beginJob();
  virtual void analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup);
  virtual void endJob();
  virtual ~DijetTreeProducer();
    

 private:  
  void initialize();
  // For JECs
  bool redoJECs_;
  //edm::FileInPath L1corrAK4_, L2corrAK4_, L3corrAK4_, ResCorrAK4_, L1corrAK8_, L2corrAK8_, L3corrAK8_, ResCorrAK8_;
  edm::FileInPath L1corrAK4_DATA_, L2corrAK4_DATA_, L3corrAK4_DATA_, ResCorrAK4_DATA_, L1corrAK8_DATA_, L2corrAK8_DATA_, L3corrAK8_DATA_, ResCorrAK8_DATA_;
  edm::FileInPath L1corrAK4_MC_, L2corrAK4_MC_, L3corrAK4_MC_, L1corrAK8_MC_, L2corrAK8_MC_, L3corrAK8_MC_;
  JetCorrectorParameters *L1ParAK4_DATA;
  JetCorrectorParameters *L2ParAK4_DATA;
  JetCorrectorParameters *L3ParAK4_DATA;
  JetCorrectorParameters *L2L3ResAK4_DATA;
  FactorizedJetCorrector *JetCorrectorAK4_DATA;
  JetCorrectorParameters *L1ParAK4_MC;
  JetCorrectorParameters *L2ParAK4_MC;
  JetCorrectorParameters *L3ParAK4_MC;
  FactorizedJetCorrector *JetCorrectorAK4_MC;
  JetCorrectorParameters *L1ParAK8_DATA;
  JetCorrectorParameters *L2ParAK8_DATA;
  JetCorrectorParameters *L3ParAK8_DATA;
  JetCorrectorParameters *L2L3ResAK8_DATA;
  FactorizedJetCorrector *JetCorrectorAK8_DATA;
  JetCorrectorParameters *L1ParAK8_MC;
  JetCorrectorParameters *L2ParAK8_MC;
  JetCorrectorParameters *L3ParAK8_MC;
  FactorizedJetCorrector *JetCorrectorAK8_MC;
  //---- configurable parameters --------   
  double ptMinAK4_,ptMinAK8_;
  bool isData_;
    
  // Migrate to consumes-system for running in 80X
    
  edm::EDGetTokenT<pat::JetCollection> srcJetsAK4_;
  edm::EDGetTokenT<pat::JetCollection> srcJetsAK8_;

  edm::EDGetTokenT<double> srcRho_;
  edm::EDGetTokenT<std::vector<pat::MET> > srcMET_;
  edm::EDGetTokenT<reco::VertexCollection> srcVrtx_;

  edm::EDGetTokenT<reco::GenJetCollection> srcGenJetsAK4_;
  edm::EDGetTokenT<reco::GenJetCollection> srcGenJetsAK8_;
  edm::EDGetTokenT<reco::GenParticleCollection> srcPrunedGenParticles_;
    
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > srcPU_;
  edm::EDGetTokenT<GenEventInfoProduct> srcGenInfo_;

  edm::Service<TFileService> fs_;
  TTree *outTree_;

  //---- TRIGGER -------------------------
  triggerExpression::Data triggerCache_;
  std::vector<triggerExpression::Evaluator*> vtriggerSelector_;
  std::vector<std::string> vtriggerAlias_,vtriggerSelection_;
  TH1F *triggerPassHisto_,*triggerNamesHisto_,*puHisto_;
  //---- output TREE variables ------
  //---- global event variables -----
  int   run_,evt_,nVtx_,lumi_;
  int   nJetsAK4_, nJetsAK8_, nGenJetsAK4_, nGenJetsAK8_;
  float rho_,met_,metSig_;
  float htAK4_,mjjAK4_,dEtajjAK4_,dPhijjAK4_;
  float htAK8_,mjjAK8_,dEtajjAK8_,dPhijjAK8_;
  //    float htCA8_,mjjCA8_,dEtajjCA8_,dPhijjCA8_;
  std::vector<bool> *triggerResult_;

  //---- NOISE FILTERS -------------------------
  triggerExpression::Data noiseFilterCache_;

  triggerExpression::Evaluator * HBHENoiseFilter_Selector_;
  triggerExpression::Evaluator * CSCHaloNoiseFilter_Selector_;
  triggerExpression::Evaluator * HCALlaserNoiseFilter_Selector_;
  triggerExpression::Evaluator * ECALDeadCellNoiseFilter_Selector_;
  triggerExpression::Evaluator * GoodVtxNoiseFilter_Selector_;
  triggerExpression::Evaluator * TrkFailureNoiseFilter_Selector_;
  triggerExpression::Evaluator * EEBadScNoiseFilter_Selector_;
  triggerExpression::Evaluator * ECALlaserNoiseFilter_Selector_;
  triggerExpression::Evaluator * TrkPOGNoiseFilter_Selector_;
  triggerExpression::Evaluator * TrkPOG_manystrip_NoiseFilter_Selector_;
  triggerExpression::Evaluator * TrkPOG_toomanystrip_NoiseFilter_Selector_;
  triggerExpression::Evaluator * TrkPOG_logError_NoiseFilter_Selector_;

  bool passFilterHBHE_;
  bool passFilterCSCHalo_;
  bool passFilterHCALlaser_;
  bool passFilterECALDeadCell_;
  bool passFilterGoodVtx_;
  bool passFilterTrkFailure_;
  bool passFilterEEBadSc_;
  bool passFilterECALlaser_;
  bool passFilterTrkPOG_;
  bool passFilterTrkPOG_manystrip_;
  bool passFilterTrkPOG_toomanystrip_;
  bool passFilterTrkPOG_logError_;

  //---- jet and genJet variables --------------
  std::vector<float> *ptAK4_,*jecAK4_,*etaAK4_,*phiAK4_,*massAK4_,*energyAK4_,*areaAK4_,*csvAK4_,*chfAK4_,*nhfAK4_,*phfAK4_,*elfAK4_,*mufAK4_,*nemfAK4_,*cemfAK4_;
  std::vector<int> *idLAK4_,*idTAK4_, *chHadMultAK4_, *chMultAK4_, *neHadMultAK4_, *neMultAK4_, *phoMultAK4_,*pFlavourAK4_,*hFlavourAK4_,*nbHadAK4_,*ncHadAK4_,*pFlavourAK8_,*hFlavourAK8_,*nbHadAK8_,*ncHadAK8_;
  std::vector<float> *hf_hfAK4_, *hf_emfAK4_, *hofAK4_;
  //std::vector<float> *cutbasedJetId_, *fullJetId_, *fullJetDiscriminant_;
  std::vector<float> *ptGenAK4_,*etaGenAK4_,*phiGenAK4_,*massGenAK4_,*energyGenAK4_;
    
  std::vector<float> *ptAK8_,*jecAK8_,*etaAK8_,*phiAK8_,*massAK8_,*energyAK8_,*areaAK8_,*csvAK8_,*chfAK8_,*nhfAK8_,*phfAK8_,*elfAK8_,*mufAK8_,*nemfAK8_,*cemfAK8_, *massPrunedAK8_, *massSoftDropAK8_, *dR_AK8_,*tau1AK8_,*tau2AK8_, *tau3AK8_ ;
  std::vector<int> *idLAK8_,*idTAK8_, *chHadMultAK8_, *chMultAK8_, *neHadMultAK8_, *neMultAK8_, *phoMultAK8_;
  std::vector<float> *hf_hfAK8_, *hf_emfAK8_, *hofAK8_;
  std::vector<float> *ptGenAK8_,*etaGenAK8_,*phiGenAK8_,*massGenAK8_,*energyGenAK8_;

  //---- MC variables ---------------
  std::vector<float> *npu_; 
  std::vector<int> *Number_interactions;
  std::vector <int> *OriginBX;
   
  //-----gen particles from hard scattering ------
  
  float ptHat_; 
  int processID_;
  float weight_;
  std::vector<float> *gen_eta, *gen_phi, *gen_p, *gen_px, *gen_py, *gen_pz, *gen_pt, *gen_energy,  *gen_vx, *gen_vy, *gen_vz;
  std::vector<int> *gen_numDaught, *gen_status, *gen_index, *gen_motherIndex, *gen_pdgId;  

};

#endif
