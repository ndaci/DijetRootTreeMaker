#include <iostream>
#include <sstream>
#include <istream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <functional>
#include <vector>
#include <cassert>
#include "TFile.h"
#include "TH1D.h"
#include "TMath.h"
#include "TLorentzVector.h"

#include "CMSDIJET/DijetRootTreeMaker/plugins/DijetTreeProducer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/JetCollection.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

using namespace std;
using namespace reco;
using namespace pat;
using namespace edm;

DijetTreeProducer::DijetTreeProducer(edm::ParameterSet const& cfg)
{

  // Extra AK8 jet collections from JetToolbox
  useJetTB_ = cfg.getParameter<bool>("useJetTB");
  if(useJetTB_) {
    srcJetsAK8_TB_ = (consumes<pat::JetCollection>(    cfg.getParameter<InputTag>("jetsAK8_TB")));
  }

  // Migrate to Consumes-system. Skip Calo-stuff
  srcJetsAK4_    = (consumes<pat::JetCollection>(    cfg.getParameter<InputTag>("jetsAK4")));
  srcJetsAK8_    = (consumes<pat::JetCollection>(    cfg.getParameter<InputTag>("jetsAK8")));
  
  srcRho_     = (consumes<double>(                cfg.getParameter<InputTag>("rho")));
  srcMET_     = (consumes<vector <pat::MET> >(    cfg.getParameter<InputTag>("met")));
  srcVrtx_    = (consumes<reco::VertexCollection>(cfg.getParameter<InputTag>("vtx")));
  
  ptMinAK4_   = cfg.getParameter<double>("ptMinAK4");
  ptMinAK8_   = cfg.getParameter<double>("ptMinAK8");
  
  srcPU_      = consumes<std::vector<PileupSummaryInfo> >(cfg.getUntrackedParameter<InputTag>("pu"));
  
  // These are now causing data run to fail. Weird it used to work with 2015 version?!
  isData_ = cfg.getParameter<bool>("isData");
  if (!isData_){
    srcGenJetsAK4_      = (consumes<GenJetCollection>(cfg.getParameter<InputTag>("genJetsAK4")));
    srcGenJetsAK8_      = (consumes<GenJetCollection>(cfg.getParameter<InputTag>("genJetsAK8")));
    srcPrunedGenParticles_ = (consumes<reco::GenParticleCollection>(cfg.getParameter<InputTag> ("genParticles")));
    srcGenInfo_           = consumes<GenEventInfoProduct>(cfg.getUntrackedParameter<InputTag>  ("ptHat"));
    //srcPU_              = cfg.getUntrackedParameter<InputTag>    ("pu",InputTag(""));
    //srcGenInfo_         = cfg.getUntrackedParameter<InputTag>  ("ptHat",InputTag());
  }

  triggerCache_       = triggerExpression::Data(cfg.getParameterSet("triggerConfiguration"),consumesCollector());
  vtriggerAlias_      = cfg.getParameter<std::vector<std::string> > ("triggerAlias");
  vtriggerSelection_  = cfg.getParameter<std::vector<std::string> > ("triggerSelection");
  noiseFilterCache_   = triggerExpression::Data(cfg.getParameterSet("noiseFilterConfiguration"),consumesCollector());

  HBHENoiseFilter_Selector_                 = triggerExpression::parse( cfg.getParameter<std::string> ("noiseFilterSelection_HBHENoiseFilter") );
  CSCHaloNoiseFilter_Selector_              = triggerExpression::parse( cfg.getParameter<std::string> ("noiseFilterSelection_CSCTightHaloFilter") );
  HCALlaserNoiseFilter_Selector_            = triggerExpression::parse( cfg.getParameter<std::string> ("noiseFilterSelection_hcalLaserEventFilter") );
  ECALDeadCellNoiseFilter_Selector_         = triggerExpression::parse( cfg.getParameter<std::string> ("noiseFilterSelection_EcalDeadCellTriggerPrimitiveFilter") );
  GoodVtxNoiseFilter_Selector_              = triggerExpression::parse( cfg.getParameter<std::string> ("noiseFilterSelection_goodVertices") );
  //TrkFailureNoiseFilter_Selector_           = triggerExpression::parse( cfg.getParameter<std::string> ("noiseFilterSelection_trackingFailureFilter") );
  EEBadScNoiseFilter_Selector_              = triggerExpression::parse( cfg.getParameter<std::string> ("noiseFilterSelection_eeBadScFilter") );
  ECALlaserNoiseFilter_Selector_            = triggerExpression::parse( cfg.getParameter<std::string> ("noiseFilterSelection_ecalLaserCorrFilter") );
  TrkPOGNoiseFilter_Selector_               = triggerExpression::parse( cfg.getParameter<std::string> ("noiseFilterSelection_trkPOGFilters") );
  TrkPOG_manystrip_NoiseFilter_Selector_    = triggerExpression::parse( cfg.getParameter<std::string> ("noiseFilterSelection_trkPOG_manystripclus53X") );
  TrkPOG_toomanystrip_NoiseFilter_Selector_ = triggerExpression::parse( cfg.getParameter<std::string> ("noiseFilterSelection_trkPOG_toomanystripclus53X") );
  TrkPOG_logError_NoiseFilter_Selector_     = triggerExpression::parse( cfg.getParameter<std::string> ("noiseFilterSelection_trkPOG_logErrorTooManyClusters") );

  if (vtriggerAlias_.size() != vtriggerSelection_.size()) {
    cout<<"ERROR: the number of trigger aliases does not match the number of trigger names !!!"<<endl;
    return;
  }
  for(unsigned i=0;i<vtriggerSelection_.size();i++) {
    vtriggerSelector_.push_back(triggerExpression::parse(vtriggerSelection_[i]));
  }
  
  // For JECs
  redoJECs_ = cfg.getParameter<bool>("redoJECs");
  // AK4 DATA
  L1corrAK4_DATA_  = cfg.getParameter<edm::FileInPath>("L1corrAK4_DATA");
  L2corrAK4_DATA_  = cfg.getParameter<edm::FileInPath>("L2corrAK4_DATA");
  L3corrAK4_DATA_  = cfg.getParameter<edm::FileInPath>("L3corrAK4_DATA");
  ResCorrAK4_DATA_ = cfg.getParameter<edm::FileInPath>("ResCorrAK4_DATA");
  // AK4 MC 
  L1corrAK4_MC_ = cfg.getParameter<edm::FileInPath>("L1corrAK4_MC");
  L2corrAK4_MC_ = cfg.getParameter<edm::FileInPath>("L2corrAK4_MC");
  L3corrAK4_MC_ = cfg.getParameter<edm::FileInPath>("L3corrAK4_MC");
  // AK8 DATA
  L1corrAK8_DATA_  = cfg.getParameter<edm::FileInPath>("L1corrAK8_DATA");
  L2corrAK8_DATA_  = cfg.getParameter<edm::FileInPath>("L2corrAK8_DATA");
  L3corrAK8_DATA_  = cfg.getParameter<edm::FileInPath>("L3corrAK8_DATA");
  ResCorrAK8_DATA_ = cfg.getParameter<edm::FileInPath>("ResCorrAK8_DATA");
  // AK8 MC
  L1corrAK8_MC_ = cfg.getParameter<edm::FileInPath>("L1corrAK8_MC");
  L2corrAK8_MC_ = cfg.getParameter<edm::FileInPath>("L2corrAK8_MC");
  L3corrAK8_MC_ = cfg.getParameter<edm::FileInPath>("L3corrAK8_MC");

  if(redoJECs_)
    {
      // AK4
      L1ParAK4_DATA = new JetCorrectorParameters(L1corrAK4_DATA_.fullPath());
      L2ParAK4_DATA = new JetCorrectorParameters(L2corrAK4_DATA_.fullPath());
      L3ParAK4_DATA = new JetCorrectorParameters(L3corrAK4_DATA_.fullPath());
      L2L3ResAK4_DATA = new JetCorrectorParameters(ResCorrAK4_DATA_.fullPath());
      L1ParAK4_MC = new JetCorrectorParameters(L1corrAK4_MC_.fullPath());
      L2ParAK4_MC = new JetCorrectorParameters(L2corrAK4_MC_.fullPath());
      L3ParAK4_MC = new JetCorrectorParameters(L3corrAK4_MC_.fullPath());

      std::vector<JetCorrectorParameters> vParAK4_DATA;
      std::vector<JetCorrectorParameters> vParAK4_MC;
      vParAK4_DATA.push_back(*L1ParAK4_DATA);
      vParAK4_DATA.push_back(*L2ParAK4_DATA);
      vParAK4_DATA.push_back(*L3ParAK4_DATA);
      vParAK4_DATA.push_back(*L2L3ResAK4_DATA);
      vParAK4_MC.push_back(*L1ParAK4_MC);
      vParAK4_MC.push_back(*L2ParAK4_MC);
      vParAK4_MC.push_back(*L3ParAK4_MC);

      JetCorrectorAK4_DATA = new FactorizedJetCorrector(vParAK4_DATA);
      JetCorrectorAK4_MC = new FactorizedJetCorrector(vParAK4_MC);

      // AK8
      L1ParAK8_DATA = new JetCorrectorParameters(L1corrAK8_DATA_.fullPath());
      L2ParAK8_DATA = new JetCorrectorParameters(L2corrAK8_DATA_.fullPath());
      L3ParAK8_DATA = new JetCorrectorParameters(L3corrAK8_DATA_.fullPath());
      L2L3ResAK8_DATA = new JetCorrectorParameters(ResCorrAK8_DATA_.fullPath());
      L1ParAK8_MC = new JetCorrectorParameters(L1corrAK8_MC_.fullPath());
      L2ParAK8_MC = new JetCorrectorParameters(L2corrAK8_MC_.fullPath());
      L3ParAK8_MC = new JetCorrectorParameters(L3corrAK8_MC_.fullPath());

      std::vector<JetCorrectorParameters> vParAK8_DATA;
      std::vector<JetCorrectorParameters> vParAK8_MC;
      vParAK8_DATA.push_back(*L1ParAK8_DATA);
      vParAK8_DATA.push_back(*L2ParAK8_DATA);
      vParAK8_DATA.push_back(*L3ParAK8_DATA);
      vParAK8_DATA.push_back(*L2L3ResAK8_DATA);
      vParAK8_MC.push_back(*L1ParAK8_MC);
      vParAK8_MC.push_back(*L2ParAK8_MC);
      vParAK8_MC.push_back(*L3ParAK8_MC);

      JetCorrectorAK8_DATA = new FactorizedJetCorrector(vParAK8_DATA);
      JetCorrectorAK8_MC = new FactorizedJetCorrector(vParAK8_MC);
    }
}

//////////////////////////////////////////////////////////////////////////////////////////
void DijetTreeProducer::beginJob() 
{

  //--- book the trigger histograms ---------
  triggerNamesHisto_ = fs_->make<TH1F>("TriggerNames","TriggerNames",1,0,1);
  //triggerNamesHisto_->SetBit(TH1::kCanRebin); // Does now work in CMSSW 806
  //triggerNamesHisto_->GetXaxis()->SetCanExtend(true);
  
  // Now the 'SetBit' procedure also could be omitted altogether, as ROOT6 change blog
  // suggests that it's not needed in this case:
  // "TAxis::kCanExtend bit is set on automatically for axis where all bins have label (i.e. when the axis is alphanumeric)."
  // https://root.cern.ch/content/main-histogram-changes-root-6
  //
  // Code compiles fine without this bit.
  
  for(unsigned i=0;i<vtriggerSelection_.size();i++) {
    triggerNamesHisto_->Fill(vtriggerSelection_[i].c_str(),1);
  }
  triggerPassHisto_ = fs_->make<TH1F>("TriggerPass","TriggerPass",1,0,1);
  //triggerPassHisto_->SetBit(TH1::kCanRebin); // Does now work in CMSSW 806
  //triggerPassHisto_->GetXaxis()->SetCanExtend(true);

  // Instantiate all vectors to define branches
  InstantiateVectorForBranches();
  
  // Book the tree and define branches
  outTree_ = fs_->make<TTree>("events","events");
  DefineBranches();

}

//////////////////////////////////////////////////////////////////////////////////////////
void DijetTreeProducer::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) 
{
  initialize();

  //edm::Handle<edm::View<pat::Jet> > jetsAK4;
  Handle<pat::JetCollection> jetsAK4;
  iEvent.getByToken(srcJetsAK4_,jetsAK4);

  //edm::Handle<edm::View<pat::Jet> > jetsAK8;
  Handle<pat::JetCollection> jetsAK8;
  iEvent.getByToken(srcJetsAK8_,jetsAK8);

  Handle<pat::JetCollection> jetsAK8_TB;
  if(useJetTB_) {
    iEvent.getByToken(srcJetsAK8_TB_,jetsAK8_TB);
  }

  //edm::Handle<edm::View<reco::GenJet> > handle_genJetsAK4;
  Handle<reco::GenJetCollection> handle_genJetsAK4;
  if (!iEvent.isRealData())
    iEvent.getByToken(srcGenJetsAK4_,handle_genJetsAK4);

  //edm::Handle<edm::View<reco::GenJet> > handle_genJetsAK8;
  Handle<reco::GenJetCollection> handle_genJetsAK8;
  if (!iEvent.isRealData())
    iEvent.getByToken(srcGenJetsAK8_,handle_genJetsAK8);
  
  Handle<double>  rho;
  iEvent.getByToken(srcRho_,rho);

  //edm::Handle<edm::View<pat::MET> >  met;'
  Handle<vector<pat::MET> > met;
  iEvent.getByToken(srcMET_,met);

  Handle<reco::VertexCollection> recVtxs;
  iEvent.getByToken(srcVrtx_,recVtxs);

  //-------------- Event Info -----------------------------------
  rho_    = *rho;
  met_    = (*met)[0].et();
  if ((*met)[0].sumEt() > 0) {
    metSig_ = (*met)[0].et()/(*met)[0].sumEt();
  }
  nVtx_   = recVtxs->size();
  run_    = iEvent.id().run();
  evt_    = iEvent.id().event();
  lumi_   = iEvent.id().luminosityBlock();

  //---------- pu -----------------------
  edm::Handle<std::vector<PileupSummaryInfo> > PupInfo;
  if (!iEvent.isRealData()) {
    iEvent.getByToken(srcPU_,PupInfo);
    
    //std::cout << "PupInfo.isValid()? : " << PupInfo.isValid() << endl;

    if(PupInfo.isValid()) {
      for( std::vector<PileupSummaryInfo>::const_iterator it = PupInfo->begin(); it != PupInfo->end(); ++it ) {
	npu_ -> push_back ( it -> getTrueNumInteractions() );
	Number_interactions -> push_back ( it->getPU_NumInteractions() ); 
	OriginBX -> push_back ( it -> getBunchCrossing());                
	
      }
    }
    else {
      //edm::LogError("DijetTreeProducer: PileUpError") << "Error! Can't get the product " << srcPU_;
      cout << "an edm::LogError call for PileUpError used to be here, but that does not work anymore -Juska" << endl;
    }
    
    // std::vector<PileupSummaryInfo>::const_iterator PUI;
    // for(PUI = PupInfo->begin(); PUI != PupInfo->end(); ++PUI) {
    //   if (PUI->getBunchCrossing() == 0) {
    //     npu_ = PUI->getTrueNumInteractions();
    //   }
    // }
  }// if MC
  
  //-------------- Gen Event Info -----------------------------------
  if (!iEvent.isRealData()) {

    edm::Handle<GenEventInfoProduct> genEvtInfo;
    iEvent.getByToken(srcGenInfo_,genEvtInfo);
    
    if( !genEvtInfo.isValid() ) {
      edm::LogInfo("GenEvtInfo") << "ERROR: genEvtInfo not valid! " << genEvtInfo;
    }

    if( genEvtInfo.isValid() ) {
      edm::LogInfo("GenEvtInfo") << "Successfully obtained " << genEvtInfo;
      ptHat_ = (genEvtInfo->hasBinningValues() ? genEvtInfo->binningValues()[0] : -999.);
      processID_ = genEvtInfo->signalProcessID();
      weight_ = genEvtInfo->weight();            
    }
    

    //------------------ Gen particles hard scattering -------------------
    //    (to be implemented)

    // to be saved only for partons that start the jet -> from genJets take the costituents -> 
    //see hypernews https://hypernews.cern.ch/HyperNews/CMS/get/csa14/49/2.html
    //and https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD#Advanced_topics_re_clustering_ev 


    edm::Handle<reco::GenParticleCollection> prunedGenParticles;
    if (!iEvent.isRealData())
      iEvent.getByToken(srcPrunedGenParticles_, prunedGenParticles);
    

    // std::cout << "-------------------------------" << endl;
    // std::cout << "   DEBUG   gen particles" << endl;
    // std::cout << "-------------------------------" << endl;
    // std::cout << "prunedGenParticles.failedToGet() = " << prunedGenParticles.isValid() << endl;
    // std::cout << "prunedGenParticles.isValid() = " << prunedGenParticles.isValid() << endl;
    
    if( prunedGenParticles.isValid() ) {
            
      for( reco::GenParticleCollection::const_iterator it = prunedGenParticles->begin(); it != prunedGenParticles->end(); ++it ) {
        // exit from loop when you reach the required number of GenParticles
        //if(eta->size() >= maxSize)
        //  break;
	
    	//save only particles from hard scattering 
	//already done from the pruner
    	//if(it->status()<21 || it->status()>29) continue; 
    	int idx = std::distance(prunedGenParticles->begin(),it);

        // fill in all the vectors
        gen_eta		->push_back( it->eta() );
        gen_phi		->push_back( it->phi() );
        gen_p		->push_back( it->p() );
        gen_px		->push_back( it->px() );
        gen_py		->push_back( it->py() );
        gen_pz		->push_back( it->pz() );
        gen_pt		->push_back( it->pt() );
        gen_energy	->push_back( it->energy() );
        gen_pdgId	->push_back( it->pdgId() );
        gen_vx		->push_back( it->vx() );
        gen_vy		->push_back( it->vy() );
        gen_vz		->push_back( it->vz() );
        gen_numDaught	->push_back( it->numberOfDaughters() );
        gen_status	->push_back( it->status() );
    	gen_index   	->push_back( idx );
  
    	int midx = -1;

	for( reco::GenParticleCollection::const_iterator mit = prunedGenParticles->begin(); mit != prunedGenParticles->end(); ++mit ) {
	
    	  if( it->mother()==&(*mit) ) {
    	    midx = std::distance(prunedGenParticles->begin(),mit);
    	    break;
    	  }
    	}
    	gen_motherIndex->push_back( midx );
	
	//cout << "id : " << idx << "   pdgId : " << it->pdgId() << "   status : " <<  it->status() << "   mother index : " << midx  << "  pt : " << it->pt() << "  pz : " << it->pz() << endl; 

      }//loop over genParticles
      //std::cout << "N gen particles saved = " << gen_index->size() << std::endl;  

    }
    
  }// if MC
  
  //-------------- Trigger Info -----------------------------------
  triggerPassHisto_->Fill("totalEvents",1);
  if (triggerCache_.setEvent(iEvent,iSetup)) {
    for(unsigned itrig=0;itrig<vtriggerSelector_.size();itrig++) {
      bool result(false);
      if (vtriggerSelector_[itrig]) {
        if (triggerCache_.configurationUpdated()) {
          vtriggerSelector_[itrig]->init(triggerCache_);
        }
        result = (*(vtriggerSelector_[itrig]))(triggerCache_);
      }
      if (result) {
        triggerPassHisto_->Fill(vtriggerAlias_[itrig].c_str(),1);
      }
      triggerResult_->push_back(result);
    }
  }

  // if (!iEvent.isRealData())
  //   {
      
  //-------------- Noise Filter Info -----------------------------------
  if (noiseFilterCache_.setEvent(iEvent,iSetup)) {
    
    if(noiseFilterCache_.configurationUpdated()) {
      HBHENoiseFilter_Selector_                -> init(noiseFilterCache_);
      CSCHaloNoiseFilter_Selector_             -> init(noiseFilterCache_);
      HCALlaserNoiseFilter_Selector_           -> init(noiseFilterCache_);
      ECALDeadCellNoiseFilter_Selector_        -> init(noiseFilterCache_);
      GoodVtxNoiseFilter_Selector_             -> init(noiseFilterCache_);
      //TrkFailureNoiseFilter_Selector_          -> init(noiseFilterCache_);
      EEBadScNoiseFilter_Selector_             -> init(noiseFilterCache_);
      ECALlaserNoiseFilter_Selector_           -> init(noiseFilterCache_);
      TrkPOGNoiseFilter_Selector_              -> init(noiseFilterCache_);
      TrkPOG_manystrip_NoiseFilter_Selector_   -> init(noiseFilterCache_);
      TrkPOG_toomanystrip_NoiseFilter_Selector_-> init(noiseFilterCache_);
      TrkPOG_logError_NoiseFilter_Selector_    -> init(noiseFilterCache_);
    }
	
    passFilterHBHE_                = (*HBHENoiseFilter_Selector_)                (noiseFilterCache_);    
    passFilterCSCHalo_             = (*CSCHaloNoiseFilter_Selector_)             (noiseFilterCache_);    
    passFilterHCALlaser_           = (*HCALlaserNoiseFilter_Selector_)           (noiseFilterCache_);    
    passFilterECALDeadCell_        = (*ECALDeadCellNoiseFilter_Selector_)        (noiseFilterCache_);    
    passFilterGoodVtx_             = (*GoodVtxNoiseFilter_Selector_)             (noiseFilterCache_);    
    //passFilterTrkFailure_          = (*TrkFailureNoiseFilter_Selector_)          (noiseFilterCache_);    
    passFilterEEBadSc_             = (*EEBadScNoiseFilter_Selector_)             (noiseFilterCache_);    
    passFilterECALlaser_           = (*ECALlaserNoiseFilter_Selector_)           (noiseFilterCache_);    
    passFilterTrkPOG_              = (*TrkPOGNoiseFilter_Selector_)              (noiseFilterCache_);    
    passFilterTrkPOG_manystrip_    = (*TrkPOG_manystrip_NoiseFilter_Selector_)   (noiseFilterCache_);    
    passFilterTrkPOG_toomanystrip_ = (*TrkPOG_toomanystrip_NoiseFilter_Selector_)(noiseFilterCache_);    
    passFilterTrkPOG_logError_     = (*TrkPOG_logError_NoiseFilter_Selector_)    (noiseFilterCache_);    
  }
  
  //----- at least one good vertex -----------
  //bool cut_vtx = (recVtxs->size() > 0);
  
  //if (cut_vtx) {

  // AK4
  std::vector<double> jecFactorsAK4;
  std::vector<unsigned> sortedAK4JetIdx;
  if(redoJECs_)
    {
      // sort AK4 jets by increasing pT
      std::multimap<double, unsigned> sortedAK4Jets;
      //for(edm::View<pat::Jet>::const_iterator ijet = jetsAK4->begin();ijet != jetsAK4->end(); ++ijet)
      for(pat::JetCollection::const_iterator ijet = jetsAK4->begin();ijet != jetsAK4->end(); ++ijet)
	{
	  double correction = 1.;
	  JetCorrectorAK4_DATA->setJetEta(ijet->eta());
	  JetCorrectorAK4_DATA->setJetPt(ijet->correctedJet(0).pt());
	  JetCorrectorAK4_DATA->setJetA(ijet->jetArea());
	  JetCorrectorAK4_DATA->setRho(rho_);
	  JetCorrectorAK4_MC->setJetEta(ijet->eta());
	  JetCorrectorAK4_MC->setJetPt(ijet->correctedJet(0).pt());
	  JetCorrectorAK4_MC->setJetA(ijet->jetArea());
	  JetCorrectorAK4_MC->setRho(rho_);
	  if (iEvent.isRealData()) 
	    correction = JetCorrectorAK4_DATA->getCorrection();
	  else
	    correction = JetCorrectorAK4_MC->getCorrection();


	  jecFactorsAK4.push_back(correction);
	  sortedAK4Jets.insert(std::make_pair(ijet->correctedJet(0).pt()*correction, ijet - jetsAK4->begin()));
	}
      // get jet indices in decreasing pT order
      for(std::multimap<double, unsigned>::const_reverse_iterator it = sortedAK4Jets.rbegin(); it != sortedAK4Jets.rend(); ++it)
        sortedAK4JetIdx.push_back(it->second);
    }
  else
    {
      for(pat::JetCollection::const_iterator ijet = jetsAK4->begin();ijet != jetsAK4->end(); ++ijet)
	{
	  jecFactorsAK4.push_back(1./ijet->jecFactor(0));
	  sortedAK4JetIdx.push_back(ijet - jetsAK4->begin());
	}
    }

  nJetsAK4_ = 0;
  float htAK4(0.0);
  vector<TLorentzVector> vP4AK4;
  for(std::vector<unsigned>::const_iterator i = sortedAK4JetIdx.begin(); i != sortedAK4JetIdx.end(); ++i) {

    pat::JetCollection::const_iterator ijet = (jetsAK4->begin() + *i);
    double chf = ijet->chargedHadronEnergyFraction();
    double nhf = ijet->neutralHadronEnergyFraction(); // + ijet->HFHadronEnergyFraction();
    double phf = ijet->photonEnergy()/(ijet->jecFactor(0) * ijet->energy());
    double elf = ijet->electronEnergy()/(ijet->jecFactor(0) * ijet->energy());
    //double muf = ijet->muonEnergy()/(ijet->jecFactor(0) * ijet->energy());
    double muf = ijet->muonEnergyFraction();

    double hf_hf = ijet->HFHadronEnergyFraction();
    double hf_emf= ijet->HFEMEnergyFraction();
    double hof   = ijet->hoEnergyFraction();

    int chm    = ijet->chargedHadronMultiplicity();
      
    int chMult = ijet->chargedMultiplicity();
    int neMult = ijet->neutralMultiplicity();
    int npr    = chMult + neMult;

    int chHadMult = chm; //ijet->chargedHadronMultiplicity();
    int neHadMult = ijet->neutralHadronMultiplicity();
    int phoMult = ijet->photonMultiplicity();
      
    // Juska's added fractions for identical JetID with recommendations
    double nemf = ijet->neutralEmEnergyFraction();
    double cemf = ijet->chargedEmEnergyFraction();
    int NumConst = npr;

    float eta  = ijet->eta(); // removed fabs() -Juska
    float pt   = ijet->correctedJet(0).pt()*jecFactorsAK4.at(*i); // Is this OK? Correct corrected? -Juska

    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID
    int idL = (nhf<0.99 && nemf<0.99 && NumConst>1 && muf < 0.8) && ((fabs(eta) <= 2.4 && chf>0 && chMult>0 && cemf<0.99) || fabs(eta)>2.4);
    int idT = (nhf<0.90 && nemf<0.90 && NumConst>1 && muf<0.8) && ((fabs(eta)<=2.4 && chf>0 && chMult>0 && cemf<0.90) || fabs(eta)>2.4);

       
      
    if (pt > ptMinAK4_) {
      htAK4 += pt;
      nJetsAK4_++;

      vP4AK4.push_back(TLorentzVector(ijet->correctedJet(0).px()*jecFactorsAK4.at(*i),ijet->correctedJet(0).py()*jecFactorsAK4.at(*i),ijet->correctedJet(0).pz()*jecFactorsAK4.at(*i),ijet->correctedJet(0).energy()*jecFactorsAK4.at(*i)));
      chfAK4_           ->push_back(chf);
      nhfAK4_           ->push_back(nhf);
      phfAK4_           ->push_back(phf);
      elfAK4_           ->push_back(elf);
      mufAK4_           ->push_back(muf);
      nemfAK4_          ->push_back(nemf);
      cemfAK4_          ->push_back(cemf);
      hf_hfAK4_         ->push_back(hf_hf);
      hf_emfAK4_        ->push_back(hf_emf);
      hofAK4_           ->push_back(hof);
      jecAK4_           ->push_back(jecFactorsAK4.at(*i));
      ptAK4_            ->push_back(pt);
      phiAK4_           ->push_back(ijet->phi());
      etaAK4_           ->push_back(ijet->eta());
      massAK4_          ->push_back(ijet->correctedJet(0).mass()*jecFactorsAK4.at(*i));
      energyAK4_        ->push_back(ijet->correctedJet(0).energy()*jecFactorsAK4.at(*i));
      areaAK4_          ->push_back(ijet->jetArea());
      csvAK4_           ->push_back(ijet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
      pFlavourAK4_      ->push_back(ijet->partonFlavour());
      hFlavourAK4_      ->push_back(ijet->hadronFlavour());
      nbHadAK4_         ->push_back(ijet->jetFlavourInfo().getbHadrons().size());
      ncHadAK4_         ->push_back(ijet->jetFlavourInfo().getcHadrons().size());
      idLAK4_           ->push_back(idL);
      idTAK4_           ->push_back(idT);
      chHadMultAK4_     ->push_back(chHadMult);
      chMultAK4_        ->push_back(chMult);
      neHadMultAK4_     ->push_back(neHadMult);  
      neMultAK4_        ->push_back(neMult);
      phoMultAK4_       ->push_back(phoMult); 

    }

  }// jet loop  

  htAK4_     = htAK4;
  if (nJetsAK4_ > 1) { //assuming jets are ordered by pt in the pat collection
    mjjAK4_    = (vP4AK4[0]+vP4AK4[1]).M();
    dEtajjAK4_ = fabs((*etaAK4_)[0]-(*etaAK4_)[1]); 
    dPhijjAK4_ = fabs(deltaPhi((*phiAK4_)[0],(*phiAK4_)[1]));
  }

  // Fill reco AK8 Jets from MINIAOD contents
  FillJetsAK8(iEvent, jetsAK8, 0); // idxColl=0 <-> original MINIAOD collection
  if(useJetTB_) {
    FillJetsAK8(iEvent, jetsAK8_TB, 1); // idxColl=1 <-> jets from JetToolbox
  }    

  //-------------- Gen Jets Info -----------------------------------

  if (!iEvent.isRealData()) {

    //AK4
    nGenJetsAK4_ = 0;
    vector<TLorentzVector> vP4GenAK4;
    reco::GenJetCollection genJetsAK4 = *handle_genJetsAK4; 
    for(reco::GenJetCollection::const_iterator ijet = genJetsAK4.begin();ijet != genJetsAK4.end(); ++ijet) { 	
      //float eta  = fabs(ijet->eta());
      float pt   = ijet->pt();
      if (pt > ptMinAK4_) {
	nGenJetsAK4_++;
	vP4GenAK4.push_back(TLorentzVector(ijet->px(),ijet->py(),ijet->pz(),ijet->energy()));
	ptGenAK4_            ->push_back(pt);
	phiGenAK4_           ->push_back(ijet->phi());
	etaGenAK4_           ->push_back(ijet->eta());
	massGenAK4_          ->push_back(ijet->mass());
	energyGenAK4_        ->push_back(ijet->energy());
      }
    }// jet loop  
      
    //AK8
    nGenJetsAK8_ = 0;
    vector<TLorentzVector> vP4GenAK8;      
    reco::GenJetCollection genJetsAK8 = *handle_genJetsAK8;
    for(reco::GenJetCollection::const_iterator ijet = genJetsAK8.begin();ijet != genJetsAK8.end(); ++ijet) { 	
      //float eta  = fabs(ijet->eta());
      float pt   = ijet->pt();
      if (pt > ptMinAK8_) {
	nGenJetsAK8_++;
	vP4GenAK8.push_back(TLorentzVector(ijet->px(),ijet->py(),ijet->pz(),ijet->energy()));
	ptGenAK8_            ->push_back(pt);
	phiGenAK8_           ->push_back(ijet->phi());
	etaGenAK8_           ->push_back(ijet->eta());
	massGenAK8_          ->push_back(ijet->mass());
	energyGenAK8_        ->push_back(ijet->energy());
      }
    }// jet loop  
      
  }//if MC 

  //  }// if vtx
  
  
  //---- Fill Tree ---
  //if (mjjAK4_ > mjjMin_ && dEtajjAK4_ < dEtaMax_) {
  outTree_->Fill();     
  //}
  //------------------
  
  
}//end analyze for each event

int DijetTreeProducer::FillJetsAK8(edm::Event const& iEvent, const edm::Handle<pat::JetCollection> &jetsAK8, const u_int idxCollAK8)
{

  if(idxCollAK8 >= nCollAK8) {
    cout << "DijetTreeProducer::Error FillJetsAK8() called with idxCollAK8=" 
	 << idxCollAK8 << " >= " << nCollAK8 << endl;
    return -1;
  }

  // Name of jettiness depending on input collection (0: standard miniaod; 1: JetToolbox)
  TString nameJettiness[nCollAK8] = {"NjettinessAK8", "NjettinessAK8CHS"};

  // AK8
  std::vector<double> jecFactorsAK8;
  std::vector<unsigned> sortedAK8JetIdx;
  if(redoJECs_)
    {
      // sort AK8 jets by increasing pT
      std::multimap<double, unsigned> sortedAK8Jets;
      for(pat::JetCollection::const_iterator ijet = jetsAK8->begin();ijet != jetsAK8->end(); ++ijet)
	{
	  double correction = 1.;

	  JetCorrectorAK8_DATA->setJetEta(ijet->eta());
	  JetCorrectorAK8_DATA->setJetPt(ijet->correctedJet(0).pt());
	  JetCorrectorAK8_DATA->setJetA(ijet->jetArea());
	  JetCorrectorAK8_DATA->setRho(rho_);
	  JetCorrectorAK8_MC->setJetEta(ijet->eta());
	  JetCorrectorAK8_MC->setJetPt(ijet->correctedJet(0).pt());
	  JetCorrectorAK8_MC->setJetA(ijet->jetArea());
	  JetCorrectorAK8_MC->setRho(rho_);

	  if (iEvent.isRealData()) 
	    correction = JetCorrectorAK8_DATA->getCorrection();
	  else
	    correction = JetCorrectorAK8_MC->getCorrection();

	  jecFactorsAK8.push_back(correction);
	  sortedAK8Jets.insert(std::make_pair(ijet->correctedJet(0).pt()*correction, ijet - jetsAK8->begin()));
	}
      // get jet indices in decreasing pT order
      for(std::multimap<double, unsigned>::const_reverse_iterator it = sortedAK8Jets.rbegin(); it != sortedAK8Jets.rend(); ++it)
        sortedAK8JetIdx.push_back(it->second);
    }
  else
    {
      for(pat::JetCollection::const_iterator ijet = jetsAK8->begin();ijet != jetsAK8->end(); ++ijet)
	{
	  jecFactorsAK8.push_back(1./ijet->jecFactor(0));
	  sortedAK8JetIdx.push_back(ijet - jetsAK8->begin());
	}
    }

  nJetsAK8_ = 0;
  float htAK8(0.0);
  vector<TLorentzVector> vP4AK8;
  TLorentzVector tempP4;
  TLorentzVector puppi_softdrop, puppi_softdrop_subjet;

  for(std::vector<unsigned>::const_iterator i = sortedAK8JetIdx.begin(); i != sortedAK8JetIdx.end(); ++i) {

    pat::JetCollection::const_iterator ijet = (jetsAK8->begin() + *i);

    double chf = ijet->chargedHadronEnergyFraction();
    double nhf = ijet->neutralHadronEnergyFraction(); // + ijet->HFHadronEnergyFraction();
    double phf = ijet->photonEnergy()/(ijet->jecFactor(0) * ijet->energy());
    double elf = ijet->electronEnergy()/(ijet->jecFactor(0) * ijet->energy());
    //double muf = ijet->muonEnergy()/(ijet->jecFactor(0) * ijet->energy());
    double muf = ijet->muonEnergyFraction();

    double hf_hf = ijet->HFHadronEnergyFraction();
    double hf_emf= ijet->HFEMEnergyFraction();
    double hof    = ijet->hoEnergyFraction();

    int chm    = ijet->chargedHadronMultiplicity();
      
    int chMult = ijet->chargedMultiplicity();
    int neMult = ijet->neutralMultiplicity();
    int npr    = chMult + neMult;

    int chHadMult = chm; //ijet->chargedHadronMultiplicity();
    int neHadMult = ijet->neutralHadronMultiplicity();
    int phoMult   = ijet->photonMultiplicity();
      
    // Juska's added fractions for identical JetID with recommendations
    double nemf  = ijet->neutralEmEnergyFraction();
    double cemf  = ijet->chargedEmEnergyFraction();
    int NumConst = npr;

    float eta  = ijet->eta(); // removed fabs() -Juska
    float pt   = ijet->correctedJet(0).pt()*jecFactorsAK8.at(*i); // Is this OK? Correct corrected? -Juska

    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID
    int idL = (nhf<0.99 && nemf<0.99 && NumConst>1 && muf < 0.8) && ((fabs(eta) <= 2.4 && chf>0 && chMult>0 && cemf<0.99) || fabs(eta)>2.4);
    int idT = (nhf<0.90 && nemf<0.90 && NumConst>1 && muf<0.8) && ((fabs(eta)<=2.4 && chf>0 && chMult>0 && cemf<0.90) || fabs(eta)>2.4);

    // Sub-structure
    float jchs_tau1, jchs_tau2, jchs_tau3, jchs_m_pruned, jchs_m_softdrop;
    float jpuppi_pt, jpuppi_eta, jpuppi_phi, jpuppi_mass, jpuppi_tau1, jpuppi_tau2, jpuppi_tau3;
    ///
    jchs_tau1 = GetUserFloat( *ijet , nameJettiness[idxCollAK8]+":tau1");
    jchs_tau2 = GetUserFloat( *ijet , nameJettiness[idxCollAK8]+":tau2");
    jchs_tau3 = GetUserFloat( *ijet , nameJettiness[idxCollAK8]+":tau3");
    //
    jchs_m_pruned   = GetUserFloat( *ijet , "ak8PFJetsCHSPrunedMass" );
    jchs_m_softdrop = GetUserFloat( *ijet , "ak8PFJetsCHSSoftDropMass" );
    //
    jpuppi_pt   = GetUserFloat( *ijet , "ak8PFJetsPuppiValueMap:pt" );
    jpuppi_eta  = GetUserFloat( *ijet , "ak8PFJetsPuppiValueMap:eta" );
    jpuppi_phi  = GetUserFloat( *ijet , "ak8PFJetsPuppiValueMap:phi" );
    jpuppi_mass = GetUserFloat( *ijet , "ak8PFJetsPuppiValueMap:mass" );
    jpuppi_tau1 = GetUserFloat( *ijet , "ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau1" );
    jpuppi_tau2 = GetUserFloat( *ijet , "ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau2" );
    jpuppi_tau3 = GetUserFloat( *ijet , "ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau3" );
      
    if (pt > ptMinAK8_) {

      // Counting jets
      nJetsAK8_++;

      // HT
      htAK8 += pt;

      // 4-momentum
      tempP4 = TLorentzVector( ijet->correctedJet(0).px()*jecFactorsAK8.at(*i),
			       ijet->correctedJet(0).py()*jecFactorsAK8.at(*i),
			       ijet->correctedJet(0).pz()*jecFactorsAK8.at(*i),
			       ijet->correctedJet(0).energy()*jecFactorsAK8.at(*i) 
			       );
      vP4AK8.push_back(tempP4);

      // Kinematics
      ptAK8_            [idxCollAK8]->push_back(pt);
      phiAK8_           [idxCollAK8]->push_back(ijet->phi());
      etaAK8_           [idxCollAK8]->push_back(ijet->eta());

      // Energy
      jecAK8_           [idxCollAK8]->push_back(jecFactorsAK8.at(*i));
      massAK8_          [idxCollAK8]->push_back(ijet->correctedJet(0).mass()*jecFactorsAK8.at(*i));
      energyAK8_        [idxCollAK8]->push_back(ijet->correctedJet(0).energy()*jecFactorsAK8.at(*i));
      areaAK8_          [idxCollAK8]->push_back(ijet->jetArea());

      // Energy fractions
      chfAK8_           [idxCollAK8]->push_back(chf);
      nhfAK8_           [idxCollAK8]->push_back(nhf);
      phfAK8_           [idxCollAK8]->push_back(phf);
      elfAK8_           [idxCollAK8]->push_back(elf);
      mufAK8_           [idxCollAK8]->push_back(muf);
      nemfAK8_          [idxCollAK8]->push_back(nemf);
      cemfAK8_          [idxCollAK8]->push_back(cemf);
      hf_hfAK8_         [idxCollAK8]->push_back(hf_hf);
      hf_emfAK8_        [idxCollAK8]->push_back(hf_emf);
      hofAK8_           [idxCollAK8]->push_back(hof);

      // Multiplicities
      chHadMultAK8_     [idxCollAK8]->push_back(chHadMult);
      chMultAK8_        [idxCollAK8]->push_back(chMult);
      neHadMultAK8_     [idxCollAK8]->push_back(neHadMult);  
      neMultAK8_        [idxCollAK8]->push_back(neMult);
      phoMultAK8_       [idxCollAK8]->push_back(phoMult); 

      // ID variables
      idLAK8_           [idxCollAK8]->push_back(idL);
      idTAK8_           [idxCollAK8]->push_back(idT);

      // Flavour
      csvAK8_           [idxCollAK8]->push_back(ijet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
      pFlavourAK8_      [idxCollAK8]->push_back(ijet->partonFlavour());
      hFlavourAK8_      [idxCollAK8]->push_back(ijet->hadronFlavour());
      nbHadAK8_         [idxCollAK8]->push_back(ijet->jetFlavourInfo().getbHadrons().size());
      ncHadAK8_         [idxCollAK8]->push_back(ijet->jetFlavourInfo().getcHadrons().size());

      // Sub-structure
      tau1AK8_          [idxCollAK8]->push_back(jchs_tau1);
      tau2AK8_          [idxCollAK8]->push_back(jchs_tau2);
      tau3AK8_          [idxCollAK8]->push_back(jchs_tau3);
      massPrunedAK8_    [idxCollAK8]->push_back(jchs_m_pruned);
      massSoftDropAK8_  [idxCollAK8]->push_back(jchs_m_softdrop);

      // Associated PUPPI jets //
      ///
      /// kinematics
      ptAK8_Puppi_   [idxCollAK8]->push_back(jpuppi_pt);
      etaAK8_Puppi_  [idxCollAK8]->push_back(jpuppi_eta);
      phiAK8_Puppi_  [idxCollAK8]->push_back(jpuppi_phi);
      massAK8_Puppi_ [idxCollAK8]->push_back(jpuppi_mass);
      ///
      /// sub-structure
      tau1AK8_Puppi_ [idxCollAK8]->push_back(jpuppi_tau1);
      tau2AK8_Puppi_ [idxCollAK8]->push_back(jpuppi_tau2);
      tau3AK8_Puppi_ [idxCollAK8]->push_back(jpuppi_tau3);
      ///

      /*
      /// soft-drop mass
      puppi_softdrop        = TLorentzVector();
      puppi_softdrop_subjet = TLorentzVector();
      auto const & sdSubjetsPuppi = ijet->subjets("SoftDropPuppi");
      for ( auto const & it : sdSubjetsPuppi ) {
	puppi_softdrop_subjet.SetPtEtaPhiM(it->correctedP4(0).pt(),it->correctedP4(0).eta(),it->correctedP4(0).phi(),it->correctedP4(0).mass());
	puppi_softdrop+=puppi_softdrop_subjet;
      }
      massSoftDropAK8_Puppi_ [idxCollAK8]->push_back( puppi_softdrop.M() );
      */

      //---- match with the pruned jet collection -----
      // double dRmin(1000);
      // double auxm(0.0);
      // for(edm::View<pat::Jet>::const_iterator ijetpr = jetsAK8->begin();ijetpr != jetsAK8->end(); ++ijetpr) { 
      //   float dR = deltaR(ijet->eta(),ijet->phi(),ijetpr->eta(),ijetpr->phi());
      //   if (dR < dRmin) {
      //     auxm = ijetpr->mass();
      //     dRmin = dR;
      //   } 
      // } 
      // massPruned_->push_back(auxm);
      // dR_->push_back(dRmin);
	
    }
  }// jet loop  

  htAK8_     = htAK8;
  if (nJetsAK8_ > 1) { //assuming jets are ordered by pt in the pat collection
    mjjAK8_    = (vP4AK8[0]+vP4AK8[1]).M();
    dEtajjAK8_ = fabs( ( *etaAK8_[idxCollAK8] )[0] - ( *etaAK8_[idxCollAK8] )[1] ); 
    dPhijjAK8_ = fabs(deltaPhi( (*phiAK8_[idxCollAK8])[0] , (*phiAK8_[idxCollAK8])[1] ) );
  }

  return 0;

}


//////////////////////////////////////////////////////////////////////////////////////////
DijetTreeProducer::~DijetTreeProducer() 
{
}

int DijetTreeProducer::InstantiateVectorForBranches()
{

  triggerResult_ = new std::vector<bool>;

  gen_eta          = new std::vector<float>;
  gen_phi          = new std::vector<float>;
  gen_p            = new std::vector<float>;
  gen_px           = new std::vector<float>;
  gen_py           = new std::vector<float>;
  gen_pz           = new std::vector<float>;
  gen_pt           = new std::vector<float>;
  gen_energy	   = new std::vector<float>; 
  gen_pdgId	   = new std::vector<int>; 
  gen_vx	   = new std::vector<float>; 
  gen_vy	   = new std::vector<float>; 
  gen_vz	   = new std::vector<float>;	 
  gen_numDaught    = new std::vector<int>;      
  gen_status	   = new std::vector<int>; 
  gen_index	   = new std::vector<int>; 
  gen_motherIndex  = new std::vector<int>; 

  ptAK4_             = new std::vector<float>;
  jecAK4_            = new std::vector<float>;
  etaAK4_            = new std::vector<float>;
  phiAK4_            = new std::vector<float>;
  massAK4_           = new std::vector<float>;
  energyAK4_         = new std::vector<float>;
  areaAK4_           = new std::vector<float>;
  csvAK4_            = new std::vector<float>;
  pFlavourAK4_       = new std::vector<int>;
  hFlavourAK4_       = new std::vector<int>;
  nbHadAK4_          = new std::vector<int>;
  ncHadAK4_          = new std::vector<int>;
  chfAK4_            = new std::vector<float>;
  nhfAK4_            = new std::vector<float>;
  phfAK4_            = new std::vector<float>;
  mufAK4_            = new std::vector<float>;
  elfAK4_            = new std::vector<float>;
  nemfAK4_           = new std::vector<float>;
  cemfAK4_           = new std::vector<float>;
  // Hadronic forward hadrons
  hf_hfAK4_          = new std::vector<float>;
  // Hadronic forward electromagnetic fraction
  hf_emfAK4_         = new std::vector<float>;
  hofAK4_            = new std::vector<float>;
  idLAK4_            = new std::vector<int>;
  idTAK4_            = new std::vector<int>;
  chHadMultAK4_     = new std::vector<int>;
  chMultAK4_         = new std::vector<int>;
  neHadMultAK4_      = new std::vector<int>;
  neMultAK4_         = new std::vector<int>;
  phoMultAK4_        = new std::vector<int>;


  // instantiate vectors for all AK8 collections available
  for(u_int i=0 ; i<nCollAK8 ; i++) {

    ptAK8_[i]             = new std::vector<float>;
    jecAK8_[i]            = new std::vector<float>;
    etaAK8_[i]            = new std::vector<float>;
    phiAK8_[i]            = new std::vector<float>;
    massAK8_[i]           = new std::vector<float>;
    energyAK8_[i]         = new std::vector<float>;
    areaAK8_[i]           = new std::vector<float>;
    csvAK8_[i]            = new std::vector<float>;
    pFlavourAK8_[i]       = new std::vector<int>;
    hFlavourAK8_[i]       = new std::vector<int>;
    nbHadAK8_[i]          = new std::vector<int>;
    ncHadAK8_[i]          = new std::vector<int>;
    chfAK8_[i]            = new std::vector<float>;
    nhfAK8_[i]            = new std::vector<float>;
    phfAK8_[i]            = new std::vector<float>;
    mufAK8_[i]            = new std::vector<float>;
    elfAK8_[i]            = new std::vector<float>;
    nemfAK8_[i]           = new std::vector<float>;
    cemfAK8_[i]           = new std::vector<float>;
    // Hadronic forward hadrons
    hf_hfAK8_[i]          = new std::vector<float>;
    // Hadronic forward photons
    hf_emfAK8_[i]         = new std::vector<float>;
    hofAK8_[i]            = new std::vector<float>;
    idLAK8_[i]            = new std::vector<int>;
    idTAK8_[i]            = new std::vector<int>;
    massPrunedAK8_[i]     = new std::vector<float>;
    massSoftDropAK8_[i]   = new std::vector<float>;
    tau1AK8_[i]           = new std::vector<float>;
    tau2AK8_[i]           = new std::vector<float>;
    tau3AK8_[i]           = new std::vector<float>;
    chHadMultAK8_[i]      = new std::vector<int>;   
    chMultAK8_[i]         = new std::vector<int>;
    neHadMultAK8_[i]      = new std::vector<int>; 
    neMultAK8_[i]         = new std::vector<int>;
    phoMultAK8_[i]        = new std::vector<int>;
    //
    ptAK8_Puppi_[i]       = new std::vector<float>;
    etaAK8_Puppi_[i]      = new std::vector<float>;
    phiAK8_Puppi_[i]      = new std::vector<float>;
    massAK8_Puppi_[i]     = new std::vector<float>;
    tau1AK8_Puppi_[i]     = new std::vector<float>;
    tau2AK8_Puppi_[i]     = new std::vector<float>;
    tau3AK8_Puppi_[i]     = new std::vector<float>;
    massSoftDropAK8_Puppi_[i] = new std::vector<float>;

  }

  npu_                = new std::vector<float>;  
  Number_interactions = new std::vector<int>;
  OriginBX            = new std::vector<int>;

  ptGenAK4_             = new std::vector<float>;
  etaGenAK4_            = new std::vector<float>;
  phiGenAK4_            = new std::vector<float>;
  massGenAK4_           = new std::vector<float>;
  energyGenAK4_         = new std::vector<float>;
  ptGenAK8_             = new std::vector<float>;
  etaGenAK8_            = new std::vector<float>;
  phiGenAK8_            = new std::vector<float>;
  massGenAK8_           = new std::vector<float>;
  energyGenAK8_         = new std::vector<float>;

  return 0;
}

int DijetTreeProducer::DefineBranches()
{

  outTree_->Branch("runNo"                ,&run_               ,"run_/I");
  outTree_->Branch("evtNo"                ,&evt_               ,"evt_/I");
  outTree_->Branch("lumi"                 ,&lumi_              ,"lumi_/I");
  outTree_->Branch("nvtx"                 ,&nVtx_              ,"nVtx_/I");
  outTree_->Branch("rho"                  ,&rho_               ,"rho_/F");
  outTree_->Branch("met"                  ,&met_               ,"met_/F");
  outTree_->Branch("metSig"               ,&metSig_            ,"metSig_/F");

  outTree_->Branch("gen_eta"		,"vector<float>" , &gen_eta      	);
  outTree_->Branch("gen_phi"		,"vector<float>" , &gen_phi 		);
  outTree_->Branch("gen_p"		,"vector<float>",   &gen_p 		);
  outTree_->Branch("gen_px"		,"vector<float>",  &gen_px 		);
  outTree_->Branch("gen_py"		,"vector<float>",  &gen_py 		);
  outTree_->Branch("gen_pz"		,"vector<float>",  &gen_pz 		);
  outTree_->Branch("gen_pt"		,"vector<float>",  &gen_pt 		);
  outTree_->Branch("gen_energy"	    	,"vector<float>",  &gen_energy		);
  outTree_->Branch("gen_pdgId"	    	,"vector<int>",   &gen_pdgId  		);
  outTree_->Branch("gen_vx"	    	,"vector<float>", &gen_vx		);
  outTree_->Branch("gen_vy"	    	,"vector<float>", &gen_vy		);
  outTree_->Branch("gen_vz"	    	,"vector<float>", &gen_vz		);
  outTree_->Branch("gen_numDaught"  	,"vector<int>",    &gen_numDaught       );
  outTree_->Branch("gen_status"	    	,"vector<int>",  &gen_status     	);
  outTree_->Branch("gen_index"	    	,"vector<int>",  &gen_index      	);
  outTree_->Branch("gen_motherIndex"	,"vector<int>", &gen_motherIndex 	);

  outTree_->Branch("nJetsAK4"           ,&nJetsAK4_          ,"nJetsAK4_/I"		);
  outTree_->Branch("htAK4"              ,&htAK4_             ,"htAK4_/F"		);
  outTree_->Branch("mjjAK4"             ,&mjjAK4_            ,"mjjAK4_/F"		);
  outTree_->Branch("dEtajjAK4"          ,&dEtajjAK4_         ,"dEtajjAK4_/F"		);
  outTree_->Branch("dPhijjAK4"          ,&dPhijjAK4_         ,"dPhijjAK4_/F"		); 
  outTree_->Branch("nJetsAK8"           ,&nJetsAK8_          ,"nJetsAK8_/I"		);
  outTree_->Branch("htAK8"              ,&htAK8_             ,"htAK8_/F"		);
  outTree_->Branch("mjjAK8"             ,&mjjAK8_            ,"mjjAK8_/F"		);
  outTree_->Branch("dEtajjAK8"          ,&dEtajjAK8_         ,"dEtajjAK8_/F"		);
  outTree_->Branch("dPhijjAK8"          ,&dPhijjAK8_         ,"dPhijjAK8_/F"		); 
  // outTree_->Branch("nJetsCA8"           ,&nJetsCA8_          ,"nJetsCA8_/I"		);
  // outTree_->Branch("htCA8"           ,&htCA8_             ,"htCA8_/F"		);
  // outTree_->Branch("mjjCA8"          ,&mjjCA8_            ,"mjjCA8_/F"		);
  // outTree_->Branch("dEtajjCA8"       ,&dEtajjCA8_         ,"dEtajjCA8_/F"	);
  // outTree_->Branch("dPhijjCA8"       ,&dPhijjCA8_         ,"dPhijjCA8_/F"	); 

  //------------------------------------------------------------------

  outTree_->Branch("jetPtAK4"                ,"vector<float>"     ,&ptAK4_);
  outTree_->Branch("jetJecAK4"               ,"vector<float>"     ,&jecAK4_);
  outTree_->Branch("jetEtaAK4"               ,"vector<float>"     ,&etaAK4_);
  outTree_->Branch("jetPhiAK4"               ,"vector<float>"     ,&phiAK4_);
  outTree_->Branch("jetMassAK4"              ,"vector<float>"     ,&massAK4_);
  outTree_->Branch("jetEnergyAK4"            ,"vector<float>"     ,&energyAK4_);
  outTree_->Branch("jetAreaAK4"              ,"vector<float>"     ,&areaAK4_);
  outTree_->Branch("jetCSVAK4"               ,"vector<float>"     ,&csvAK4_);
  outTree_->Branch("pFlavourAK4"             ,"vector<int>"       ,&pFlavourAK4_);
  outTree_->Branch("hFlavourAK4"             ,"vector<int>"       ,&hFlavourAK4_);
  outTree_->Branch("nbHadAK4"                ,"vector<int>"       ,&nbHadAK4_);
  outTree_->Branch("ncHadAK4"                ,"vector<int>"       ,&ncHadAK4_);
  outTree_->Branch("jetChfAK4"               ,"vector<float>"     ,&chfAK4_);
  outTree_->Branch("jetNhfAK4"               ,"vector<float>"     ,&nhfAK4_);
  outTree_->Branch("jetPhfAK4"               ,"vector<float>"     ,&phfAK4_);
  outTree_->Branch("jetMufAK4"               ,"vector<float>"     ,&mufAK4_);
  outTree_->Branch("jetElfAK4"               ,"vector<float>"     ,&elfAK4_);
  outTree_->Branch("jetNemfAK4"              ,"vector<float>"     ,&nemfAK4_);
  outTree_->Branch("jetCemfAK4"              ,"vector<float>"     ,&cemfAK4_);
  outTree_->Branch("jetHf_hfAK4"             ,"vector<float>"     ,&hf_hfAK4_);
  outTree_->Branch("jetHf_emfAK4"            ,"vector<float>"    ,&hf_emfAK4_);
  outTree_->Branch("jetHofAK4"               ,"vector<float>"    ,&hofAK4_);
  outTree_->Branch("idLAK4"                  ,"vector<int>"      ,&idLAK4_);   
  outTree_->Branch("idTAK4"                  ,"vector<int>"      ,&idTAK4_);   
  outTree_->Branch("chHadMultAK4"          ,"vector<int>"      ,&chHadMultAK4_);   
  outTree_->Branch("chMultAK4"              ,"vector<int>"      ,&chMultAK4_);   
  outTree_->Branch("neHadMultAK4"           ,"vector<int>"      ,&neHadMultAK4_);   
  outTree_->Branch("neMultAK4"              ,"vector<int>"      ,&neMultAK4_);   
  outTree_->Branch("phoMultAK4"             ,"vector<int>"      ,&phoMultAK4_);   


  // Fill AK8 branches for all AK8 collections available
  TString nameAK8[nCollAK8] = {"miniaod", "toolbox"};

  for(u_int i=0 ; i<nCollAK8 ; i++) {

    outTree_->Branch("jetPtAK8_"+nameAK8[i]                ,"vector<float>"     ,&ptAK8_[i]);
    outTree_->Branch("jetJecAK8_"+nameAK8[i]               ,"vector<float>"     ,&jecAK8_[i]);
    outTree_->Branch("jetEtaAK8_"+nameAK8[i]               ,"vector<float>"     ,&etaAK8_[i]);
    outTree_->Branch("jetPhiAK8_"+nameAK8[i]               ,"vector<float>"     ,&phiAK8_[i]);
    outTree_->Branch("jetMassAK8_"+nameAK8[i]              ,"vector<float>"     ,&massAK8_[i]);
    outTree_->Branch("jetEnergyAK8_"+nameAK8[i]            ,"vector<float>"     ,&energyAK8_[i]);
    outTree_->Branch("jetAreaAK8_"+nameAK8[i]              ,"vector<float>"     ,&areaAK8_[i]);
    outTree_->Branch("jetCSVAK8_"+nameAK8[i]               ,"vector<float>"     ,&csvAK8_[i]);
    outTree_->Branch("pFlavourAK8_"+nameAK8[i]             ,"vector<int>"       ,&pFlavourAK8_[i]);
    outTree_->Branch("hFlavourAK8_"+nameAK8[i]             ,"vector<int>"       ,&hFlavourAK8_[i]);
    outTree_->Branch("nbHadAK8_"+nameAK8[i]                ,"vector<int>"       ,&nbHadAK8_[i]);
    outTree_->Branch("ncHadAK8_"+nameAK8[i]                ,"vector<int>"       ,&ncHadAK8_[i]);
    outTree_->Branch("jetChfAK8_"+nameAK8[i]               ,"vector<float>"     ,&chfAK8_[i]);
    outTree_->Branch("jetNhfAK8_"+nameAK8[i]               ,"vector<float>"     ,&nhfAK8_[i]);
    outTree_->Branch("jetPhfAK8_"+nameAK8[i]               ,"vector<float>"     ,&phfAK8_[i]);
    outTree_->Branch("jetMufAK8_"+nameAK8[i]               ,"vector<float>"     ,&mufAK8_[i]);
    outTree_->Branch("jetElfAK8_"+nameAK8[i]               ,"vector<float>"     ,&elfAK8_[i]); 
    outTree_->Branch("jetNemfAK8_"+nameAK8[i]              ,"vector<float>"     ,&nemfAK8_[i]);
    outTree_->Branch("jetCemfAK8_"+nameAK8[i]              ,"vector<float>"     ,&cemfAK8_[i]);
    outTree_->Branch("jetHf_hfAK8_"+nameAK8[i]             ,"vector<float>"     ,&hf_hfAK8_[i]);
    outTree_->Branch("jetHf_emfAK8_"+nameAK8[i]            ,"vector<float>"     ,&hf_emfAK8_[i]);
    outTree_->Branch("jetHofAK8_"+nameAK8[i]               ,"vector<float>"     ,&hofAK8_[i]);
    outTree_->Branch("idLAK8_"+nameAK8[i]                  ,"vector<int>"      ,&idLAK8_[i]);   
    outTree_->Branch("idTAK8_"+nameAK8[i]                  ,"vector<int>"      ,&idTAK8_[i]);   
    //
    // Sub-structure
    outTree_->Branch("jetMassPrunedAK8_"+nameAK8[i]        ,"vector<float>"     ,&massPrunedAK8_[i]);
    outTree_->Branch("jetMassSoftDropAK8_"+nameAK8[i]      ,"vector<float>"     ,&massSoftDropAK8_[i]);
    outTree_->Branch("jetTau1AK8_"+nameAK8[i]              ,"vector<float>"     ,&tau1AK8_[i]);
    outTree_->Branch("jetTau2AK8_"+nameAK8[i]              ,"vector<float>"     ,&tau2AK8_[i]);
    outTree_->Branch("jetTau3AK8_"+nameAK8[i]              ,"vector<float>"     ,&tau3AK8_[i]); 
    //
    // Associated PUPPI jets
    outTree_->Branch("jetPtAK8_Puppi_"+nameAK8[i]          , "vector<float>", &ptAK8_Puppi_[i]);
    outTree_->Branch("jetEtaAK8_Puppi_"+nameAK8[i]         , "vector<float>", &etaAK8_Puppi_[i]);
    outTree_->Branch("jetPhiAK8_Puppi_"+nameAK8[i]         , "vector<float>", &phiAK8_Puppi_[i]);
    outTree_->Branch("jetMassAK8_Puppi_"+nameAK8[i]        , "vector<float>", &massAK8_Puppi_[i]);
    outTree_->Branch("jetTau1AK8_Puppi_"+nameAK8[i]        , "vector<float>", &tau1AK8_Puppi_[i]);
    outTree_->Branch("jetTau2AK8_Puppi_"+nameAK8[i]        , "vector<float>", &tau2AK8_Puppi_[i]);
    outTree_->Branch("jetTau3AK8_Puppi_"+nameAK8[i]        , "vector<float>", &tau3AK8_Puppi_[i]);
    outTree_->Branch("jetMassSoftDropAK8_Puppi_"+nameAK8[i], "vector<float>", &massSoftDropAK8_Puppi_[i]);
    //
    //outTree_->Branch("jetDRAK8_"+nameAK8[i]                ,"vector<float>"     ,&dRAK8_[i]); 
    outTree_->Branch("chHadMultAK8_"+nameAK8[i]          ,"vector<int>"      ,&chHadMultAK8_[i]);   
    outTree_->Branch("chMultAK8_"+nameAK8[i]              ,"vector<int>"      ,&chMultAK8_[i]);   
    outTree_->Branch("neHadMultAK8_"+nameAK8[i]           ,"vector<int>"      ,&neHadMultAK8_[i]);   
    outTree_->Branch("neMultAK8_"+nameAK8[i]              ,"vector<int>"      ,&neMultAK8_[i]);   
    outTree_->Branch("phoMultAK8_"+nameAK8[i]             ,"vector<int>"      ,&phoMultAK8_[i]);   

  }
 
  //------------------------------------------------------------------
  outTree_->Branch("triggerResult","vector<bool>",&triggerResult_);

  //------------------------------------------------------------------
  outTree_->Branch("passFilterHBHE"                 ,&passFilterHBHE_                ,"passFilterHBHE_/O");
  outTree_->Branch("passFilterCSCHalo"              ,&passFilterCSCHalo_             ,"passFilterCSCHalo_/O");
  outTree_->Branch("passFilterHCALlaser"            ,&passFilterHCALlaser_           ,"passFilterHCALlaser_/O");
  outTree_->Branch("passFilterECALDeadCell"         ,&passFilterECALDeadCell_        ,"passFilterECALDeadCell_/O");
  outTree_->Branch("passFilterGoodVtx"              ,&passFilterGoodVtx_             ,"passFilterGoodVtx_/O");
  //outTree_->Branch("passFilterTrkFailure"           ,&passFilterTrkFailure_          ,"passFilterTrkFailure_/O");
  outTree_->Branch("passFilterEEBadSc"              ,&passFilterEEBadSc_             ,"passFilterEEBadSc_/O");
  outTree_->Branch("passFilterECALlaser"            ,&passFilterECALlaser_           ,"passFilterECALlaser_/O");
  outTree_->Branch("passFilterTrkPOG"               ,&passFilterTrkPOG_              ,"passFilterTrkPOG_/O");
  outTree_->Branch("passFilterTrkPOG_manystrip"     ,&passFilterTrkPOG_manystrip_    ,"passFilterTrkPOG_manystrip_/O");
  outTree_->Branch("passFilterTrkPOG_toomanystrip"  ,&passFilterTrkPOG_toomanystrip_ ,"passFilterTrkPOG_toomanystrip_/O");
  outTree_->Branch("passFilterTrkPOG_logError"      ,&passFilterTrkPOG_logError_     ,"passFilterTrkPOG_logError_/O");

  //------------------- MC --------------------------------- 
  outTree_->Branch("npu"                  ,"vector<float>"       , &npu_ );
  outTree_->Branch("PileupInteractions"   ,"vector<int>"       , &Number_interactions );
  outTree_->Branch("PileupOriginBX"       ,"vector<int>"       , &OriginBX );
  outTree_->Branch("ptHat"                ,&ptHat_             ,"ptHat_/F");
  outTree_->Branch("processID"            ,&processID_         ,"processID_/I");
  outTree_->Branch("weight"               ,&weight_            ,"weight_/F");

  outTree_->Branch("nGenJetsAK4"             ,&nGenJetsAK4_          ,"nGenJetsAK4_/I");
  outTree_->Branch("nGenJetsAK8"             ,&nGenJetsAK8_          ,"nGenJetsAK8_/I");

  outTree_->Branch("jetPtGenAK4"                ,"vector<float>"     ,&ptGenAK4_);
  outTree_->Branch("jetEtaGenAK4"               ,"vector<float>"     ,&etaGenAK4_);
  outTree_->Branch("jetPhiGenAK4"               ,"vector<float>"     ,&phiGenAK4_);
  outTree_->Branch("jetMassGenAK4"              ,"vector<float>"     ,&massGenAK4_);
  outTree_->Branch("jetEnergyGenAK4"            ,"vector<float>"     ,&energyGenAK4_);
  outTree_->Branch("jetPtGenAK8"                ,"vector<float>"     ,&ptGenAK8_);
  outTree_->Branch("jetEtaGenAK8"               ,"vector<float>"     ,&etaGenAK8_);
  outTree_->Branch("jetPhiGenAK8"               ,"vector<float>"     ,&phiGenAK8_);
  outTree_->Branch("jetMassGenAK8"              ,"vector<float>"     ,&massGenAK8_);
  outTree_->Branch("jetEnergyGenAK8"            ,"vector<float>"     ,&energyGenAK8_);

  return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////
void DijetTreeProducer::initialize()
{
  run_            = -999;
  evt_            = -999;
  lumi_           = -999;
  nVtx_           = -999;
  rho_            = -999;
  met_            = -999;
  metSig_         = -999;
  nJetsAK4_       = -999;
  htAK4_          = -999;
  mjjAK4_         = -999; 
  dEtajjAK4_      = -999; 
  dPhijjAK4_      = -999;

  ptAK4_             ->clear();
  etaAK4_            ->clear();
  phiAK4_            ->clear();
  massAK4_           ->clear();
  energyAK4_         ->clear();
  areaAK4_           ->clear();
  csvAK4_            ->clear();
  pFlavourAK4_       ->clear();
  hFlavourAK4_       ->clear();
  nbHadAK4_          ->clear();
  ncHadAK4_          ->clear();
  chfAK4_            ->clear();
  nhfAK4_            ->clear();
  phfAK4_            ->clear();
  elfAK4_            ->clear();
  mufAK4_            ->clear();
  nemfAK4_           ->clear();
  cemfAK4_           ->clear();
  hf_hfAK4_          ->clear();
  hf_emfAK4_         ->clear();
  hofAK4_            ->clear();
  jecAK4_            ->clear();
  jecAK4_            ->clear();
  idLAK4_            ->clear();
  idTAK4_            ->clear();
  // Juska's fix
  chHadMultAK4_     ->clear();
  chMultAK4_        ->clear();
  neHadMultAK4_     ->clear();
  neMultAK4_        ->clear();
  phoMultAK4_       ->clear();
  //massPrunedAK4_     ->clear();
  //tau1AK4_           ->clear();
  //tau2AK4_           ->clear();
  //dRAK4_             ->clear();
  //cutbasedJetId_      ->clear();
  //fullJetId_          ->clear();
  //fullJetDiscriminant_->clear();
  //ptAK4matchCaloJet_  ->clear();
  //emfAK4matchCaloJet_ ->clear(); 
  
  nJetsAK8_          = -999;
  htAK8_             = -999;
  mjjAK8_            = -999; 
  dEtajjAK8_         = -999; 
  dPhijjAK8_         = -999;
  //
  for(u_int i=0 ; i<nCollAK8 ; i++) {

    ptAK8_[i]             ->clear();
    etaAK8_[i]            ->clear();
    phiAK8_[i]            ->clear();
    massAK8_[i]           ->clear();
    energyAK8_[i]         ->clear();
    areaAK8_[i]           ->clear();
    csvAK8_[i]            ->clear();
    pFlavourAK8_[i]       ->clear();
    hFlavourAK8_[i]       ->clear();
    nbHadAK8_[i]          ->clear();
    ncHadAK8_[i]          ->clear();
    chfAK8_[i]            ->clear();
    nhfAK8_[i]            ->clear();
    phfAK8_[i]            ->clear();
    elfAK8_[i]            ->clear();
    mufAK8_[i]            ->clear();
    nemfAK8_[i]           ->clear();
    cemfAK8_[i]           ->clear();
    hf_hfAK8_[i]          ->clear();
    hf_emfAK8_[i]         ->clear();
    hofAK8_[i]            ->clear();
    jecAK8_[i]            ->clear();
    jecAK8_[i]            ->clear();
    idLAK8_[i]            ->clear();
    idTAK8_[i]            ->clear();
    massPrunedAK8_[i]     ->clear();
    massSoftDropAK8_[i]   ->clear();
    tau1AK8_[i]           ->clear();
    tau2AK8_[i]           ->clear();
    tau3AK8_[i]           ->clear();

    // Juska's fix
    chHadMultAK8_[i]     ->clear();
    chMultAK8_[i]        ->clear();
    neHadMultAK8_[i]     ->clear();
    neMultAK8_[i]        ->clear();
    phoMultAK8_[i]        ->clear();
    //dRAK8_[i]             ->clear();

    ptAK8_Puppi_[i] ->clear();   
    etaAK8_Puppi_[i] ->clear();  
    phiAK8_Puppi_[i] ->clear();  
    massAK8_Puppi_[i] ->clear(); 
    tau1AK8_Puppi_[i] ->clear(); 
    tau2AK8_Puppi_[i] ->clear(); 
    tau3AK8_Puppi_[i] ->clear(); 
    massSoftDropAK8_Puppi_[i] ->clear(); 

  }
  
  triggerResult_     ->clear();
  
  passFilterHBHE_                  = false;
  passFilterCSCHalo_               = false;
  passFilterHCALlaser_             = false;
  passFilterECALDeadCell_          = false;
  passFilterGoodVtx_               = false;
  //passFilterTrkFailure_            = false;
  passFilterEEBadSc_               = false;
  passFilterECALlaser_             = false;
  passFilterTrkPOG_                = false;
  passFilterTrkPOG_manystrip_      = false;
  passFilterTrkPOG_toomanystrip_   = false;
  passFilterTrkPOG_logError_       = false;

  //----- MC -------
  npu_ ->clear();
  Number_interactions ->clear();
  OriginBX            -> clear();
  
  ptHat_ = -999; 
  processID_ = -999; 
  weight_ = -999;

  nGenJetsAK4_ = -999;
  nGenJetsAK8_ = -999;
  
  ptGenAK4_    ->clear();
  phiGenAK4_   ->clear();
  etaGenAK4_   ->clear();
  massGenAK4_  ->clear();
  energyGenAK4_->clear();
  ptGenAK8_    ->clear();
  phiGenAK8_   ->clear();
  etaGenAK8_   ->clear();
  massGenAK8_  ->clear();
  energyGenAK8_->clear();
  
  gen_eta		->clear();
  gen_phi		->clear();
  gen_p		        ->clear();
  gen_px		->clear();
  gen_py		->clear();
  gen_pz		->clear();
  gen_pt		->clear();
  gen_energy    	->clear();
  gen_pdgId	        ->clear();
  gen_vx		->clear();
  gen_vy		->clear();
  gen_vz		->clear();
  gen_numDaught	        ->clear();
  gen_status	        ->clear();
  gen_index   	        ->clear();
  gen_motherIndex       ->clear();  
  
}

//////////////////////////////////////////////////////////////////////////////////////////
void DijetTreeProducer::endJob() 
{  
  delete triggerResult_;


  delete gen_eta	;
  delete gen_phi	;
  delete gen_p		;
  delete gen_px	;
  delete gen_py	;
  delete gen_pz	;
  delete gen_pt	;
  delete gen_energy    ;
  delete gen_pdgId	;
  delete gen_vx	;
  delete gen_vy	;
  delete gen_vz	;
  delete gen_numDaught	;
  delete gen_status	;
  delete gen_index   	;
  delete gen_motherIndex;

  delete ptAK4_;
  delete jecAK4_;
  delete etaAK4_;
  delete phiAK4_;
  delete massAK4_;
  delete energyAK4_;
  delete areaAK4_;
  delete csvAK4_;
  delete pFlavourAK4_;
  delete hFlavourAK4_;
  delete nbHadAK4_;
  delete ncHadAK4_;
  delete chfAK4_;
  delete nhfAK4_;
  delete phfAK4_;
  delete mufAK4_;
  delete elfAK4_;
  delete nemfAK4_;
  delete cemfAK4_;
  delete hf_hfAK4_;
  delete hf_emfAK4_;
  delete hofAK4_;
  delete idLAK4_;
  delete idTAK4_;
  delete chHadMultAK4_ ;
  delete chMultAK4_    ;
  delete neHadMultAK4_ ;
  delete neMultAK4_    ;
  delete phoMultAK4_   ;

  for(u_int i=0 ; i<nCollAK8 ; i++) {

    delete ptAK8_[i];
    delete jecAK8_[i];
    delete etaAK8_[i];
    delete phiAK8_[i];
    delete massAK8_[i];
    delete energyAK8_[i];
    delete areaAK8_[i];
    delete csvAK8_[i];
    delete pFlavourAK8_[i];
    delete hFlavourAK8_[i];
    delete nbHadAK8_[i];
    delete ncHadAK8_[i];
    delete chfAK8_[i];
    delete nhfAK8_[i];
    delete phfAK8_[i];
    delete mufAK8_[i];
    delete elfAK8_[i];
    delete nemfAK8_[i];
    delete cemfAK8_[i];
    delete hf_hfAK8_[i];
    delete hf_emfAK8_[i];
    delete hofAK8_[i];
    delete idLAK8_[i];
    delete idTAK8_[i];
    delete massPrunedAK8_[i];
    delete massSoftDropAK8_[i];
    delete tau1AK8_[i];
    delete tau2AK8_[i];
    delete tau3AK8_[i];
    delete chHadMultAK8_[i];
    delete chMultAK8_[i]    ;
    delete neHadMultAK8_[i] ;
    delete neMultAK8_[i]    ;
    delete phoMultAK8_[i]   ;

    delete ptAK8_Puppi_[i] ;   
    delete etaAK8_Puppi_[i] ;  
    delete phiAK8_Puppi_[i] ;  
    delete massAK8_Puppi_[i] ; 
    delete tau1AK8_Puppi_[i] ; 
    delete tau2AK8_Puppi_[i] ; 
    delete tau3AK8_Puppi_[i] ; 
    delete massSoftDropAK8_Puppi_[i] ; 

  }

  for(unsigned i=0;i<vtriggerSelector_.size();i++) {
    delete vtriggerSelector_[i];
  }

  delete HBHENoiseFilter_Selector_;
  delete CSCHaloNoiseFilter_Selector_;
  delete HCALlaserNoiseFilter_Selector_;
  delete ECALDeadCellNoiseFilter_Selector_;
  delete GoodVtxNoiseFilter_Selector_;
  //delete TrkFailureNoiseFilter_Selector_;
  delete EEBadScNoiseFilter_Selector_;
  delete ECALlaserNoiseFilter_Selector_;
  delete TrkPOGNoiseFilter_Selector_;
  delete TrkPOG_manystrip_NoiseFilter_Selector_;
  delete TrkPOG_toomanystrip_NoiseFilter_Selector_;
  delete TrkPOG_logError_NoiseFilter_Selector_;

}

template<typename T> float DijetTreeProducer::GetUserFloat(const pat::PATObject< T > &obj,  const TString tag)
{
  float result = -8888888;
  if(obj.hasUserFloat(tag)) {
    result = obj.userFloat(tag);
  }
  return result;
}

DEFINE_FWK_MODULE(DijetTreeProducer);
