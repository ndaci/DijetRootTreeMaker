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
  
  // Migrate to Consumes-system. Skip Calo-stuff
  srcJetsAK4_ = (consumes<pat::JetCollection>(    cfg.getParameter<InputTag>("jetsAK4")));
  srcJetsAK8_ = (consumes<pat::JetCollection>(    cfg.getParameter<InputTag>("jetsAK8")));
  
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
  TrkFailureNoiseFilter_Selector_           = triggerExpression::parse( cfg.getParameter<std::string> ("noiseFilterSelection_trackingFailureFilter") );
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
  
  //--- book the tree -----------------------
  outTree_ = fs_->make<TTree>("events","events");
  outTree_->Branch("runNo"                ,&run_               ,"run_/I");
  outTree_->Branch("evtNo"                ,&evt_               ,"evt_/I");
  outTree_->Branch("lumi"                 ,&lumi_              ,"lumi_/I");
  outTree_->Branch("nvtx"                 ,&nVtx_              ,"nVtx_/I");
  outTree_->Branch("rho"                  ,&rho_               ,"rho_/F");
  outTree_->Branch("met"                  ,&met_               ,"met_/F");
  outTree_->Branch("metSig"               ,&metSig_            ,"metSig_/F");

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

  ptAK8_             = new std::vector<float>;
  jecAK8_            = new std::vector<float>;
  etaAK8_            = new std::vector<float>;
  phiAK8_            = new std::vector<float>;
  massAK8_           = new std::vector<float>;
  energyAK8_         = new std::vector<float>;
  areaAK8_           = new std::vector<float>;
  csvAK8_            = new std::vector<float>;
  pFlavourAK8_       = new std::vector<int>;
  hFlavourAK8_       = new std::vector<int>;
  nbHadAK8_          = new std::vector<int>;
  ncHadAK8_          = new std::vector<int>;
  chfAK8_            = new std::vector<float>;
  nhfAK8_            = new std::vector<float>;
  phfAK8_            = new std::vector<float>;
  mufAK8_            = new std::vector<float>;
  elfAK8_            = new std::vector<float>;
  nemfAK8_           = new std::vector<float>;
  cemfAK8_           = new std::vector<float>;
  // Hadronic forward hadrons
  hf_hfAK8_          = new std::vector<float>;
  // Hadronic forward photons
  hf_emfAK8_         = new std::vector<float>;
  hofAK8_            = new std::vector<float>;
  idLAK8_            = new std::vector<int>;
  idTAK8_            = new std::vector<int>;
  massPrunedAK8_     = new std::vector<float>;
  massSoftDropAK8_   = new std::vector<float>;
  tau1AK8_           = new std::vector<float>;
  tau2AK8_           = new std::vector<float>;
  tau3AK8_           = new std::vector<float>;
  chHadMultAK8_      = new std::vector<int>;   
  chMultAK8_         = new std::vector<int>;
  neHadMultAK8_      = new std::vector<int>; 
  neMultAK8_         = new std::vector<int>;
  phoMultAK8_        = new std::vector<int>;
 
  //dRAK8_             = new std::vector<float>;
  outTree_->Branch("jetPtAK8"                ,"vector<float>"     ,&ptAK8_);
  outTree_->Branch("jetJecAK8"               ,"vector<float>"     ,&jecAK8_);
  outTree_->Branch("jetEtaAK8"               ,"vector<float>"     ,&etaAK8_);
  outTree_->Branch("jetPhiAK8"               ,"vector<float>"     ,&phiAK8_);
  outTree_->Branch("jetMassAK8"              ,"vector<float>"     ,&massAK8_);
  outTree_->Branch("jetEnergyAK8"            ,"vector<float>"     ,&energyAK8_);
  outTree_->Branch("jetAreaAK8"              ,"vector<float>"     ,&areaAK8_);
  outTree_->Branch("jetCSVAK8"               ,"vector<float>"     ,&csvAK8_);
  outTree_->Branch("pFlavourAK8"             ,"vector<int>"       ,&pFlavourAK8_);
  outTree_->Branch("hFlavourAK8"             ,"vector<int>"       ,&hFlavourAK8_);
  outTree_->Branch("nbHadAK8"                ,"vector<int>"       ,&nbHadAK8_);
  outTree_->Branch("ncHadAK8"                ,"vector<int>"       ,&ncHadAK8_);
  outTree_->Branch("jetChfAK8"               ,"vector<float>"     ,&chfAK8_);
  outTree_->Branch("jetNhfAK8"               ,"vector<float>"     ,&nhfAK8_);
  outTree_->Branch("jetPhfAK8"               ,"vector<float>"     ,&phfAK8_);
  outTree_->Branch("jetMufAK8"               ,"vector<float>"     ,&mufAK8_);
  outTree_->Branch("jetElfAK8"               ,"vector<float>"     ,&elfAK8_); 
  outTree_->Branch("jetNemfAK8"              ,"vector<float>"     ,&nemfAK8_);
  outTree_->Branch("jetCemfAK8"              ,"vector<float>"     ,&cemfAK8_);
  outTree_->Branch("jetHf_hfAK8"             ,"vector<float>"     ,&hf_hfAK8_);
  outTree_->Branch("jetHf_emfAK8"            ,"vector<float>"     ,&hf_emfAK8_);
  outTree_->Branch("jetHofAK8"               ,"vector<float>"     ,&hofAK8_);
  outTree_->Branch("idLAK8"                  ,"vector<int>"      ,&idLAK8_);   
  outTree_->Branch("idTAK8"                  ,"vector<int>"      ,&idTAK8_);   
  outTree_->Branch("jetMassPrunedAK8"        ,"vector<float>"     ,&massPrunedAK8_);
  outTree_->Branch("jetMassSoftDropAK8"      ,"vector<float>"     ,&massSoftDropAK8_);
  outTree_->Branch("jetTau1AK8"              ,"vector<float>"     ,&tau1AK8_);
  outTree_->Branch("jetTau2AK8"              ,"vector<float>"     ,&tau2AK8_);
  outTree_->Branch("jetTau3AK8"              ,"vector<float>"     ,&tau3AK8_); 
  //outTree_->Branch("jetDRAK8"                ,"vector<float>"     ,&dRAK8_); 
  outTree_->Branch("chHadMultAK8"          ,"vector<int>"      ,&chHadMultAK8_);   
  outTree_->Branch("chMultAK8"              ,"vector<int>"      ,&chMultAK8_);   
  outTree_->Branch("neHadMultAK8"           ,"vector<int>"      ,&neHadMultAK8_);   
  outTree_->Branch("neMultAK8"              ,"vector<int>"      ,&neMultAK8_);   
  outTree_->Branch("phoMultAK8"             ,"vector<int>"      ,&phoMultAK8_);   
 
  //------------------------------------------------------------------
  triggerResult_ = new std::vector<bool>;
  outTree_->Branch("triggerResult","vector<bool>",&triggerResult_);

  //------------------------------------------------------------------
  outTree_->Branch("passFilterHBHE"                 ,&passFilterHBHE_                ,"passFilterHBHE_/O");
  outTree_->Branch("passFilterCSCHalo"              ,&passFilterCSCHalo_             ,"passFilterCSCHalo_/O");
  outTree_->Branch("passFilterHCALlaser"            ,&passFilterHCALlaser_           ,"passFilterHCALlaser_/O");
  outTree_->Branch("passFilterECALDeadCell"         ,&passFilterECALDeadCell_        ,"passFilterECALDeadCell_/O");
  outTree_->Branch("passFilterGoodVtx"              ,&passFilterGoodVtx_             ,"passFilterGoodVtx_/O");
  outTree_->Branch("passFilterTrkFailure"           ,&passFilterTrkFailure_          ,"passFilterTrkFailure_/O");
  outTree_->Branch("passFilterEEBadSc"              ,&passFilterEEBadSc_             ,"passFilterEEBadSc_/O");
  outTree_->Branch("passFilterECALlaser"            ,&passFilterECALlaser_           ,"passFilterECALlaser_/O");
  outTree_->Branch("passFilterTrkPOG"               ,&passFilterTrkPOG_              ,"passFilterTrkPOG_/O");
  outTree_->Branch("passFilterTrkPOG_manystrip"     ,&passFilterTrkPOG_manystrip_    ,"passFilterTrkPOG_manystrip_/O");
  outTree_->Branch("passFilterTrkPOG_toomanystrip"  ,&passFilterTrkPOG_toomanystrip_ ,"passFilterTrkPOG_toomanystrip_/O");
  outTree_->Branch("passFilterTrkPOG_logError"      ,&passFilterTrkPOG_logError_     ,"passFilterTrkPOG_logError_/O");


  //------------------- MC ---------------------------------
  npu_                = new std::vector<float>;  
  Number_interactions = new std::vector<int>;
  OriginBX            = new std::vector<int>;
 
  outTree_->Branch("npu"                  ,"vector<float>"       , &npu_ );
  outTree_->Branch("PileupInteractions"   ,"vector<int>"       , &Number_interactions );
  outTree_->Branch("PileupOriginBX"       ,"vector<int>"       , &OriginBX );
  outTree_->Branch("ptHat"                ,&ptHat_             ,"ptHat_/F");
  outTree_->Branch("processID"            ,&processID_         ,"processID_/I");
  outTree_->Branch("weight"               ,&weight_            ,"weight_/F");

  outTree_->Branch("nGenJetsAK4"             ,&nGenJetsAK4_          ,"nGenJetsAK4_/I");
  outTree_->Branch("nGenJetsAK8"             ,&nGenJetsAK8_          ,"nGenJetsAK8_/I");

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

  delete ptAK8_;
  delete jecAK8_;
  delete etaAK8_;
  delete phiAK8_;
  delete massAK8_;
  delete energyAK8_;
  delete areaAK8_;
  delete csvAK8_;
  delete pFlavourAK8_;
  delete hFlavourAK8_;
  delete nbHadAK8_;
  delete ncHadAK8_;
  delete chfAK8_;
  delete nhfAK8_;
  delete phfAK8_;
  delete mufAK8_;
  delete elfAK8_;
  delete nemfAK8_;
  delete cemfAK8_;
  delete hf_hfAK8_;
  delete hf_emfAK8_;
  delete hofAK8_;
  delete idLAK8_;
  delete idTAK8_;
  delete massPrunedAK8_;
  delete massSoftDropAK8_;
  delete tau1AK8_;
  delete tau2AK8_;
  delete tau3AK8_;
  delete chHadMultAK8_;
  delete chMultAK8_    ;
  delete neHadMultAK8_ ;
  delete neMultAK8_    ;
  delete phoMultAK8_   ;

  for(unsigned i=0;i<vtriggerSelector_.size();i++) {
    delete vtriggerSelector_[i];
  }

  delete HBHENoiseFilter_Selector_;
  delete CSCHaloNoiseFilter_Selector_;
  delete HCALlaserNoiseFilter_Selector_;
  delete ECALDeadCellNoiseFilter_Selector_;
  delete GoodVtxNoiseFilter_Selector_;
  delete TrkFailureNoiseFilter_Selector_;
  delete EEBadScNoiseFilter_Selector_;
  delete ECALlaserNoiseFilter_Selector_;
  delete TrkPOGNoiseFilter_Selector_;
  delete TrkPOG_manystrip_NoiseFilter_Selector_;
  delete TrkPOG_toomanystrip_NoiseFilter_Selector_;
  delete TrkPOG_logError_NoiseFilter_Selector_;

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
      TrkFailureNoiseFilter_Selector_          -> init(noiseFilterCache_);
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
    passFilterTrkFailure_          = (*TrkFailureNoiseFilter_Selector_)          (noiseFilterCache_);    
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

      //tau1AK4_          ->push_back(ijet->userFloat("NjettinessAK4:tau1"));
      //tau2AK4_          ->push_back(ijet->userFloat("NjettinessAK4:tau2"));
      //cutbasedJetId_      ->push_back(ijet->userInt("pileupJetIdEvaluator:cutbasedId"));
      //fullJetId_          ->push_back(ijet->userFloat("pileupJetIdEvaluator:fullDiscriminant"));
      //fullJetDiscriminant_->push_back(ijet->userInt("pileupJetIdEvaluator:fullId"));

    }

  }// jet loop  

  htAK4_     = htAK4;
  if (nJetsAK4_ > 1) { //assuming jets are ordered by pt in the pat collection
    mjjAK4_    = (vP4AK4[0]+vP4AK4[1]).M();
    dEtajjAK4_ = fabs((*etaAK4_)[0]-(*etaAK4_)[1]); 
    dPhijjAK4_ = fabs(deltaPhi((*phiAK4_)[0],(*phiAK4_)[1]));
  }

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
    int phoMult = ijet->photonMultiplicity();
      
    // Juska's added fractions for identical JetID with recommendations
    double nemf = ijet->neutralEmEnergyFraction();
    double cemf = ijet->chargedEmEnergyFraction();
    int NumConst = npr;

    float eta  = ijet->eta(); // removed fabs() -Juska
    float pt   = ijet->correctedJet(0).pt()*jecFactorsAK8.at(*i); // Is this OK? Correct corrected? -Juska

    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID
    int idL = (nhf<0.99 && nemf<0.99 && NumConst>1 && muf < 0.8) && ((fabs(eta) <= 2.4 && chf>0 && chMult>0 && cemf<0.99) || fabs(eta)>2.4);
    int idT = (nhf<0.90 && nemf<0.90 && NumConst>1 && muf<0.8) && ((fabs(eta)<=2.4 && chf>0 && chMult>0 && cemf<0.90) || fabs(eta)>2.4);
      
      
    if (pt > ptMinAK8_) {
      htAK8 += pt;
      nJetsAK8_++;

      vP4AK8.push_back(TLorentzVector(ijet->correctedJet(0).px()*jecFactorsAK8.at(*i),ijet->correctedJet(0).py()*jecFactorsAK8.at(*i),ijet->correctedJet(0).pz()*jecFactorsAK8.at(*i),ijet->correctedJet(0).energy()*jecFactorsAK8.at(*i)));
      chfAK8_           ->push_back(chf);
      nhfAK8_           ->push_back(nhf);
      phfAK8_           ->push_back(phf);
      elfAK8_           ->push_back(elf);
      mufAK8_           ->push_back(muf);
      nemfAK8_          ->push_back(nemf);
      cemfAK8_          ->push_back(cemf);
      hf_hfAK8_         ->push_back(hf_hf);
      hf_emfAK8_        ->push_back(hf_emf);
      hofAK8_           ->push_back(hof);
      jecAK8_           ->push_back(jecFactorsAK8.at(*i));
      ptAK8_            ->push_back(pt);
      phiAK8_           ->push_back(ijet->phi());
      etaAK8_           ->push_back(ijet->eta());
      massAK8_          ->push_back(ijet->correctedJet(0).mass()*jecFactorsAK8.at(*i));
      energyAK8_        ->push_back(ijet->correctedJet(0).energy()*jecFactorsAK8.at(*i));
      areaAK8_          ->push_back(ijet->jetArea());
      csvAK8_           ->push_back(ijet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
      pFlavourAK8_      ->push_back(ijet->partonFlavour());
      hFlavourAK8_      ->push_back(ijet->hadronFlavour());
      nbHadAK8_         ->push_back(ijet->jetFlavourInfo().getbHadrons().size());
      ncHadAK8_         ->push_back(ijet->jetFlavourInfo().getcHadrons().size());
      idLAK8_           ->push_back(idL);
      idTAK8_           ->push_back(idT);
      
      // Disable, causes crash. See the interesting error below. Juska
      /*
	Exception Message:
	Requested UserFloat NjettinessAK8:tau1 is not available! Possible UserFloats are: 
	ak8PFJetsCHSPrunedMass ak8PFJetsCHSSoftDropMass ak8PFJetsCHSTrimmedMass ak8PFJetsCHSFilteredMass NjettinessAK8:tau1 NjettinessAK8:tau2 NjettinessAK8:tau3 
	tau1AK8_          ->push_back(ijet->userFloat("NjettinessAK8:tau1"));
	tau2AK8_          ->push_back(ijet->userFloat("NjettinessAK8:tau2"));
	tau3AK8_          ->push_back(ijet->userFloat("NjettinessAK8:tau3"));
	massPrunedAK8_    ->push_back(ijet->userFloat("ak8PFJetsCHSPrunedMass"));
	massSoftDropAK8_  ->push_back(ijet->userFloat("ak8PFJetsCHSSoftDropMass"));
      */
      
      chHadMultAK8_     ->push_back(chHadMult);
      chMultAK8_        ->push_back(chMult);
      neHadMultAK8_     ->push_back(neHadMult);  
      neMultAK8_        ->push_back(neMult);
      phoMultAK8_       ->push_back(phoMult); 
	
	
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
    dEtajjAK8_ = fabs((*etaAK8_)[0]-(*etaAK8_)[1]); 
    dPhijjAK8_ = fabs(deltaPhi((*phiAK8_)[0],(*phiAK8_)[1]));
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
  nJetsAK4_          = -999;
  htAK4_             = -999;
  mjjAK4_            = -999; 
  dEtajjAK4_         = -999; 
  dPhijjAK4_         = -999;
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
  hf_hfAK4_             ->clear();
  hf_emfAK4_            ->clear();
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
  phoMultAK4_        ->clear();
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
  ptAK8_             ->clear();
  etaAK8_            ->clear();
  phiAK8_            ->clear();
  massAK8_           ->clear();
  energyAK8_         ->clear();
  areaAK8_           ->clear();
  csvAK8_            ->clear();
  pFlavourAK8_       ->clear();
  hFlavourAK8_       ->clear();
  nbHadAK8_          ->clear();
  ncHadAK8_          ->clear();
  chfAK8_            ->clear();
  nhfAK8_            ->clear();
  phfAK8_            ->clear();
  elfAK8_            ->clear();
  mufAK8_            ->clear();
  nemfAK8_           ->clear();
  cemfAK8_           ->clear();
  hf_hfAK8_          ->clear();
  hf_emfAK8_         ->clear();
  hofAK8_            ->clear();
  jecAK8_            ->clear();
  jecAK8_            ->clear();
  idLAK8_            ->clear();
  idTAK8_            ->clear();
  massPrunedAK8_     ->clear();
  massSoftDropAK8_   ->clear();
  tau1AK8_           ->clear();
  tau2AK8_           ->clear();
  tau3AK8_           ->clear();
  // Juska's fix
  chHadMultAK8_     ->clear();
  chMultAK8_        ->clear();
  neHadMultAK8_     ->clear();
  neMultAK8_        ->clear();
  phoMultAK8_        ->clear();
  //dRAK8_             ->clear();
  
  triggerResult_     ->clear();
  
  passFilterHBHE_                  = false;
  passFilterCSCHalo_               = false;
  passFilterHCALlaser_             = false;
  passFilterECALDeadCell_          = false;
  passFilterGoodVtx_               = false;
  passFilterTrkFailure_            = false;
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
DijetTreeProducer::~DijetTreeProducer() 
{
}

DEFINE_FWK_MODULE(DijetTreeProducer);
