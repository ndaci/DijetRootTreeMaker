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
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

using namespace std;
using namespace reco;

DijetTreeProducer::DijetTreeProducer(edm::ParameterSet const& cfg) 
{
  srcJetsAK4_         = cfg.getParameter<edm::InputTag>             ("jetsAK4");
  srcJetsAK8_         = cfg.getParameter<edm::InputTag>             ("jetsAK8");
  //srcJetsCA8_         = cfg.getParameter<edm::InputTag>             ("jetsCA8");
  srcGenJetsAK4_      = cfg.getParameter<edm::InputTag>             ("genJetsAK4");
  srcGenJetsAK8_      = cfg.getParameter<edm::InputTag>             ("genJetsAK8");
  //srcGenJetsCA8_      = cfg.getParameter<edm::InputTag>             ("genJetsCA8");
  srcRho_             = cfg.getParameter<edm::InputTag>             ("rho");
  srcMET_             = cfg.getParameter<edm::InputTag>             ("met");
  srcVrtx_            = cfg.getParameter<edm::InputTag>             ("vtx");
  srcPU_              = cfg.getUntrackedParameter<edm::InputTag>    ("pu",edm::InputTag(""));
  srcGenInfo_           = cfg.getUntrackedParameter<edm::InputTag>  ("ptHat",edm::InputTag());
  srcPrunedGenParticles_ = cfg.getParameter<edm::InputTag>          ("genParticles");

  ptMinAK4_           = cfg.getParameter<double>                    ("ptMinAK4");
  ptMinAK8_           = cfg.getParameter<double>                    ("ptMinAK8");
  //ptMinCA8_           = cfg.getParameter<double>                    ("ptMinCA8");
  //mjjMin_             = cfg.getParameter<double>                    ("mjjMin");
  //dEtaMax_            = cfg.getParameter<double>                    ("dEtaMax");
  triggerCache_       = triggerExpression::Data(cfg.getParameterSet("triggerConfiguration"),consumesCollector());
  vtriggerAlias_      = cfg.getParameter<std::vector<std::string> > ("triggerAlias");
  vtriggerSelection_  = cfg.getParameter<std::vector<std::string> > ("triggerSelection");

  if (vtriggerAlias_.size() != vtriggerSelection_.size()) {
    cout<<"ERROR: the number of trigger aliases does not match the number of trigger names !!!"<<endl;
    return;
  }

  for(unsigned i=0;i<vtriggerSelection_.size();i++) {
    vtriggerSelector_.push_back(triggerExpression::parse(vtriggerSelection_[i]));
  }
  
  // For JECs
  redoJECs_ = cfg.getParameter<bool>("redoJECs");
  // AK4
  L1corrAK4_ = cfg.getParameter<edm::FileInPath>("L1corrAK4");
  L2corrAK4_ = cfg.getParameter<edm::FileInPath>("L2corrAK4");
  L3corrAK4_ = cfg.getParameter<edm::FileInPath>("L3corrAK4");
  // AK8
  L1corrAK8_ = cfg.getParameter<edm::FileInPath>("L1corrAK8");
  L2corrAK8_ = cfg.getParameter<edm::FileInPath>("L2corrAK8");
  L3corrAK8_ = cfg.getParameter<edm::FileInPath>("L3corrAK8");

  if(redoJECs_)
  {
    // AK4
    L1ParAK4 = new JetCorrectorParameters(L1corrAK4_.fullPath());
    L2ParAK4 = new JetCorrectorParameters(L2corrAK4_.fullPath());
    L3ParAK4 = new JetCorrectorParameters(L3corrAK4_.fullPath());

    std::vector<JetCorrectorParameters> vParAK4;
    vParAK4.push_back(*L1ParAK4);
    vParAK4.push_back(*L2ParAK4);
    vParAK4.push_back(*L3ParAK4);

    JetCorrectorAK4 = new FactorizedJetCorrector(vParAK4);

    // AK8
    L1ParAK8 = new JetCorrectorParameters(L1corrAK8_.fullPath());
    L2ParAK8 = new JetCorrectorParameters(L2corrAK8_.fullPath());
    L3ParAK8 = new JetCorrectorParameters(L3corrAK8_.fullPath());

    std::vector<JetCorrectorParameters> vParAK8;
    vParAK8.push_back(*L1ParAK8);
    vParAK8.push_back(*L2ParAK8);
    vParAK8.push_back(*L3ParAK8);

    JetCorrectorAK8 = new FactorizedJetCorrector(vParAK8);
  }
}
//////////////////////////////////////////////////////////////////////////////////////////
void DijetTreeProducer::beginJob() 
{
  //--- book the trigger histograms ---------
  triggerNamesHisto_ = fs_->make<TH1F>("TriggerNames","TriggerNames",1,0,1);
  triggerNamesHisto_->SetBit(TH1::kCanRebin);
  for(unsigned i=0;i<vtriggerSelection_.size();i++) {
    triggerNamesHisto_->Fill(vtriggerSelection_[i].c_str(),1);
  }
  triggerPassHisto_ = fs_->make<TH1F>("TriggerPass","TriggerPass",1,0,1);
  triggerPassHisto_->SetBit(TH1::kCanRebin);
  
  //--- book the tree -----------------------
  outTree_ = fs_->make<TTree>("events","events");
  outTree_->Branch("runNo"                ,&run_               ,"run_/I");
  outTree_->Branch("evtNo"                ,&evt_               ,"evt_/I");
  outTree_->Branch("lumi"                 ,&lumi_              ,"lumi_/I");
  outTree_->Branch("nvtx"                 ,&nVtx_              ,"nVtx_/I");
  outTree_->Branch("rho"                  ,&rho_               ,"rho_/F");
  outTree_->Branch("met"                  ,&met_               ,"met_/F");
  outTree_->Branch("metSig"               ,&metSig_            ,"metSig_/F");

  gen_eta           = new std::vector<float>;
  gen_phi           = new std::vector<float>;
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
  chfAK4_            = new std::vector<float>;
  nhfAK4_            = new std::vector<float>;
  phfAK4_            = new std::vector<float>;
  mufAK4_            = new std::vector<float>;
  elfAK4_            = new std::vector<float>;
  // Hadronic forward hadrons
  hf_hfAK4_          = new std::vector<float>;
  // Hadronic forward electromagnetic fraction
  hf_emfAK4_         = new std::vector<float>;
  idLAK4_            = new std::vector<int>;
  idTAK4_            = new std::vector<int>;
  chHadMultAK4_     = new std::vector<int>;
  chMultAK4_         = new std::vector<int>;
  neHadMultAK4_      = new std::vector<int>;
  neMultAK4_         = new std::vector<int>;
  phoMultAK4_        = new std::vector<int>;
  
  
  //massPrunedAK4_     = new std::vector<float>;
  //tau1AK4_           = new std::vector<float>;
  //tau2AK4_           = new std::vector<float>;
  //dRAK4_             = new std::vector<float>;

  //cutbasedJetId_       = new std::vector<float>;
  //fullJetId_           = new std::vector<float>;
  //fullJetDiscriminant_ = new std::vector<float>;

  outTree_->Branch("jetPtAK4"                ,"vector<float>"     ,&ptAK4_);
  outTree_->Branch("jetJecAK4"               ,"vector<float>"     ,&jecAK4_);
  outTree_->Branch("jetEtaAK4"               ,"vector<float>"     ,&etaAK4_);
  outTree_->Branch("jetPhiAK4"               ,"vector<float>"     ,&phiAK4_);
  outTree_->Branch("jetMassAK4"              ,"vector<float>"     ,&massAK4_);
  outTree_->Branch("jetEnergyAK4"            ,"vector<float>"     ,&energyAK4_);
  outTree_->Branch("jetAreaAK4"              ,"vector<float>"     ,&areaAK4_);
  outTree_->Branch("jetChfAK4"               ,"vector<float>"     ,&chfAK4_);
  outTree_->Branch("jetNhfAK4"               ,"vector<float>"     ,&nhfAK4_);
  outTree_->Branch("jetPhfAK4"               ,"vector<float>"     ,&phfAK4_);
  outTree_->Branch("jetMufAK4"               ,"vector<float>"     ,&mufAK4_);
  outTree_->Branch("jetElfAK4"               ,"vector<float>"     ,&elfAK4_);
  outTree_->Branch("jetHf_hfAK4"             ,"vector<float>"     ,&hf_hfAK4_);
  outTree_->Branch("jetHf_emfAK4"            ,"vector<float>"    ,&hf_emfAK4_);   
  outTree_->Branch("idLAK4"                  ,"vector<int>"      ,&idLAK4_);   
  outTree_->Branch("idTAK4"                  ,"vector<int>"      ,&idTAK4_);   
  outTree_->Branch("chHadMultAK4"          ,"vector<int>"      ,&chHadMultAK4_);   
  outTree_->Branch("chMultAK4"              ,"vector<int>"      ,&chMultAK4_);   
  outTree_->Branch("neHadMultAK4"           ,"vector<int>"      ,&neHadMultAK4_);   
  outTree_->Branch("neMultAK4"              ,"vector<int>"      ,&neMultAK4_);   
  outTree_->Branch("phoMultAK4"             ,"vector<int>"      ,&phoMultAK4_);   
  //outTree_->Branch("jetMassPrunedAK4"        ,"vector<float>"     ,&massPrunedAK4_);
  //outTree_->Branch("jetTau1AK4"              ,"vector<float>"     ,&tau1AK4_);
  //outTree_->Branch("jetTau2AK4"              ,"vector<float>"     ,&tau2AK4_);
  //outTree_->Branch("jetDRAK4"                ,"vector<float>"     ,&dRAK4_); 
  //outTree_->Branch("cutbasedJetId"             ,"vector<float>"     ,&cutbasedJetId_);
  //outTree_->Branch("fullJetId"                 ,"vector<float>"     ,&fullJetId_);
  //outTree_->Branch("fullJetDiscriminant"       ,"vector<float>"     ,&fullJetDiscriminant_);

  ptAK8_             = new std::vector<float>;
  jecAK8_            = new std::vector<float>;
  etaAK8_            = new std::vector<float>;
  phiAK8_            = new std::vector<float>;
  massAK8_           = new std::vector<float>;
  energyAK8_         = new std::vector<float>;
  areaAK8_           = new std::vector<float>;
  chfAK8_            = new std::vector<float>;
  nhfAK8_            = new std::vector<float>;
  phfAK8_            = new std::vector<float>;
  mufAK8_            = new std::vector<float>;
  elfAK8_            = new std::vector<float>;
  // Hadronic forward hadrons
  hf_hfAK8_          = new std::vector<float>;
  // Hadronic forward photons
  hf_emfAK8_         = new std::vector<float>;
  idLAK8_            = new std::vector<int>;
  idTAK8_            = new std::vector<int>;
  massPrunedAK8_     = new std::vector<float>;
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
  outTree_->Branch("jetChfAK8"               ,"vector<float>"     ,&chfAK8_);
  outTree_->Branch("jetNhfAK8"               ,"vector<float>"     ,&nhfAK8_);
  outTree_->Branch("jetPhfAK8"               ,"vector<float>"     ,&phfAK8_);
  outTree_->Branch("jetMufAK8"               ,"vector<float>"     ,&mufAK8_);
  outTree_->Branch("jetElfAK8"               ,"vector<float>"     ,&elfAK8_); 
  outTree_->Branch("jetElfAK8"               ,"vector<float>"     ,&elfAK8_);
  outTree_->Branch("jetHf_hfAK8"             ,"vector<float>"     ,&hf_hfAK8_);
  outTree_->Branch("jetHf_emfAK8"            ,"vector<float>"     ,&hf_emfAK8_);  
  outTree_->Branch("idLAK8"                  ,"vector<int>"      ,&idLAK8_);   
  outTree_->Branch("idTAK8"                  ,"vector<int>"      ,&idTAK8_);   
  outTree_->Branch("jetMassPrunedAK8"        ,"vector<float>"     ,&massPrunedAK8_);
  outTree_->Branch("jetTau1AK8"              ,"vector<float>"     ,&tau1AK8_);
  outTree_->Branch("jetTau2AK8"              ,"vector<float>"     ,&tau2AK8_);
  outTree_->Branch("jetTau3AK8"              ,"vector<float>"     ,&tau3AK8_); 
  //outTree_->Branch("jetDRAK8"                ,"vector<float>"     ,&dRAK8_); 
  outTree_->Branch("chHadMultAK8"          ,"vector<int>"      ,&chHadMultAK8_);   
  outTree_->Branch("chMultAK8"              ,"vector<int>"      ,&chMultAK8_);   
  outTree_->Branch("neHadMultAK8"           ,"vector<int>"      ,&neHadMultAK8_);   
  outTree_->Branch("neMultAK8"              ,"vector<int>"      ,&neMultAK8_);   
  outTree_->Branch("phoMultAK8"             ,"vector<int>"      ,&phoMultAK8_);   
 
  // ptCA8_             = new std::vector<float>;
  // jecCA8_            = new std::vector<float>;
  // etaCA8_            = new std::vector<float>;
  // phiCA8_            = new std::vector<float>;
  // massCA8_           = new std::vector<float>;
  // energyCA8_         = new std::vector<float>;
  // chfCA8_            = new std::vector<float>;
  // nhfCA8_            = new std::vector<float>;
  // phfCA8_            = new std::vector<float>;
  // mufCA8_            = new std::vector<float>;
  // elfCA8_            = new std::vector<float>;
  // idLCA8_            = new std::vector<int>;
  // idTCA8_            = new std::vector<int>;
  // massPrunedCA8_     = new std::vector<float>;
  // tau1CA8_           = new std::vector<float>;
  // tau2CA8_           = new std::vector<float>;
  // tau3CA8_           = new std::vector<float>;
  // //dRCA8_             = new std::vector<float>;
  // outTree_->Branch("jetPtCA8"                ,"vector<float>"     ,&ptCA8_);
  // outTree_->Branch("jetJecCA8"               ,"vector<float>"     ,&jecCA8_);
  // outTree_->Branch("jetEtaCA8"               ,"vector<float>"     ,&etaCA8_);
  // outTree_->Branch("jetPhiCA8"               ,"vector<float>"     ,&phiCA8_);
  // outTree_->Branch("jetMassCA8"              ,"vector<float>"     ,&massCA8_);
  // outTree_->Branch("jetEnergyCA8"            ,"vector<float>"     ,&energyCA8_);
  // outTree_->Branch("jetChfCA8"               ,"vector<float>"     ,&chfCA8_);
  // outTree_->Branch("jetNhfCA8"               ,"vector<float>"     ,&nhfCA8_);
  // outTree_->Branch("jetPhfCA8"               ,"vector<float>"     ,&phfCA8_);
  // outTree_->Branch("jetMufCA8"               ,"vector<float>"     ,&mufCA8_);
  // outTree_->Branch("jetElfCA8"               ,"vector<float>"     ,&elfCA8_);   
  // outTree_->Branch("idLCA8"                  ,"vector<int>"      ,&idLCA8_);   
  // outTree_->Branch("idTCA8"                  ,"vector<int>"      ,&idTCA8_);   
  // outTree_->Branch("jetMassPrunedCA8"        ,"vector<float>"     ,&massPrunedCA8_);
  // outTree_->Branch("jetTau1CA8"              ,"vector<float>"     ,&tau1CA8_);
  // outTree_->Branch("jetTau2CA8"              ,"vector<float>"     ,&tau2CA8_);
  // outTree_->Branch("jetTau3CA8"              ,"vector<float>"     ,&tau3CA8_);
  // //outTree_->Branch("jetDRCA8"                ,"vector<float>"     ,&dRCA8_); 



  //------------------------------------------------------------------
  triggerResult_ = new std::vector<bool>;
  outTree_->Branch("triggerResult","vector<bool>",&triggerResult_);

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
  //  outTree_->Branch("nGenJetsCA8"             ,&nGenJetsCA8_          ,"nGenJetsCA8_/I");

  // ptGenCA8_             = new std::vector<float>;
  // etaGenCA8_            = new std::vector<float>;
  // phiGenCA8_            = new std::vector<float>;
  // massGenCA8_           = new std::vector<float>;
  // energyGenCA8_         = new std::vector<float>;
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
  // outTree_->Branch("jetPtGenCA8"                ,"vector<float>"     ,&ptGenCA8_);
  // outTree_->Branch("jetEtaGenCA8"               ,"vector<float>"     ,&etaGenCA8_);
  // outTree_->Branch("jetPhiGenCA8"               ,"vector<float>"     ,&phiGenCA8_);
  // outTree_->Branch("jetMassGenCA8"              ,"vector<float>"     ,&massGenCA8_);
  // outTree_->Branch("jetEnergyGenCA8"            ,"vector<float>"     ,&energyGenCA8_);




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
  delete chfAK4_;
  delete nhfAK4_;
  delete phfAK4_;
  delete mufAK4_;
  delete elfAK4_;
  delete hf_hfAK4_;
  delete hf_emfAK4_;
  delete idLAK4_;
  delete idTAK4_;
  delete chHadMultAK4_ ;
  delete chMultAK4_    ;
  delete neHadMultAK4_ ;
  delete neMultAK4_    ;
  delete phoMultAK4_   ;
  //delete massPrunedAK4_;
  //delete tau1AK4_;
  //delete tau2AK4_;
  //delete dRAK4_;
  //delete cutbasedJetId_      ;
  //delete fullJetId_          ;
  //delete fullJetDiscriminant_;
  delete ptAK8_;
  delete jecAK8_;
  delete etaAK8_;
  delete phiAK8_;
  delete massAK8_;
  delete energyAK8_;
  delete areaAK8_;
  delete chfAK8_;
  delete nhfAK8_;
  delete phfAK8_;
  delete mufAK8_;
  delete elfAK8_;
  delete hf_hfAK8_;
  delete hf_emfAK8_;
  delete idLAK8_;
  delete idTAK8_;
  delete massPrunedAK8_;
  delete tau1AK8_;
  delete tau2AK8_;
  delete tau3AK8_;
  delete chHadMultAK8_;
  delete chMultAK8_    ;
  delete neHadMultAK8_ ;
  delete neMultAK8_    ;
  delete phoMultAK8_   ;

  //delete dRAK8_;
  
  // delete ptCA8_;
  // delete jecCA8_;
  // delete etaCA8_;
  // delete phiCA8_;
  // delete massCA8_;
  // delete energyCA8_;
  // delete chfCA8_;
  // delete nhfCA8_;
  // delete phfCA8_;
  // delete mufCA8_;
  // delete elfCA8_;
  // delete idLCA8_;
  // delete idTCA8_;
  // delete massPrunedCA8_;
  // delete tau1CA8_;
  // delete tau2CA8_;
  // delete tau3CA8_;
  // //delete dRCA8_;
  
  // delete ptGenCA8_      ;
  // delete etaGenCA8_     ;
  // delete phiGenCA8_     ;
  // delete massGenCA8_    ;
  // delete energyGenCA8_  ;
  // delete ptGenAK4_      ;
  // delete etaGenAK4_     ;
  // delete phiGenAK4_     ;
  // delete massGenAK4_    ;
  // delete energyGenAK4_  ;
  // delete ptGenAK8_      ;
  // delete etaGenAK8_     ;
  // delete phiGenAK8_     ;
  // delete massGenAK8_    ;
  // delete energyGenAK8_  ;


  
  for(unsigned i=0;i<vtriggerSelector_.size();i++) {
    delete vtriggerSelector_[i];
  }
}
//////////////////////////////////////////////////////////////////////////////////////////
void DijetTreeProducer::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) 
{
  initialize();

  edm::Handle<edm::View<pat::Jet> > jetsAK4;
  iEvent.getByLabel(srcJetsAK4_,jetsAK4);

  edm::Handle<edm::View<pat::Jet> > jetsAK8;
  iEvent.getByLabel(srcJetsAK8_,jetsAK8);

  // edm::Handle<edm::View<pat::Jet> > jetsCA8;
  // iEvent.getByLabel(srcJetsCA8_,jetsCA8);
  // edm::View<pat::Jet> pat_jetsCA8 = *jetsCA8;

  edm::Handle<edm::View<reco::GenJet> > handle_genJetsAK4;
  iEvent.getByLabel(srcGenJetsAK4_,handle_genJetsAK4);
  edm::View<reco::GenJet> genJetsAK4 = *handle_genJetsAK4;
  
  edm::Handle<edm::View<reco::GenJet> > handle_genJetsAK8;
  iEvent.getByLabel(srcGenJetsAK8_,handle_genJetsAK8);
  edm::View<reco::GenJet> genJetsAK8 = *handle_genJetsAK8;
  
  // edm::Handle<edm::View<reco::GenJet> > handle_genJetsCA8;
  // iEvent.getByLabel(srcGenJetsCA8_,handle_genJetsCA8);
  // edm::View<reco::GenJet> genJetsCA8 = *handle_genJetsCA8;

  edm::Handle<double>  rho;
  iEvent.getByLabel(srcRho_,rho);

  edm::Handle<edm::View<MET> >  met;
  iEvent.getByLabel(srcMET_,met);

  edm::Handle<reco::VertexCollection> recVtxs;
  iEvent.getByLabel(srcVrtx_,recVtxs);
  
  //---------- pu -----------------------
  edm::Handle<std::vector<PileupSummaryInfo> > PupInfo;
  if (!iEvent.isRealData()) {
    iEvent.getByLabel(srcPU_,PupInfo);
    
    //std::cout << "PupInfo.isValid()? : " << PupInfo.isValid() << endl;

    if(PupInfo.isValid()) {
      for( std::vector<PileupSummaryInfo>::const_iterator it = PupInfo->begin(); it != PupInfo->end(); ++it ) {
	npu_ -> push_back ( it -> getTrueNumInteractions() );
	Number_interactions -> push_back ( it->getPU_NumInteractions() ); 
	OriginBX -> push_back ( it -> getBunchCrossing());                
	
      }
    }
    else {
      edm::LogError("DijetTreeProducer: PileUpError") << "Error! Can't get the product " << srcPU_;
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
    iEvent.getByLabel(srcGenInfo_,genEvtInfo);
    
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
    iEvent.getByLabel(srcPrunedGenParticles_, prunedGenParticles);
    

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
     
  //----- at least one good vertex -----------
  bool cut_vtx = (recVtxs->size() > 0);
  
  if (cut_vtx) {
    
    // Event
    rho_    = *rho;
    met_    = (*met)[0].et();
    if ((*met)[0].sumEt() > 0) {
      metSig_ = (*met)[0].et()/(*met)[0].sumEt();
    }
    nVtx_   = recVtxs->size();
    run_    = iEvent.id().run();
    evt_    = iEvent.id().event();
    lumi_   = iEvent.id().luminosityBlock();
    
    // AK4
    std::vector<double> jecFactorsAK4;
    std::vector<unsigned> sortedAK4JetIdx;
    if(redoJECs_)
    {
      // sort AK4 jets by increasing pT
      std::multimap<double, unsigned> sortedAK4Jets;
      for(edm::View<pat::Jet>::const_iterator ijet = jetsAK4->begin();ijet != jetsAK4->end(); ++ijet)
      {
        JetCorrectorAK4->setJetEta(ijet->eta());
        JetCorrectorAK4->setJetPt(ijet->correctedJet(0).pt());
        JetCorrectorAK4->setJetA(ijet->jetArea());
        JetCorrectorAK4->setRho(rho_);

        double correction = JetCorrectorAK4->getCorrection();

        jecFactorsAK4.push_back(correction);
        sortedAK4Jets.insert(std::make_pair(ijet->correctedJet(0).pt()*correction, ijet - jetsAK4->begin()));
      }
      // get jet indices in decreasing pT order
      for(std::multimap<double, unsigned>::const_reverse_iterator it = sortedAK4Jets.rbegin(); it != sortedAK4Jets.rend(); ++it)
        sortedAK4JetIdx.push_back(it->second);
    }
    else
    {
      for(edm::View<pat::Jet>::const_iterator ijet = jetsAK4->begin();ijet != jetsAK4->end(); ++ijet)
      {
        jecFactorsAK4.push_back(1./ijet->jecFactor(0));
        sortedAK4JetIdx.push_back(ijet - jetsAK4->begin());
      }
    }

    nJetsAK4_ = 0;
    float htAK4(0.0);
    vector<TLorentzVector> vP4AK4;
    for(std::vector<unsigned>::const_iterator i = sortedAK4JetIdx.begin(); i != sortedAK4JetIdx.end(); ++i) {

      edm::View<pat::Jet>::const_iterator ijet = (jetsAK4->begin() + *i);
      double chf = ijet->chargedHadronEnergyFraction();
      double nhf = ijet->neutralHadronEnergyFraction(); // + ijet->HFHadronEnergyFraction();
      double phf = ijet->photonEnergy()/(ijet->jecFactor(0) * ijet->energy());
      double elf = ijet->electronEnergy()/(ijet->jecFactor(0) * ijet->energy());
      double muf = ijet->muonEnergy()/(ijet->jecFactor(0) * ijet->energy());
      
      double hf_hf = ijet->HFHadronEnergyFraction();
      double hf_emf= ijet->HFEMEnergyFraction();
      
      int chm    = ijet->chargedHadronMultiplicity();
      int npr    = ijet->chargedMultiplicity() + ijet->neutralMultiplicity(); 
      float eta  = fabs(ijet->eta());
      float pt   = ijet->correctedJet(0).pt()*jecFactorsAK4.at(*i);
      int idL   = (npr>1 && phf<0.99 && nhf<0.99);
      int idT   = (idL && ((eta<=2.4 && nhf<0.9 && phf<0.9 && elf<0.99 && muf<0.99 && chf>0 && chm>0) || eta>2.4));
      int chHadMult = ijet->chargedHadronMultiplicity();
      int chMult = ijet->chargedMultiplicity();
      int neHadMult = ijet->neutralHadronMultiplicity();
      int neMult = ijet->neutralMultiplicity();
      int phoMult = ijet->photonMultiplicity();
       
      
      if (pt > ptMinAK4_) {
        htAK4 += pt;
        nJetsAK4_++;

        vP4AK4.push_back(TLorentzVector(ijet->correctedJet(0).px()*jecFactorsAK4.at(*i),ijet->correctedJet(0).py()*jecFactorsAK4.at(*i),ijet->correctedJet(0).pz()*jecFactorsAK4.at(*i),ijet->correctedJet(0).energy()*jecFactorsAK4.at(*i)));
        chfAK4_           ->push_back(chf);
        nhfAK4_           ->push_back(nhf);
        phfAK4_           ->push_back(phf);
        elfAK4_           ->push_back(elf);
        mufAK4_           ->push_back(muf);
        hf_hfAK4_         ->push_back(hf_hf);
        hf_emfAK4_        ->push_back(hf_emf);
        jecAK4_           ->push_back(jecFactorsAK4.at(*i));
        ptAK4_            ->push_back(pt);
        phiAK4_           ->push_back(ijet->phi());
        etaAK4_           ->push_back(ijet->eta());
        massAK4_          ->push_back(ijet->correctedJet(0).mass()*jecFactorsAK4.at(*i));
        energyAK4_        ->push_back(ijet->correctedJet(0).energy()*jecFactorsAK4.at(*i));
        areaAK4_          ->push_back(ijet->jetArea());
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

	//---- match with the pruned jet collection -----
        // already in jettoolbok (?)

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

      }// matching with pruned jets


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
      for(edm::View<pat::Jet>::const_iterator ijet = jetsAK8->begin();ijet != jetsAK8->end(); ++ijet)
      {
        JetCorrectorAK8->setJetEta(ijet->eta());
        JetCorrectorAK8->setJetPt(ijet->correctedJet(0).pt());
        JetCorrectorAK8->setJetA(ijet->jetArea());
        JetCorrectorAK8->setRho(rho_);

        double correction = JetCorrectorAK8->getCorrection();

        jecFactorsAK8.push_back(correction);
        sortedAK8Jets.insert(std::make_pair(ijet->correctedJet(0).pt()*correction, ijet - jetsAK8->begin()));
      }
      // get jet indices in decreasing pT order
      for(std::multimap<double, unsigned>::const_reverse_iterator it = sortedAK8Jets.rbegin(); it != sortedAK8Jets.rend(); ++it)
        sortedAK8JetIdx.push_back(it->second);
    }
    else
    {
      for(edm::View<pat::Jet>::const_iterator ijet = jetsAK8->begin();ijet != jetsAK8->end(); ++ijet)
      {
        jecFactorsAK8.push_back(1./ijet->jecFactor(0));
        sortedAK8JetIdx.push_back(ijet - jetsAK8->begin());
      }
    }

    nJetsAK8_ = 0;
    float htAK8(0.0);
    vector<TLorentzVector> vP4AK8;
    for(std::vector<unsigned>::const_iterator i = sortedAK8JetIdx.begin(); i != sortedAK8JetIdx.end(); ++i) {

      edm::View<pat::Jet>::const_iterator ijet = (jetsAK8->begin() + *i);
      double chf = ijet->chargedHadronEnergyFraction();
      double nhf = ijet->neutralHadronEnergyFraction();// + ijet->HFHadronEnergyFraction();
      double phf = ijet->photonEnergy()/(ijet->jecFactor(0) * ijet->energy());
      double elf = ijet->electronEnergy()/(ijet->jecFactor(0) * ijet->energy());
      double muf = ijet->muonEnergy()/(ijet->jecFactor(0) * ijet->energy());
      
      double hf_hf = ijet->HFHadronEnergyFraction();
      double hf_emf= ijet->HFEMEnergyFraction();
      
      int chm    = ijet->chargedHadronMultiplicity();
      int npr    = ijet->chargedMultiplicity() + ijet->neutralMultiplicity(); 
      float eta  = fabs(ijet->eta());
      float pt   = ijet->correctedJet(0).pt()*jecFactorsAK8.at(*i);
      int idL   = (npr>1 && phf<0.99 && nhf<0.99);
      int idT   = (idL && ((eta<=2.4 && nhf<0.9 && phf<0.9 && elf<0.99 && muf<0.99 && chf>0 && chm>0) || eta>2.4));
      int chHadMult = ijet->chargedHadronMultiplicity();
      int chMult = ijet->chargedMultiplicity();
      int neHadMult = ijet->neutralHadronMultiplicity();
      int neMult = ijet->neutralMultiplicity();
      int phoMult = ijet->photonMultiplicity();
      if (pt > ptMinAK8_) {
	htAK8 += pt;
	nJetsAK8_++;

        vP4AK8.push_back(TLorentzVector(ijet->correctedJet(0).px()*jecFactorsAK8.at(*i),ijet->correctedJet(0).py()*jecFactorsAK8.at(*i),ijet->correctedJet(0).pz()*jecFactorsAK8.at(*i),ijet->correctedJet(0).energy()*jecFactorsAK8.at(*i)));
        chfAK8_           ->push_back(chf);
        nhfAK8_           ->push_back(nhf);
        phfAK8_           ->push_back(phf);
        elfAK8_           ->push_back(elf);
        mufAK8_           ->push_back(muf);
        hf_hfAK8_         ->push_back(hf_hf);
        hf_emfAK8_        ->push_back(hf_emf);
        jecAK8_           ->push_back(jecFactorsAK8.at(*i));
        ptAK8_            ->push_back(pt);
        phiAK8_           ->push_back(ijet->phi());
        etaAK8_           ->push_back(ijet->eta());
        massAK8_          ->push_back(ijet->correctedJet(0).mass()*jecFactorsAK8.at(*i));
        energyAK8_        ->push_back(ijet->correctedJet(0).energy()*jecFactorsAK8.at(*i));
        areaAK8_          ->push_back(ijet->jetArea());
	idLAK8_           ->push_back(idL);
	idTAK8_           ->push_back(idT);
        tau1AK8_          ->push_back(ijet->userFloat("NjettinessAK8:tau1"));
        tau2AK8_          ->push_back(ijet->userFloat("NjettinessAK8:tau2"));
        tau3AK8_          ->push_back(ijet->userFloat("NjettinessAK8:tau3"));
	massPrunedAK8_    ->push_back(ijet->userFloat("ak8PFJetsCHSPrunedLinks"));
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

    
    // // CA8
    // nJetsCA8_ = 0;
    // float htCA8(0.0);
    // vector<TLorentzVector> vP4CA8;
    // for(edm::View<pat::Jet>::const_iterator ijet = pat_jetsCA8.begin();ijet != pat_jetsCA8.end(); ++ijet) { 
    //   double chf = ijet->chargedHadronEnergyFraction();
    //   double nhf = ijet->neutralHadronEnergyFraction() + ijet->HFHadronEnergyFraction();
    //   double phf = ijet->photonEnergy()/(ijet->jecFactor(0) * ijet->energy());
    //   double elf = ijet->electronEnergy()/(ijet->jecFactor(0) * ijet->energy());
    //   double muf = ijet->muonEnergy()/(ijet->jecFactor(0) * ijet->energy());
    //   int chm    = ijet->chargedHadronMultiplicity();
    //   int npr    = ijet->chargedMultiplicity() + ijet->neutralMultiplicity(); 
    //   float eta  = fabs(ijet->eta());
    //   float pt   = ijet->pt();
    //   int idL   = (npr>1 && phf<0.99 && nhf<0.99);
    //   int idT   = (idL && ((eta<=2.4 && nhf<0.9 && phf<0.9 && elf<0.99 && muf<0.99 && chf>0 && chm>0) || eta>2.4));
    //   if (pt > ptMinCA8_) {
    //     htCA8 += pt;
    //     nJetsCA8_++;
	
    //     vP4CA8.push_back(TLorentzVector(ijet->px(),ijet->py(),ijet->pz(),ijet->energy()));
    //     chfCA8_           ->push_back(chf);
    //     nhfCA8_           ->push_back(nhf);
    //     phfCA8_           ->push_back(phf);
    //     elfCA8_           ->push_back(elf);
    //     mufCA8_           ->push_back(muf);
    //     jecCA8_           ->push_back(1./ijet->jecFactor(0));
    //     ptCA8_            ->push_back(pt);
    //     phiCA8_           ->push_back(ijet->phi());
    //     etaCA8_           ->push_back(ijet->eta());
    //     massCA8_          ->push_back(ijet->mass());
    //     energyCA8_        ->push_back(ijet->energy());
    // 	idLCA8_           ->push_back(idL);
    // 	idTCA8_           ->push_back(idT);
    //     tau1CA8_          ->push_back(ijet->userFloat("NjettinessCA8:tau1"));
    //     tau2CA8_          ->push_back(ijet->userFloat("NjettinessCA8:tau2"));
    //     tau3CA8_          ->push_back(ijet->userFloat("NjettinessCA8:tau3"));
    // 	massPrunedCA8_    ->push_back(ijet->userFloat("ca8PFJetsCHSPrunedLinks"));
	
    //   } 
    // }// jet loop  
    // htCA8_     = htCA8;
    // if (nJetsCA8_ > 1) { //assuming jets are ordered by pt in the pat collection
    //   mjjCA8_    = (vP4CA8[0]+vP4CA8[1]).M();
    //   dEtajjCA8_ = fabs((*etaCA8_)[0]-(*etaCA8_)[1]); 
    //   dPhijjCA8_ = fabs(deltaPhi((*phiCA8_)[0],(*phiCA8_)[1]));
    // }
  
    
    //-------------- Gen Jets Info -----------------------------------

    
    nGenJetsAK4_ = 0;
    vector<TLorentzVector> vP4GenAK4;
    if (!iEvent.isRealData()) {
      for(edm::View<reco::GenJet>::const_iterator ijet = genJetsAK4.begin();ijet != genJetsAK4.end(); ++ijet) { 
	
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
      
      for(edm::View<reco::GenJet>::const_iterator ijet = genJetsAK8.begin();ijet != genJetsAK8.end(); ++ijet) { 
	
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
      
      // nGenJetsCA8_ = 0;
      // vector<TLorentzVector> vP4GenCA8;
      
      // for(edm::View<pat::Jet>::const_iterator ijet = pat_jetsCA8.begin();ijet != pat_jetsCA8.end(); ++ijet) { 
	
	
      // 	//float eta  = fabs(ijet->eta());
      // 	float pt   = ijet->pt();
      // 	if (pt > ptMinCA8_) {
      // 	  nGenJetsCA8_++;
      // 	  vP4GenCA8.push_back(TLorentzVector(ijet->px(),ijet->py(),ijet->pz(),ijet->energy()));
      // 	  ptGenCA8_            ->push_back(pt);
      // 	  phiGenCA8_           ->push_back(ijet->phi());
      // 	  etaGenCA8_           ->push_back(ijet->eta());
      // 	  massGenCA8_          ->push_back(ijet->mass());
      // 	  energyGenCA8_        ->push_back(ijet->energy());
      // 	}
      // }// jet loop  
    }//if MC 
  }// if vtx
  
  
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
  chfAK4_            ->clear();
  nhfAK4_            ->clear();
  phfAK4_            ->clear();
  elfAK4_            ->clear();
  mufAK4_            ->clear();
  hf_hfAK4_             ->clear();
  hf_emfAK4_            ->clear();
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
  chfAK8_            ->clear();
  nhfAK8_            ->clear();
  phfAK8_            ->clear();
  elfAK8_            ->clear();
  mufAK8_            ->clear();
  hf_hfAK8_             ->clear();
  hf_emfAK8_            ->clear();
  jecAK8_            ->clear();
  jecAK8_            ->clear();
  idLAK8_            ->clear();
  idTAK8_            ->clear();
  massPrunedAK8_     ->clear();
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
  
  
  // nJetsCA8_          = -999;
  // htCA8_             = -999;
  // mjjCA8_            = -999; 
  // dEtajjCA8_         = -999; 
  // dPhijjCA8_         = -999;
  // ptCA8_             ->clear();
  // etaCA8_            ->clear();
  // phiCA8_            ->clear();
  // massCA8_           ->clear();
  // energyCA8_         ->clear();
  // chfCA8_            ->clear();
  // nhfCA8_            ->clear();
  // phfCA8_            ->clear();
  // elfCA8_            ->clear();
  // mufCA8_            ->clear();
  // jecCA8_            ->clear();
  // jecCA8_            ->clear();
  // idLCA8_            ->clear();
  // idTCA8_            ->clear();
  // massPrunedCA8_     ->clear();
  // tau1CA8_           ->clear();
  // tau2CA8_           ->clear();
  // //dRCA8_             ->clear();
  
  
  triggerResult_     ->clear();
  
  //----- MC -------
  npu_ ->clear();
  Number_interactions ->clear();
  OriginBX            -> clear();
  
  ptHat_ = -999; 
  processID_ = -999; 
  weight_ = -999;
  
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
  // ptGenCA8_    ->clear();
  // phiGenCA8_   ->clear();
  // etaGenCA8_   ->clear();
  // massGenCA8_  ->clear();
  // energyGenCA8_->clear();
  
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



