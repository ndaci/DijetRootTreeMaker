#include <iostream>
#include <sstream>
#include "TChain.h"
#include "TString.h"
#include "TH1F.h"
#include "TCanvas.h"

using namespace std;

Int_t plots(TString process="mc_zz_qqqq")
{

  // Define chain
  TChain* ch = new TChain("dijets/events");
  ch->Add("tree_ZZ4Q_20k.root");
  
  // Define jet collections
  const UInt_t nColl=2; // 0: MINIAOD ; 1: JetToolbox
  TString collName [nColl] = {"miniaod"  , "toolbox"};
  TString collTitle[nColl] = {"(miniaod)", "(toolbox)"};
  TString mode     [nColl] = {"H", "HSAME"};
  Int_t   color    [nColl] = {kBlue, kRed};

  // Define variables
  const UInt_t nVar = 14;
  TString varName [nVar] = {"jetPtAK8"  , "jetEtaAK8",                               // AK8 CHS   Jet kinematics
			    "jetMassAK8", "jetMassPrunedAK8", "jetMassSoftDropAK8",  // AK8 CHS   Jet masses
			    "jetTau1AK8", "jetTau2AK8",                              // AK8 CHS   Jet substructure
			    "jetPtAK8_Puppi"  , "jetEtaAK8_Puppi",                   // AK8 PUPPI Jet kinematics
			    "jetMassAK8_Puppi", "jetMassPrunedAK8_Puppi", "jetMassSoftDropAK8_Puppi", // AK8 PUPPI Jet masses
			    "jetTau1AK8_Puppi", "jetTau2AK8_Puppi" };                // AK8 PUPPI Jet substructure

  TString varTitle[nVar] = {"AK8 CHS Jet p_{T}", "AK8 CHS Jet #eta", 
			    "AK8 CHS Jet Mass" , "AK8 CHS Jet Pruned Mass", "AK8 CHS Jet SoftDrop Mass", 
			    "AK8 CHS #tau_{1}" , "AK8 CHS #tau_{2}",
			    "AK8 PUPPI Jet p_{T}", "AK8 PUPPI Jet #eta", 
			    "AK8 PUPPI Jet Mass" , "AK8 PUPPI Jet Pruned Mass", "AK8 PUPPI Jet SoftDrop Mass", 
			    "AK8 PUPPI #tau_{1}" , "AK8 PUPPI #tau_{2}" };

  UInt_t  nBins   [nVar] = {100  ,   50,  40,  40,  40, 20, 20,
			    100  ,   50,  40,  40,  40, 20, 20};

  Float_t firstbin[nVar] = {160  , -3.0,   0,   0,   0,  0,  0,
			    160  , -3.0,   0,   0,   0,  0,  0};

  Float_t lastbin [nVar] = {2160 ,  3.0, 200, 200, 200,  1,  1,
			    2160 ,  3.0, 200, 200, 200,  1,  1};
  //
  UInt_t  nSplit  [nVar];
  for(UInt_t iV=0; iV<nVar; iV++) {
    nSplit[iV] = 3;
  }

  // Split variables (first, second, third jet...)
  TString Order[4]={"First ", "Second ", "Third ", "Fourth "};
  vector<TString> varIdx     [nVar];
  vector<TString> varIdxName [nVar];
  vector<TString> varIdxTitle[nVar];
  ostringstream oss;
  TString toss;
  //
  for(UInt_t iV=0; iV<nVar; iV++) {
    if(nSplit[iV]==0) {
      varIdx     [iV].push_back("");
      varIdxName [iV].push_back("");
      varIdxTitle[iV].push_back("");
    }
    else {
      for(UInt_t iS=0; iS<nSplit[iV]; iS++) {
	toss="";
	oss.str("");
	oss << iS;
	toss = TString(oss.str().c_str());
	//
	varIdx     [iV].push_back("["+toss+"]");
	varIdxName [iV].push_back("J"+toss);
	varIdxTitle[iV].push_back(Order[iS]);
      }
    }
  }

  // Declare histograms
  map<TString, TH1F*> mapHistos;
  TH1F* hTemp;
  vector<TString> bname[nVar][nColl];
  vector<TString> bcut [nVar][nColl];
  vector<TString> name [nVar][nColl];
  vector<TString> title[nVar][nColl];
  //
  for(UInt_t iV=0; iV<nVar; iV++) {
    for(UInt_t iC=0; iC<nColl; iC++) {
      for(UInt_t iS=0; iS<varIdx[iV].size(); iS++) {

	// histogram title & name, branch name
	title[iV][iC].push_back(varIdxTitle[iV][iS]+varTitle[iV]+" ("+collName[iC]+")");
	name [iV][iC].push_back(varName [iV]+"_"   +collName[iC]+"_" +varIdxName [iV][iS]);
	bname[iV][iC].push_back(varName [iV]+"_"   +collName[iC]     +varIdx     [iV][iS]);

	// selection cut
	bcut[iV][iC].push_back("jetPtAK8_"+collName[iC]+varIdx[iV][iS]+">=170");
	//cout << bcut[iV][iC][iS] << endl;

	// declare histogram
	hTemp = new TH1F(name[iV][iC][iS], title[iV][iC][iS], nBins[iV], firstbin[iV], lastbin[iV]);
	hTemp -> SetXTitle(title[iV][iC][iS]);
	hTemp -> SetLineColor(color[iC]);
	hTemp -> SetLineWidth(2);
	hTemp -> SetYTitle("Events");

	// map histogram
	mapHistos[name[iV][iC][iS]] = hTemp;

      }
    }
  }

  // Create TCanvas
  TString plotName="";
  TCanvas *c = new TCanvas("c","c",100,100,600,600);
  gPad->SetLogy();

  // Loop over variables and indices iS (first/second... jet)
  for(UInt_t iV=0; iV<nVar; iV++) {
    for(UInt_t iS=0; iS<varIdx[iV].size(); iS++) {

      // avoid un-existing collections
      if(varName[iV]=="jetMassPrunedAK8_Puppi") continue; //fixme

      // make 1 plot with both collections
      for(UInt_t iC=0; iC<nColl; iC++) {
	hTemp = mapHistos[name[iV][iC][iS]];
	ch->Draw(bname[iV][iC][iS]+">>"+hTemp->GetName(), bcut[iV][iC][iS], mode[iC]);	
      }

      // print plot
      plotName = "plot_"+varName[iV]+"_"+varIdxName[iV][iS]+".pdf";
      c->Print("plots/"+plotName,"pdf");

    }
  }

  // THE END //
  return 0;
}
