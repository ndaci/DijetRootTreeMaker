#!usr/bin/python
import sys, getopt
import ROOT
from ROOT import TFile, TH1F, TCanvas, TPad, TMath
from ROOT import gROOT, gPad 
from ROOT import RooRealVar, RooDataHist, RooPlot, RooArgList, RooArgSet,  RooBernstein, RooCBShape, RooAddPdf, RooFit, RooGenericPdf, RooWorkspace, RooMsgService, RooHistPdf
from setTDRStyle import setTDRStyle

import optparse
usage = "usage: %prog [options]"
parser = optparse.OptionParser(usage)
parser.add_option("-s","--fitSig",action="store_true",default=False,dest="fitSig")
parser.add_option("-H","--histpdfSig",action="store_true",default=False,dest="histpdfSig")
parser.add_option("-d","--fitDat",action="store_true",default=False,dest="fitDat")
parser.add_option("-m","--mass",action="store",type="int",dest="mass",default=5000)
parser.add_option("-o","--path",action="store",type="string",dest="PATH")


parser.add_option("--lumi",action="store",type="float",dest="lumi",default=1000.)
parser.add_option("--sigEff",action="store",type="float",dest="sigEff",default=0.629)
parser.add_option("--sigXS",action="store",type="float",dest="sigXS",default=0.0182)


(options, args) = parser.parse_args()

mass = options.mass
fitSig = options.fitSig
histpdfSig = options.histpdfSig
fitDat = options.fitDat
PATH = options.PATH


gROOT.Reset()
setTDRStyle()
gROOT.ForceStyle()
gROOT.SetStyle('tdrStyle')

RooMsgService.instance().setSilentMode(ROOT.kTRUE)
RooMsgService.instance().setStreamStatus(0,ROOT.kFALSE)
RooMsgService.instance().setStreamStatus(1,ROOT.kFALSE)

# -----------------------------------------
# get histograms
## ---- CERN -------
PATH = '/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_2_1_DiJet/src/CMSDIJET/DijetRootTreeAnalyzer/output/'

#filenameSig = PATH+'Resonance_Shapes_qg_PU20_13TeV_newJEC.root'
filenameSig = PATH+'rootfile_QstarToJJ_M_'+str(mass)+'_Tune4C_13TeV_pythia8__Phys14DR-PU20bx25_PHYS14_25_V1-v1__MINIAODSIM_reduced_skim.root'
#filenameDat = PATH+'../test_fit/dijetFitResults_FuncType0_nParFit4_MC_1fb-1.root'
filenameDat = PATH+'../scripts/histo_bkg_mjj.root'

#inputHistName = 'hist_mass_1GeV'
inputHistName = 'hist_allCutsQCD'


infSig = TFile.Open(filenameSig)
#hSig = filenameSig.Get('h_qg_'+mass)

tree_sig = infSig.Get("rootTupleTree/tree")
hSig = TH1F("hist_mass_1GeV","",13999,1,14000)
tree_sig.Project("hist_mass_1GeV","mjj","deltaETAjj<1.3")
hSig.Print()
#hSig.Draw()

#hSig   = infSig.Get(inputHistName)
#hSig.Rebin(20)

infDat = TFile.Open(filenameDat)
hDat   = infDat.Get(inputHistName)
#hDat.Rebin(20)

# -----------------------------------------
# define observable
x = RooRealVar('mjj','mjj',1118,6099)

if fitSig: 

    # define parameters for signal fit
    m = RooRealVar('mean','mean',float(mass),float(mass)-200,float(mass)+200)
    s = RooRealVar('sigma','sigma',0.1*float(mass),0,10000)
    a = RooRealVar('alpha','alpha',1,-10,10)
    n = RooRealVar('n','n',1,0,100)
    sig = RooCBShape('sig','sig',x,m,s,a,n)        

    p  = RooRealVar('p','p',1,0,5)
    x0 = RooRealVar('x0','x0',1000,100,5000)

    bkg = RooGenericPdf('bkg','1/(exp(pow(@0/@1,@2))+1)',RooArgList(x,x0,p))

    fsig= RooRealVar('fsig','fsig',0.5,0.,1.)
    signal = RooAddPdf('signal','signal',sig,bkg,fsig)

    # -----------------------------------------
    # fit signal
    canSname = 'can_Mjj'+str(mass)
    canS = TCanvas(canSname,canSname,900,600)
    gPad.SetLogy() 

    roohistSig = RooDataHist('roohist','roohist',RooArgList(x),hSig)

    roohistSig.Print() 
    res_sig = signal.fitTo(roohistSig, RooFit.Save(ROOT.kTRUE))
    res_sig.Print()
    frame = x.frame()
    roohistSig.plotOn(frame,RooFit.Binning(166))
    signal.plotOn(frame)
    signal.plotOn(frame,RooFit.Components('bkg'),RooFit.LineColor(ROOT.kRed),RooFit.LineWidth(2),RooFit.LineStyle(ROOT.kDashed))
    frame.GetXaxis().SetRangeUser(1118,6099)
    frame.GetXaxis().SetTitle('m_{jj} (GeV)')
    frame.Draw()

    parsSig = signal.getParameters(roohistSig)
    parsSig.setAttribAll('Constant', True)


if histpdfSig:
    
    # -----------------------------------------
    # hist pdf signal
    canSname = 'can_Mjj'+str(mass)
    canS = TCanvas(canSname,canSname,900,600)
    gPad.SetLogy() 

    roohistSig = RooDataHist('roohist','roohist',RooArgList(x),hSig)
    roohistSig.Print()
    signal = RooHistPdf('signal','signal',RooArgSet(x),roohistSig)
    signal.Print()
    frame = x.frame()
    roohistSig.plotOn(frame,RooFit.Binning(166))
    signal.plotOn(frame,RooFit.Binning(166),RooFit.LineColor(ROOT.kRed),RooFit.LineWidth(2),RooFit.LineStyle(ROOT.kDashed))

    frame.GetXaxis().SetRangeUser(1118,6099)
    frame.GetXaxis().SetTitle('m_{jj} (GeV)')
    frame.Draw()

    parsSig = signal.getParameters(roohistSig)
    parsSig.setAttribAll('Constant', True)



if fitDat: 

    # -----------------------------------------
    # define parameters for background
    # function 0 (standard parametrization)
    NBINS = 166
    #if fitModel==0:
    p1 = RooRealVar('p1','p1',6.21862535247,0.,100)
    p2 = RooRealVar('p2','p2',6.48308946408,0,60)
    p3 = RooRealVar('p3','p3',0.217160577769,0,10)

    background = RooGenericPdf('background','(pow(1-@0/13000,@1)/pow(@0/13000,@2+@3*log(@0/13000)))',RooArgList(x,p1,p2,p3))
    
    ##variation 1, with one more parameter 
    #if fitModel==1:
    #  TMath::Power(1-x/8000,[1])*(1+[4]*x/8000) ) / ( TMath::Power(x/8000,[2]+[3]*log(x/8000))
    #
    
    roohistBkg = RooDataHist('roohist','roohist',RooArgList(x),hDat)
    roohistBkg.Print()
    res = background.fitTo(roohistBkg, RooFit.Save(ROOT.kTRUE))
    res.Print()

    # -----------------------------------------
    # plot background
    canBname = 'can_Mjj_Data'
    canB = TCanvas(canBname,canBname,900,600)
    gPad.SetLogy() 
    canB.cd(1).SetBottomMargin(0.4)

    frame1 = x.frame()
    frame2 = x.frame()
    roohistBkg.plotOn(frame1,RooFit.Binning(NBINS))
    background.plotOn(frame1)
    hpull = frame1.pullHist()
    frame2.addPlotable(hpull,'p')

    frame1.SetMinimum(0.5)
    frame1.GetXaxis().SetTitle('')
    frame1.GetXaxis().SetLabelSize(0.0)
    frame1.GetYaxis().SetTickLength(0.06)
    frame1.Draw()

    pad = TPad('pad','pad',0.,0.,1.,1.)
    pad.SetTopMargin(0.6)
    pad.SetFillColor(0)
    pad.SetFillStyle(0)
    pad.Draw()
    pad.cd(0)
    frame2.SetMinimum(-5)
    frame2.SetMaximum(5)
    frame2.GetYaxis().SetNdivisions(505)
    frame2.GetXaxis().SetTitleOffset(0.9)
    frame2.GetYaxis().SetTitleOffset(0.8)
    frame2.GetYaxis().SetTickLength(0.06)
    frame2.GetYaxis().SetTitleSize(0.05)
    frame2.GetYaxis().SetLabelSize(0.03)
    frame2.GetYaxis().SetTitle('(Data-Fit)/Error')
    frame2.GetXaxis().SetTitle('m_{jj} (GeV)')
    frame2.Draw();

    parsBkg = background.getParameters(roohistBkg)
    #parsBkg.setAttribAll('Constant', True)

if ((histpdfSig or fitSig) and fitDat):
    
    # -----------------------------------------
    # write everything to a workspace to make a datacard
    dcFN = 'Qstar'+str(mass)+'_datacard.txt'
    wsFN = 'Qstar'+str(mass)+'_workspace.root'
    

    nObs = roohistBkg.sumEntries();
    
    w = RooWorkspace('w','workspace')
    getattr(w,'import')(signal)
    getattr(w,'import')(background)
    getattr(w,'import')(roohistBkg,RooFit.Rename("data_obs"))  
    w.Print()
    w.writeToFile(wsFN)
    
    # -----------------------------------------
    # write a datacard
    LUMI = options.lumi
    signalCrossSection = options.sigXS
    signalEfficiency = options.sigEff
    ExpectedSignalRate = signalCrossSection*LUMI*signalEfficiency

    datacard = open(dcFN,'w')
    datacard.write('imax 1\n')
    datacard.write('jmax 1\n')
    datacard.write('kmax *\n')
    datacard.write('---------------\n')
    datacard.write('shapes * * '+wsFN+' w:$PROCESS\n')
    datacard.write('---------------\n')
    datacard.write('bin 1\n')    
    datacard.write('observation '+str(nObs)+'\n')
    datacard.write('------------------------------\n')
    datacard.write('bin          1          1\n')          
    datacard.write('process      signal     background\n')
    datacard.write('process      0          1\n')          
    datacard.write('rate         '+str(ExpectedSignalRate)+'         '+str(nObs)+'\n')
    datacard.write('------------------------------\n')      
    #nuisance parameters --- gaussian prior
    #datacard.write('p1  param    '+str(p1.getValV())+'   '+str(p1.getError())+'\n')
    #datacard.write('p2  param    '+str(p2.getValV())+'   '+str(p2.getError())+'\n')
    #datacard.write('p3  param    '+str(p3.getValV())+'   '+str(p3.getError())+'\n')
    #flat parameters --- flat prior
    datacard.write('p1  flatParam \n')
    datacard.write('p2  flatParam \n')
    datacard.write('p3  flatParam \n')


 
##----- keep the GUI alive ------------
#if __name__ == '__main__':
#  rep = ''
#  while not rep in ['q','Q']:
#    rep = raw_input('enter "q" to quit: ')
#    if 1 < len(rep):
#      rep = rep[0]
