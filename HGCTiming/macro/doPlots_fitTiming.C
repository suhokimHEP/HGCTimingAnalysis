//example to run 
//photon (PDG 22) - kaon 0 long (PDG 130)
//optionType is th einput folder with the produced .root
//
//root -l doPlots_fitTiming.C'("22", "TDRset_0PU_allEta")

#include "TLegend.h"
#include "TLatex.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <iostream>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include "TTree.h"
#include "TChain.h"
#include <vector>
#include <fstream>
#include <string>
#include "TROOT.h"
#include "TSystem.h"
#include <iostream>
#include <fstream>


#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooCBShape.h"
#include "RooArgusBG.h"
#include "RooFFTConvPdf.h"
#include "RooAddPdf.h"
#include "RooWorkspace.h"
#include "RooFitResult.h"
#include "TText.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"


using namespace RooFit;

float getXmax(TH1F* histo, float& YMax){

  float yVal = 0.;
  int xBin = 1;
  for(int iB=1; iB<histo->GetNbinsX(); ++iB){
    if(histo->GetBinContent(iB) > yVal){
      xBin = iB;
      yVal = histo->GetBinContent(iB);
      YMax = yVal;
      if(yVal > 0 && histo->GetBinContent(iB) < yVal) break;
    }
  }

  std::cout << histo->GetName() << " " << histo->GetBinCenter(xBin) << std::endl; 
  return histo->GetBinCenter(xBin);
}


void doPlots_fitTiming(){
  gROOT->Macro("/afs/cern.ch/user/a/amartell/public/setStyle.C");

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  
  gROOT->Reset();
  gROOT->Macro("~/public/setStyle.C");
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1);

  gROOT->SetBatch(kTRUE);

  std::cout << " inizio ci sono " << std::endl; 


  bool doAllTheFits = true;
  
  //  int iColors[16] = {kRed, kOrange+4, kOrange-3, kOrange-2, kBlue, kBlue-9, kAzure-9, kAzure+10, kCyan, kGreen+1, kCyan-2, kYellow+2}; //kGray+1};
  int iColors[7] = {kOrange-3, kRed, kMagenta, kBlue, kCyan, kGreen+1, kGray+2};
  int iStyle[4] = {20, 21, 22, 23}; //kGray+1};
  //  int iColors[6] = {kGreen+1, kBlue, kRed};
  
  int nBinsSide = 2;
  
  float Emin = 0.;
  float Emax = 20.;
  int numberOfBins = 2000;
  float EbinWidth = (Emax - Emin) / numberOfBins;

  int nFiles = 1;

  TH1F* hDummy[7];
  TH1F* hDummySt[2];


  std::vector<float> energyThreshold;
  /*
  energyThreshold.push_back(0.5);
  energyThreshold.push_back(0.8);
  energyThreshold.push_back(1.);
  energyThreshold.push_back(1.2);
  energyThreshold.push_back(1.4);
  energyThreshold.push_back(1.7);
  */
  energyThreshold.push_back(0.05);
  energyThreshold.push_back(0.1);
  energyThreshold.push_back(0.5);
  energyThreshold.push_back(1.);
  energyThreshold.push_back(1.5);
  energyThreshold.push_back(2.);
  energyThreshold.push_back(2.5);
  
  TGraphErrors* tg[7][2];
  TGraphErrors* tgM[7][2];
  TGraphErrors* tgMall[7][2];
  TGraphErrors* tgFrEvt[7][2];
  for(int iT=0; iT<energyThreshold.size(); ++iT){
    hDummy[iT] = new TH1F(Form("hDummy%d", iT), "", 1, 0., 1);
    hDummy[iT]->SetLineColor(iColors[iT]);

    for(int iR=0; iR<nBinsSide; ++iR){
      if(iT == 0){
	hDummySt[iR] = new TH1F(Form("hDummySt%d", iR), "", 1, 0., 1);
	hDummySt[iR]->SetMarkerStyle(iStyle[iR]);
      }

      tg[iT][iR] = new TGraphErrors();
      tg[iT][iR]->SetName(Form("resolution_%d_side%d", iT, iR));
      tg[iT][iR]->SetPoint(0, -1, -1);
      
      tg[iT][iR]->SetLineColor(iColors[iT]);   
      tg[iT][iR]->SetMarkerColor(iColors[iT]); 
      tg[iT][iR]->SetMarkerSize(1.2); 
      tg[iT][iR]->SetLineWidth(2);
      tg[iT][iR]->SetMarkerStyle(iStyle[iR]);
      ////////      
      tgM[iT][iR] = new TGraphErrors();
      tgM[iT][iR]->SetName(Form("mean_%d_side%d", iT, iR));
      tgM[iT][iR]->SetPoint(0, -1, -1);
      
      tgM[iT][iR]->SetLineColor(iColors[iT]);
      tgM[iT][iR]->SetLineWidth(2);
      tgM[iT][iR]->SetMarkerSize(1.2); 
      tgM[iT][iR]->SetMarkerColor(iColors[iT]);
      tgM[iT][iR]->SetMarkerStyle(iStyle[iR]);
      ///////////
      tgMall[iT][iR] = new TGraphErrors();
      tgMall[iT][iR]->SetName(Form("meanAll_%d_side%d", iT, iR));
      tgMall[iT][iR]->SetPoint(0, -1, -1);
      
      tgMall[iT][iR]->SetLineColor(iColors[iT]);
      tgMall[iT][iR]->SetLineWidth(2);
      tgMall[iT][iR]->SetMarkerSize(1.2); 
      tgMall[iT][iR]->SetMarkerColor(iColors[iT]);
      tgMall[iT][iR]->SetMarkerStyle(iStyle[iR]);

      //////////
           
      tgFrEvt[iT][iR] = new TGraphErrors();
      tgFrEvt[iT][iR]->SetName(Form("fractionEvt_%d_side%d", iT, iR));
      tgFrEvt[iT][iR]->SetPoint(0, -1, -1);
      
      tgFrEvt[iT][iR]->SetLineColor(iColors[iT]);
      tgFrEvt[iT][iR]->SetLineWidth(2);
      tgFrEvt[iT][iR]->SetMarkerSize(1.2); 
      tgFrEvt[iT][iR]->SetMarkerColor(iColors[iT]);
      tgFrEvt[iT][iR]->SetMarkerStyle(iStyle[iR]);
  
    }
  }  

  std::cout << " >>> now load files " << std::endl;

  TFile* inF[1];
  for(int ij=0; ij<nFiles; ++ij){
    //inF[ij] = TFile::Open("../test/testTimeCalib1p2Side/JOB_testTimeCalib1p2Side/OutTimeHGC_RecHitsCalib_testTimeCalib1p2Side.root");
    //with randomVtx
    //inF[ij] = TFile::Open("../test/testTimeCalib1p2Side_NuGun/JOB_testTimeCalib1p2Side_NuGun/OutTimeHGC_RecHitsCalib_testTimeCalib1p2Side_NuGun.root");
    //inF[ij] = TFile::Open("../test/OutTimeHGC_RecHitsCalib_testTimeCalib1p2Side_NuGun_Nuvtx.root");
    inF[ij] = TFile::Open("../test/OutTimeHGC_RecHitsCalib_testTimeCalib1p2Side_NuGun_highestPtRECOVtx.root");
  }
  
  std::cout << " >>> files loaded " << std::endl;


  TH2F* h2_timeHitsEnergy_layer;
  TCanvas* c1;

  for(int iL=0; iL<41; ++iL){

    for(int ij=0; ij<energyThreshold.size(); ++ij){

      int ieta = (energyThreshold.at(ij) - Emin) / EbinWidth;
      //    for(int ieta=0; ieta<numberOfBins; ++ieta){
      for(int iRad=0; iRad<nBinsSide; ++iRad){
	
	if(iRad < nBinsSide){
	  h2_timeHitsEnergy_layer = (TH2F*)(inF[0]->Get(Form("ana/h2_timeHitsEnergy_layer%d_iS%d", iL, iRad)));
	  
	  TH1F* slice = (TH1F*)h2_timeHitsEnergy_layer->ProjectionY(Form("timeD_iL%d_iS%d", iL, iRad), ieta, 2000.);

	  // Declare observable x
	  RooRealVar x("x","x",-1,2);
	  
	  // Create a binned dataset that imports contents of TH1 and associates its contents to observable 'x'
	  RooDataHist dT("dT","dT", x, Import(*slice));
	  
	  RooPlot* frameS = x.frame();
	  dT.plotOn(frameS) ; 
	  RooPlot* frame = x.frame();
	  dT.plotOn(frame) ; 
	  
	  // Fit a Gaussian p.d.f to the data
	  RooRealVar meanCBS("meanCBS", "", 0., -1, 2);
	  RooRealVar sigmaCBS("sigmaCBS", "", 0.4, 0.2, 0.6);
	  RooRealVar alphaCBS("alphaCBS", "", -1., -3., 0.);
	  RooRealVar nCBS("nCBS", "", 1, 0., 3);
	  RooCBShape modelCBS("modelCBS", "", x, meanCBS, sigmaCBS, alphaCBS, nCBS);

	  // RooRealVar meanS("meanS","",0,-1.,2.);
          // RooRealVar sigmaS("sigmaS","", 0.4, 0.2, 0.5);
          // RooGaussian gaussS("gaussS","", x, meanS, sigmaS);

	  RooFitResult *rS = modelCBS.fitTo(dT, Minimizer("Minuit2"), Save(true));
	  modelCBS.plotOn(frameS);
	  float chi2_S = frameS->chiSquare(6);

	  RooRealVar mean("mean","",0,-1.,2.);
	  RooRealVar sigma("sigma","", 0.04, 0.02, 0.05);
	  RooGaussian gauss("gauss","", x, mean, sigma);

	  RooRealVar meanB("meanB","",0,-1.,2.);
	  RooRealVar sigmaB("sigmaB","", 0.4, 0.2, 0.5);
	  RooGaussian gaussB("gaussB","", x, meanB, sigmaB);

	  RooRealVar meanCB("meanCB", "", 0., -1, 2);
	  RooRealVar sigmaCB("sigmaCB", "", 0.4, 0.2, 0.6);
	  RooRealVar alphaCB("alphaCB", "", -1., -3., 0.);
	  RooRealVar nCB("nCB", "", 1, 0., 3);
	  RooCBShape modelCB("modelCB", "", x, meanCB, sigmaCB, alphaCB, nCB);

	  // RooRealVar tau("tau", "", -0.8, -0.3, -1.5);
	  // RooExponential expB("expB", "", x, tau);

	  // RooFFTConvPdf lxg("lxg", "expo (X) gauss", x, gaussB, expB);


	  RooRealVar nsig("nsig","#signal events", slice->Integral() / 0.001, 0., slice->Integral() * 0.1);
	  RooRealVar nbkg("nbkg","#background events", slice->Integral(), 0., slice->Integral() * 1.1);
	  //RooAddPdf model("model","g+a",RooArgList(gauss,lxg),RooArgList(nsig,nbkg));
	  //RooAddPdf model("model","",RooArgList(gauss,gaussB),RooArgList(nsig,nbkg));


	  RooAddPdf model("model","",RooArgList(gauss,modelCB),RooArgList(nsig,nbkg));


	  RooFitResult *rJ = model.fitTo(dT, Minimizer("Minuit2"), Save(true));
	  //	  RooFitResult *rJ = gaussB.fitTo(dT, Minimizer("Minuit2"), Save(true));
	  model.plotOn(frame);

	  float chi2_J = frame->chiSquare(6);
	  std::cout << ">>> chi2_S = " << chi2_S << " chi2_J  = " << chi2_J << std::endl;

	  bool bkgOnly = false;
	  if(chi2_S < 1.21 || chi2_S < chi2_J) bkgOnly = true;
	  bkgOnly = true;

	  if(!bkgOnly){
	    model.plotOn(frame, Components("gauss"),LineColor(kRed+1));
	    model.paramOn(frame, RooFit::Layout(0.6,0.8,0.9),RooFit::Format("NEA",AutoPrecision(1)));
	    frame->getAttLine()->SetLineColorAlpha(kWhite, 0.2);
	    frame->getAttText()->SetTextSize(0.03);
	    frame->getAttText()->SetTextFont(42);
	  }
	  else{
	    modelCBS.plotOn(frameS);
	    //modelCBS.plotOn(frameS, RooFit::Layout(0.6,0.8,0.9),RooFit::Format("NEA",AutoPrecision(1)));
	    // frameS->getAttLine()->SetLineColorAlpha(kWhite, 0.2);
            // frameS->getAttText()->SetTextSize(0.03);
            // frameS->getAttText()->SetTextFont(42);
	  }
	  std::cout << " >>> qui passato " << std::endl;

	  if(!bkgOnly){
	    RooRealVar* parS_J = (RooRealVar*) rJ->floatParsFinal().find("nsig");
	    auto nEv_Sig = parS_J->getValV();
	    auto nEvError_Sig = parS_J->getError();
	    
	    RooRealVar* parMeanJ = (RooRealVar*) rJ->floatParsFinal().find("mean");
	    RooRealVar* parSigmaJ = (RooRealVar*) rJ->floatParsFinal().find("sigma");
	    float meanValJ = parMeanJ->getValV();
	    float sigmaValJ = parSigmaJ->getValV();
	    
	    
	    std::cout << "\n  iL = " << iL << " nEvt = " << nEv_Sig << "+/-" << nEvError_Sig 
		      << " meanSig = " << meanValJ << "+/-" << parMeanJ->getError() 
		      << " sigmaSig = " << sigmaValJ << "+/-" << parSigmaJ->getError() << std::endl;

	    // if(iL == 1){
	    //   tgM[ij][iRad]->SetPoint(0, 0, -1);
	    //   tgFrEvt[ij][iRad]->SetPoint(0, 0, 0);
	    // }
	    tgM[ij][iRad]->SetPoint(iL+1, iL+1 + 0.1*(ij - 3), meanValJ);
	    tgM[ij][iRad]->SetPointError(iL+1, 0, sigmaValJ);

	    tgFrEvt[ij][iRad]->SetPoint(iL+1, iL+1 + 0.1*(ij - 3), slice->GetEntries());
	    //tgFrEvt[ij][iRad]->SetPoint(iL, iL + 0.1*(ij - 3), meanValJ);
	    //	    tgFrEvt[ij][iRad]->SetPointError(iL, 0, 0);
	    std::cout << " >>> slice->GetEntries() = " << slice->GetEntries() 
		      << " iL = " << iL << " ij = " << ij << " \n "<< std::endl;
	  }
	  
	  TLatex tL1;
	  tL1.SetNDC();
	  tL1.SetTextSize(0.03);
	  tL1.SetTextFont(132);
	  TLatex tL2;
	  tL2.SetNDC();
	  tL2.SetTextSize(0.03);
	  tL2.SetTextFont(132);

	  c1 = new TCanvas();
	  c1->cd();
	  if(!bkgOnly) frame->Draw();
	  else frameS->Draw();
	  tL1.DrawLatex(0.2, 0.85, Form("chi2bkg = %.2f", chi2_S));
	  tL2.DrawLatex(0.2, 0.80, Form("chi2sig = %.2f", chi2_J));
	  c1->Print(Form("plots/time_iL%d_iS%d_iE%.2f.png", iL+1, iRad, energyThreshold.at(ij)), "png");
	  c1->Delete();
	  h2_timeHitsEnergy_layer->Delete();

	}//side
      }//side
    }//E threshold
  }//loop over layers



  TLegend *legTGM = new TLegend(0.65,0.70,0.85,0.95,NULL,"brNDC");
  legTGM->SetTextFont(42);
  legTGM->SetTextSize(0.04);
  legTGM->SetFillColor(kWhite);
  legTGM->SetLineColor(kWhite);
  legTGM->SetShadowColor(kWhite);
  for(int iF=0; iF<energyThreshold.size(); ++iF){
    legTGM->AddEntry(hDummy[iF], Form("E > %.1f GeV", energyThreshold.at(iF)), "l");
  }


  TLegend *legTGM2s = new TLegend(0.40,0.70,0.60,0.85,NULL,"brNDC");
  legTGM2s->SetTextFont(42);
  legTGM2s->SetTextSize(0.05);
  legTGM2s->SetFillColor(kWhite);
  legTGM2s->SetLineColor(kWhite);
  legTGM2s->SetShadowColor(kWhite);
  legTGM2s->AddEntry(hDummySt[0], " z < 0", "p");
  legTGM2s->AddEntry(hDummySt[1], " z > 0", "p");


  TCanvas* ch_M_zpos = new TCanvas();
  ch_M_zpos->cd();
  tgM[0][0]->GetXaxis()->SetTitle("n layer (z > 0)");
  tgM[0][0]->GetYaxis()->SetTitle("mean(t) ns");
  tgM[0][0]->GetYaxis()->SetRangeUser(-0.5, 2.);
  tgM[0][0]->GetXaxis()->SetRangeUser(0., 50.);
  tgM[0][0]->Draw("ap");
  //  tgM[0][1]->Draw("p, same");
  for(int iF=1; iF<energyThreshold.size(); ++iF){
    tgM[iF][0]->Draw("p, same");
    //    tgM[iF][1]->Draw("p, same");
  }
  legTGM->Draw("same");
  //  legTGM2s->Draw("same");
  ch_M_zpos->Print("plotsTrend/h_timeMeanVsLayer_zpos_zoomOut.png", "png");
  ch_M_zpos->Print("plotsTrend/h_timeMeanVsLayer_zpos.pdf", "pdf");  
  ch_M_zpos->Print("plotsTrend/h_timeMeanVsLayer_zpos.root", "root");
  tgM[0][0]->GetYaxis()->SetRangeUser(-0.1, 0.1);
  ch_M_zpos->Print("plotsTrend/h_timeMeanVsLayer_zpos_zoomOut.png", "png");
  ch_M_zpos->Print("plotsTrend/h_timeMeanVsLayer_zpos_zoomOut.pdf", "pdf");  
  ch_M_zpos->Print("plotsTrend/h_timeMeanVsLayer_zpos_zoomOut.root", "root");
  tgM[0][0]->GetYaxis()->SetRangeUser(-0.1, 0.4);
  ch_M_zpos->Print("plotsTrend/h_timeMeanVsLayer_zpos.png", "png");
  ch_M_zpos->Print("plotsTrend/h_timeMeanVsLayer_zpos.pdf", "pdf");  
  ch_M_zpos->Print("plotsTrend/h_timeMeanVsLayer_zpos.root", "root");
  tgM[0][0]->GetXaxis()->SetRangeUser(0.5, 28.5);
  ch_M_zpos->Print("plotsTrend/h_timeMeanVsLayer_zposEE.png", "png");
  ch_M_zpos->Print("plotsTrend/h_timeMeanVsLayer_zposEE.pdf", "pdf");  
  ch_M_zpos->Print("plotsTrend/h_timeMeanVsLayer_zposEE.root", "root");
  tgM[0][0]->GetXaxis()->SetRangeUser(28.5, 41.5);
  ch_M_zpos->Print("plotsTrend/h_timeMeanVsLayer_zposFH.png", "png");
  ch_M_zpos->Print("plotsTrend/h_timeMeanVsLayer_zposFH.pdf", "pdf");  
  ch_M_zpos->Print("plotsTrend/h_timeMeanVsLayer_zposFH.root", "root");


  TCanvas* ch_M_zneg = new TCanvas();
  ch_M_zneg->cd();
  tgM[0][1]->GetXaxis()->SetTitle("n layer (z > 0)");
  tgM[0][1]->GetYaxis()->SetTitle("mean(t) ns");
  tgM[0][1]->GetYaxis()->SetRangeUser(-0.5, 2.);
  tgM[0][1]->GetXaxis()->SetRangeUser(0., 50.);
  tgM[0][1]->Draw("ap");
  //  tgM[0][1]->Draw("p, same");
  for(int iF=1; iF<energyThreshold.size(); ++iF){
    tgM[iF][1]->Draw("p, same");
    //    tgM[iF][1]->Draw("p, same");
  }
  legTGM->Draw("same");
  //  legTGM2s->Draw("same");
  ch_M_zneg->Print("plotsTrend/h_timeMeanVsLayer_zneg_zoomOut.png", "png");
  ch_M_zneg->Print("plotsTrend/h_timeMeanVsLayer_zneg_zoomOut.pdf", "pdf");  
  ch_M_zneg->Print("plotsTrend/h_timeMeanVsLayer_zneg_zoomOut.root", "root");
  tgM[0][1]->GetYaxis()->SetRangeUser(-0.1, 0.1);
  ch_M_zneg->Print("plotsTrend/h_timeMeanVsLayer_zneg_zoomIn.png", "png");
  ch_M_zneg->Print("plotsTrend/h_timeMeanVsLayer_zneg_zoomIn.pdf", "pdf");  
  ch_M_zneg->Print("plotsTrend/h_timeMeanVsLayer_zneg_zoomIn.root", "root");
  tgM[0][1]->GetYaxis()->SetRangeUser(-0.1, 0.4);
  ch_M_zneg->Print("plotsTrend/h_timeMeanVsLayer_zneg.png", "png");
  ch_M_zneg->Print("plotsTrend/h_timeMeanVsLayer_zneg.pdf", "pdf");  
  ch_M_zneg->Print("plotsTrend/h_timeMeanVsLayer_zneg.root", "root");
  tgM[0][1]->GetXaxis()->SetRangeUser(0.5, 28.5);
  ch_M_zneg->Print("plotsTrend/h_timeMeanVsLayer_znegEE.png", "png");
  ch_M_zneg->Print("plotsTrend/h_timeMeanVsLayer_znegEE.pdf", "pdf");  
  ch_M_zneg->Print("plotsTrend/h_timeMeanVsLayer_znegEE.root", "root");
  tgM[0][1]->GetXaxis()->SetRangeUser(28.5, 41.5);
  ch_M_zneg->Print("plotsTrend/h_timeMeanVsLayer_znegFH.png", "png");
  ch_M_zneg->Print("plotsTrend/h_timeMeanVsLayer_znegFH.pdf", "pdf");  
  ch_M_zneg->Print("plotsTrend/h_timeMeanVsLayer_znegFH.root", "root");


  /////////

  TCanvas* ch_fr_zpos = new TCanvas();
  gPad->SetLogy();
  ch_fr_zpos->cd();
  tgFrEvt[0][0]->GetXaxis()->SetTitle("n layer (z > 0)");
  tgFrEvt[0][0]->GetYaxis()->SetTitle("nHits");
  tgFrEvt[0][0]->GetYaxis()->SetRangeUser(1.e3, 1.e8);
  tgFrEvt[0][0]->GetXaxis()->SetRangeUser(0., 50.);
  tgFrEvt[0][0]->Draw("ap");
  //  tgM[0][1]->Draw("p, same");
  for(int iF=1; iF<energyThreshold.size(); ++iF){
    tgFrEvt[iF][0]->Draw("p, same");
    //    tgM[iF][1]->Draw("p, same");
  }
  legTGM->Draw("same");
  //  legTGM2s->Draw("same");
  ch_fr_zpos->Print("plotsTrend/h_fractionHits_zpos.png", "png");
  ch_fr_zpos->Print("plotsTrend/h_fractionHits_zpos.pdf", "pdf");  
  ch_fr_zpos->Print("plotsTrend/h_fractionHits_zpos.root", "root");



  /////
  TCanvas* ch_fr_zneg = new TCanvas();
  gPad->SetLogy();
  ch_fr_zneg->cd();
  tgFrEvt[0][0]->GetXaxis()->SetTitle("n layer (z < 0)");
  tgFrEvt[0][0]->GetYaxis()->SetTitle("nHits");
  tgFrEvt[0][0]->GetYaxis()->SetRangeUser(1.e3, 1.e8);
  tgFrEvt[0][0]->GetXaxis()->SetRangeUser(0., 50.);
  tgFrEvt[0][0]->Draw("ap");
  //  tgM[0][1]->Draw("p, same");
  for(int iF=1; iF<energyThreshold.size(); ++iF){
    tgFrEvt[iF][0]->Draw("p, same");
    //    tgM[iF][1]->Draw("p, same");
  }
  legTGM->Draw("same");
  //  legTGM2s->Draw("same");
  ch_fr_zneg->Print("plotsTrend/h_fractionHits_zneg.png", "png");
  ch_fr_zneg->Print("plotsTrend/h_fractionHits_zneg.pdf", "pdf");  
  ch_fr_zneg->Print("plotsTrend/h_fractionHits_zneg.root", "root");




  return;
}
