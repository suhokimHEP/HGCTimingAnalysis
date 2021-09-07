#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>

#include "TROOT.h"
#include "TSystem.h"
#include "TKey.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TString.h"
#include "TCut.h"
#include "TMath.h"
#include "TApplication.h"
#include "TError.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TPad.h"
#include "TMultiGraph.h"
#include "TStyle.h"
#include "TChain.h"

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
using namespace ROOT::VecOps;

void fitMIP(){

  gROOT->Reset();
  gROOT->Macro("~/public/setStyle.C");

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptTitle(0);
 
  RooWorkspace w("w");    
  RooRealVar x("x", "", 0.5, 10.5);
  //  RooDataSet data ("data", "data", RooArgSet(x));

  std::string inputFileList = "outMIP_SoN4.root";
  ROOT::RDataFrame d("MIPtree", inputFileList.c_str());
    //ROOT::RDataFrame d("MIPtree", "outMIP_SoN4.root");
  auto allMIP = d.Histo1D( {"allMIP", "", 120, 0.5, 10.5}, "MIP_val");

  // TCanvas* allMIP_tc = new TCanvas();
  // allMIP_tc->cd();
  // allMIP->GetXaxis()->SetTitle("MIP values");
  // allMIP->Draw();

  RooDataHist dataHist("dataHist","dataHist", x, Import(*(allMIP)) ); 
  
  w.import(x);
  w.factory("nNoise[10, 0, 1.e4]"); 
  w.factory("nSig[1.e5, 0.0, 1.e6]");   

  w.factory("RooLandau::landau(x, m_landau[1.1, 0.5, 5.], s_landau[0.2, 0.01, 0.2])");
  w.factory("RooGaussian::gauss(x, m_gauss[0.0, 0.0, 0.], s_gauss[0.1, 0.01, 0.5])");
  w.factory("RooFFTConvPdf::lxg(x, landau, gauss)");

  w.factory("RooGaussian::gaussB(x, m_gaussNoise[0., 0., 0.], s_gaussNoise[0.01, 0.001, 0.1])");
  w.factory("SUM::model(nNoise * gaussB, nSig* lxg)");
  //w.factory("SUM::model(nBkg * gaussB, nSig* landau)");
  RooAbsPdf * model = w.pdf("model");


  RooFitResult * rJ = model->fitTo(dataHist, Minimizer("Minuit2"),Save(true));

  RooPlot* frame = x.frame();
  frame->SetXTitle("MIP values");
  dataHist.plotOn(frame, Binning(120));

  model->plotOn(frame, LineColor(kBlue));
  model->plotOn(frame, Components("lxg"), LineColor(kRed+1));
  model->paramOn(frame, RooFit::Layout(0.6,0.8,0.9),RooFit::Format("NEA",AutoPrecision(1)));
  frame->getAttLine()->SetLineColorAlpha(kWhite, 0.2);
  frame->getAttText()->SetTextSize(0.03);
  frame->getAttText()->SetTextFont(42);

  TCanvas* c = new TCanvas();
  c->cd();
  frame->Draw();
  c->Print("mip_allDetector.png", "png");

  float chi2_J = frame->chiSquare();

  std::cout << " chi2_J  = " << chi2_J << std::endl;

  ///single tile fits
  int firstLayer = 37;

  RooRealVar xL("xL", "", 0.5, 10.5);
  w.import(xL);
  w.factory("nNoise_l[10, 0, 1.e4]"); 
  w.factory("nSig_l[1.e5, 0.0, 1.e6]");   

  w.factory("RooLandau::landauL(xL, m_landau_l[1.1, 0.5, 5.], s_landau_l[0.2, 0.01, 0.2])");
  //  w.var("m_gauss")->setConstant();
  //  w.var("s_gauss")->setConstant();
  w.var("m_gaussNoise")->setConstant();
  w.var("s_gaussNoise")->setConstant();
  w.factory("RooGaussian::gaussL(xL, m_gauss_l[0.0, 0.0, 0.], s_gauss_l[0.1, 0.01, 0.5])");
  //  w.factory("RooGaussian::gaussL(xL, m_gauss, s_gauss)");
  w.factory("RooFFTConvPdf::lxgL(xL, landauL, gaussL)");
  w.factory("RooGaussian::gaussBL(xL, m_gaussNoise, s_gaussNoise)");
  w.factory("SUM::modelL(nNoise_l * gaussBL, nSig_l* lxgL)");
  RooAbsPdf * modelL = w.pdf("modelL");

  TH1F* MIP_values_layer = new TH1F("MIP_values_layer", "", 500, 0., 5.);
  TH1F* MIP_values_layer_ring = new TH1F("MIP_values_layer_ring", "", 500, 0., 5.);
  TGraph* MIPerror_vs_nEntries_layer_ring = new TGraph();
  TH1F* MIP_values_layer_ring_phi = new TH1F("MIP_values_layer_ring_phi", "", 500, 0., 5.);
  TH2F* iR_vs_layer = new TH2F("iR_vs_layer", "", 14, 0., 13, 75, -37, 37.);

  auto min_iL = d.Min("MIP_layer").GetValue();
  auto max_iL = d.Max("MIP_layer").GetValue();
  auto min_iR_all = d.Min("MIP_iR").GetValue();
  auto max_iR_all = d.Max("MIP_iR").GetValue();
  auto min_iPhi_all = d.Min("MIP_iPhi").GetValue();
  auto max_iPhi_all = d.Max("MIP_iPhi").GetValue();

  std::cout << " min_iR = " << min_iR_all << " max_iR = " << max_iR_all 
	    << " min_iL = " << min_iL << " max_iL = " << max_iL 
	    << " min_iPhi = " << min_iPhi_all << " max_iPhi = " << max_iPhi_all << std::endl;



  for(int ij_L=min_iL; ij_L < max_iL; ++ij_L){
    std::cout << " fitting for layer = " << ij_L << std::endl;

    auto filterMinMAx = [ij_L](int lval) { return (ij_L == lval); };
    auto radiiForLayer = d.Filter(filterMinMAx, {"MIP_layer"});
    auto min_iR = radiiForLayer.Min("MIP_iR").GetValue();
    auto max_iR = radiiForLayer.Max("MIP_iR").GetValue();

    for(int ij_R=min_iR; ij_R < max_iR; ++ij_R){
      std::cout << " fitting for layer = " << ij_L << " ring = " << ij_R << std::endl;
      //auto filterSel = [ij_L, ij_R](int lval, int rval) { return (ij_L == lval && ij_R == rval); };                                                
      //auto mipH = d.Filter(filterSel, {"MIP_layer", "MIP_iR"}).Histo1D({"mipH", "", 120, 0.5, 10.5}, "MIP_val"); 
      // iR_vs_layer->Fill(ij_L, ij_R, mipH->GetEntries());


      auto filterMinMAx = [ij_L, ij_R](int lval, int rval) { return (ij_L == lval && ij_R == rval); };
      auto phiForLayerRadii = d.Filter(filterMinMAx, {"MIP_layer", "MIP_iR"});
      auto min_iPhi = phiForLayerRadii.Min("MIP_iPhi").GetValue();
      auto max_iPhi = phiForLayerRadii.Max("MIP_iPhi").GetValue();

      /*
      for(int ij_P=min_iPhi; ij_P < max_iPhi; ++ij_P){	
	std::cout << " fitting for layer = " << ij_L << " ring = " << ij_R << " phi = " << ij_P << std::endl;

	//RA FIXME LAMBDA
	auto filterSel = [ij_P, ij_L, ij_R](int pval, int lval, int rval) { return (ij_P == pval && ij_L == lval && ij_R == rval); }; 
	auto mipH = d.Filter(filterSel, {"MIP_iPhi", "MIP_layer", "MIP_iR"}).Histo1D({"mipH", "", 120, 0.5, 10.5}, "MIP_val");
	if(mipH->GetEntries() == 0) continue;

	RooDataHist locData("locData","locData", xL, Import(*(mipH)) );

	//fit and save
	RooFitResult * rJL = modelL->fitTo(locData, Minimizer("Minuit2"),Save(true));
	auto fittedVal = w.var("m_landau_l")->getVal();
	MIP_values_layer_ring_phi->Fill(fittedVal);


	RooPlot* frameL = xL.frame();
	frameL->SetXTitle("MIP values");
	locData.plotOn(frameL, Binning(120));
	modelL->plotOn(frameL, LineColor(kBlue));
	modelL->plotOn(frameL, Components("lxgL"), LineColor(kRed+1));
	modelL->paramOn(frameL, RooFit::Layout(0.6,0.8,0.9),RooFit::Format("NEA",AutoPrecision(1)));
	frameL->getAttLine()->SetLineColorAlpha(kWhite, 0.2);
	frameL->getAttText()->SetTextSize(0.03);
	frameL->getAttText()->SetTextFont(42);

	TCanvas* cL = new TCanvas();
	cL->cd();
	frameL->Draw();
	cL->Print("Form(singleFits/mip_L%d_R%d_P%d.png, ij_L+firstLayer, ij_R, ij_P)", "png");

	float chi2_JL = frameL->chiSquare();

	std::cout << " chi2_JL  = " << chi2_JL << std::endl;

	std::cout << "\n ij_L = " << ij_L << " ij_R = " << ij_R << " ij_P = " << ij_P << std::endl;
      }
      */
      
      auto filterSel = [ij_L, ij_R](int lval, int rval) { return (ij_L == lval && ij_R == rval); };                            
      auto mipH = d.Filter(filterSel, {"MIP_layer", "MIP_iR"}).Histo1D({"mipH", "", 120, 0.5, 10.5}, "MIP_val");
      float nEntries = mipH->GetEntries();
      if(nEntries == 0) continue;     
      RooDataHist locData("locData","locData", xL, Import(*(mipH)) );

      //fit and save
      RooFitResult * rJL = modelL->fitTo(locData, Minimizer("Minuit2"),Save(true));
      auto fittedVal = w.var("m_landau_l")->getVal();
      MIP_values_layer_ring->Fill(fittedVal);
      auto fittedValError = w.var("s_landau_l")->getVal();
      MIPerror_vs_nEntries_layer_ring->SetPoint(MIPerror_vs_nEntries_layer_ring->GetN(), nEntries, fittedValError);

      RooPlot* frameL = xL.frame();                                                                                                                            
      frameL->SetXTitle("MIP values");                                                                                                                         
      locData.plotOn(frameL, Binning(120));                                                                                                                    
      modelL->plotOn(frameL, LineColor(kBlue));                                                                                                                
      modelL->plotOn(frameL, Components("lxgL"), LineColor(kRed+1));                                                                                           
      modelL->paramOn(frameL, RooFit::Layout(0.6,0.8,0.9),RooFit::Format("NEA",AutoPrecision(1)));                                                             
      frameL->getAttLine()->SetLineColorAlpha(kWhite, 0.2);                                                                                                    
      frameL->getAttText()->SetTextSize(0.03);                                                                                                                 
      frameL->getAttText()->SetTextFont(42);                                                                                                                   

      TCanvas* cL = new TCanvas();                                                                                                                             
      cL->cd();                                                                                                                                                
      frameL->Draw();                                                                                    
      cL->Print(Form("singleFits/mip_L%d_R%d.png", ij_L+firstLayer, ij_R), ".png"); 
      
      delete cL;
      delete frameL;
      //delete mipH;
    }
    /*
    auto filterSel = [ij_L](int lval) { return (ij_L == lval); };  
    auto mipH = d.Filter(filterSel, {"MIP_layer"}).Histo1D({"mipH", "", 120, 0.5, 10.5}, "MIP_val"); 
    RooDataHist locData("locData","locData", xL, Import(*(mipH)) ); 
    
 
    //fit and save  
    RooFitResult * rJL = modelL->fitTo(locData, Minimizer("Minuit2"),Save(true)); 
    auto fittedVal = w.var("m_landau_l")->getVal();
    std::cout << " fittedVal = " << fittedVal << std::endl;
    MIP_values_layer->Fill(fittedVal);
    
    RooPlot* frameL = xL.frame();
    frameL->SetXTitle("MIP values");
    locData.plotOn(frameL, Binning(120));
    modelL->plotOn(frameL, LineColor(kBlue)); 
    modelL->plotOn(frameL, Components("lxgL"), LineColor(kRed+1)); 
    modelL->paramOn(frameL, RooFit::Layout(0.6,0.8,0.9),RooFit::Format("NEA",AutoPrecision(1))); 
    frameL->getAttLine()->SetLineColorAlpha(kWhite, 0.2); 
    frameL->getAttText()->SetTextSize(0.03); 
    frameL->getAttText()->SetTextFont(42); 
    
    TCanvas* cL = new TCanvas(); 
    cL->cd();
    frameL->Draw();
    cL->Print(Form("singleFits/mip_L%d.png", ij_L+firstLayer), ".png");
    */
  }

  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1);
  TF1* hfit_layer = new TF1("hfit_layer", "gaus", 0., 2.);
  TCanvas* tc_fitVal = new TCanvas();
  tc_fitVal->cd();
  /*
  hfit_layer->SetParameters(MIP_values_layer->GetEntries()/2., 1, 0.02);
  MIP_values_layer->GetXaxis()->SetTitle("MPV values for layers");
  MIP_values_layer->Draw();
  MIP_values_layer->Fit("hfit_layer", "R");
  tc_fitVal->Print("MIP_values_layer.png", "png");
  */

  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1);
  
  hfit_layer->SetParameters(MIP_values_layer_ring->GetEntries()/2., 1, 0.02);
  MIP_values_layer_ring->GetXaxis()->SetTitle("MPV values for layers and rings");
  MIP_values_layer_ring->Draw();
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1);
  MIP_values_layer_ring->Fit("hfit_layer", "R");
  tc_fitVal->Print("MIP_values_layer_ring.png", "png");
  /*
  hfit_layer->SetParameters(MIP_values_layer_ring_phi->GetEntries()/2., 1, 0.02);
  MIP_values_layer_ring_phi->GetXaxis()->SetTitle("MPV values for tile");
  MIP_values_layer_ring_phi->Draw();
  MIP_values_layer_ring_phi->Fit("hfit_layer", "R");
  tc_fitVal->Print("MIP_values_layer_ring_phi.png", "png");
  */

  TCanvas* tc_scatter = new TCanvas();
  tc_scatter->cd();
  MIPerror_vs_nEntries_layer_ring->Sort();
  MIPerror_vs_nEntries_layer_ring->GetXaxis()->SetTilte("nEntries");
  MIPerror_vs_nEntries_layer_ring->GetXaxis()->SetTilte("precision");
  MIPerror_vs_nEntries_layer_ring->Draw("ap");
  tc_scatter->Print("MIPerror_vs_nEntries_layer_ring.ring", "png");

  /*
  TCanvas* tc_occu = new TCanvas();
  tc_occu->cd();
  iR_vs_layer->Draw("colz");
  */

  /*
  TFile* inF = TFile::Open(inputFileList.c_str());
  TTree* tree = (TTree*)inF->Get("MIPtree");

  Float_t MIP_val;
  Int_t MIP_layer;
  Int_t MIP_iR;
  Int_t MIP_iPhi;
  Float_t MIP_SoN;
  tree->SetBranchAddress("MIP_layer", &MIP_layer);
  tree->SetBranchAddress("MIP_iR", &MIP_iR);
  tree->SetBranchAddress("MIP_iPhi", &MIP_iPhi);
  tree->SetBranchAddress("MIP_SoN", &MIP_SoN);
  tree->SetBranchAddress("MIP_val", &MIP_val);
  auto totN = tree->GetEntriesFast();
  */
}
