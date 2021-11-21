//g++  -o fitMIP  fitMIP.cpp `root-config --cflags --glibs`

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
#include "TProfile.h"
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
#include "TPaveStats.h"

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <iostream>
#include <string>
#include "TStopwatch.h"
#include <vector>
#include <ROOT/RDF/RInterface.hxx>
#include <ROOT/RDF/InterfaceUtils.hxx>

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
//using namespace ROOT::VecOps;

void radfitMIP(){
  //int main(){
  gROOT->Reset();
  gROOT->Macro("./setStyle.C");

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptTitle(0);
 
  RooWorkspace w("w");    
  RooRealVar x("x", "", 0.5, 10.5);
  //  RooDataSet data ("data", "data", RooArgSet(x));

  std::string inputFileList = "new.root";
 TFile* inF = TFile::Open(TString(inputFileList));
 TTree* newT = (TTree*)inF->Get("MIPtree");
  std::vector<int> Lentries;
Lentries.clear();
  for(int Rn=16; Rn <= 38; ++Rn){
        int entries = newT->GetEntries(Form("MIP_iR==%d",Rn));
	Lentries.push_back(entries);
	printf("%d:\n",entries);
	}
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

  w.factory("RooLandau::landau(x, m_landau[1.1, 0.5, 5.], s_landau[0.2, 0.01, 0.5])");
  w.factory("RooGaussian::gauss(x, m_gauss[0.0, 0.0, 0.], s_gauss[0.1, 0.01, 1.])");
  w.factory("RooFFTConvPdf::lxg(x, landau, gauss)");

  w.factory("RooGaussian::gaussB(x, m_gaussNoise[0., 0., 0.], s_gaussNoise[0.5, 0.01, 2.])");
  w.factory("SUM::model(nNoise * gaussB, nSig* lxg)");
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

  w.factory("RooLandau::landauL(xL, m_landau_l[1.1, 0.5, 5.], s_landau_l[0.2, 0.01, 0.5])");
  w.var("m_gaussNoise")->setConstant();
  w.var("s_gaussNoise")->setConstant();
  w.factory("RooGaussian::gaussL(xL, m_gauss_l[0.0, 0.0, 0.], s_gauss_l[0.1, 0.01, 1.])");
  w.factory("RooFFTConvPdf::lxgL(xL, landauL, gaussL)");
  w.factory("RooGaussian::gaussBL(xL, m_gaussNoise, s_gaussNoise)");
  w.factory("SUM::modelL(nNoise_l * gaussBL, nSig_l* lxgL)");
  RooAbsPdf * modelL = w.pdf("modelL");

  TH1F* MIP_values_layer = new TH1F("MIP_values_layer", "", 500, 0., 5.);
  TH1F* MIP_values_layer_ring = new TH1F("MIP_values_layer_ring", "", 500, 0., 5.);
  TGraph* MIPerror_vs_nEntries_layer_ring = new TGraph();
  TProfile* tp_MIPerror_vs_nEntries_layer_ring = new TProfile("tp_MIPerror_vs_nEntries_layer_ring", "", 400, 0., 400);
  TH1F* MIP_values_layer_ring_phi = new TH1F("MIP_values_layer_ring_phi", "", 500, 0., 5.);
  TH2F* iR_vs_layer = new TH2F("iR_vs_layer", "", 14, 0., 13, 75, -37, 37.);

  auto min_iL = d.Min("MIP_layer").GetValue();
  auto max_iL = d.Max("MIP_layer").GetValue();
   auto min_iR_all = d.Min("MIP_iR").GetValue();
   auto max_iR_all = d.Max("MIP_iR").GetValue();
  // auto min_iPhi_all = d.Min("MIP_iPhi").GetValue();
  // auto max_iPhi_all = d.Max("MIP_iPhi").GetValue();

  // std::cout << " min_iR = " << min_iR_all << " max_iR = " << max_iR_all 
  // 	    << " min_iL = " << min_iL << " max_iL = " << max_iL 
  // 	    << " min_iPhi = " << min_iPhi_all << " max_iPhi = " << max_iPhi_all << std::endl;


  //std::vector<int> minEvents;
  //minEvents.push_back(2500);
  //minEvents.push_back(1000);
  //int minESize = minEvents.size();

 // for(int ij_L=min_iL; ij_L <= max_iL; ++ij_L){
 // //  for(int ij_L=13; ij_L < max_iL; ++ij_L){
 //   std::cout << " fitting for layer = " << ij_L << std::endl;

 //   auto filterMinMAx = [ij_L](int lval) { return (ij_L == lval); };
 //   auto radiiForLayer = d.Filter(filterMinMAx, {"MIP_layer"});
//    auto absRadiiForLayer = radiiForLayer.Define("abs_MIP_iR", "std::abs(MIP_iR)");
    auto absRadiiForLayer = d.Define("abs_MIP_iR", "std::abs(MIP_iR)");
    auto min_iR = absRadiiForLayer.Min("abs_MIP_iR").GetValue();
    auto max_iR = absRadiiForLayer.Max("abs_MIP_iR").GetValue();

    for(int ij_R=16; ij_R <= 38; ++ij_R){
    //  std::cout << "\n  ==>> Fitting for layer = " << ij_L << " ring = " << ij_R << std::endl;


      auto filterSel = [ij_R](int rval) { return (ij_R == rval); };                            
    //  auto filterSel = [ij_L, ij_R](int lval, int rval) { return (ij_L == lval && ij_R == rval); };                            
      auto filterTree = d.Filter(filterSel, {"MIP_iR"});
    //  auto filterTree = d.Filter(filterSel, {"MIP_layer", "MIP_iR"});

      //for(int minE = 0; minE < minESize; ++minE){
	int numberOfE = Lentries.at(ij_R-16);
	//auto mipH = radiiForLayer.Range(0, numberOfE).Histo1D({"mipH", "", 120, 0.5, 10.5}, "MIP_val");
	auto mipH = filterTree.Range(0, numberOfE).Histo1D({"mipH", "", 120, 0.5, 10.5}, "MIP_val");

      float nEntries = mipH->GetEntries();
      //if(nEntries == 0) continue;     
      if(nEntries < numberOfE) continue;     
      RooDataHist locData("locData","locData", xL, Import(*(mipH)) );

      //fit and save
      RooFitResult * rJL = modelL->fitTo(locData, Minimizer("Minuit2"),Save(true));
      auto fittedVal = w.var("m_landau_l")->getVal();
      MIP_values_layer_ring->Fill(fittedVal);
      auto fittedValError = w.var("m_landau_l")->getError();
      MIPerror_vs_nEntries_layer_ring->SetPoint(MIPerror_vs_nEntries_layer_ring->GetN(), nEntries, fittedValError*100.);
      tp_MIPerror_vs_nEntries_layer_ring->Fill(numberOfE, fittedValError*100.);

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
      cL->Print(Form("singleFits/mip_R%d_nEvts%d.png", ij_R, numberOfE), ".png"); 
      //cL->Print(Form("singleFits/mip_L%d_R%d_nEvts%d.png", ij_L+firstLayer, ij_R, numberOfE), ".png"); 
      
      delete cL;
      delete frameL;
      //delete mipH;
      //}
      //}
  }


  TLatex tL;
  tL.SetNDC();
  tL.SetTextSize(0.05);
  tL.SetTextFont(132);

  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1);
  TF1* hfit_layer = new TF1("hfit_layer", "gaus", 0., 2.);

  TCanvas* tc_fitVal = new TCanvas();
  tc_fitVal->cd();
  hfit_layer->SetParameters(MIP_values_layer_ring->GetEntries()/2., 1, 0.02);
  MIP_values_layer_ring->GetXaxis()->SetTitle("MPV values for layers and rings");
  MIP_values_layer_ring->Draw();
  MIP_values_layer_ring->Fit("hfit_layer", "q");
  //MIP_values_layer_ring->Fit("hfit_layer", "R");
  gPad->Update();
  tL.DrawLatex(0.5, 0.8, Form("mean = %.2f +/- %.2f", hfit_layer->GetParameter(1), hfit_layer->GetParError(1)));
  tL.DrawLatex(0.5, 0.7, Form("sigma = %.2f +/- %.2f", hfit_layer->GetParameter(2), hfit_layer->GetParError(2)));
  tc_fitVal->Print("MIP_values_layer_ring.png", "png");

  TF1* hfit_precision = new TF1("hfit_precision", "[0]/sqrt(x) + [1]", 0., 500.);
  hfit_precision->SetParameter(1, 0.01);
  hfit_precision->SetParameter(0, 0.5);

  TCanvas* tc_scatter = new TCanvas();
  tc_scatter->cd();
  MIPerror_vs_nEntries_layer_ring->Sort();
  MIPerror_vs_nEntries_layer_ring->GetXaxis()->SetTitle("nEntries");
  MIPerror_vs_nEntries_layer_ring->GetYaxis()->SetTitle("precision (%)");
  MIPerror_vs_nEntries_layer_ring->GetYaxis()->SetRangeUser(0., 100.);
  MIPerror_vs_nEntries_layer_ring->Draw("ap");
  MIPerror_vs_nEntries_layer_ring->Fit("hfit_precision","q");
  //MIPerror_vs_nEntries_layer_ring->Fit("hfit_precision");
  gPad->Update();
  tL.DrawLatex(0.5, 0.8, "a / sqrt(x) + b");
  tL.DrawLatex(0.5, 0.7, Form("a = %.2f +/- %.2f", hfit_precision->GetParameter(0), hfit_precision->GetParError(0)));
  tL.DrawLatex(0.5, 0.6, Form("b = %.2f +/- %.2f", hfit_precision->GetParameter(1), hfit_precision->GetParError(1)));
  tc_scatter->Print("MIPerror_vs_nEntries_layer_ring.png", "png");


  tc_scatter->cd();
  tp_MIPerror_vs_nEntries_layer_ring->GetXaxis()->SetTitle("nEntries");
  tp_MIPerror_vs_nEntries_layer_ring->GetYaxis()->SetTitle("precision (%)");
  tp_MIPerror_vs_nEntries_layer_ring->GetYaxis()->SetRangeUser(0., 100.);
  tp_MIPerror_vs_nEntries_layer_ring->SetMarkerStyle(20);
  tp_MIPerror_vs_nEntries_layer_ring->SetLineWidth(2);
  tp_MIPerror_vs_nEntries_layer_ring->Draw("e");
  hfit_precision->SetParameter(1, 0.01);
  hfit_precision->SetParameter(0, 0.5);
  tp_MIPerror_vs_nEntries_layer_ring->Fit("hfit_precision","q");
  //tp_MIPerror_vs_nEntries_layer_ring->Fit("hfit_precision");
  gPad->Update();
  tL.DrawLatex(0.5, 0.8, "a / sqrt(x) + b");
  tL.DrawLatex(0.5, 0.7, Form("a = %.2f +/- %.2f", hfit_precision->GetParameter(0), hfit_precision->GetParError(0)));
  tL.DrawLatex(0.5, 0.6, Form("b = %.2f +/- %.2f", hfit_precision->GetParameter(1), hfit_precision->GetParError(1)));
  tL.DrawLatex(0.5, 0.5, Form("#Chi^{2} = %.2f", hfit_precision->GetChisquare()/hfit_precision->GetNDF()));
  tc_scatter->Print("tp_MIPerror_vs_nEntries_layer_ring.png", "png");

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
