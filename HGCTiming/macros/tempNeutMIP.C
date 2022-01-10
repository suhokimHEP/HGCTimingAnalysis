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
  struct spread {
    spread()
        : a(0.), b(0.), c(0.), d(0.), e(0.), f(0.), g(0.), h(0.), i(0.), j(0.), k(0.),l(0.),m(0.),n(0.),o(0.),p(0.),q(0.),r(0.),s(0.),t(0.),u(0.),v(0.)  {}
    float a, b, c, d, e, f, g, h, i, j, k,l,m,n,o,p,q,r,s,t,u,v;
  };

void tempNeutMIP(){
  //int main(){
  gROOT->Reset();
  gROOT->Macro("./setStyle.C");

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptTitle(0);
 
 TText* title = new TText(1,1,"") ;
 title->SetTextSize(0.04);
 title->SetTextColor(kBlack);
 title->SetTextAlign(11);
 title->SetTextFont(62);
 
 TText* extra = new TText(1,1,"") ;
 extra->SetTextSize(0.03);
 extra->SetTextColor(kBlack);
 extra->SetTextAlign(11);
 extra->SetTextFont(52);
   
 TText* extra2 = new TText(1,1,"") ;
 extra2->SetTextSize(0.025);
 extra2->SetTextColor(kBlack);
 extra2->SetTextAlign(11);
 extra2->SetTextFont(62);
 TText* LRInfo = new TText(1,1,"") ;
 LRInfo->SetTextSize(0.025);
 LRInfo->SetTextColor(kBlack);
 LRInfo->SetTextAlign(11);
 LRInfo->SetTextFont(62);

    char lumistring [50];
     int dummy; 


  //  RooDataSet data ("data", "data", RooArgSet(x));

  //std::string inputFileList = "Gun50_singlelayer.root";
  //std::string inputFileList = "Rand40mu_startup_sn2.0.root";
  //std::string nametag = "startup";
  std::string nametag = "eol";
 FILE * outfulltable;
 FILE * outfulltable2;
 FILE * outfulltable3;
      outfulltable = fopen(Form("tex/mlandau_%s.tex",nametag.c_str()),"w");
      outfulltable2 = fopen(Form("tex/mlandauerr_%s.tex",nametag.c_str()),"w");
      outfulltable3 = fopen(Form("tex/mlandaurat_%s.tex",nametag.c_str()),"w");
  std::string spreadFileList = "outMIP_Neutfit_"+nametag+".tex";
  std::ifstream infile(spreadFileList);
  if (!infile.is_open()) {
    std::cout << "Unable to open '" << spreadFileList << "'" << std::endl;
  }
  std::string line;
int linecount =0;
std::pair<int, vector<float>> fitspread;
std::vector<pair<int,vector<float>>> spreadmap;
vector<float> vecspr;
spread spr;
  while (getline(infile, line)) {
    std::stringstream linestream(line);
    linestream >> spr.a >> spr.b>> spr.c >>spr.d >>spr.e >>spr.f >>spr.g >>spr.h >>spr.i >>spr.j >>spr.k >>spr.l >>spr.m >>spr.n >>spr.o >>spr.p >>spr.q >>spr.r >>spr.s >>spr.t >>spr.u >>spr.v;
    vecspr.push_back(spr.a);
    vecspr.push_back(spr.b);
    vecspr.push_back(spr.c);
    vecspr.push_back(spr.d);
    vecspr.push_back(spr.e);
    vecspr.push_back(spr.f);
    vecspr.push_back(spr.g);
    vecspr.push_back(spr.h);
    vecspr.push_back(spr.i);
    vecspr.push_back(spr.j);
    vecspr.push_back(spr.k);
    vecspr.push_back(spr.l);
    vecspr.push_back(spr.m);
    vecspr.push_back(spr.n);
    vecspr.push_back(spr.o);
    vecspr.push_back(spr.p);
    vecspr.push_back(spr.q);
    vecspr.push_back(spr.r);
    vecspr.push_back(spr.s);
    vecspr.push_back(spr.t);
    vecspr.push_back(spr.u);
    vecspr.push_back(spr.v);
    fitspread.second = vecspr;
    fitspread.first = linecount;
    spreadmap.push_back(fitspread);
    vecspr.clear();
    linecount++;
   } 
  std::string inputFileList = "Rand40mu_"+nametag+".root";
  std::string outtag = inputFileList.substr(0,inputFileList.find(".root"));
  printf("%s",outtag.c_str());
  ROOT::RDataFrame d("MIPtree", inputFileList.c_str());
  auto allMIP = d.Histo1D( {"allMIP", "", 120, 0.5, 10.5}, "MIP_val");
 TFile* inF = TFile::Open(TString(inputFileList));
 TTree* newT = (TTree*)inF->Get("MIPtree");


  auto min_iL = d.Min("MIP_layer").GetValue();
  auto max_iL = d.Max("MIP_layer").GetValue();
  // auto min_iR_all = d.Min("MIP_iR").GetValue();
  // auto max_iR_all = d.Max("MIP_iR").GetValue();
  // auto min_iPhi_all = d.Min("MIP_iPhi").GetValue();
  // auto max_iPhi_all = d.Max("MIP_iPhi").GetValue();
  
  ///single tile fits
  int firstLayer = 37;
  for(int ij_L=min_iL; ij_L <= 8; ++ij_L){
  //for(int ij_L=min_iL; ij_L <= max_iL; ++ij_L){
    std::cout << " fitting for layer = " << ij_L << std::endl;

    auto filterMinMAx = [ij_L](int lval) { return (ij_L == lval); };
    auto radiiForLayer = d.Filter(filterMinMAx, {"MIP_layer"});
    auto absRadiiForLayer = radiiForLayer.Define("abs_MIP_iR", "std::abs(MIP_iR)");
    auto min_iR = absRadiiForLayer.Min("abs_MIP_iR").GetValue();
    auto max_iR = absRadiiForLayer.Max("abs_MIP_iR").GetValue();
        //fprintf (outfulltable, "%f  ", min_iR);		
        //fprintf (outfulltable2, "%f  ", min_iR);		
        //fprintf (outfulltable, "%f  ", max_iR);		
        //fprintf (outfulltable2, "%f  ", max_iR);		

 bool minbool = true;
  //for(int ij_R=23; ij_R <= 23; ++ij_R){
  for(int ij_R=min_iR; ij_R <= max_iR; ++ij_R){
      std::cout << "\n  ==>> Fitting for layer = " << ij_L << " ring = " << ij_R << std::endl;

        int entries = newT->GetEntries(Form("MIP_iR==%d&&MIP_layer==%d",ij_R,ij_L));
        int Mentries = newT->GetEntries(Form("MIP_iR==%d&&MIP_layer==%d&&MIP_val<0.58",ij_R,ij_L));
        float tilespr = spreadmap[ij_L].second[ij_R-min_iR];
	if(tilespr==0.) tilespr=.3168;
  printf("%f\n",tilespr);
  RooWorkspace w("w");    
  RooRealVar xL("xL", "", 0.5, 10.5);
  w.import(xL);
  w.factory("nSig_l[6.e3, 3.e3, 2.e4]");   
  //w.factory("nSig_l[2.e3, 0.0, 2.e4]");   
  //w.factory(Form("nNoise_l[1.e2, 50, 3.e2]")); 
  w.factory(Form("nNoise_l[4.e2, %d, 1.e4]",Mentries)); 

  w.factory("RooLandau::landauL(xL, m_landau_l[1.1, 0.5, 5.], s_landau_l[0.2, 0.01, 0.5])");
  w.factory("RooGaussian::gaussL(xL, m_gauss_l[0.0, 0.0, 0.], s_gauss_l[0.1, 0.01, 1.])");
  w.factory("RooFFTConvPdf::lxgL(xL, landauL, gaussL)");
  //w.factory(Form("RooGaussian::gaussBL(xL, m_gaussNoise[0.0], s_gaussNoise[%f])",tilespr));
  w.factory("RooGaussian::gaussBL(xL, m_gaussNoise[0.0], s_gaussNoise[0.318])");
  w.factory("SUM::modelL(nNoise_l * gaussBL, nSig_l* lxgL)");
  RooAbsPdf * modelL = w.pdf("modelL");

  TH1F* MIP_values_layer = new TH1F("MIP_values_layer", "", 500, 0., 5.);
  TH1F* MIP_values_layer_ring = new TH1F("MIP_values_layer_ring", "", 500, 0., 5.);
  TGraph* MIPerror_vs_nEntries_layer_ring = new TGraph();
  TProfile* tp_MIPerror_vs_nEntries_layer_ring = new TProfile("tp_MIPerror_vs_nEntries_layer_ring", "", 400, 0., 400);
  TH1F* MIP_values_layer_ring_phi = new TH1F("MIP_values_layer_ring_phi", "", 500, 0., 5.);
  TH2F* iR_vs_layer = new TH2F("iR_vs_layer", "", 14, 0., 13, 75, -37, 37.);






      auto filterSel = [ij_L, ij_R](int lval, int rval) { return (ij_L == lval && ij_R == rval); };                            
      auto filterTree = d.Filter(filterSel, {"MIP_layer", "MIP_iR"});
      //      auto mipH = d.Filter(filterSel, {"MIP_layer", "MIP_iR"}).Histo1D({"mipH", "", 120, 0.5, 10.5}, "MIP_val");
	//int numberOfE = entries;
	int numberOfE = 6000;
	auto mipH = filterTree.Range(0, numberOfE).Histo1D({"mipH", "", 120, 0.5, 10.5}, "MIP_val");

      float nEntries = mipH->GetEntries();
      if(nEntries < numberOfE) {
	continue;
	}     
	if(minbool) {
        fprintf (outfulltable, "%d  ", ij_R);		
        fprintf (outfulltable2, "%d  ", ij_R);		
        fprintf (outfulltable3, "%d  ", ij_R);		
	minbool=false;
	}  
    RooDataHist locData("locData","locData", xL, Import(*(mipH)) );

      //fit and save
      RooFitResult * rJL = modelL->fitTo(locData, Minimizer("Minuit2"),Save(true));
      auto fittedVal = w.var("m_landau_l")->getVal();
      MIP_values_layer_ring->Fill(fittedVal);
      auto fittedValError = w.var("m_landau_l")->getError();
	//printf("thevals:%f,%f",fittedVal,fittedValError);
        fprintf (outfulltable, "%f   ", fittedVal);		
        fprintf (outfulltable2, "%f   ", fittedValError);		
        fprintf (outfulltable3, "%f   ", fittedValError/fittedVal);		
      MIPerror_vs_nEntries_layer_ring->SetPoint(MIPerror_vs_nEntries_layer_ring->GetN(), nEntries, fittedValError*100.);
      tp_MIPerror_vs_nEntries_layer_ring->Fill(numberOfE, fittedValError*100.);

      RooPlot* frameL = xL.frame();                                                                                                                            
      frameL->SetXTitle("MIP values");                                                                                                                         
      locData.plotOn(frameL, Binning(120));                                                                                                                    
      modelL->plotOn(frameL, LineColor(kBlue));                                                                                                                
      modelL->plotOn(frameL, Components("lxgL"), LineColor(kRed+1));                                                                                           
      modelL->paramOn(frameL, RooFit::Layout(0.6,0.35,0.9),RooFit::Format("NEA",AutoPrecision(4)));                                                             
      frameL->getAttLine()->SetLineColorAlpha(kWhite, 0.2);                                                                                                    
      frameL->getAttText()->SetTextSize(0.03);                                                                                                                 
      frameL->getAttText()->SetTextFont(42);                                                                                                                   

      TCanvas* cL = new TCanvas();                                                                                                                             

      cL->cd();                                                                                                                                                
      frameL->Draw();                                                                                    
     title->DrawTextNDC(0.2,0.96,"CMS");
     extra->DrawTextNDC(0.3,0.96,"Work In Progress");
     extra2->DrawTextNDC(0.51,0.96,"Phase2");
     LRInfo->DrawTextNDC(0.58,0.96,Form("Layer=%d, Ring=%d, EOL",ij_L+firstLayer,ij_R));
      cL->Print(Form("FixedEvt/Uncer_wNeutfit_%s_mip_L%d_R%d_nEvts%d.png", outtag.c_str(),ij_L+firstLayer, ij_R, numberOfE), ".png"); 
      
      delete cL;
      delete frameL;
      //delete mipH;
	}
        fprintf (outfulltable, "\n");		
        fprintf (outfulltable2, "\n");		
        fprintf (outfulltable3, "\n");		
  }



//  TLatex tL;
//  tL.SetNDC();
//  tL.SetTextSize(0.05);
//  tL.SetTextFont(132);
//
//  gStyle->SetOptStat(1);
//  gStyle->SetOptFit(1);
//  TF1* hfit_layer = new TF1("hfit_layer", "gaus", 0., 2.);
//
//  TCanvas* tc_fitVal = new TCanvas();
//  tc_fitVal->cd();
//  hfit_layer->SetParameters(MIP_values_layer_ring->GetEntries()/2., 1, 0.02);
//  MIP_values_layer_ring->GetXaxis()->SetTitle("MPV values for layers and rings");
//  MIP_values_layer_ring->Draw();
//  MIP_values_layer_ring->Fit("hfit_layer", "q");
//  //MIP_values_layer_ring->Fit("hfit_layer", "R");
//  gPad->Update();
//  tL.DrawLatex(0.5, 0.8, Form("mean = %.2f +/- %.2f", hfit_layer->GetParameter(1), hfit_layer->GetParError(1)));
//  tL.DrawLatex(0.5, 0.7, Form("sigma = %.2f +/- %.2f", hfit_layer->GetParameter(2), hfit_layer->GetParError(2)));
//  tc_fitVal->Print("MIP_values_layer_ring.png", "png");
//
//  TF1* hfit_precision = new TF1("hfit_precision", "[0]/sqrt(x) + [1]", 0., 500.);
//  hfit_precision->SetParameter(1, 0.01);
//  hfit_precision->SetParameter(0, 0.5);
//
//  TCanvas* tc_scatter = new TCanvas();
//  tc_scatter->cd();
//  MIPerror_vs_nEntries_layer_ring->Sort();
//  MIPerror_vs_nEntries_layer_ring->GetXaxis()->SetTitle("nEntries");
//  MIPerror_vs_nEntries_layer_ring->GetYaxis()->SetTitle("precision (%)");
//  MIPerror_vs_nEntries_layer_ring->GetYaxis()->SetRangeUser(0., 100.);
//  MIPerror_vs_nEntries_layer_ring->Draw("ap");
//  MIPerror_vs_nEntries_layer_ring->Fit("hfit_precision","q");
//  //MIPerror_vs_nEntries_layer_ring->Fit("hfit_precision");
//  gPad->Update();
//  tL.DrawLatex(0.5, 0.8, "a / sqrt(x) + b");
//  tL.DrawLatex(0.5, 0.7, Form("a = %.2f +/- %.2f", hfit_precision->GetParameter(0), hfit_precision->GetParError(0)));
//  tL.DrawLatex(0.5, 0.6, Form("b = %.2f +/- %.2f", hfit_precision->GetParameter(1), hfit_precision->GetParError(1)));
//  tc_scatter->Print("MIPerror_vs_nEntries_layer_ring.png", "png");
//
//
//  tc_scatter->cd();
//  tp_MIPerror_vs_nEntries_layer_ring->GetXaxis()->SetTitle("nEntries");
//  tp_MIPerror_vs_nEntries_layer_ring->GetYaxis()->SetTitle("precision (%)");
//  tp_MIPerror_vs_nEntries_layer_ring->GetYaxis()->SetRangeUser(0., 100.);
//  tp_MIPerror_vs_nEntries_layer_ring->SetMarkerStyle(20);
//  tp_MIPerror_vs_nEntries_layer_ring->SetLineWidth(2);
//  tp_MIPerror_vs_nEntries_layer_ring->Draw("e");
//  hfit_precision->SetParameter(1, 0.01);
//  hfit_precision->SetParameter(0, 0.5);
//  tp_MIPerror_vs_nEntries_layer_ring->Fit("hfit_precision","q");
//  //tp_MIPerror_vs_nEntries_layer_ring->Fit("hfit_precision");
//  gPad->Update();
//  tL.DrawLatex(0.5, 0.8, "a / sqrt(x) + b");
//  tL.DrawLatex(0.5, 0.7, Form("a = %.2f +/- %.2f", hfit_precision->GetParameter(0), hfit_precision->GetParError(0)));
//  tL.DrawLatex(0.5, 0.6, Form("b = %.2f +/- %.2f", hfit_precision->GetParameter(1), hfit_precision->GetParError(1)));
//  tL.DrawLatex(0.5, 0.5, Form("#Chi^{2} = %.2f", hfit_precision->GetChisquare()/hfit_precision->GetNDF()));
//  tc_scatter->Print("tp_MIPerror_vs_nEntries_layer_ring.png", "png");

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
