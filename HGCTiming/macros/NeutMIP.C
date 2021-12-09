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
void NeutMIP(){
  //int main(){
  gROOT->Reset();
  gROOT->Macro("./setStyle.C");

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptTitle(0);
 
  //  RooDataSet data ("data", "data", RooArgSet(x));

  //std::string inputFileList = "Gun50_singlelayer.root";
  //std::string inputFileList = "Rand40mu_startup_sn2.0.root";
  std::string inputFileList = "outMIP_Neutfit_eol.root";
  std::string outtag = inputFileList.substr(0,inputFileList.find(".root"));
  printf("%s",outtag.c_str());
  ROOT::RDataFrame d("MIPtree", inputFileList.c_str());
    //ROOT::RDataFrame d("MIPtree", "outMIP_SoN4.root");
  auto allMIP = d.Histo1D( {"allMIP", "", 60, -3, 3}, "MIP_val");
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
  for(int ij_L=min_iL; ij_L <= max_iL; ++ij_L){
  //for(int ij_L=min_iL; ij_L <= 3; ++ij_L){
    std::cout << " fitting for layer = " << ij_L << std::endl;

    auto filterMinMAx = [ij_L](int lval) { return (ij_L == lval); };
    auto radiiForLayer = d.Filter(filterMinMAx, {"MIP_layer"});
    auto absRadiiForLayer = radiiForLayer.Define("abs_MIP_iR", "std::abs(MIP_iR)");
    auto min_iR = absRadiiForLayer.Min("abs_MIP_iR").GetValue();
    auto max_iR = absRadiiForLayer.Max("abs_MIP_iR").GetValue();


  for(int ij_R=min_iR; ij_R <= max_iR; ++ij_R){
      std::cout << "\n  ==>> Fitting for layer = " << ij_L << " ring = " << ij_R << std::endl;

        int entries = newT->GetEntries(Form("MIP_iR==%d&&MIP_layer==%d",ij_R,ij_L));
        int Mentries = newT->GetEntries(Form("MIP_iR==%d&&MIP_layer==%d&&MIP_val<0.543",ij_R,ij_L));
  RooWorkspace w("w");    
  RooRealVar xL("xL", "", -3, 3);
  w.import(xL);
  w.factory("nNoise_l[5.e7, 5.e5, 5.e8]"); 

  w.factory("RooGaussian::gaussBL(xL, m_gauss_l[0.0, 0.0, 0.], s_gauss_l[0.1, 0.01, 1.])");
  w.factory("SUM::modelL(nNoise_l * gaussBL)");
  RooAbsPdf * modelL = w.pdf("modelL");






      auto filterSel = [ij_L, ij_R](int lval, int rval) { return (ij_L == lval && ij_R == rval); };                            
      auto filterTree = d.Filter(filterSel, {"MIP_layer", "MIP_iR"});
      //      auto mipH = d.Filter(filterSel, {"MIP_layer", "MIP_iR"}).Histo1D({"mipH", "", 60, 0.5, 10.5}, "MIP_val");
	int numberOfE = entries;
	//int numberOfE = 200;
	auto mipH = filterTree.Range(0, numberOfE).Histo1D({"mipH", "", 60, -3, 3}, "MIP_val");

      float nEntries = mipH->GetEntries();
      if(nEntries < numberOfE) continue;     
      RooDataHist locData("locData","locData", xL, Import(*(mipH)) );

      //fit and save
      RooFitResult * rJL = modelL->fitTo(locData, Minimizer("Minuit2"),Save(true));
      auto fittedVal = w.var("m_gauss_l")->getVal();
      auto fittedValError = w.var("m_gauss_l")->getError();

      RooPlot* frameL = xL.frame();                                                                                                                            
      frameL->SetXTitle("MIP values");                                                                                                                         
      locData.plotOn(frameL, Binning(60));                                                                                                                    
      modelL->plotOn(frameL, LineColor(kBlue));                                                                                                                
      //modelL->plotOn(frameL, Components("lxgL"), LineColor(kRed+1));                                                                                           
      modelL->paramOn(frameL, RooFit::Layout(0.6,0.8,0.9),RooFit::Format("NEA",AutoPrecision(1)));                                                             
      frameL->getAttLine()->SetLineColorAlpha(kWhite, 0.2);                                                                                                    
      frameL->getAttText()->SetTextSize(0.03);                                                                                                                 
      frameL->getAttText()->SetTextFont(42);                                                                                                                   

      TCanvas* cL = new TCanvas();                                                                                                                             
      cL->cd();                                                                                                                                                
      frameL->Draw();                                                                                    
     
      cL->Print(Form("singleFits/Neutrino%s_mip_L%d_R%d_nEvts%d.png", outtag.c_str(),ij_L+firstLayer, ij_R, numberOfE), ".png"); 
      
      delete cL;
      delete frameL;
      //delete mipH;
	}
  }


}
