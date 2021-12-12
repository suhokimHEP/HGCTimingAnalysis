//g++  -o lookAtMIP_correctRate  lookAtMIP_correctRate.cpp `root-config --cflags --libs`
//#include <boost/filesystem.hpp>
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
#include "TFitResult.h"
#include "TFile.h"
#include "TProfile2D.h"
using namespace std;
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

  //  std::cout << histo->GetName() << " " << histo->GetBinCenter(xBin) << std::endl; 
  return histo->GetBinCenter(xBin);
}



int Neutrino_correctRate(){
	TString filename="Neutfit_startup_sn2.0_40";
	TString outname = "outMIP_"+filename;
	  gROOT->Macro("./setStyle.C");

	  gStyle->SetOptStat(0);
	  gStyle->SetOptTitle(0);
	  
	  gROOT->Reset();
	  gStyle->SetOptStat(1);
	  gStyle->SetOptFit(1);

	  std::cout << " inizio ci sono " << std::endl; 
	  TFile* inF = TFile::Open("../test/"+filename+".root");
	  //TFile* inF = TFile::Open("../test/gitignore/"+aversion+"/"+filename+".root");
	  //TFile* inF = TFile::Open("root://cmsxrootd.fnal.gov//store/user/skim2/"+filename);

	  //config
	  int firstLayer = 37;
	  TTree* newT = (TTree*)inF->Get("ana/newT");
	  std::vector<float> *muonP = 0;
	  std::vector<float> *muonPt = 0;
	  //only with pat muons if enabled
	  //  std::vector<float> *muonSegC = 0;
	  std::vector<float> *muonEta = 0;
	  std::vector<float> *muonPhi = 0;
	  std::vector<float> *crossX = 0;
	  std::vector<float> *crossY = 0;
	  std::vector<float> *crossEX = 0;
	  std::vector<float> *crossEY = 0;
	  std::vector<float> *crossZ = 0;
	  std::vector<float> *crossL = 0;
	  std::vector<float> *crossM = 0;
	  std::vector<float> *crossEta = 0;
	  std::vector<float> *crossPhi = 0;
	  std::vector<float> *crossCellEta = 0;
	  std::vector<float> *crossCellPhi = 0;
	  std::vector<float> *crossP = 0;
	  std::vector<float> *crossPt = 0;
	  std::vector<int> *crossiR = 0;
	  std::vector<int> *crossiP = 0;
	  std::vector<float> *muonChi = 0;
	  std::vector<short int> *muonTrkQ = 0;
	  std::vector<float> *recHitX = 0;
	  std::vector<float> *recHitY = 0;
	  std::vector<float> *recHitZ = 0;
	  std::vector<float> *recHitSize = 0;
	  std::vector<int> *recHitiPhi = 0;
	  std::vector<int> *recHitiR = 0;
	  std::vector<int> *recHitL = 0;
	  std::vector<float> *recHitEta = 0;
	  std::vector<float> *recHitPhi = 0;
	  std::vector<float> *recHitEne = 0;
	  std::vector<float> *recHitMip = 0;
	  std::vector<float> *recHitNoise = 0;


	  newT->SetBranchAddress("muonP", &muonP);
	  newT->SetBranchAddress("muonPt", &muonPt);
	  //  newT->SetBranchAddress("muonSegC", &muonSegC);
	  newT->SetBranchAddress("muonEta", &muonEta);
	  newT->SetBranchAddress("muonPhi", &muonPhi);
	  newT->SetBranchAddress("crossX", &crossX);
	  newT->SetBranchAddress("crossY", &crossY);
	  newT->SetBranchAddress("crossEX", &crossEX);
	  newT->SetBranchAddress("crossEY", &crossEY);
	  newT->SetBranchAddress("crossZ", &crossZ);
	  newT->SetBranchAddress("crossL", &crossL);
	  newT->SetBranchAddress("crossM", &crossM);
	  newT->SetBranchAddress("crossEta", &crossEta);
	  newT->SetBranchAddress("crossPhi", &crossPhi);
	  newT->SetBranchAddress("crossCellEta", &crossCellEta);
	  newT->SetBranchAddress("crossCellPhi", &crossCellPhi);
	  newT->SetBranchAddress("crossP", &crossP);
	  newT->SetBranchAddress("crossPt", &crossPt);
	  newT->SetBranchAddress("crossiR", &crossiR);
	  newT->SetBranchAddress("crossiP", &crossiP);
	  newT->SetBranchAddress("muonChi", &muonChi);
	  newT->SetBranchAddress("muonTrkQ", &muonTrkQ);
	  newT->SetBranchAddress("recHitX", &recHitX);
	  newT->SetBranchAddress("recHitY", &recHitY);
	  newT->SetBranchAddress("recHitZ", &recHitZ);
	  newT->SetBranchAddress("recHitSize", &recHitSize);
	  newT->SetBranchAddress("recHitiR", &recHitiR);
	  newT->SetBranchAddress("recHitiPhi", &recHitiPhi);
	  newT->SetBranchAddress("recHitL", &recHitL);
	  newT->SetBranchAddress("recHitEta", &recHitEta);
	  newT->SetBranchAddress("recHitPhi", &recHitPhi);
	  newT->SetBranchAddress("recHitEne", &recHitEne);
	  newT->SetBranchAddress("recHitMip", &recHitMip);
	  newT->SetBranchAddress("recHitNoise", &recHitNoise);

	   int nEvents = newT->GetEntries();
	   std::cout << " nEvents = " << nEvents << std::endl;
	  TFile outMIP(outname+".root", "recreate");
	  //TFile outMIP("outMIP.root", "recreate");
	  TTree* tt = new TTree("MIPtree", "MIPtree");
	  Float_t MIP_val;
	  Int_t MIP_layer; 
	  Int_t MIP_iR;
	  Int_t MIP_iPhi;
	  Float_t MIP_SoN;
	  tt->Branch("MIP_layer", &MIP_layer, "MIP_layer/I");
	  tt->Branch("MIP_iR", &MIP_iR, "MIP_iR/I");
	  tt->Branch("MIP_iPhi", &MIP_iPhi, "MIP_iPhi/I");
	  tt->Branch("MIP_SoN", &MIP_SoN, "MIP_SoN/F");
	  tt->Branch("MIP_val", &MIP_val, "MIP_val/F");
	   for(int ij=0; ij<nEvents; ++ij){
	  if (ij%100 == 0){ std::cout << " entry " << ij << std::endl; }
	     //for(int ij=14; ij<16; ++ij){
	     newT->GetEntry(ij);

	     for(int iR=0; iR<recHitX->size(); ++iR){
	       int recL = recHitL->at(iR) - firstLayer;
	       MIP_layer=recL;
	       MIP_iR=recHitiR->at(iR);
	       MIP_iPhi=recHitPhi->at(iR);
	       MIP_val=recHitMip->at(iR);
	       MIP_SoN=2.5;
		tt->Fill();	
	     }
	}

	  tt->Write();
	  outMIP.Close();

  std::cout << " all done now print results " << std::endl;


  std::cout << " all done - ciao " << std::endl;

  return 0;



}
