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



//int main(int argc, char** argv){
//int lookAtMIP_correctRate(char** argv){
//int lookAtMIP_correctRate(TString argv){
//int lookAtMIP_correctRate(int argc, char** argv){
int lookAtMIP_correctRate(TString filename, TString aversion){
	//TString filename(argv[1]);
	TString outname = "outMIP_"+filename;
	  gROOT->Macro("./setStyle.C");

	  gStyle->SetOptStat(0);
	  gStyle->SetOptTitle(0);
	  
	  gROOT->Reset();
	  gStyle->SetOptStat(1);
	  gStyle->SetOptFit(1);

	  std::cout << " inizio ci sono " << std::endl; 

	  
	  int iColors[14] = {kBlack, kGray+1, kGreen+1, kSpring-2, kYellow+1, kOrange-2, kRed+1, kPink-2, kMagenta+1, kViolet-2, kBlue+1, kAzure-2, kCyan+1, kTeal-2};
	  int iStyle[4] = {20, 21, 22, 23}; //kGray+1};


	  TFile* inF = TFile::Open(filename+".root");
	  //TFile* inF = TFile::Open("../test/gitignore/"+aversion+"/"+filename+".root");
	  //TFile* inF = TFile::Open("root://cmsxrootd.fnal.gov//store/user/skim2/"+filename);

	  //config
	  int firstLayer = 37;
	  // float Nevt = 2.09206e+06;
	  //  float Nevt = 12492501;
	  float Nevt = 4981488;
	  float timeEqui = Nevt / 40.e6; //s

	  float weights[51] = {0, 
			       8.89454, 10.9379, 10.9379, 10.9379, 10.9379, 10.9379, 10.9379, 10.9379, 10.9379, 10.9329, 
			       10.9329, 10.9379, 10.9379, 10.9382, 10.9382, 10.9382, 10.9382, 10.9382, 10.9382, 10.9382, 
			       10.9382, 10.9382, 10.9382, 10.9382, 10.9382, 10.9382, 10.9382, 32.3321, 51.5743, 51.4442, 
			       51.4442, 51.4442, 51.4442, 51.4442, 51.4442, 51.4442, 51.4442, 51.4442, 51.4442, 69.5131, 
	 /*da 41 a 50*/	       87.582, 87.582, 87.582, 87.582, 87.582, 87.2146, 86.8883, 86.9295, 86.9295, 86.9295}; 

	  float zVal[14] = {411.29, 416.739, 422.187, 427.636, 436.172, 444.722, 453.263, 461.817, 470.371, 478.925, 487.47, 496.024, 504.577, 513.127};
	  float pigreco = 3.1415926535;

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

	  std::map<std::pair<int, std::pair<int, int> >, float > okChannels;
	  std::map<std::pair<int, std::pair<int, int> >, std::vector<float> > enePerRh;
	  std::map<std::pair<int, std::pair<int, int> >, std::vector<float> > mipPerRh;
	 
	  std::map<std::pair<int, int>, float > okChannelsPhi;
	  std::map<std::pair<int, int>, std::vector<float> > enePerPhi;
	  std::map<std::pair<int, int>, std::vector<float> > mipPerPhi;
	  std::map<std::pair<int, int>, float > okChannelsPhi_3H;
	  std::map<std::pair<int, int>, std::vector<float> > mipPerPhi_3H;
	  std::map<std::pair<int, std::pair<int, int> >, float > okChannels_3H;
	  std::map<std::pair<int, std::pair<int, int> >, std::vector<float> > mipPerRh_3H;

	  TH1F* mipAll = new TH1F("mipAll", "", 80, -1., 19.);
	  TH2F* h2_YvsX[14];
	  TH2F* h2_iRvsiPhi[14];
	  TH1F* h_dX[14];
	  TH1F* h_dY[14];
	  TH1F* h_dR[14];
	  TH1F* h_dRep[14];
	  TH1F* h_dRerr[14];
	  TH1F* h_dXerr[14];
	  TH1F* h_dYerr[14];
	  TH1F* h_d_iR[14];
	  TH1F* h_d_iPhi[14];
	  TH2F* h_d_iR_vs_iPhi[14];

	  TH2F* h2_dR_vsChi[14];
	  TH2F* h2_dRep_vsChi[14];
	  TH2F* h2_dR_vsTrkQ[14];

	  TH2F* h_dRvsP[14];
	  
	  TProfile2D* h2_iR_vsEta_size[14];
	  TH2F* h2_iR_vsEta[14];
	  TH2F* h2_iR_vsR[14];

	  //  TH1F* h_etaMuon = new TH1F("h_etaMuon", "", 300, 1., 4.);
	  TH1F* h_etaMuon = new TH1F("h_etaMuon", "", 400, -10., 10.);
	  TH1F* h_ptMuon = new TH1F("h_ptMuon", "", 100, 0., 100.);

	  TH1F* h_TrkChi_all = new TH1F("h_TrkChi_all", "", 1000, 0., 20.);
	  TH1F* h_TrkQ_all = new TH1F("h_TrkQ_all", "", 10, 0., 10.);
	  TH1F* h_TrkChi_passed = new TH1F("h_TrkChi_passed", "", 1000, 0., 20.);
	  TH1F* h_TrkQ_passed = new TH1F("h_TrkQ_passed", "", 10, 0., 10.);

	  TH1F* h_trkP[14];
	  TH1F* h_trkPt[14];

	  TH1F* h_dtrkP[14];

	  //////// analyze rate
	  TH1F* nMuons_vsEta[14];
	  TH1F* nMuonsPercm2_vsEta[14];
	  TH1F* ratePercm2_vsEta[14];
	  TH1F* cm2_vsEta[14];

	  TH1F* nMuons_3H_vsEta[14];
	  TH1F* nMuonsPercm2_3H_vsEta[14];
	  TH1F* ratePercm2_3H_vsEta[14];

	  TH1F* nMuons_vsEta_All[14];
	  TH1F* nMuonsPercm2_vsEta_All[14];
	  TH1F* ratePercm2_vsEta_All[14];
	  TH1F* cm2_vsEta_All[14];


	  TH1F* h_muonSeg_ok = new TH1F("h_muonSeg_ok", "", 1100, 0., 1.1);
	  TH1F* h_muonSeg_bad = new TH1F("h_muonSeg_bad", "", 1100, 0., 1.1);

	  //check for tracklets

	  TH1F* h_MuonCross[14];
	  TH1F* h_MuonNHit[14];

	  TH1F* h_MuonCrossK[14][14];
	  TH1F* h_MuonNHitK[14][14];

	  TH1F* hDen_etaRange[14];
	  TH1F* hNum_etaRange[14];
	  TH1F* hDen_etaRange_3H[14];
	  TH1F* hNum_etaRange_3H[14];

	  TH1F* numHits_etaRange[14];
	//boost::filesystem::create_directory("plotdir");
	//boost::filesystem::create_directory("plotdir/cells");
	//boost::filesystem::create_directory("plotdir/cellsPhi");
	   for(int ij=0; ij<14; ++ij){
	     h2_YvsX[ij] = new TH2F(Form("h2_YvsX_L%d", ij), "", 600, -300., 300., 600, -300., 300.);
	     h2_iRvsiPhi[ij] = new TH2F(Form("h2_iRvsiPhi_L%d", ij), "", 100, -50., 50., 300, 0., 300.);
	     h_dX[ij] = new TH1F(Form("h_dX_L%d", ij), "", 100, -10., 10.);
	     h_dY[ij] = new TH1F(Form("h_dY_L%d", ij), "", 100, -10., 10.);
	     h_dR[ij] = new TH1F(Form("h_dR_L%d", ij), "", 100, 0., 10.);
	     h_dRep[ij] = new TH1F(Form("h_dRep_L%d", ij), "", 500, 0., 0.5);
	     h_dRerr[ij] = new TH1F(Form("h_dRerr_L%d", ij), "", 100, 0., 1.);
	     h_dXerr[ij] = new TH1F(Form("h_dXerr_L%d", ij), "", 100, 0., 1.);
	     h_dYerr[ij] = new TH1F(Form("h_dYerr_L%d", ij), "", 100, 0., 1.);

	     h_d_iR[ij] = new TH1F(Form("h_d_iR_L%d", ij), "", 20, -10., 10.);
	     h_d_iPhi[ij] = new TH1F(Form("h_d_iPhi_L%d", ij), "", 20, -10., 10.);
	     h_d_iR_vs_iPhi[ij] = new TH2F(Form("h_d_iR_vs_iPhi_L%d", ij), "", 20, -10., 10., 20, -10., 10.);

	     h2_dR_vsChi[ij] = new TH2F(Form("h2_dR_vsChi_L%d", ij), "", 100, 0., 10., 100, 0., 10.);
	     h2_dRep_vsChi[ij] = new TH2F(Form("h2_dRep_vsChi_L%d", ij), "", 500, 0., 0.5, 100, 0., 10.);
	     h2_dR_vsTrkQ[ij] = new TH2F(Form("h2_dR_vsTrkQ_L%d", ij), "", 100, 0., 10., 10, 0., 10.);
	     h_dRvsP[ij] = new TH2F(Form("h_dRvsP_L%d", ij), "", 100, 0., 10., 100, 0., 10.);
	     h2_iR_vsEta_size[ij] = new TProfile2D(Form("h2_iR_vsEta_size_L%d", ij), "", 300, 1., 4., 50, 0., 50.);
	     h2_iR_vsEta[ij] = new TH2F(Form("h2_iR_vsEta_L%d", ij), "", 300, 1., 4., 50, 0., 50.);
	     h2_iR_vsR[ij] = new TH2F(Form("h2_iR_vsR_L%d", ij), "", 500, 100., 300., 50, 0., 50.);

	     nMuons_vsEta[ij] = new TH1F(Form("nMuons_vsEta_L%d", ij+firstLayer), "", 400, 0., 4.);
	     nMuonsPercm2_vsEta[ij] = new TH1F(Form("nMuonsPercm2_vsEta_L%d", ij+firstLayer), "", 400, 0., 4.);
	     ratePercm2_vsEta[ij] = new TH1F(Form("ratePercm2_vsEta_L%d", ij+firstLayer), "", 400, 0., 4.);
	     cm2_vsEta[ij] = new TH1F(Form("cm2_vsEta_L%d", ij+firstLayer), "", 400, 0., 4.);

	     nMuons_3H_vsEta[ij] = new TH1F(Form("nMuons_3H_vsEta_L%d", ij+firstLayer), "", 400, 0., 4.);
	     nMuonsPercm2_3H_vsEta[ij] = new TH1F(Form("nMuonsPercm2_3H_vsEta_L%d", ij+firstLayer), "", 400, 0., 4.);
	     ratePercm2_3H_vsEta[ij] = new TH1F(Form("ratePercm2_3H_vsEta_L%d", ij+firstLayer), "", 400, 0., 4.);

	     nMuons_vsEta_All[ij] = new TH1F(Form("nMuons_vsEta_All_L%d", ij+firstLayer), "", 400, 0., 4.);
	     nMuonsPercm2_vsEta_All[ij] = new TH1F(Form("nMuonsPercm2_vsEta_All_L%d", ij+firstLayer), "", 400, 0., 4.);
	     ratePercm2_vsEta_All[ij] = new TH1F(Form("ratePercm2_vsEta_All_L%d", ij+firstLayer), "", 400, 0., 4.);
	  
	     h_MuonCross[ij] = new TH1F(Form("h_MuonCross_L%d", ij+firstLayer), "", 300, 1., 4.);
	     h_MuonNHit[ij] = new TH1F(Form("h_MuonNHit_L%d_N%d", ij+firstLayer, 1), "", 300, 1., 4.);

	     hDen_etaRange[ij] = new TH1F(Form("hDen_etaRange_L%d", ij+firstLayer), "", 1, 1.55, 1.71);
	     hNum_etaRange[ij] = new TH1F(Form("hNum_etaRange_L%d", ij+firstLayer), "", 1, 1.55, 1.71);

	     hDen_etaRange_3H[ij] = new TH1F(Form("hDen_etaRange_3H_L%d", ij+firstLayer), "", 1, 1.55, 1.71);
	     hNum_etaRange_3H[ij] = new TH1F(Form("hNum_etaRange_3H_L%d", ij+firstLayer), "", 1, 1.55, 1.71);

	     numHits_etaRange[ij] = new TH1F(Form("numHits_etaRange_L%d", ij+firstLayer), "", 20, 0., 20.);

	     for(int kl=0; kl<14; ++kl){
	       h_MuonNHitK[kl][ij] = new TH1F(Form("h_MuonNHitK_L%d_N%d", ij+firstLayer, kl), "", 300, 1., 4.);
	       h_MuonCrossK[kl][ij] = new TH1F(Form("h_MuonCrossK_L%d_N%d", ij+firstLayer, kl), "", 300, 1., 4.);
	     }

	     h_dtrkP[ij] = new TH1F(Form("h_dtrkP_L%d", ij+firstLayer), "", 1000, 0., 500.);
	     h_trkP[ij] = new TH1F(Form("h_trkP_L%d", ij+firstLayer), "", 1000, 0., 100.);
	     h_trkPt[ij] = new TH1F(Form("h_trkPt_L%d", ij+firstLayer), "", 1000, 0., 100.);


	   }

	   bool matchByID = true;
	   bool matchByID_p009 = true;
	   float edgeID_value = 0.009;
	   bool matchByID_error = true;
	   bool matchBydR = false;


	   float errorMax[14] = {0.3, 0.31, 0.32, 0.38, 0.4, 0.42, 0.45, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85};

	   int nEvents = newT->GetEntries();
	   std::cout << " nEvents = " << nEvents << std::endl;
	   //for(int ij=0; ij<nEvents; ++ij){
	   for(int ij=0; ij<40; ++ij){
	  if (ij%10 == 0){ std::cout << " entry " << ij << std::endl; }
	     //for(int ij=14; ij<16; ++ij){
	     newT->GetEntry(ij);

	     //     std::cout << "\n  evento = " << ij << std::endl;
	     std::map<int, std::vector<float>> rhX;
	     std::map<int, std::vector<float>> rhY;
	     std::map<int, std::vector<float>> rhZ;
	     std::map<int, std::vector<float>> rhA;
	     std::map<int, std::vector<int>> rhiPhi;
	     std::map<int, std::vector<int>> rhiR;
	     std::map<int, std::vector<float>> rhEta;
	     std::map<int, std::vector<float>> rhPhi;
	     std::map<int, std::vector<float>> rhMip;
	     std::map<int, std::vector<float>> rhEne;
	     std::map<int, std::vector<int>> trkN;
	     std::map<int, std::vector<float>> trkX;
	     std::map<int, std::vector<float>> trkY;
	     std::map<int, std::vector<float>> trkXerr;
	     std::map<int, std::vector<float>> trkYerr;
	     std::map<int, std::vector<float>> trkEta;
	     std::map<int, std::vector<float>> trkPhi;
	     std::map<int, std::vector<float>> trkCEta;
	     std::map<int, std::vector<float>> trkCPhi;
	     std::map<int, std::vector<int>> trkiR;
	     std::map<int, std::vector<int>> trkiP;
	     std::map<int, std::vector<float>> trkSegC;
	     std::map<int, std::vector<float>> trkPt;
	     std::map<int, std::vector<float>> trkP;
	     std::map<int, std::vector<float>> trkChi;
	     std::map<int, std::vector<float>> trkGeneralPt;
	     std::map<int, std::vector<short int>> trkQ;
	     std::map<std::pair<int,int>, float> seltrkP;
	     //    std::cout << " muonP->size() = " << muonP->size() << std::endl;
	     for(int iM=0; iM<muonP->size(); ++iM){
	       int trkL = crossL->at(iM) - firstLayer;
	       trkN[trkL].push_back(int(crossM->at(iM)));
	       std::pair<int, int> selTrkIdx (trkL, int(crossM->at(iM)));
	       seltrkP[selTrkIdx] = crossP->at(iM);
	       trkX[trkL].push_back(crossX->at(iM));
	       trkY[trkL].push_back(crossY->at(iM));
	       trkXerr[trkL].push_back(crossEX->at(iM));
	       trkYerr[trkL].push_back(crossEY->at(iM));
	       trkEta[trkL].push_back(crossEta->at(iM));
	       trkPhi[trkL].push_back(crossPhi->at(iM));
	       trkCEta[trkL].push_back(crossCellEta->at(iM));
	       trkCPhi[trkL].push_back(crossCellPhi->at(iM));
	       trkP[trkL].push_back(crossP->at(iM));
	       trkPt[trkL].push_back(crossPt->at(iM));
	       trkiR[trkL].push_back(crossiR->at(iM));
	       trkiP[trkL].push_back(crossiP->at(iM));
	       //only  in latest pT 5
	       //       trkSegC[trkL].push_back(muonSegC->at(iM));
	       trkGeneralPt[trkL].push_back(muonPt->at(iM));
	       trkChi[trkL].push_back(muonChi->at(iM));
	       trkQ[trkL].push_back(muonTrkQ->at(iM));
	       if(trkL == 0){
		 h_etaMuon->Fill(std::abs(muonEta->at(iM)));
		 h_ptMuon->Fill(muonPt->at(iM));

		 h_TrkChi_all->Fill(muonChi->at(iM));
		 h_TrkQ_all->Fill(muonTrkQ->at(iM));
	       }
	     }

	     for(int iR=0; iR<recHitX->size(); ++iR){
	       int recL = recHitL->at(iR) - firstLayer;
	       rhX[recL].push_back(recHitX->at(iR));
	       rhY[recL].push_back(recHitY->at(iR));
	       rhZ[recL].push_back(recHitZ->at(iR));
	       rhA[recL].push_back(recHitSize->at(iR));
	       rhiPhi[recL].push_back(recHitiPhi->at(iR));
	       rhiR[recL].push_back(recHitiR->at(iR));
	       rhEta[recL].push_back(recHitEta->at(iR));
	       rhPhi[recL].push_back(recHitPhi->at(iR));
	       rhMip[recL].push_back(recHitMip->at(iR));
	       rhEne[recL].push_back(recHitEne->at(iR));
	     }
	     //  std::cout << " end trk and rh loops " << std::endl;

	     std::map<int, bool> foundFirstHit;
	     std::map<int, float> x1stHit;
	     std::map<int, float> y1stHit;
	     std::map<int, float> xPrevHit;
	     std::map<int, float> yPrevHit;
	     std::map<int, int> matchedTrack;
	     std::map<int, int> matchedTrackLast8L;
	     std::map<int, int> matchedTrackL[14];
	     std::map<int, int> startLMatchedTrack;

	     for(int iL=0; iL < 14; ++iL){
	       int iTc = -1;

	       for(auto ij : seltrkP){
		 int layer = ij.first.first;
		 int idx = ij.first.second;
		 float momentum = ij.second;

		 std::pair<int, int> prevTrkIdx(layer-1, idx);
		 //rafix
		 if(seltrkP.find(prevTrkIdx) != seltrkP.end()){
		   h_dtrkP[layer]->Fill(std::abs(momentum - seltrkP[prevTrkIdx])*1.e3);
		   //std::cout << " diff = " << std::abs(momentum - seltrkP[prevTrkIdx])*1.e3 << std::endl;
		 }
	       }

	       bool foundMatchInLayer = false;

	       for(auto iT : trkX[iL]){
		 ++iTc;

		 int trackIdx = trkN[iL][iTc];
		 h_trkP[iL]->Fill(trkP[iL][iTc]);
		 h_trkPt[iL]->Fill(trkPt[iL][iTc]);

		 if(trkGeneralPt[iL][iTc] < 10.) continue;
		 if(trkChi[iL][iTc] > 1.5) continue;
		 
		 //std::cout << " trackIdx = " << trackIdx << std::endl;
		 if(foundFirstHit.find(trackIdx) == foundFirstHit.end()) foundFirstHit[trackIdx] = false;

		 float trackEta = trkEta[iL][iTc];
		 float trackPhi = trkPhi[iL][iTc];
		 float trackY = trkY[iL][iTc];
		 float trkR = sqrt(iT*iT + trackY*trackY);
		 float dXerr = trkXerr[iL][iTc];
		 float dYerr = trkYerr[iL][iTc];
		 float trkEtaOnLayer = asinh(zVal[iL]/trkR);                            
		 nMuons_vsEta_All[iL]->Fill(std::abs(trkEtaOnLayer));  

		 if(matchByID_error){
		   if(dXerr > errorMax[iL] || dYerr > errorMax[iL])
		     continue;
		 }
		 if(matchByID_p009){
		   float dPhi = std::abs(trkCPhi[iL][iTc] - trackPhi);
		   if(dPhi > pigreco) dPhi -= 2* pigreco;
		   if(std::abs(trkCEta[iL][iTc] - trackEta) > edgeID_value || std::abs(dPhi) > edgeID_value ){
		     continue;
		   }
		 }
		 //std::cout << " fill den iL = " << iL << " trackIdx = " << trackIdx << std::endl;
		 h_MuonCross[iL]->Fill(trkEtaOnLayer); 
		 if(trkEtaOnLayer > 1.55 && trkEtaOnLayer < 1.7) hDen_etaRange[iL]->Fill(trkEtaOnLayer); 

		 bool alreadyFilled = false;
		 int iRc = -1;
		 int nearestHit = -1;
		 float nearestdR = 99.;

		 for(auto iR : rhX[iL]){
		   ++iRc;

		   if(trkEta[iL][iTc] * rhZ[iL][iRc] < 0.) { 
		     continue;
		   }

		   float dX = iR - iT;
		   float recEta = rhEta[iL][iRc];
		   float recPhi = rhPhi[iL][iRc];
		   float recY = rhY[iL][iRc];
		   float dY = recY - trackY;
		   float dR2 = dX*dX + dY*dY;
		   float dR = sqrt(dR2);
		   float dRerr = sqrt(dXerr*dXerr + dYerr*dYerr);

		   float reciR = rhiR[iL][iRc];
		   auto dp = std::abs(recPhi - trackPhi);
		   if (dp > pigreco) dp -= (2 * pigreco);
		   float dRep2 = (recEta - trackEta) * (recEta - trackEta) + dp * dp; 
		   float dRep = sqrt(dRep2);

		   h_dX[iL]->Fill(dX);
		   h_dY[iL]->Fill(dY);
		   h_dR[iL]->Fill(dR);
		   h_dRep[iL]->Fill(dRep);
		   //h_dRvsP[iL]->Fill(dR, trkP[iL][iTc]);
		   //	   std::cout << " dRerr = " << dRerr << " dXerr = " << dXerr << " dYerr = " << dYerr  << std::endl;
		   h_dRerr[iL]->Fill(dRerr);
		   h_dXerr[iL]->Fill(dXerr);
		   h_dYerr[iL]->Fill(dYerr);

		   h2_dR_vsChi[iL]->Fill(dR, trkChi[iL][iTc]);
		   h2_dRep_vsChi[iL]->Fill(dRep, trkChi[iL][iTc]);
		   h2_dR_vsTrkQ[iL]->Fill(dR, trkQ[iL][iTc]);

		   std::pair<int, int> channelPhi = std::pair<int, int>(iL, std::abs(reciR));

		   if(trkChi[iL][iTc] > 1.5){
		     continue;
		   }

		   if(matchByID){
		     if(rhiR[iL][iRc] != trkiR[iL][iTc] || rhiPhi[iL][iRc] != trkiP[iL][iTc]) continue;

		     if(matchByID_error){
		       if(dXerr > errorMax[iL] || dYerr > errorMax[iL])
			 continue;
		     }

		     if(matchByID_p009){
		       float dPhi = std::abs(trkPhi[iL][iTc] - rhPhi[iL][iRc]);
		       if(dPhi > pigreco) dPhi -= 2* pigreco;
		       if(std::abs(trkEta[iL][iTc] - rhEta[iL][iRc]) > edgeID_value || std::abs(dPhi) > edgeID_value ){
			 continue;
		       }
		     }

		     if(okChannelsPhi.find(channelPhi) != okChannelsPhi.end()) {
		       okChannelsPhi[channelPhi] += 1;
		     }
		     else {  
		       okChannelsPhi[channelPhi] = 1;
		     }
		     
		     h_MuonNHit[iL]->Fill(trkEtaOnLayer);
		     if(trkEtaOnLayer >= 1.55 && trkEtaOnLayer < 1.7){
		       hNum_etaRange[iL]->Fill(trkEtaOnLayer); 
		       //	     std::cout << " hNum_etaRange fillato trkEtaOnLayer = " << trkEtaOnLayer  << std::endl;
		     }
		     // else 	     std::cout << " trkEtaOnLayer = " << trkEtaOnLayer  << std::endl;
		     // //  std::cout << " fill num iL = " << iL << " trackIdx = " << trackIdx << std::endl;
		     

		     if(foundFirstHit[trackIdx] == false){
		       matchedTrack[trackIdx] = 1;
		       startLMatchedTrack[trackIdx] = iL;
		       if(iL > 5)
			 matchedTrackLast8L[trackIdx] = 1;
		     }
		     else {
		       matchedTrack[trackIdx] += 1;
		       if(iL > 5)
			 matchedTrackLast8L[trackIdx] += 1;
		     }
		     foundFirstHit[trackIdx] = true;

		     matchedTrackL[iL][trackIdx] = 1;
		     

		     int localBin = h_MuonNHit[iL]->GetXaxis()->FindBin(trkEtaOnLayer);
		     if(h_MuonNHit[iL]->GetBinContent(localBin) > h_MuonCross[iL]->GetBinContent(localBin)){
		       std::cout << " event a = " << ij << std::endl;
		       return 100;
		     }
		     
		   
		     //here rate for trk
		     if(!alreadyFilled){
		       nMuons_vsEta[iL]->Fill(std::abs(trkEtaOnLayer));
		       alreadyFilled = true;
		     }
		     
		     mipPerPhi[channelPhi].push_back(rhMip[iL][iRc]);
		     continue;
		   }
		   else if(matchBydR){
		     float hitTrack_dx = iR - iT;
		     float hitTrack_dy = recY - trackY;
		     float hitTrack_dR = sqrt(hitTrack_dx*hitTrack_dx + hitTrack_dy*hitTrack_dy);
		     float cellSize = sqrt(rhA[iL][iRc]) / 2.;

		     if(std::abs(hitTrack_dx) > cellSize+dXerr || std::abs(hitTrack_dy) > cellSize+dYerr || hitTrack_dR > nearestdR) continue;
		     nearestHit = iRc;
		     nearestdR = hitTrack_dR;
		   }
		 }//loop over recHits
		 if(matchBydR && nearestHit != -1){
		   float reciR = rhiR[iL][nearestHit];
		   std::pair<int, int> channelPhi = std::pair<int, int>(iL, std::abs(reciR));
		   if(okChannelsPhi.find(channelPhi) != okChannelsPhi.end()) {
		     okChannelsPhi[channelPhi] += 1;
		   }
		   else {  
		     okChannelsPhi[channelPhi] = 1;
		   }

		   h_MuonNHit[iL]->Fill(trkEtaOnLayer);
		   if(trkEtaOnLayer >= 1.55 && trkEtaOnLayer < 1.7){
		     hNum_etaRange[iL]->Fill(trkEtaOnLayer); 
		     //	     std::cout << " hNum_etaRange fillato trkEtaOnLayer = " << trkEtaOnLayer  << std::endl;
		   }
		   // else 	     std::cout << " trkEtaOnLayer = " << trkEtaOnLayer  << std::endl;
		   // //  std::cout << " fill num iL = " << iL << " trackIdx = " << trackIdx << std::endl;


		   if(foundFirstHit[trackIdx] == false){
		     matchedTrack[trackIdx] = 1;
		     startLMatchedTrack[trackIdx] = iL;
		   }
		   else matchedTrack[trackIdx] += 1;
		   foundFirstHit[trackIdx] = true;

		   int localBin = h_MuonNHit[iL]->GetXaxis()->FindBin(trkEtaOnLayer);
		   if(h_MuonNHit[iL]->GetBinContent(localBin) > h_MuonCross[iL]->GetBinContent(localBin)){
		     std::cout << " event a = " << ij << std::endl;
		     return 100;
		   }
		 
		   
		   //here rate for trk
		   if(!alreadyFilled){
		     nMuons_vsEta[iL]->Fill(std::abs(trkEtaOnLayer));
		     alreadyFilled = true;
		   }

		   mipPerPhi[channelPhi].push_back(rhMip[iL][nearestHit]);
		   continue;
		   /////
		 }// if matchBydR is ok

		 if(matchedTrackL[iL][trackIdx] == 0){
		   int iRc = -1;
		   for(auto iR : rhX[iL]){
		     ++iRc;
		     if(trkEta[iL][iTc] * rhZ[iL][iRc] < 0.) {
		       continue;
		     }
		     h_d_iR[iL]->Fill(rhiR[iL][iRc] - trkiR[iL][iTc]);
		     h_d_iPhi[iL]->Fill(rhiPhi[iL][iRc] - trkiP[iL][iTc]);
		     float  lastBinCenter = h_dRep[iL]->GetBinCenter(h_dRep[iL]->GetNbinsX()-2);
		     //	     h_dRep[iL]->Fill(lastBinCenter);
		     h_d_iR_vs_iPhi[iL]->Fill(rhiPhi[iL][iRc] - trkiP[iL][iTc], rhiR[iL][iRc] - trkiR[iL][iTc]);
		   }
		 }


	       }//track
	     }// layer
	 
 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 


	     //     continue;
	     int thevar =3;
	     for(int nHits = 1; nHits<5; ++nHits){
	       if (nHits != thevar) continue;
	     foundFirstHit.clear();
	     x1stHit.clear();
	     y1stHit.clear();
	     xPrevHit.clear();
	     yPrevHit.clear();
	     for(int iL=0; iL < 14; ++iL){
	       int iTc = -1;

	       for(auto iT : trkX[iL]){
		 ++iTc;
		 int trackIdx = trkN[iL][iTc];

		 float trackY = trkY[iL][iTc];
		 float trkR = sqrt(iT*iT + trackY*trackY);
		 float trkEtaOnLayer = asinh(zVal[iL]/trkR);

		 if(trkEtaOnLayer > 1.55 && trkEtaOnLayer < 1.7 && matchedTrack[trackIdx] != 0) numHits_etaRange[iL]->Fill(matchedTrack[trackIdx]);

		 //std::cout << " nHits = " << nHits << " matchedTrack[trackIdx] = " << matchedTrack[trackIdx] << " startLMatchedTrack[trackIdx] = " << startLMatchedTrack[trackIdx] << " iL = " << iL << std::endl;

		 //questo ok ma conta ultimi layer per max 1, 2, 3 hits...
		 if( (nHits > matchedTrack[trackIdx] &&  (matchedTrack[trackIdx] != (14 - startLMatchedTrack[trackIdx]) || iL < startLMatchedTrack[trackIdx])) || matchedTrack[trackIdx] == 0){
		   continue;
		 }
		 

		 //if(nHits > matchedTrack[trackIdx]) continue;
		 if(matchedTrack[trackIdx] < thevar) continue;
		 //	 std::cout << " ok " << std::endl;

		 //just using the back 8 layers
		 if(nHits > matchedTrackLast8L[trackIdx]) continue;

		 //if(trkGeneralPt[iL][iTc] < 10.) continue;
		 if(trkChi[iL][iTc] > 1.5) continue;

		 float dXerr = trkXerr[iL][iTc];
		 float dYerr = trkYerr[iL][iTc];	     
		  
		 if(matchByID_error){
		   if(dXerr > errorMax[iL] || dYerr > errorMax[iL])
		     continue;
		 }

		 if(matchByID_p009){
		   float dPhi = std::abs(trkCPhi[iL][iTc] - trkPhi[iL][iTc]);
		   if(dPhi > pigreco) dPhi -= 2* pigreco;
		   if(std::abs(trkCEta[iL][iTc] - trkEta[iL][iTc]) > edgeID_value || std::abs(dPhi) > edgeID_value ){
		     continue;
		   }
		 }
		 if(foundFirstHit.find(trackIdx) == foundFirstHit.end()) foundFirstHit[trackIdx] = false;

		 h_MuonCrossK[nHits][iL]->Fill(trkEtaOnLayer); 
		 if(trkEtaOnLayer > 1.55 && trkEtaOnLayer < 1.7) hDen_etaRange_3H[iL]->Fill(trkEtaOnLayer);
		 //std::cout << " fill den K iL = " << iL << " trackIdx = " << trackIdx << std::endl;

		 //if(iL < startLMatchedTrack[trackIdx]) continue;

		 bool alreadyFilled = false;
		 bool alreadyFilled_3H = false;
		 int iRc = -1;

		 int nearestHit = -1;
		 float nearestdR = 99.;
		 for(auto iR : rhX[iL]){
		   ++iRc;
		   if(trkEta[iL][iTc] * rhZ[iL][iRc] < 0.) {
		     //std::cout << " opposite side " << std::endl;
		     continue;
		   }

		   float dX = iR - iT;
		   float recY = rhY[iL][iRc];
		   float dY = recY - trackY;
		   float dR2 = dX*dX + dY*dY;
		   float dR = sqrt(dR2);

		   if(matchByID){
		     if(rhiR[iL][iRc] != trkiR[iL][iTc] || rhiPhi[iL][iRc] != trkiP[iL][iTc]) continue;
		   
		     if(matchByID_error){
		       if(dXerr > errorMax[iL] || dYerr > errorMax[iL])
			 continue;
		     }

		     if(matchByID_p009){
		       float dPhi = std::abs(trkPhi[iL][iTc] - rhPhi[iL][iRc]);
		       if(dPhi > pigreco) dPhi -= 2* pigreco;
		       if(std::abs(trkEta[iL][iTc] - rhEta[iL][iRc]) > edgeID_value || std::abs(dPhi) > edgeID_value ){
			 continue;
		       }
		     }

		     float reciR = rhiR[iL][iRc];
		     std::pair<int, int> channelPhi = std::pair<int, int>(iL, std::abs(reciR));
		     if(nHits == thevar){
		       if(okChannelsPhi_3H.find(channelPhi) != okChannelsPhi.end()) {
			 okChannelsPhi_3H[channelPhi] += 1;
		       }
		       else {
			 okChannelsPhi_3H[channelPhi] = 1;
		       }
		       std::pair<int, int> channelRh = std::pair<int, int>(reciR, rhiPhi[iL][iRc]);
		       std::pair<int, std::pair<int, int>> channelRh_L = std::pair<int, std::pair<int, int>>(iL, channelRh);
		       if(okChannels_3H.find(channelRh_L) != okChannels_3H.end()) {
			 okChannels_3H[channelRh_L] += 1;
			 mipPerRh_3H[channelRh_L].push_back(rhMip[iL][iRc]); 
		       }
		       else {
			 okChannels_3H[channelRh_L] = 1;
			 mipPerRh_3H[channelRh_L].push_back(rhMip[iL][iRc]); 
		       }
		     }

		     h_MuonNHitK[nHits][iL]->Fill(trkEtaOnLayer);
		     //std::cout << " fill num K iL = " << iL << " trackIdx = " << trackIdx << std::endl;
		     
		     if(trkEtaOnLayer >= 1.55 && trkEtaOnLayer < 1.7){
		       hNum_etaRange_3H[iL]->Fill(trkEtaOnLayer);
		     }
		     

		     int localBin = h_MuonNHitK[nHits][iL]->GetXaxis()->FindBin(trkEtaOnLayer);
		     if(h_MuonNHitK[nHits][iL]->GetBinContent(localBin) > h_MuonCrossK[nHits][iL]->GetBinContent(localBin)){
		       std::cout << " event a = " << ij << std::endl;
		       return 100;
		     }
		     if(iL == 0 && nHits == 13 && h_MuonNHitK[nHits][iL]->GetBinContent(localBin) == h_MuonNHitK[1][iL]->GetBinContent(localBin) &&
			h_MuonNHitK[nHits][iL]->GetBinContent(localBin) > h_MuonNHitK[nHits-1][iL]->GetBinContent(localBin)){
		       std::cout << " event b = " << ij << std::endl;
		       return 100;
		     }

		     if(nHits == thevar){
		       //here rate for trk with 3 Hits                                                                                                                                                          
		       if(!alreadyFilled_3H){
			 nMuons_3H_vsEta[iL]->Fill(std::abs(trkEtaOnLayer));
			 alreadyFilled_3H = true;
		       }
		       mipPerPhi_3H[channelPhi].push_back(rhMip[iL][iRc]);
		       //	       h_muonSeg_ok->Fill(trkSegC[iL][iTc]);
		     }
		     //else if (nHits < 3) h_muonSeg_bad->Fill(trkSegC[iL][iTc]);
		   }
		   else if(matchBydR){
		     float hitTrack_dx = iR - iT;
		     float hitTrack_dy = recY - trackY;
		     float hitTrack_dR = sqrt(hitTrack_dx*hitTrack_dx + hitTrack_dy*hitTrack_dy);
		     float cellSize = sqrt(rhA[iL][iRc]) / 2.;

		     if(std::abs(hitTrack_dx) > cellSize+dXerr || std::abs(hitTrack_dy) > cellSize+dYerr || hitTrack_dR > nearestdR) continue;
		     nearestHit = iRc;
		     nearestdR = hitTrack_dR;
		   }
		 }//loop over recHits
		 if(matchBydR && nearestHit != -1){
		   h_MuonNHitK[nHits][iL]->Fill(trkEtaOnLayer);
		   if(trkEtaOnLayer >= 1.55 && trkEtaOnLayer < 1.7){
		     hNum_etaRange_3H[iL]->Fill(trkEtaOnLayer);
		   }
		 }
	       }// trk loop
	     }// layer loop

	     }// loop on nhits
	   }// events;

	   std::cout << " okChannels.size = " << okChannels.size() << " okChannelPhi.size() = " << okChannelsPhi.size() << std::endl;

	   TFile pippo("pippo.root", "recreate");
	   pippo.cd();
	   h_muonSeg_ok->Write();
	   h_muonSeg_bad->Write();
	   for(int ij =0 ; ij<14; ++ij) numHits_etaRange[ij]->Write(Form("numHits_etaRange_L%d", ij));

	   pippo.Close();


	   //rates
	   for(int ij=0; ij<14; ++ij){
	     float binW = (4. - 0. ) / 400.;
	     for(int iB=2; iB<nMuonsPercm2_vsEta[ij]->GetNbinsX()-1; ++iB){
	      float rMax = zVal[ij]/sinh(0. + binW * iB);     
	      float rMin = zVal[ij]/sinh(0. + binW * (iB+1));     
	      float area = 2. * (rMax*rMax - rMin*rMin) * pigreco;
	      // std::cout << " rMin = " << rMin << " rMax = " << rMax << " etaL = " << 1 + binW * (iB) << " etaH  = " << (1 + binW * (iB+1)) << " area = " << area << std::endl;
	      nMuonsPercm2_vsEta[ij]->SetBinContent(iB, nMuons_vsEta[ij]->GetBinContent(iB)/area);
	      ratePercm2_vsEta[ij]->SetBinContent(iB, nMuons_vsEta[ij]->GetBinContent(iB)/area/timeEqui);
	      cm2_vsEta[ij]->SetBinContent(iB, area);

	      nMuonsPercm2_vsEta_All[ij]->SetBinContent(iB, nMuons_vsEta_All[ij]->GetBinContent(iB)/area);
	      ratePercm2_vsEta_All[ij]->SetBinContent(iB, nMuons_vsEta_All[ij]->GetBinContent(iB)/area/timeEqui);

	      nMuonsPercm2_3H_vsEta[ij]->SetBinContent(iB, nMuons_3H_vsEta[ij]->GetBinContent(iB)/area);
	      ratePercm2_3H_vsEta[ij]->SetBinContent(iB, nMuons_3H_vsEta[ij]->GetBinContent(iB)/area/timeEqui);
	     }
	     
	     TCanvas* tmuo = new TCanvas();
	     tmuo->cd();
	     nMuonsPercm2_vsEta[ij]->GetXaxis()->SetTitle(Form("#eta layer %d", ij+firstLayer));
	     nMuonsPercm2_vsEta[ij]->GetYaxis()->SetTitle("n Muons / cm2");
	     nMuonsPercm2_vsEta[ij]->Draw();
	     tmuo->Print(Form("plotdir/nMuonsPercm2_vsEta_L%d.png", ij+firstLayer), "png");
	     //    tmuo->Print(Form("plotdir/nMuonsPercm2_vsEta_L%d.root", ij+firstLayer), "root");
	     
	     nMuons_vsEta[ij]->GetXaxis()->SetTitle(Form("#eta layer %d", ij+firstLayer));
	     nMuons_vsEta[ij]->Draw();
	     tmuo->Print(Form("plotdir/nMuons_vsEta_L%d.png", ij+firstLayer), "png");

	     ratePercm2_vsEta[ij]->GetYaxis()->SetRangeUser(0., 0.5);
	     ratePercm2_vsEta[ij]->GetXaxis()->SetTitle(Form("#eta layer %d", ij+firstLayer));
	     ratePercm2_vsEta[ij]->GetYaxis()->SetTitle("n Muons / cm2 / s");
	     ratePercm2_vsEta[ij]->Draw();
	     tmuo->Print(Form("plotdir/ratePercm2_vsEta_L%d.png", ij+firstLayer), "png");  

	     gPad->SetLogy();
	     cm2_vsEta[ij]->GetYaxis()->SetRangeUser(0.5, 5.e4);
	     cm2_vsEta[ij]->GetXaxis()->SetTitle(Form("#eta layer %d", ij+firstLayer));
	     cm2_vsEta[ij]->GetYaxis()->SetTitle("cm2");
	     cm2_vsEta[ij]->Draw();
	     tmuo->Print(Form("plotdir/cm2_vsEta_L%d.png", ij+firstLayer), "png");

	     TCanvas* tmuo_3H = new TCanvas();
	     tmuo_3H->cd();
	     nMuonsPercm2_3H_vsEta[ij]->GetXaxis()->SetTitle(Form("#eta layer %d", ij+firstLayer));
	     nMuonsPercm2_3H_vsEta[ij]->GetYaxis()->SetTitle("n Muons / cm2");
	     nMuonsPercm2_3H_vsEta[ij]->Draw();
	     tmuo_3H->Print(Form("plotdir/nMuonsPercm2_3H_vsEta_L%d.png", ij+firstLayer), "png");
	     //    tmuo->Print(Form("plotdir/nMuonsPercm2_vsEta_L%d.root", ij+firstLayer), "root");
	     
	     nMuons_3H_vsEta[ij]->GetXaxis()->SetTitle(Form("#eta layer %d", ij+firstLayer));
	     nMuons_3H_vsEta[ij]->Draw();
	     tmuo_3H->Print(Form("plotdir/nMuons_3H_vsEta_L%d.png", ij+firstLayer), "png");

	     ratePercm2_3H_vsEta[ij]->GetYaxis()->SetRangeUser(0., 0.5);
	     ratePercm2_3H_vsEta[ij]->GetXaxis()->SetTitle(Form("#eta layer %d", ij+firstLayer));
	     ratePercm2_3H_vsEta[ij]->GetYaxis()->SetTitle("n Muons / cm2 / s");
	     ratePercm2_3H_vsEta[ij]->Draw();
	     tmuo_3H->Print(Form("plotdir/ratePercm2_3H_vsEta_L%d.png", ij+firstLayer), "png");  


	     
	     TCanvas* tmuoAll = new TCanvas();
	     tmuoAll->cd();
	     nMuonsPercm2_vsEta_All[ij]->GetXaxis()->SetTitle(Form("#eta layer %d", ij+firstLayer));
	     nMuonsPercm2_vsEta_All[ij]->GetYaxis()->SetTitle("n Muons / cm2");
	     nMuonsPercm2_vsEta_All[ij]->Draw();
	     tmuoAll->Print(Form("plotdir/nMuonsPercm2_vsEta_All_L%d.png", ij+firstLayer), "png");
	     //    tmuo->Print(Form("plotdir/nMuonsPercm2_vsEta_L%d.root", ij+firstLayer), "root");
	     
	     nMuons_vsEta_All[ij]->GetXaxis()->SetTitle(Form("#eta layer %d", ij+firstLayer));
	     nMuons_vsEta_All[ij]->Draw();
	     tmuoAll->Print(Form("plotdir/nMuons_vsEta_All_L%d.png", ij+firstLayer), "png");

	     ratePercm2_vsEta_All[ij]->GetYaxis()->SetRangeUser(0., 20.);
	     ratePercm2_vsEta_All[ij]->GetXaxis()->SetTitle(Form("#eta layer %d", ij+firstLayer));
	     ratePercm2_vsEta_All[ij]->GetYaxis()->SetTitle("n Muons / cm2 / s");
	     ratePercm2_vsEta_All[ij]->Draw();
	     tmuoAll->Print(Form("plotdir/ratePercm2_vsEta_All_L%d.png", ij+firstLayer), "png");  
	   }
	   
	   
	   //   return 10;

	   TCanvas* tChi = new TCanvas();
	   tChi->cd();
	   h_TrkChi_all->GetXaxis()->SetTitle("trk chi all"); 
	   h_TrkChi_all->Draw();
	   tChi->Print("plotdir/h_TrkChi_all.png", "png");
	   tChi->Print("plotdir/h_TrkChi_all.root", "root");

	   h_TrkQ_all->GetXaxis()->SetTitle("trk Q all"); 
	   h_TrkQ_all->Draw();
	   tChi->Print("plotdir/h_TrkQ_all.png", "png");
	   tChi->Print("plotdir/h_TrkQ_all.root", "root");

	   h_TrkChi_passed->GetXaxis()->SetTitle("trk chi all"); 
	   h_TrkChi_passed->Draw();
	   tChi->Print("plotdir/h_TrkChi_passed.png", "png");
	   tChi->Print("plotdir/h_TrkChi_passed.root", "root");

	   h_TrkQ_passed->GetXaxis()->SetTitle("trk Q all"); 
	   h_TrkQ_passed->Draw();
	   tChi->Print("plotdir/h_TrkQ_passed.png", "png");
	   tChi->Print("plotdir/h_TrkQ_passed.root", "root");

	   for(int ij=0; ij<14; ++ij){
	     tChi->cd();
	     h_trkP[ij]->GetXaxis()->SetTitle(Form("track p on layer %d", ij+firstLayer));
	     h_trkP[ij]->Draw();
	     tChi->Print(Form("plotdir/h_trkP_L%d.png", ij+firstLayer), "png");

	     gPad->SetLogy();
	     tChi->cd();
	     h_trkPt[ij]->GetXaxis()->SetTitle(Form("track pt on layer %d", ij+firstLayer));
	     h_trkPt[ij]->Draw();
	     tChi->Print(Form("plotdir/h_trkPt_L%d.png", ij+firstLayer), "png");

	     tChi->cd();
	     h_dtrkP[ij]->GetXaxis()->SetTitle(Form("dp wrt previous layer (layer %d) (MeV)", ij+firstLayer));
	     h_dtrkP[ij]->Draw();
	     tChi->Print(Form("plotdir/h_dtrkP_L%d.png", ij+firstLayer), "png");

	   }
	   //   return;
	   //   if(1 == 2){
	   TCanvas* tMuo = new TCanvas();
	   tMuo->cd();
	   h_etaMuon->GetXaxis()->SetTitle("muon #eta");
	   h_etaMuon->Draw();
	   tMuo->Print("plotdir/etaMuon_all_L37.png", "png");
	   tMuo->Print("plotdir/etaMuon_all_L37.root", "root");

	   gPad->SetLogy();
	   h_ptMuon->GetXaxis()->SetTitle("muon pT");
	   h_ptMuon->Draw();
	   tMuo->Print("plotdir/ptMuon_all_L37.png", "png");
	   tMuo->Print("plotdir/ptMuon_all_L37.root", "root");
	   //   }



	   //about efficiencies
	   TGraph* etaRange_eff = new TGraph();
	   TGraph* etaRange_3H_eff = new TGraph();

	   /*
	   TH1F* MuonNHit_eff[14];
	   TH1F* MuonNHitK_eff[14][14];

	   TGraph* MuonNHit_eff_tg[14];
	   TGraph* MuonNHitK_eff_tg[14][14];
	   */
	   for(int ij=0; ij<14; ++ij){
	     TCanvas* tc = new TCanvas();
	     tc->cd();

	     if(ij == 0) {
	       etaRange_eff->SetPoint(0, 0, 0);
	       etaRange_3H_eff->SetPoint(0, 0, 0);
	     }
	     etaRange_eff->SetPoint(etaRange_eff->GetN()+1, ij+firstLayer, hNum_etaRange[ij]->GetBinContent(1)/hDen_etaRange[ij]->GetBinContent(1));
	     hNum_etaRange[ij]->Draw();
	     tc->Print(Form("plotdir/hNum_etaRange_L%d.png", ij+firstLayer), "png");
	     hDen_etaRange[ij]->Draw();
	     tc->Print(Form("plotdir/hDen_etaRange_L%d.png", ij+firstLayer), "png");
	     etaRange_3H_eff->SetPoint(etaRange_3H_eff->GetN()+1, ij+firstLayer, hNum_etaRange_3H[ij]->GetBinContent(1)/hDen_etaRange_3H[ij]->GetBinContent(1));
	   }
	   etaRange_eff->SetPoint(etaRange_eff->GetN()+1, 60, 1.5);
	   TCanvas* tgEta = new TCanvas();
	   tgEta->cd();
	   etaRange_eff->GetXaxis()->SetTitle("layer");
	   etaRange_eff->GetXaxis()->SetRangeUser(30., 55);
	   etaRange_eff->Draw("ap");
	   tgEta->Print("plotdir/tgEta.png", "png");

	   etaRange_3H_eff->SetPoint(etaRange_3H_eff->GetN()+1, 60, 1.5);
	   tgEta->cd();
	   etaRange_3H_eff->GetXaxis()->SetTitle("layer");
	   etaRange_3H_eff->GetXaxis()->SetRangeUser(30., 55);
	   etaRange_3H_eff->Draw("ap");
	   tgEta->Print("plotdir/tgEta_3H.png", "png");


	   //topological plots
	   gStyle->SetOptStat(1);
	   for(int ij=0; ij<14; ++ij){
	     TCanvas* tc = new TCanvas();
	     tc->cd();
	     h2_YvsX[ij]->GetXaxis()->SetTitle(Form("layer %d", ij));
	     h2_YvsX[ij]->Draw("colz");
	     tc->Print(Form("plotdir/h2_YvsX_L%d.png", ij), "png");

	     h2_iRvsiPhi[ij]->GetXaxis()->SetTitle(Form("|iR| layer %d", ij));
	     h2_iRvsiPhi[ij]->GetYaxis()->SetTitle("iPhi");
	     h2_iRvsiPhi[ij]->Draw("colz");
	     tc->Print(Form("plotdir/h2_iRvsiPhi_L%d.png", ij), "png");
	     //   continue;
	     h_dX[ij]->GetXaxis()->SetTitle(Form("dX layer %d", ij));
	     h_dX[ij]->Draw();
	     tc->Print(Form("plotdir/h_dX_L%d.png", ij), "png");

	     h_dY[ij]->GetXaxis()->SetTitle(Form("dY layer %d", ij));
	     h_dY[ij]->Draw();
	     tc->Print(Form("plotdir/h_dY_L%d.png", ij), "png");

	     h_dR[ij]->GetXaxis()->SetTitle(Form("dR layer %d", ij+firstLayer));
	     h_dR[ij]->Draw();
	     tc->Print(Form("plotdir/h_dR_L%d.png", ij+firstLayer), "png");
	     tc->Print(Form("plotdir/h_dR_L%d.root", ij+firstLayer), "root");

	     h_dRep[ij]->GetXaxis()->SetTitle(Form("dR(#eta, #phi) layer %d", ij+firstLayer));
	     h_dRep[ij]->Draw();
	     tc->Print(Form("plotdir/h_dRep_L%d.png", ij+firstLayer), "png");
	     tc->Print(Form("plotdir/h_dRep_L%d.root", ij+firstLayer), "root");
	     h_dRep[ij]->GetXaxis()->SetRangeUser(0., 0.05);
	     tc->Print(Form("plotdir/h_dRep_L%d_zoomIn.png", ij+firstLayer), "png");

	     h_d_iR[ij]->GetXaxis()->SetTitle(Form("d_iR(recHit - track) layer %d", ij+firstLayer));
	     h_d_iR[ij]->Draw();
	     tc->Print(Form("plotdir/h_d_iR_L%d.png", ij+firstLayer), "png");
	     tc->Print(Form("plotdir/h_d_iR_L%d.root", ij+firstLayer), "root");

	     h_d_iPhi[ij]->GetXaxis()->SetTitle(Form("d_iPhi(recHit - track) layer %d", ij+firstLayer));
	     h_d_iPhi[ij]->Draw();
	     tc->Print(Form("plotdir/h_d_iPhi_L%d.png", ij+firstLayer), "png");
	     tc->Print(Form("plotdir/h_d_iPhi_L%d.root", ij+firstLayer), "root");

	     h_d_iR_vs_iPhi[ij]->GetXaxis()->SetTitle(Form("d_iPhi(recHit - track) layer %d", ij+firstLayer));
	     h_d_iR_vs_iPhi[ij]->GetYaxis()->SetTitle(Form("d_iR(recHit - track) layer %d", ij+firstLayer));
	     h_d_iR_vs_iPhi[ij]->Draw("colz");
	     tc->Print(Form("plotdir/h_d_iR_vs_iPhi_L%d.png", ij+firstLayer), "png");
	     tc->Print(Form("plotdir/h_d_iR_vs_iPhi_L%d.root", ij+firstLayer), "root");


	     h_dRerr[ij]->GetXaxis()->SetTitle(Form("dRerror(x, y) layer %d", ij+firstLayer));
	     h_dRerr[ij]->Draw();
	     tc->Print(Form("plotdir/h_dRerr_L%d.png", ij+firstLayer), "png");
	     tc->Print(Form("plotdir/h_dRerr_L%d.root", ij+firstLayer), "root");
	     //     h_dRerr[ij]->GetXaxis()->SetRangeUser(0., 0.05);
	     tc->Print(Form("plotdir/h_dRerr_L%d.png", ij+firstLayer), "png");
	     tc->Print(Form("plotdir/h_dRerr_L%d.root", ij+firstLayer), "root");

	     h_dXerr[ij]->GetXaxis()->SetTitle(Form("dXerror (cm) layer %d", ij+firstLayer));
	     h_dXerr[ij]->Draw();
	     tc->Print(Form("plotdir/h_dXerr_L%d.png", ij+firstLayer), "png");
	     tc->Print(Form("plotdir/h_dXerr_L%d.root", ij+firstLayer), "root");
	     //     h_dRerr[ij]->GetXaxis()->SetRangeUser(0., 0.05);
	     tc->Print(Form("plotdir/h_dXerr_L%d.png", ij+firstLayer), "png");
	     tc->Print(Form("plotdir/h_dXerr_L%d.root", ij+firstLayer), "root");

	     h_dYerr[ij]->GetXaxis()->SetTitle(Form("dYerror (cm) layer %d", ij+firstLayer));
	     h_dYerr[ij]->Draw();
	     tc->Print(Form("plotdir/h_dYerr_L%d.png", ij+firstLayer), "png");
	     tc->Print(Form("plotdir/h_dYerr_L%d.root", ij+firstLayer), "root");
	     //     h_dRerr[ij]->GetXaxis()->SetRangeUser(0., 0.05);
	     tc->Print(Form("plotdir/h_dYerr_L%d.png", ij+firstLayer), "png");
	     tc->Print(Form("plotdir/h_dYerr_L%d.root", ij+firstLayer), "root");

	     h2_dR_vsChi[ij]->GetXaxis()->SetTitle(Form("dR layer %d", ij));
	     h2_dR_vsChi[ij]->GetYaxis()->SetTitle("muon trk chi2");
	     h2_dR_vsChi[ij]->Draw("colz");
	     tc->Print(Form("plotdir/h2_dR_vsChi_L%d.png", ij), "png");

	     h2_dRep_vsChi[ij]->GetXaxis()->SetTitle(Form("dR(#eta, #phi) layer %d", ij));
	     h2_dRep_vsChi[ij]->GetYaxis()->SetTitle("muon trk chi2");
	     h2_dRep_vsChi[ij]->Draw("colz");
	     tc->Print(Form("plotdir/h2_dRep_vsChi_L%d.png", ij), "png");


	     TCanvas* tc2 = new TCanvas();
	     tc2->cd();
	     h2_dR_vsTrkQ[ij]->GetXaxis()->SetTitle(Form("dR layer %d", ij));
	     h2_dR_vsTrkQ[ij]->GetYaxis()->SetTitle("muon trk quality");
	     h2_dR_vsTrkQ[ij]->Draw("colz");
	     tc2->Print(Form("plotdir/h2_dR_vsTrkQ_L%d.png", ij), "png");

	     tc2->cd();
	     h2_iR_vsEta_size[ij]->GetXaxis()->SetTitle(Form("eta layer %d", ij+firstLayer));
	     h2_iR_vsEta_size[ij]->GetYaxis()->SetTitle(Form("|iR| layer %d", ij+firstLayer));
	     h2_iR_vsEta_size[ij]->Draw("colz");
	     h2_iR_vsEta_size[ij]->SetMarkerStyle(7);
	     tc2->Print(Form("plotdir/h2_iR_vsEta_size_L%d.png", ij+firstLayer), "png");

	     tc2->cd();
	     h2_iR_vsEta[ij]->GetXaxis()->SetTitle(Form("eta layer %d", ij+firstLayer));
	     h2_iR_vsEta[ij]->GetYaxis()->SetTitle(Form("|iR| layer %d", ij+firstLayer));
	     h2_iR_vsEta[ij]->Draw(""); //"colz");
	     h2_iR_vsEta[ij]->SetMarkerStyle(7);
	     tc2->Print(Form("plotdir/h2_iR_vsEta_L%d.png", ij+firstLayer), "png");

	     tc2->cd();
	     h2_iR_vsR[ij]->GetXaxis()->SetTitle(Form("R layer %d", ij+firstLayer));
	     h2_iR_vsR[ij]->GetYaxis()->SetTitle(Form("|iR| layer %d", ij+firstLayer));
	     h2_iR_vsR[ij]->Draw(""); //"colz");
	     h2_iR_vsR[ij]->SetMarkerStyle(7);
	     tc2->Print(Form("plotdir/h2_iR_vsR_L%d.png", ij+firstLayer), "png");

	   }


	   //now fitting the MIP peaks

	   TLatex tL;
	   tL.SetNDC();
	   tL.SetTextSize(0.05);
	   tL.SetTextFont(132);

	  std::cout << " now plotting " << std::endl;
	  //  return;
	  // now integrated vs Phi
	  //  TH1F* h_Ene_Phi[400];
	  TH1F* h_Mip_Phi[400];
	  // TH1F* h_SoN_Phi[400];
	  // TH2F* h_MipvsP_Phi[400];
	  TH1F* h_Count_Phi = new TH1F("h_Count_Phi", "", 500, 0., 2000.);

	  TH2F* h2_iRvsLayer_MIP = new TH2F("h2_iRvsLayer_MIP", "", 15, 36, 51, 50, 0., 50.);
	  TH2F* h2_iRvsLayer_MIPerr = new TH2F("h2_iRvsLayer_MIPerr", "", 15, 36, 51, 50, 0., 50.);
	  TH1F* h_minNcounts = new TH1F("h_minNcounts", "", 400, 0., 800.);
	  TH1F* h_MPV_values = new TH1F("h_MPV_values", "", 500., 0., 5.);

	  TH1F* h2_MIP_vs_Layer = new TH1F("h2_MIP_vs_Layer", "", 15, 36.5, 50.5); //, 500, 0., 5.);
	  TH1F* h_MPV_layer_values = new TH1F("h_MPV_layer_values", "", 500., 0., 5.);
	  TH1F* h_Mip_Layer[14];
	  for(int ij=0; ij<14; ++ij){
	    h_Mip_Layer[ij] = new TH1F(Form("h_Mip_Layer_%d", ij+firstLayer), "", 110, -2., 20.);
	  }

	  TH2F* h2_iRvsLayer_MIP_3H = new TH2F("h2_iRvsLayer_MIP_3H", "", 15, 36, 51, 50, 0., 50.);
	  TH1F* h_Mip_Phi_3H[400];
	  TH1F* h_MPV_values_3H = new TH1F("h_MPV_values_3H", "", 500., 0., 5.);
	  TH1F* h2_MIP_vs_Layer_3H = new TH1F("h2_MIP_vs_Layer_3H", "", 15, 36.5, 50.5); //, 500, 0., 5.);
	  TH1F* h_MPV_layer_values_3H = new TH1F("h_MPV_layer_values_3H", "", 500., 0., 5.);
	  TH1F* h_Mip_Layer_3H[14];
	  for(int ij=0; ij<14; ++ij){
	    h_Mip_Layer_3H[ij] = new TH1F(Form("h_Mip_Layer_3H_%d", ij+firstLayer), "", 110, -2., 20.);
	  }


	  int savedCout = 0;
	  int iC = -1;

	  TF1* hfithisto = new TF1("hfithisto", "landaun",0., 10);
	  hfithisto->SetLineColor(kBlue);
	  for(auto ic : okChannelsPhi){
	    ++iC;
	    
	    int iL = ic.first.first;
	    int iR = ic.first.second;

	    int iVal = ic.second;
	    std::cout << " iVal = " << iVal << std::endl;
	    //    if(iVal > 500) std::cout << " iL = " << iL << " iR = " << iR << " iPhi = " << iPhi << " iVal = " << iVal << std::endl;

	    int countFilled = 0;
	    int rejectedBy1stdR = 0;
	    int rejectedByPrevdR = 0;

	    TH1F* dummyMIP = new TH1F("dummyMIP", "", 110, -2., 20.);
	    for(int ij=0; ij<iVal; ++ij){
	      if(mipPerPhi[ic.first].size() != iVal) std::cout << " PROBLEM!!!  " << std::endl;
	      float valueMIP = mipPerPhi[ic.first][ij];
	      if(valueMIP == -1 ) ++rejectedByPrevdR;
	      else if(valueMIP == -2 ) ++rejectedBy1stdR;
	      else {
		dummyMIP->Fill(mipPerPhi[ic.first][ij]);
		++countFilled;
		//	std::cout << " Filling layer = " << iL << std::endl;
		h_Mip_Layer[iL]->Fill(mipPerPhi[ic.first][ij]);
	      }
	    }
	    h_Count_Phi->Fill(countFilled);

	    float yMax;
	    getXmax(dummyMIP, yMax);
	    hfithisto->SetParameters(yMax, 1., 0.1);
	    hfithisto->SetParLimits(0, yMax/2., yMax * 1.2);
	    hfithisto->SetParLimits(1, 0., 10.);
	    hfithisto->SetParLimits(2, 0., 0.3);
	    //std::cout << " histo entries = " << dummyMIP->GetEntries() << std::endl;
	    //TFitResultPtr r = dummyMIP->Fit("hfithisto", "RBS");
	    TFitResultPtr r = dummyMIP->Fit("hfithisto", "q");
	    if(r != 0) continue;
	    float meanV = hfithisto->GetParameter(1);
	    delete dummyMIP;
	    if(meanV < 0) continue;
	    h_minNcounts->Fill(iVal);
	    h2_iRvsLayer_MIP->Fill(iL+firstLayer, iR, hfithisto->GetParameter(1));
	    h2_iRvsLayer_MIPerr->Fill(iL+firstLayer, iR, hfithisto->GetParError(1));
	    h_MPV_values->Fill(hfithisto->GetParameter(1));

	    //    continue;
	    //    if(iVal < 500) continue;

	    std::cout << " savedCout = " << savedCout << std::endl;
	    if(savedCout < 400){
	      //      h_Ene_Phi[savedCout] = new TH1F(Form("h_Ene_Phi_%d_%d", iL, iR), "", 1000, 0., 1.);
	      h_Mip_Phi[savedCout] = new TH1F(Form("h_Mip_Phi_%d_%d", iL, iR), "", 110, -2., 20.);
	      // h_SoN_Phi[savedCout] = new TH1F(Form("h_SoN_Phi_%d_%d", iL, iR), "", 100, 0., 1.e3);
	      // h_MipvsP_Phi[savedCout] = new TH2F(Form("h_MipvsP_Phi_%d_%d", iL, iR), "", 100, 0., 1.e3, 100, 0., 100.);

	      for(int ij=0; ij<iVal; ++ij){
		float valueMIP = mipPerPhi[ic.first][ij];
		//	h_Ene_Phi[savedCout]->Fill(enePerPhi[ic.first][ij]);
		if(valueMIP >= 0) h_Mip_Phi[savedCout]->Fill(mipPerPhi[ic.first][ij]);
		//std::cout << " energy = " << enePerRh[ic.first][ij]  << std::endl;
		// h_SoN_Phi[savedCout]->Fill(sonPerPhi[ic.first][ij]);
		// h_MipvsP_Phi[savedCout]->Fill(mipPerPhi[ic.first][ij], momPerPhi[ic.first][ij]);
	      }

	      TCanvas* tcM = new TCanvas();
	      tcM->cd();

	      h_Mip_Phi[savedCout]->GetXaxis()->SetTitle("MIP");
	      h_Mip_Phi[savedCout]->Draw();
	      hfithisto->SetParameters(yMax, 1., 0.1);
	      hfithisto->SetParLimits(0, yMax/2., yMax * 1.2);
	      hfithisto->SetParLimits(1, 0., 10.);
	      hfithisto->SetParLimits(2, 0., 0.3);
	      //h_Mip_Phi[savedCout]->Fit("hfithisto", "RB");
	      h_Mip_Phi[savedCout]->Fit("hfithisto", "q");

	      //      tL.DrawLatex(0.55, 0.70, Form("by Chi2 = %d", rejectedByChi));
	      // tL.DrawLatex(0.55, 0.60, Form("by 1st dR = %d", rejectedBy1stdR));
	      // tL.DrawLatex(0.55, 0.50, Form("by prev dR = %d", rejectedByPrevdR));
	      // tL.DrawLatex(0.55, 0.40, Form("found = %d", countFilled));

	      tcM->Print(Form("plotdir/cellsPhi/h_Mip_Phi_%d_%d.png", iL+firstLayer, iR), "png");
	      if(iL == 0) tcM->Print(Form("plotdir/cells/h_Mip_Phi_%d_%d.root", iL+firstLayer, iR), "root");

	      ++savedCout;
	    }
	  }

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

	  int iCounts = -1;
	  for(auto ic : okChannels_3H){
	    ++iCounts;

	    int iL = ic.first.first;
	    int iR = ic.first.second.first;
	    int iP = ic.first.second.second;
	    int iVal = ic.second;
	    for(int ij=0; ij<iVal; ++ij){
	      if(mipPerRh_3H[ic.first].size() != iVal) std::cout << " PROBLEM_MIP_3H!!!  " << std::endl;
	      std::cout << " iL = " << iL << " iR = " << iR << " iPhi = " << iP << " MIP = " << mipPerRh_3H[ic.first][ij] << std::endl;
	      MIP_layer = iL;
	      MIP_iR = iR;
	      MIP_iPhi = iP;
	      MIP_val = mipPerRh_3H[ic.first][ij];
	      MIP_SoN = 2.5;
	      tt->Fill();
	    }
	  }

	  tt->Write();
	  outMIP.Close();



  //3hits
  savedCout = 0;
  iC = -1;
  for(auto ic : okChannelsPhi_3H){
    ++iC;
    
    int iL = ic.first.first;
    int iR = ic.first.second;

    int iVal = ic.second;
    std::cout << " iVal = " << iVal << std::endl;
    //    if(iVal > 500) std::cout << " iL = " << iL << " iR = " << iR << " iPhi = " << iPhi << " iVal = " << iVal << std::endl;

    int rejectedBy1stdR = 0;
    int rejectedByPrevdR = 0;

    TH1F* dummyMIP = new TH1F("dummyMIP", "", 110, -2., 20.);
    for(int ij=0; ij<iVal; ++ij){
      if(mipPerPhi_3H[ic.first].size() != iVal) std::cout << " PROBLEM_3H!!!  " << std::endl;
      float valueMIP = mipPerPhi_3H[ic.first][ij];
      if(valueMIP == -1 ) ++rejectedByPrevdR;
      else if(valueMIP == -2 ) ++rejectedBy1stdR;
      else {
	dummyMIP->Fill(mipPerPhi_3H[ic.first][ij]);
	//	std::cout << " Filling layer = " << iL << std::endl;
	h_Mip_Layer_3H[iL]->Fill(mipPerPhi_3H[ic.first][ij]);
      }
    }

    float yMax;
    getXmax(dummyMIP, yMax);
    hfithisto->SetParameters(yMax, 1., 0.1);
    hfithisto->SetParLimits(0, yMax/2., yMax * 1.2);
    hfithisto->SetParLimits(1, 0., 10.);
    hfithisto->SetParLimits(2, 0., 0.3);
    //std::cout << " histo entries = " << dummyMIP->GetEntries() << std::endl;
    //TFitResultPtr r = dummyMIP->Fit("hfithisto", "RBS");
    TFitResultPtr r = dummyMIP->Fit("hfithisto", "q");
    if(r != 0) continue;
    float meanV = hfithisto->GetParameter(1);
    delete dummyMIP;
    if(meanV < 0) continue;
    h_minNcounts->Fill(iVal);
    h2_iRvsLayer_MIP_3H->Fill(iL+firstLayer, iR, hfithisto->GetParameter(1));
    h_MPV_values_3H->Fill(hfithisto->GetParameter(1));

    //    continue;
    //    if(iVal < 500) continue;

    std::cout << " savedCout = " << savedCout << std::endl;
    if(savedCout < 400){
      h_Mip_Phi_3H[savedCout] = new TH1F(Form("h_Mip_Phi_3H_%d_%d", iL, iR), "", 110, -2., 20.);
      for(int ij=0; ij<iVal; ++ij){
	float valueMIP = mipPerPhi_3H[ic.first][ij];
	if(valueMIP >= 0) h_Mip_Phi_3H[savedCout]->Fill(mipPerPhi_3H[ic.first][ij]);
      }
      
      TCanvas* tcM = new TCanvas();
      tcM->cd();

      h_Mip_Phi_3H[savedCout]->GetXaxis()->SetTitle("MIP");
      h_Mip_Phi_3H[savedCout]->Draw();
      hfithisto->SetParameters(yMax, 1., 0.1);
      hfithisto->SetParLimits(0, yMax/2., yMax * 1.2);
      hfithisto->SetParLimits(1, 0., 10.);
      hfithisto->SetParLimits(2, 0., 0.3);
      //h_Mip_Phi_3H[savedCout]->Fit("hfithisto", "RB");
      h_Mip_Phi_3H[savedCout]->Fit("hfithisto", "q");

      tcM->Print(Form("plotdir/cellsPhi/h_Mip_Phi_3H_%d_%d.png", iL+firstLayer, iR), "png");
      if(iL == 0) tcM->Print(Form("plotdir/cells/h_Mip_Phi_3H_%d_%d.root", iL+firstLayer, iR), "root");

      ++savedCout;
    }
  }



  //now fit h_Mip_Layer
  for(int ij=0; ij<14; ++ij){
    float yMax;
    getXmax(h_Mip_Layer[ij], yMax);
    hfithisto->SetParameters(yMax, 1., 0.1);
    hfithisto->SetParLimits(0, yMax/2., yMax * 1.2);
    hfithisto->SetParLimits(1, 0., 10.);
    hfithisto->SetParLimits(2, 0., 0.3);
    //std::cout << " histo per layer entries = " << h_Mip_Layer[ij]->GetEntries() << std::endl;
    //TFitResultPtr r = h_Mip_Layer[ij]->Fit("hfithisto", "RBS");
    TFitResultPtr r = h_Mip_Layer[ij]->Fit("hfithisto", "q");
    if(r != 0) continue;
    float meanV = hfithisto->GetParameter(1);
    if(meanV < 0) continue;
    h2_MIP_vs_Layer->SetBinContent(ij+1, hfithisto->GetParameter(1));
    h2_MIP_vs_Layer->SetBinError(ij+1, hfithisto->GetParError(1));
    std::cout << " >> layer = " << ij+firstLayer << " val = " << hfithisto->GetParameter(1) << std::endl;
    h_MPV_layer_values->Fill(hfithisto->GetParameter(1));

    TCanvas* tcM = new TCanvas();
    tcM->cd();
    h_Mip_Layer[ij]->GetXaxis()->SetTitle(Form("MIP layer %d", ij+firstLayer) );
    h_Mip_Layer[ij]->Draw();
    tcM->Print(Form("plotdir/cellsPhi/h_Mip_Layer_%d.png", ij+firstLayer), "png");

    //3hits
    getXmax(h_Mip_Layer_3H[ij], yMax);
    hfithisto->SetParameters(yMax, 1., 0.1);
    hfithisto->SetParLimits(0, yMax/2., yMax * 1.2);
    hfithisto->SetParLimits(1, 0., 10.);
    hfithisto->SetParLimits(2, 0., 0.3);
    //std::cout << " histo per layer entries = " << h_Mip_Layer[ij]->GetEntries() << std::endl;
    //TFitResultPtr r3H = h_Mip_Layer_3H[ij]->Fit("hfithisto", "RBS");
    TFitResultPtr r3H = h_Mip_Layer_3H[ij]->Fit("hfithisto", "q");
    if(r3H != 0) continue;
    meanV = hfithisto->GetParameter(1);
    if(meanV < 0) continue;
    h2_MIP_vs_Layer_3H->SetBinContent(ij+1, hfithisto->GetParameter(1));
    h2_MIP_vs_Layer_3H->SetBinError(ij+1, hfithisto->GetParError(1));
    std::cout << " >> layer = " << ij+firstLayer << " val = " << hfithisto->GetParameter(1) << std::endl;
    h_MPV_layer_values_3H->Fill(hfithisto->GetParameter(1));

    //    TCanvas* tcM = new TCanvas();
    tcM->cd();
    h_Mip_Layer_3H[ij]->GetXaxis()->SetTitle(Form("MIP layer %d", ij+firstLayer) );
    h_Mip_Layer_3H[ij]->Draw();
    tcM->Print(Form("plotdir/cellsPhi/h_Mip_Layer_3H_%d.png", ij+firstLayer), "png");



  }



  TF1* hfit = new TF1("hfit", "gaus", 0., 2.);
  TCanvas* tcC_Phi = new TCanvas();
  //  h_MPV_values->Rebin(2);
  hfit->SetParameters(h_MPV_values->GetEntries()/2., 1, 0.02);
  h_MPV_values->GetXaxis()->SetTitle("MPV values");
  h_MPV_values->Draw();
  //h_MPV_values->Fit("hfit", "R");
  h_MPV_values->Fit("hfit", "q");
  tcC_Phi->Print("plotdir/cellsPhi/h_MPV_values.png", "png");

  hfit->SetParameters(h_MPV_layer_values->GetEntries()/2., 1, 0.02);
  h_MPV_layer_values->GetXaxis()->SetTitle("MPV values for layers");
  h_MPV_layer_values->Draw();
  //h_MPV_layer_values->Fit("hfit", "R");
  h_MPV_layer_values->Fit("hfit", "q");
  tcC_Phi->Print("plotdir/cellsPhi/h_MPV_layer_values.png", "png");

  hfit->SetParameters(h_MPV_values_3H->GetEntries()/2., 1, 0.02);
  h_MPV_values_3H->GetXaxis()->SetTitle("MPV values (3hits)");
  h_MPV_values_3H->Draw();
  //h_MPV_values_3H->Fit("hfit", "R");
  h_MPV_values_3H->Fit("hfit", "q");
  tcC_Phi->Print("plotdir/cellsPhi/h_MPV_values_3H.png", "png");

  hfit->SetParameters(h_MPV_layer_values_3H->GetEntries()/2., 1, 0.02);
  h_MPV_layer_values_3H->GetXaxis()->SetTitle("MPV values for layers (3hits)");
  h_MPV_layer_values_3H->Draw();
  //h_MPV_layer_values_3H->Fit("hfit", "R");
  h_MPV_layer_values_3H->Fit("hfit", "q");
  tcC_Phi->Print("plotdir/cellsPhi/h_MPV_layer_values_3H.png", "png");


  gPad->SetLogy();
  tcC_Phi->cd();
  h_Count_Phi->GetXaxis()->SetTitle("n counts in rings over #phi");
  h_Count_Phi->Draw();
  tcC_Phi->Print("plotdir/cellsPhi/h_Count_Phi.png", "png");
  tcC_Phi->Print("plotdir/cellsPhi/h_Count_Phi.root", "root");

  h_minNcounts->GetXaxis()->SetTitle("minimal N. counts");
  h_minNcounts->Draw();
  tcC_Phi->Print("plotdir/cellsPhi/h_minNcounts.png", "png");

  gStyle->SetOptStat(0);
  TCanvas* tcC2_Phi = new TCanvas();
  h2_iRvsLayer_MIP->GetXaxis()->SetTitle("layer");
  h2_iRvsLayer_MIP->GetYaxis()->SetTitle("iR");
  h2_iRvsLayer_MIP->GetZaxis()->SetRangeUser(0., 3.);
  h2_iRvsLayer_MIP->Draw("colz");
  tcC2_Phi->Print("plotdir/cellsPhi/h2_iRvsLayer_MIP.png", "png");
  tcC2_Phi->Print("plotdir/cellsPhi/h2_iRvsLayer_MIP.root", "root");

  tcC2_Phi->cd();
  h2_MIP_vs_Layer->GetXaxis()->SetTitle("layer");
  h2_MIP_vs_Layer->GetYaxis()->SetTitle("MIP");
  h2_MIP_vs_Layer->GetYaxis()->SetRangeUser(0.8, 1.2);
  h2_MIP_vs_Layer->Draw();
  tcC2_Phi->Print("plotdir/cellsPhi/h2_MIP_vs_Layer.png", "png");
  tcC2_Phi->Print("plotdir/cellsPhi/h2_MIP_vs_Layer.root", "root");

  tcC2_Phi->cd();
  h2_iRvsLayer_MIP_3H->GetXaxis()->SetTitle("layer");
  h2_iRvsLayer_MIP_3H->GetYaxis()->SetTitle("iR");
  h2_iRvsLayer_MIP_3H->GetZaxis()->SetRangeUser(0., 3.);
  h2_iRvsLayer_MIP_3H->Draw("colz");
  tcC2_Phi->Print("plotdir/cellsPhi/h2_iRvsLayer_MIP_3H.png", "png");
  tcC2_Phi->Print("plotdir/cellsPhi/h2_iRvsLayer_MIP_3H.root", "root");


  tcC2_Phi->cd();
  h2_MIP_vs_Layer_3H->GetXaxis()->SetTitle("layer");
  h2_MIP_vs_Layer_3H->GetYaxis()->SetTitle("MIP");
  h2_MIP_vs_Layer_3H->GetYaxis()->SetRangeUser(0.8, 1.2);
  h2_MIP_vs_Layer_3H->Draw();
  tcC2_Phi->Print("plotdir/cellsPhi/h2_MIP_vs_Layer_3H.png", "png");
  tcC2_Phi->Print("plotdir/cellsPhi/h2_MIP_vs_Layer_3H.root", "root");



  std::cout << " all done now print results " << std::endl;


  std::cout << " all done - ciao " << std::endl;

  return 0;



}
