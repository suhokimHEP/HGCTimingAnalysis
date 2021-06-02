//g++  -o lookAtMIP_correctRate  lookAtMIP_correctRate.cpp `root-config --cflags --libs`
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
#include "TFitResult.h"
#include "TFile.h"
#include "TProfile2D.h"

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



int main(){
  gROOT->Macro("/afs/cern.ch/user/a/amartell/public/setStyle.C");

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  
  gROOT->Reset();
  gROOT->Macro("~/public/setStyle.C");
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1);

  std::cout << " inizio ci sono " << std::endl; 

  
  //  int iColors[16] = {kRed, kOrange+4, kOrange-3, kOrange-2, kBlue, kBlue-9, kAzure-9, kAzure+10, kCyan, kGreen+1, kCyan-2, kYellow+2}; //kGray+1};
  int iColors[14] = {kBlack, kGray+1, kGreen+1, kSpring-2, kYellow+1, kOrange-2, kRed+1, kPink-2, kMagenta+1, kViolet-2, kBlue+1, kAzure-2, kCyan+1, kTeal-2};
  int iStyle[4] = {20, 21, 22, 23}; //kGray+1};
  //  int iColors[6] = {kGreen+1, kBlue, kRed};
  
  //  TFile* inF = TFile::Open("../test/singleMuon_newGun.root");
  //  TFile* inF = TFile::Open("../test/testMinBias.root");

  // TFile* inF = TFile::Open("../test/MinBias140PU_0.root");

  //TFile* inF = TFile::Open("/tmp/amartell/MinBias_PU140_withNewMB_pt4Cut_vArea.root");
  //  TFile* inF = TFile::Open("/tmp/amartell/MinBias_PU140_withNewMB_pt4Cut_vGoodM.root");


  //  TFile* inF = TFile::Open("/tmp/amartell/MinBias_PU140_withNewMB_pt4Cut_vGlbMt.root");
  //  TFile* inF = TFile::Open("/tmp/amartell/MinBias_PU140_withNewMB_pt4Cut_vdetID.root");
  TFile* inF = TFile::Open("/tmp/amartell/MinBias_PU140_withNewMB_pt4Cut_vdetID_vMat.root");
  //  TFile* inF = TFile::Open("../test/singleMuon_newGun.root");

  //config
  int firstLayer = 37;
  float Nevt = 2.09206e+06 / 2.;
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

  TH1F* mipAll = new TH1F("mipAll", "", 80, -1., 19.);
  TH2F* h2_YvsX[14];
  TH2F* h2_iRvsiPhi[14];
  TH1F* h_dX[14];
  TH1F* h_dY[14];
  TH1F* h_dR[14];
  TH1F* h_dRep[14];
  TH1F* h_dRerr[14];

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

  //////// analyze rate
  TH1F* nMuons_vsEta[14];
  TH1F* nMuonsPercm2_vsEta[14];
  TH1F* ratePercm2_vsEta[14];
  TH1F* cm2_vsEta[14];
  TH1F* nMuons_vsEta_All[14];
  TH1F* nMuonsPercm2_vsEta_All[14];
  TH1F* ratePercm2_vsEta_All[14];
  TH1F* cm2_vsEta_All[14];


  //check for tracklets

  TH1F* h_MuonCross[14];
  TH1F* h_MuonNHit[14];

  TH1F* h_MuonCrossK[14][14];
  TH1F* h_MuonNHitK[14][14];

  TH1F* hDen_etaRange[14];
  TH1F* hNum_etaRange[14];

   for(int ij=0; ij<14; ++ij){
     h2_YvsX[ij] = new TH2F(Form("h2_YvsX_L%d", ij), "", 600, -300., 300., 600, -300., 300.);
     h2_iRvsiPhi[ij] = new TH2F(Form("h2_iRvsiPhi_L%d", ij), "", 100, -50., 50., 300, 0., 300.);
     h_dX[ij] = new TH1F(Form("h_dX_L%d", ij), "", 100, -10., 10.);
     h_dY[ij] = new TH1F(Form("h_dY_L%d", ij), "", 100, -10., 10.);
     h_dR[ij] = new TH1F(Form("h_dR_L%d", ij), "", 100, 0., 10.);
     h_dRep[ij] = new TH1F(Form("h_dRep_L%d", ij), "", 500, 0., 0.5);
     h_dRerr[ij] = new TH1F(Form("h_dRerr_L%d", ij), "", 100, 0., 1.);

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

     nMuons_vsEta_All[ij] = new TH1F(Form("nMuons_vsEta_All_L%d", ij+firstLayer), "", 400, 0., 4.);
     nMuonsPercm2_vsEta_All[ij] = new TH1F(Form("nMuonsPercm2_vsEta_All_L%d", ij+firstLayer), "", 400, 0., 4.);
     ratePercm2_vsEta_All[ij] = new TH1F(Form("ratePercm2_vsEta_All_L%d", ij+firstLayer), "", 400, 0., 4.);
  
     h_MuonCross[ij] = new TH1F(Form("h_MuonCross_L%d", ij+firstLayer), "", 300, 1., 4.);
     h_MuonNHit[ij] = new TH1F(Form("h_MuonNHit_L%d_N%d", ij+firstLayer, 1), "", 300, 1., 4.);

     hDen_etaRange[ij] = new TH1F(Form("hDen_etaRange_L%d", ij+firstLayer), "", 1, 1.55, 1.71);
     hNum_etaRange[ij] = new TH1F(Form("hNum_etaRange_L%d", ij+firstLayer), "", 1, 1.55, 1.71);

     for(int kl=0; kl<14; ++kl){
       h_MuonNHitK[kl][ij] = new TH1F(Form("h_MuonNHitK_L%d_N%d", ij+firstLayer, kl), "", 300, 1., 4.);
       h_MuonCrossK[kl][ij] = new TH1F(Form("h_MuonCrossK_L%d_N%d", ij+firstLayer, kl), "", 300, 1., 4.);
     }
   }

   bool matchByID = true;
   bool matchBydR = false;



   int nEvents = newT->GetEntries();
   std::cout << " nEvents = " << nEvents << std::endl;
   for(int ij=0; ij<nEvents; ++ij){
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
     std::map<int, std::vector<float>> trkPt;
     std::map<int, std::vector<float>> trkP;
     std::map<int, std::vector<float>> trkChi;
     std::map<int, std::vector<short int>> trkQ;

     //    std::cout << " muonP->size() = " << muonP->size() << std::endl;
     for(int iM=0; iM<muonP->size(); ++iM){
       int trkL = crossL->at(iM) - firstLayer;
       trkN[trkL].push_back(int(crossM->at(iM)));
       trkX[trkL].push_back(crossX->at(iM));
       trkY[trkL].push_back(crossY->at(iM));
       trkXerr[trkL].push_back(crossEX->at(iM));
       trkYerr[trkL].push_back(crossEY->at(iM));
       trkEta[trkL].push_back(crossEta->at(iM));
       trkPhi[trkL].push_back(crossPhi->at(iM));
       trkCEta[trkL].push_back(crossCellEta->at(iM));
       trkCPhi[trkL].push_back(crossCellPhi->at(iM));
       trkiR[trkL].push_back(crossiR->at(iM));
       trkiP[trkL].push_back(crossiP->at(iM));
       trkPt[trkL].push_back(muonPt->at(iM));
       trkP[trkL].push_back(muonP->at(iM));
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
     std::map<int, int> startLMatchedTrack;

     for(int iL=0; iL < 14; ++iL){
       int iTc = -1;

       bool foundMatchInLayer = false;

       for(auto iT : trkX[iL]){
	 ++iTc;
	 if(trkChi[iL][iTc] > 1.5) continue;

	 int trackIdx = trkN[iL][iTc];
	 //std::cout << " trackIdx = " << trackIdx << std::endl;
	 if(foundFirstHit.find(trackIdx) == foundFirstHit.end()) foundFirstHit[trackIdx] = false;

	 float trackEta = trkEta[iL][iTc];
	 float trackPhi = trkPhi[iL][iTc];
	 float trackY = trkY[iL][iTc];
	 float trkR = sqrt(iT*iT + trackY*trackY);
	 float trkEtaOnLayer = asinh(zVal[iL]/trkR);                            
	 nMuons_vsEta_All[iL]->Fill(std::abs(trkEtaOnLayer));  

	 // float dPhi = std::abs(trkCPhi[iL][iTc] - trackPhi);
	 // if(dPhi > pigreco) dPhi -= 2* pigreco;
	 // if(std::abs(trkCEta[iL][iTc] - trackEta) > 0.009 || std::abs(dPhi) > 0.009 ){
	 //   continue;
	 // }
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
	   float dXerr = trkXerr[iL][iTc];
	   float dYerr = trkYerr[iL][iTc];
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

	   h2_dR_vsChi[iL]->Fill(dR, trkChi[iL][iTc]);
	   h2_dRep_vsChi[iL]->Fill(dRep, trkChi[iL][iTc]);
	   h2_dR_vsTrkQ[iL]->Fill(dR, trkQ[iL][iTc]);

	   std::pair<int, int> channelPhi = std::pair<int, int>(iL, std::abs(reciR));

	   if(trkChi[iL][iTc] > 1.5){
	     continue;
	   }

	   if(matchByID){
	     if(rhiR[iL][iRc] != trkiR[iL][iTc] || rhiPhi[iL][iRc] != trkiP[iL][iTc]) continue;

	     // float dPhi = std::abs(trkPhi[iL][iTc] - rhPhi[iL][iRc]);
	     // if(dPhi > pigreco) dPhi -= 2* pigreco;
	     // if(std::abs(trkEta[iL][iTc] - rhEta[iL][iRc]) > 0.009 || std::abs(dPhi) > 0.009 ){
	     //   continue;
	     // }

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
       }//track
     }// layer
  
   
     //     continue;
     for(int nHits = 1; nHits<14; ++nHits){
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

	 //std::cout << " nHits = " << nHits << " matchedTrack[trackIdx] = " << matchedTrack[trackIdx] << " startLMatchedTrack[trackIdx] = " << startLMatchedTrack[trackIdx] << " iL = " << iL << std::endl;

	 if( (nHits > matchedTrack[trackIdx] &&  (matchedTrack[trackIdx] != (14 - startLMatchedTrack[trackIdx]) || iL < startLMatchedTrack[trackIdx])) || 
	     matchedTrack[trackIdx] == 0){
	   //	   if(ij == 148) std::cout << " missing hits " << std::endl; 
	   continue;
	 }
	 //	 std::cout << " ok " << std::endl;


	 if(trkChi[iL][iTc] > 1.5) continue;
	 // float dPhi = std::abs(trkCPhi[iL][iTc] - trkPhi[iL][iTc]);
         // if(dPhi > pigreco) dPhi -= 2* pigreco;
         // if(std::abs(trkCEta[iL][iTc] - trkEta[iL][iTc]) > 0.009 || std::abs(dPhi) > 0.009 ){
         //   continue;
         // }
         if(foundFirstHit.find(trackIdx) == foundFirstHit.end()) foundFirstHit[trackIdx] = false;
	 //	 if(ij == 148) std::cout << " foundFirstHit.find(trackIdx) = " << foundFirstHit[trackIdx]  << std::endl;

	 float trackY = trkY[iL][iTc];
         float trkR = sqrt(iT*iT + trackY*trackY);

         float trkEtaOnLayer = asinh(zVal[iL]/trkR);

	 h_MuonCrossK[nHits][iL]->Fill(trkEtaOnLayer); 
	 //std::cout << " fill den K iL = " << iL << " trackIdx = " << trackIdx << std::endl;

	 //if(iL < startLMatchedTrack[trackIdx]) continue;

         bool alreadyFilled = false;
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
	   
	     // float dPhi = std::abs(trkPhi[iL][iTc] - rhPhi[iL][iRc]);
	     // if(dPhi > pigreco) dPhi -= 2* pigreco;
	     // if(std::abs(trkEta[iL][iTc] - rhEta[iL][iRc]) > 0.009 || std::abs(dPhi) > 0.009 ){
	     //   continue;
	     // }
	     h_MuonNHitK[nHits][iL]->Fill(trkEtaOnLayer);
	     //std::cout << " fill num K iL = " << iL << " trackIdx = " << trackIdx << std::endl;
	     
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
	   }
	   else if(matchBydR){
	     float hitTrack_dx = iR - iT;
	     float hitTrack_dy = recY - trackY;
	     float hitTrack_dR = sqrt(hitTrack_dx*hitTrack_dx + hitTrack_dy*hitTrack_dy);
	     float dXerr = trkXerr[iL][iTc];
	     float dYerr = trkYerr[iL][iTc];	     
	     float cellSize = sqrt(rhA[iL][iRc]) / 2.;

	     if(std::abs(hitTrack_dx) > cellSize+dXerr || std::abs(hitTrack_dy) > cellSize+dYerr || hitTrack_dR > nearestdR) continue;
	     nearestHit = iRc;
	     nearestdR = hitTrack_dR;
	   }
	 }//loop over recHits
	 if(matchBydR && nearestHit != -1){
	   h_MuonNHitK[nHits][iL]->Fill(trkEtaOnLayer);
	 }
       }// trk loop
     }// layer loop

     //     if(ij == 148) std::cout << " looping nHits = " << nHits <<std::endl;
     }// loop on nhits
   }// events;

   //   return 5;
   std::cout << " okChannels.size = " << okChannels.size() << " okChannelPhi.size() = " << okChannelsPhi.size() << std::endl;
   //   return 200;


   /*
   TCanvas* tcM = new TCanvas();
   tcM->cd();
   mipAll->Draw();
   tcM->Print("plotsMinBias_correctRate/cellsPhi/mipAll.png", "png");
   tcM->Print("plotsMinBias_correctRate/cellsPhi/mipAll.root", "root");

   dX_wrt1st->GetXaxis()->SetTitle("dX hits wrt 1st layer");
   dX_wrt1st->Draw();
   tcM->Print("plotsMinBias_correctRate/x1stHit.png", "png");

   dY_wrt1st->GetXaxis()->SetTitle("dY hits wrt 1st layer");
   dY_wrt1st->Draw();
   tcM->Print("plotsMinBias_correctRate/y1stHit.png", "png");

   dX_wrtPrev->GetXaxis()->SetTitle("dX hits wrt previous layer");
   dX_wrtPrev->Draw();
   tcM->Print("plotsMinBias_correctRate/xPrevHit.png", "png");

   dY_wrtPrev->GetXaxis()->SetTitle("dY hits wrt previous layer");
   dY_wrtPrev->Draw();
   tcM->Print("plotsMinBias_correctRate/yPrevHit.png", "png");
   */
   
   //   return;

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
     }
     
     TCanvas* tmuo = new TCanvas();
     tmuo->cd();
     nMuonsPercm2_vsEta[ij]->GetXaxis()->SetTitle(Form("#eta layer %d", ij+firstLayer));
     nMuonsPercm2_vsEta[ij]->GetYaxis()->SetTitle("n Muons / cm2");
     nMuonsPercm2_vsEta[ij]->Draw();
     tmuo->Print(Form("plotsMuon_MinBias_correctRate/nMuonsPercm2_vsEta_L%d.png", ij+firstLayer), "png");
     //    tmuo->Print(Form("plotsMuon_MinBias_correctRate/nMuonsPercm2_vsEta_L%d.root", ij+firstLayer), "root");
     
     nMuons_vsEta[ij]->GetXaxis()->SetTitle(Form("#eta layer %d", ij+firstLayer));
     nMuons_vsEta[ij]->Draw();
     tmuo->Print(Form("plotsMuon_MinBias_correctRate/nMuons_vsEta_L%d.png", ij+firstLayer), "png");

     ratePercm2_vsEta[ij]->GetYaxis()->SetRangeUser(0., 5.);
     ratePercm2_vsEta[ij]->GetXaxis()->SetTitle(Form("#eta layer %d", ij+firstLayer));
     ratePercm2_vsEta[ij]->GetYaxis()->SetTitle("n Muons / cm2 / s");
     ratePercm2_vsEta[ij]->Draw();
     tmuo->Print(Form("plotsMuon_MinBias_correctRate/ratePercm2_vsEta_L%d.png", ij+firstLayer), "png");  

     gPad->SetLogy();
     cm2_vsEta[ij]->GetYaxis()->SetRangeUser(0.5, 5.e4);
     cm2_vsEta[ij]->GetXaxis()->SetTitle(Form("#eta layer %d", ij+firstLayer));
     cm2_vsEta[ij]->GetYaxis()->SetTitle("cm2");
     cm2_vsEta[ij]->Draw();
     tmuo->Print(Form("plotsMuon_MinBias_correctRate/cm2_vsEta_L%d.png", ij+firstLayer), "png");
     
     TCanvas* tmuoAll = new TCanvas();
     tmuoAll->cd();
     nMuonsPercm2_vsEta_All[ij]->GetXaxis()->SetTitle(Form("#eta layer %d", ij+firstLayer));
     nMuonsPercm2_vsEta_All[ij]->GetYaxis()->SetTitle("n Muons / cm2");
     nMuonsPercm2_vsEta_All[ij]->Draw();
     tmuoAll->Print(Form("plotsMuon_MinBias_correctRate/nMuonsPercm2_vsEta_All_L%d.png", ij+firstLayer), "png");
     //    tmuo->Print(Form("plotsMuon_MinBias_correctRate/nMuonsPercm2_vsEta_L%d.root", ij+firstLayer), "root");
     
     nMuons_vsEta_All[ij]->GetXaxis()->SetTitle(Form("#eta layer %d", ij+firstLayer));
     nMuons_vsEta_All[ij]->Draw();
     tmuoAll->Print(Form("plotsMuon_MinBias_correctRate/nMuons_vsEta_All_L%d.png", ij+firstLayer), "png");

     ratePercm2_vsEta_All[ij]->GetYaxis()->SetRangeUser(0., 20.);
     ratePercm2_vsEta_All[ij]->GetXaxis()->SetTitle(Form("#eta layer %d", ij+firstLayer));
     ratePercm2_vsEta_All[ij]->GetYaxis()->SetTitle("n Muons / cm2 / s");
     ratePercm2_vsEta_All[ij]->Draw();
     tmuoAll->Print(Form("plotsMuon_MinBias_correctRate/ratePercm2_vsEta_All_L%d.png", ij+firstLayer), "png");  
   }
   
   
   //   return;

   TCanvas* tChi = new TCanvas();
   tChi->cd();
   h_TrkChi_all->GetXaxis()->SetTitle("trk chi all"); 
   h_TrkChi_all->Draw();
   tChi->Print("plotsMuon_MinBias_correctRate/h_TrkChi_all.png", "png");
   tChi->Print("plotsMuon_MinBias_correctRate/h_TrkChi_all.root", "root");

   h_TrkQ_all->GetXaxis()->SetTitle("trk Q all"); 
   h_TrkQ_all->Draw();
   tChi->Print("plotsMuon_MinBias_correctRate/h_TrkQ_all.png", "png");
   tChi->Print("plotsMuon_MinBias_correctRate/h_TrkQ_all.root", "root");

   h_TrkChi_passed->GetXaxis()->SetTitle("trk chi all"); 
   h_TrkChi_passed->Draw();
   tChi->Print("plotsMuon_MinBias_correctRate/h_TrkChi_passed.png", "png");
   tChi->Print("plotsMuon_MinBias_correctRate/h_TrkChi_passed.root", "root");

   h_TrkQ_passed->GetXaxis()->SetTitle("trk Q all"); 
   h_TrkQ_passed->Draw();
   tChi->Print("plotsMuon_MinBias_correctRate/h_TrkQ_passed.png", "png");
   tChi->Print("plotsMuon_MinBias_correctRate/h_TrkQ_passed.root", "root");


   //   return;
   //   if(1 == 2){
   TCanvas* tMuo = new TCanvas();
   tMuo->cd();
   h_etaMuon->GetXaxis()->SetTitle("muon #eta");
   h_etaMuon->Draw();
   tMuo->Print("plotsMuon_MinBias_correctRate/etaMuon_all_L37.png", "png");
   tMuo->Print("plotsMuon_MinBias_correctRate/etaMuon_all_L37.root", "root");

   gPad->SetLogy();
   h_ptMuon->GetXaxis()->SetTitle("muon pT");
   h_ptMuon->Draw();
   tMuo->Print("plotsMuon_MinBias_correctRate/ptMuon_all_L37.png", "png");
   tMuo->Print("plotsMuon_MinBias_correctRate/ptMuon_all_L37.root", "root");
   //   }

   //   if(1 == 2){
   //return;

   //about efficiencies

   TGraph* etaRange_eff = new TGraph();

   //   TH1F* Muon1stHit_eff[14]; 
   TH1F* MuonNHit_eff[14];
   TH1F* MuonNHitK_eff[14][14];

   TGraph* MuonNHit_eff_tg[14];
   TGraph* MuonNHitK_eff_tg[14][14];

   for(int ij=0; ij<14; ++ij){
     TCanvas* tc = new TCanvas();
     tc->cd();

     h_MuonCross[ij]->Rebin(4);
     h_MuonNHit[ij]->Rebin(4);

     h_MuonCross[ij]->GetYaxis()->SetTitle("muon cross");
     h_MuonCross[ij]->GetXaxis()->SetTitle(Form("muon eta layer %d", ij+firstLayer));
     h_MuonCross[ij]->Draw();
     tc->Print(Form("plotsMinBias_correctRate/h_MuonCross_L%d.png", ij+firstLayer), "png");

     h_MuonNHit[ij]->GetYaxis()->SetTitle("muon matched hit");
     h_MuonNHit[ij]->GetXaxis()->SetTitle(Form("muon eta layer %d", ij+firstLayer));
     h_MuonNHit[ij]->Draw();
     tc->Print(Form("plotsMinBias_correctRate/h_MuonNHit_L%d.png", ij+firstLayer), "png");

     MuonNHit_eff[ij] = (TH1F*)h_MuonNHit[ij]->Clone(Form("MuonNHit_eff_%d", ij+firstLayer));
     MuonNHit_eff[ij]->Reset();
     MuonNHit_eff[ij]->Divide(h_MuonNHit[ij],h_MuonCross[ij]);
     MuonNHit_eff[ij]->GetXaxis()->SetTitle(Form("muon eta layer %d", ij+firstLayer));
     MuonNHit_eff[ij]->GetYaxis()->SetTitle("purity");
     MuonNHit_eff[ij]->Draw();
     tc->Print(Form("plotsMinBias_correctRate/MuonNHit_eff_L%d.png", ij+firstLayer), "png");
     
     MuonNHit_eff_tg[ij] = new TGraph();
     MuonNHit_eff_tg[ij]->SetTitle(Form("MuonNHit_eff_tg_L%d_nH%d",ij));
     MuonNHit_eff_tg[ij]->SetPoint(0, 0, 0);


     // h_MuonNHitK[5][ij]->GetYaxis()->SetTitle("muon N hits");
     // h_MuonNHitK[5][ij]->GetXaxis()->SetTitle(Form("muon eta layer %d", ij+firstLayer));
     // h_MuonNHitK[5][ij]->Draw();
     // tc->Print(Form("plotsMinBias_correctRate/h_MuonNHitK_L%d_N5.png", ij+firstLayer), "png");

     // h_MuonNHitK[13][ij]->GetYaxis()->SetTitle("muon N hits");
     // h_MuonNHitK[13][ij]->GetXaxis()->SetTitle(Form("muon eta layer %d", ij+firstLayer));
     // h_MuonNHitK[13][ij]->Draw();
     // tc->Print(Form("plotsMinBias_correctRate/h_MuonNHitK_L%d_N13.png", ij+firstLayer), "png");


     // Muon1stHit_eff[ij] = (TH1F*)h_Muon1stHit[ij]->Clone(Form("Muon1stHit_eff_%d", ij+firstLayer));
     // Muon1stHit_eff[ij]->Reset();
     // Muon1stHit_eff[ij]->Divide(h_Muon1stHit[ij], h_MuonCross[ij]);
     // Muon1stHit_eff[ij]->GetXaxis()->SetTitle(Form("muon eta layer %d", ij+firstLayer));
     // Muon1stHit_eff[ij]->GetYaxis()->SetTitle("purity (finding 1 hit)");
     // Muon1stHit_eff[ij]->Draw();
     // tc->Print(Form("plotsMinBias_correctRate/Muon1stHit_eff_L%d.png", ij+firstLayer), "png");

     for(int kl=0; kl<14; ++kl){
       h_MuonNHitK[kl][ij]->Rebin(4);
       h_MuonCrossK[kl][ij]->Rebin(4);
       //       if(ij == 4){std::cout << " num entries = " << h_MuonNHitK[kl][ij]->GetEntries() << " den entries = " << h_MuonCrossK[kl][ij]->GetEntries() << std::endl;}
       MuonNHitK_eff[kl][ij] = (TH1F*)h_MuonNHitK[kl][ij]->Clone(Form("MuonNHitK_eff_%d_N%d", ij+firstLayer, kl));
       MuonNHitK_eff[kl][ij]->Reset();
       MuonNHitK_eff[kl][ij]->Divide(h_MuonNHitK[kl][ij],h_MuonCrossK[kl][ij]);
       MuonNHitK_eff[kl][ij]->GetXaxis()->SetTitle(Form("muon eta layer %d", ij+firstLayer));
       MuonNHitK_eff[kl][ij]->GetYaxis()->SetTitle(Form("purity (finding %d hits)", kl));
       MuonNHitK_eff[kl][ij]->Draw();
       //       if(kl == 5 || kl == 13) tc->Print(Form("plotsMinBias_correctRate/MuonNHitK_eff_L%d_N%d.png", ij+firstLayer, kl), "png");

       MuonNHitK_eff_tg[kl][ij] = new TGraph();
       MuonNHitK_eff_tg[kl][ij]->SetTitle(Form("MuonNHitK_eff_tg_L%d_nH%d",ij, kl));
       MuonNHitK_eff_tg[kl][ij]->SetPoint(0, 0, 0);
       
     }
     if(ij == 0) etaRange_eff->SetPoint(0, 0, 0);
     etaRange_eff->SetPoint(etaRange_eff->GetN()+1, ij+firstLayer, hNum_etaRange[ij]->GetBinContent(1)/hDen_etaRange[ij]->GetBinContent(1));
     hNum_etaRange[ij]->Draw();
     tc->Print(Form("plotsMinBias_correctRate/hNum_etaRange_L%d.png", ij+firstLayer), "png");
     hDen_etaRange[ij]->Draw();
     tc->Print(Form("plotsMinBias_correctRate/hDen_etaRange_L%d.png", ij+firstLayer), "png");
   }
   etaRange_eff->SetPoint(etaRange_eff->GetN()+1, 60, 1.5);
   TCanvas* tgEta = new TCanvas();
   tgEta->cd();
   etaRange_eff->GetXaxis()->SetTitle("layer");
   etaRange_eff->GetXaxis()->SetRangeUser(30., 55);
   etaRange_eff->Draw("ap");
   tgEta->Print("plotsMinBias_correctRate/tgEta.png", "png");


   //tg
   for(int iB=1; iB<MuonNHit_eff[0]->GetNbinsX(); ++iB){
     int maxLayer = 0;
     float effMax = 0;
     float avgEta = MuonNHit_eff[0]->GetBinCenter(iB);
     for(int ij=0; ij<14; ++ij){
       float dummyEff = MuonNHit_eff[ij]->GetBinContent(iB);
       if(dummyEff > effMax){
	 effMax = dummyEff;
	 maxLayer = ij;
       }
       // if(ij == 0 && (kl == 13 || kl == 1)){
       // 	 std::cout << " layer = " << ij << " nH = " << kl << " avgEta = " << avgEta << " maxLayer = " << maxLayer 
       // 		   << " effMax = " << effMax << " dummyEff = " << dummyEff << std::endl;
       // }
     }
     MuonNHit_eff_tg[maxLayer]->SetPoint(MuonNHit_eff_tg[maxLayer]->GetN()+1, avgEta, effMax);
     //       if((kl == 13 || kl == 1))std::cout << " avgEta = " << avgEta << " effMax = " << effMax << " iL = " << maxLayer << std::endl;
   }
   for(int ij=0; ij<14; ++ij)
     MuonNHit_eff_tg[ij]->SetPoint(MuonNHit_eff_tg[ij]->GetN()+1, 5., 1.5);

   //now k tg
   for(int kl=0; kl<14; ++kl){
     for(int iB=1; iB<MuonNHitK_eff[kl][0]->GetNbinsX(); ++iB){
       int maxLayer = 0;
       float effMax = 0;
       float avgEta = MuonNHitK_eff[kl][0]->GetBinCenter(iB);
       for(int ij=0; ij<14; ++ij){
	 float dummyEff = MuonNHitK_eff[kl][ij]->GetBinContent(iB);
	 if(dummyEff > effMax){
	   effMax = dummyEff;
	   maxLayer = ij;
	 }
	 // if(ij == 0 && (kl == 13 || kl == 1)){
	 //   std::cout << " layer = " << ij << " nH = " << kl << " avgEta = " << avgEta << " maxLayer = " << maxLayer 
	 // 	     << " effMax = " << effMax << " dummyEff = " << dummyEff << std::endl;
	 // }
       }
       MuonNHitK_eff_tg[kl][maxLayer]->SetPoint(MuonNHitK_eff_tg[kl][maxLayer]->GetN()+1, avgEta, effMax);
       //       if((kl == 13 || kl == 1))std::cout << " avgEta = " << avgEta << " effMax = " << effMax << " iL = " << maxLayer << std::endl;
     }
     for(int ij=0; ij<14; ++ij)
       MuonNHitK_eff_tg[kl][ij]->SetPoint(MuonNHitK_eff_tg[kl][ij]->GetN()+1, 5., 1.5);
   }
   std::cout << " fine ciao "<< std::endl;




   gStyle->SetOptStat(0);
   gStyle->SetOptStat(0);
   TLegend *legTGM = new TLegend(0.7,0.30,1.,1.,NULL,"brNDC");
   legTGM->SetTextFont(42);
   legTGM->SetTextSize(0.05);
   legTGM->SetFillColor(kWhite);
   legTGM->SetLineColor(kWhite);
   legTGM->SetShadowColor(kWhite);

   TCanvas* tcE = new TCanvas();
   tcE->cd();
   // Muon1stHit_eff[0]->GetXaxis()->SetTitle("muon eta");
   // Muon1stHit_eff[0]->SetMarkerColor(iColors[0]);
   // Muon1stHit_eff[0]->SetLineColor(kWhite);
   // Muon1stHit_eff[0]->SetMarkerStyle(20);
   // Muon1stHit_eff[0]->Draw("p");
   // for(int ij=0; ij<14; ++ij){
   //   Muon1stHit_eff[ij]->SetMarkerColor(iColors[ij]);
   //   Muon1stHit_eff[ij]->SetMarkerStyle(20);
   //   Muon1stHit_eff[ij]->SetLineColor(kWhite);
   //   Muon1stHit_eff[ij]->Draw("same, p");
   //   legTGM->AddEntry(Muon1stHit_eff[ij], Form("layer %d", ij+firstLayer), "p");
   // }
   // gStyle->SetOptStat(0);
   // gStyle->SetOptTitle(0);
   // Muon1stHit_eff[0]->GetXaxis()->SetRangeUser(1.2, 3.);
   // legTGM->Draw("same");
   // tcE->Print("plotsMinBias_correctRate/Muon1stHit_eff.png", "png");
   // tcE->Print("plotsMinBias_correctRate/Muon1stHit_eff.root", "root");
   

   for(int kl=0; kl<14; ++kl){
     TLegend *legTGM2 = new TLegend(0.7,0.30,1.,1.,NULL,"brNDC");
     legTGM2->SetTextFont(42);
     legTGM2->SetTextSize(0.03);
     legTGM2->SetFillColor(kWhite);
     legTGM2->SetLineColor(kWhite);
     legTGM2->SetShadowColor(kWhite);

     TCanvas* tcNE = new TCanvas();
     MuonNHitK_eff[kl][0]->GetXaxis()->SetTitle("muon eta");
     MuonNHitK_eff[kl][0]->SetMarkerColor(iColors[0]);
     MuonNHitK_eff[kl][0]->SetLineColor(kWhite);
     MuonNHitK_eff[kl][0]->SetMarkerStyle(20);
     MuonNHitK_eff[kl][0]->Draw("p");
     for(int ij=0; ij<14; ++ij){
       MuonNHitK_eff[kl][ij]->SetMarkerColor(iColors[ij]);
       MuonNHitK_eff[kl][ij]->SetMarkerStyle(20);
       MuonNHitK_eff[kl][ij]->SetLineColor(kWhite);
       MuonNHitK_eff[kl][ij]->Draw("same, p");
       //       legTGM2->AddEntry(MuonNHitK_eff[kl][ij], Form("l. %d (nH. >= %d)", ij+firstLayer, kl), "p");
       if(kl < (14 - ij)) legTGM2->AddEntry(MuonNHitK_eff[kl][ij], Form("l. %d (nH. = %d)", ij+firstLayer, kl), "p");
       else legTGM2->AddEntry(MuonNHitK_eff[kl][ij], Form("l. %d (nH. = %d)", ij+firstLayer, 14 - ij), "p");
     }
     gStyle->SetOptStat(0);
     MuonNHitK_eff[kl][0]->GetXaxis()->SetRangeUser(1.2, 3.);
     MuonNHitK_eff[kl][0]->GetYaxis()->SetRangeUser(0., 1.05);
     legTGM2->Draw("same");
     tcNE->Print(Form("plotsMinBias_correctRate/MuonNHitK_eff_N%02d.png", kl), "png");
     tcNE->Print(Form("plotsMinBias_correctRate/MuonNHitK_eff_N%02d.root", kl), "root");

     if(kl == 1){
       legTGM2->Clear();
     MuonNHit_eff[0]->GetXaxis()->SetTitle("muon eta");
     MuonNHit_eff[0]->SetMarkerColor(iColors[0]);
     MuonNHit_eff[0]->SetLineColor(kWhite);
     MuonNHit_eff[0]->SetMarkerStyle(20);
     MuonNHit_eff[0]->Draw("p");
     for(int ij=0; ij<14; ++ij){
       MuonNHit_eff[ij]->SetMarkerColor(iColors[ij]);
       MuonNHit_eff[ij]->SetMarkerStyle(20);
       MuonNHit_eff[ij]->SetLineColor(kWhite);
       MuonNHit_eff[ij]->Draw("same, p");
       legTGM2->AddEntry(MuonNHit_eff[ij], Form("layer %d ", ij+firstLayer), "p");
     }
     gStyle->SetOptStat(0);
     MuonNHit_eff[0]->GetXaxis()->SetRangeUser(1.2, 3.);
     MuonNHit_eff[0]->GetYaxis()->SetRangeUser(0., 1.05);
     legTGM2->Draw("same");
     tcNE->Print(Form("plotsMinBias_correctRate/MuonNHit_eff.png", kl), "png");
     tcNE->Print(Form("plotsMinBias_correctRate/MuonNHit_eff.root", kl), "root");
     }
   }

   //now plotting tg

   TLegend *legTGM22 = new TLegend(0.7,0.30,1.,1.,NULL,"brNDC");
   legTGM22->Clear();
   legTGM22->SetTextFont(42);
   legTGM22->SetTextSize(0.03);
   legTGM22->SetFillColor(kWhite);
   legTGM22->SetLineColor(kWhite);
   legTGM22->SetShadowColor(kWhite);

   //TGraph
   TCanvas* tcNEtg = new TCanvas();
   tcNEtg->cd();
   MuonNHit_eff_tg[0]->GetXaxis()->SetTitle("muon eta");
   MuonNHit_eff_tg[0]->SetMarkerColor(iColors[0]);
   MuonNHit_eff_tg[0]->SetMarkerStyle(20);
   MuonNHit_eff_tg[0]->Draw("ap");
   for(int kl=0; kl<14; ++kl){
     MuonNHit_eff_tg[kl]->SetMarkerColor(iColors[kl]);
     MuonNHit_eff_tg[kl]->SetMarkerStyle(20);
     MuonNHit_eff_tg[kl]->SetLineColor(kWhite);
     MuonNHit_eff_tg[kl]->Draw("same, p");
     legTGM22->AddEntry(MuonNHit_eff_tg[kl], Form("layer %d", kl+firstLayer), "p");
   }
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   MuonNHit_eff_tg[0]->GetXaxis()->SetRangeUser(1.2, 3.);
   MuonNHit_eff_tg[0]->GetYaxis()->SetRangeUser(0., 1.05);
   legTGM22->Draw("same");
   tcNEtg->Print("plotsMinBias_correctRate/MuonNHit_eff_tg.png", "png");
   tcNEtg->Print("plotsMinBias_correctRate/MuonNHit_eff_tg.root", "root");


   //now k tg
   for(int kl=1; kl<14; ++kl){
     legTGM22->Clear();
     tcNEtg->cd();
     MuonNHitK_eff_tg[kl][0]->GetXaxis()->SetTitle("muon eta");
     MuonNHitK_eff_tg[kl][0]->SetMarkerColor(iColors[0]);
     MuonNHitK_eff_tg[kl][0]->SetMarkerStyle(20);
     MuonNHitK_eff_tg[kl][0]->Draw("ap");
     for(int iL=0; iL<14; ++iL){
     MuonNHitK_eff_tg[kl][iL]->SetMarkerColor(iColors[iL]);
     MuonNHitK_eff_tg[kl][iL]->SetMarkerStyle(20);
     MuonNHitK_eff_tg[kl][iL]->SetLineColor(kWhite);
     MuonNHitK_eff_tg[kl][iL]->Draw("same, p");
     if(kl < (14 - iL)) legTGM22->AddEntry(MuonNHitK_eff_tg[kl][iL], Form("l. %d (nH. = %d)", iL+firstLayer, kl), "p");
     else legTGM22->AddEntry(MuonNHitK_eff_tg[kl][iL], Form("l. %d (nH. = %d)", iL+firstLayer, 14 - iL), "p");
     }
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   MuonNHitK_eff_tg[kl][0]->GetXaxis()->SetRangeUser(1.2, 3.);
   MuonNHitK_eff_tg[kl][0]->GetYaxis()->SetRangeUser(0., 1.05);
   legTGM22->Draw("same");
   // tcNEtg->Print(Form("plotsMinBias_correctRate/MuonNHit_eff_tg_L%02d.png", kl+firstLayer), "png");
   // tcNEtg->Print(Form("plotsMinBias_correctRate/MuonNHit_eff_tg_L%02d.root", kl+firstLayer), "root");
   tcNEtg->Print(Form("plotsMinBias_correctRate/MuonNHit_eff_tg_N%02d.png", kl), "png");
   tcNEtg->Print(Form("plotsMinBias_correctRate/MuonNHit_eff_tg_N%02d.root", kl), "root");
   }
   
   //   return 20;


   gStyle->SetOptStat(1);
   for(int ij=0; ij<14; ++ij){
     TCanvas* tc = new TCanvas();
     tc->cd();
     h2_YvsX[ij]->GetXaxis()->SetTitle(Form("layer %d", ij));
     h2_YvsX[ij]->Draw("colz");
     tc->Print(Form("plotsMinBias_correctRate/h2_YvsX_L%d.png", ij), "png");

     h2_iRvsiPhi[ij]->GetXaxis()->SetTitle(Form("|iR| layer %d", ij));
     h2_iRvsiPhi[ij]->GetYaxis()->SetTitle("iPhi");
     h2_iRvsiPhi[ij]->Draw("colz");
     tc->Print(Form("plotsMinBias_correctRate/h2_iRvsiPhi_L%d.png", ij), "png");
     //   continue;
     h_dX[ij]->GetXaxis()->SetTitle(Form("dX layer %d", ij));
     h_dX[ij]->Draw();
     tc->Print(Form("plotsMinBias_correctRate/h_dX_L%d.png", ij), "png");

     h_dY[ij]->GetXaxis()->SetTitle(Form("dY layer %d", ij));
     h_dY[ij]->Draw();
     tc->Print(Form("plotsMinBias_correctRate/h_dY_L%d.png", ij), "png");

     h_dR[ij]->GetXaxis()->SetTitle(Form("dR layer %d", ij+firstLayer));
     h_dR[ij]->Draw();
     tc->Print(Form("plotsMinBias_correctRate/h_dR_L%d.png", ij+firstLayer), "png");
     tc->Print(Form("plotsMinBias_correctRate/h_dR_L%d.root", ij+firstLayer), "root");

     h_dRep[ij]->GetXaxis()->SetTitle(Form("dR(#eta, #phi) layer %d", ij+firstLayer));
     h_dRep[ij]->Draw();
     tc->Print(Form("plotsMinBias_correctRate/h_dRep_L%d.png", ij+firstLayer), "png");
     tc->Print(Form("plotsMinBias_correctRate/h_dRep_L%d.root", ij+firstLayer), "root");
     h_dRep[ij]->GetXaxis()->SetRangeUser(0., 0.05);
     tc->Print(Form("plotsMinBias_correctRate/h_dRep_L%d_zoomIn.png", ij+firstLayer), "png");

     h_dRerr[ij]->GetXaxis()->SetTitle(Form("dRerror(x, y) layer %d", ij+firstLayer));
     h_dRerr[ij]->Draw();
     tc->Print(Form("plotsMinBias_correctRate/h_dRerr_L%d.png", ij+firstLayer), "png");
     tc->Print(Form("plotsMinBias_correctRate/h_dRerr_L%d.root", ij+firstLayer), "root");
     //     h_dRerr[ij]->GetXaxis()->SetRangeUser(0., 0.05);
     tc->Print(Form("plotsMinBias_correctRate/h_dRerr_L%d.png", ij+firstLayer), "png");
     tc->Print(Form("plotsMinBias_correctRate/h_dRerr_L%d.root", ij+firstLayer), "root");

     // h_dRvsP[ij]->GetXaxis()->SetTitle(Form("dR layer %d", ij));
     // h_dRvsP[ij]->GetYaxis()->SetTitle("muon P");
     // h_dRvsP[ij]->Draw("colz");
     // tc->Print(Form("plotsMinBias_correctRate/h_dRvsdP_L%d.png", ij), "png");

     h2_dR_vsChi[ij]->GetXaxis()->SetTitle(Form("dR layer %d", ij));
     h2_dR_vsChi[ij]->GetYaxis()->SetTitle("muon trk chi2");
     h2_dR_vsChi[ij]->Draw("colz");
     tc->Print(Form("plotsMinBias_correctRate/h2_dR_vsChi_L%d.png", ij), "png");

     h2_dRep_vsChi[ij]->GetXaxis()->SetTitle(Form("dR(#eta, #phi) layer %d", ij));
     h2_dRep_vsChi[ij]->GetYaxis()->SetTitle("muon trk chi2");
     h2_dRep_vsChi[ij]->Draw("colz");
     tc->Print(Form("plotsMinBias_correctRate/h2_dRep_vsChi_L%d.png", ij), "png");


     TCanvas* tc2 = new TCanvas();
     tc2->cd();
     h2_dR_vsTrkQ[ij]->GetXaxis()->SetTitle(Form("dR layer %d", ij));
     h2_dR_vsTrkQ[ij]->GetYaxis()->SetTitle("muon trk quality");
     h2_dR_vsTrkQ[ij]->Draw("colz");
     tc2->Print(Form("plotsMinBias_correctRate/h2_dR_vsTrkQ_L%d.png", ij), "png");

     tc2->cd();
     h2_iR_vsEta_size[ij]->GetXaxis()->SetTitle(Form("eta layer %d", ij+firstLayer));
     h2_iR_vsEta_size[ij]->GetYaxis()->SetTitle(Form("|iR| layer %d", ij+firstLayer));
     h2_iR_vsEta_size[ij]->Draw("colz");
     h2_iR_vsEta_size[ij]->SetMarkerStyle(7);
     tc2->Print(Form("plotsMinBias_correctRate/h2_iR_vsEta_size_L%d.png", ij+firstLayer), "png");

     tc2->cd();
     h2_iR_vsEta[ij]->GetXaxis()->SetTitle(Form("eta layer %d", ij+firstLayer));
     h2_iR_vsEta[ij]->GetYaxis()->SetTitle(Form("|iR| layer %d", ij+firstLayer));
     h2_iR_vsEta[ij]->Draw(""); //"colz");
     h2_iR_vsEta[ij]->SetMarkerStyle(7);
     tc2->Print(Form("plotsMinBias_correctRate/h2_iR_vsEta_L%d.png", ij+firstLayer), "png");

     tc2->cd();
     h2_iR_vsR[ij]->GetXaxis()->SetTitle(Form("R layer %d", ij+firstLayer));
     h2_iR_vsR[ij]->GetYaxis()->SetTitle(Form("|iR| layer %d", ij+firstLayer));
     h2_iR_vsR[ij]->Draw(""); //"colz");
     h2_iR_vsR[ij]->SetMarkerStyle(7);
     tc2->Print(Form("plotsMinBias_correctRate/h2_iR_vsR_L%d.png", ij+firstLayer), "png");

   }
   //}//1==2


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
      }
    }
    h_Count_Phi->Fill(countFilled);

    float yMax;
    getXmax(dummyMIP, yMax);
    hfithisto->SetParameters(yMax, 1., 0.1);
    hfithisto->SetParLimits(0, yMax/2., yMax * 1.2);
    hfithisto->SetParLimits(1, 0., 10.);
    hfithisto->SetParLimits(2, 0., 0.3);
    std::cout << " histo entries = " << dummyMIP->GetEntries() << std::endl;
    TFitResultPtr r = dummyMIP->Fit("hfithisto", "RBS");
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
      h_Mip_Phi[savedCout]->Fit("hfithisto", "RB");

      //      tL.DrawLatex(0.55, 0.70, Form("by Chi2 = %d", rejectedByChi));
      // tL.DrawLatex(0.55, 0.60, Form("by 1st dR = %d", rejectedBy1stdR));
      // tL.DrawLatex(0.55, 0.50, Form("by prev dR = %d", rejectedByPrevdR));
      // tL.DrawLatex(0.55, 0.40, Form("found = %d", countFilled));

      tcM->Print(Form("plotsMinBias_correctRate/cellsPhi/h_Mip_Phi_%d_%d.png", iL+firstLayer, iR), "png");
      if(iL == 0) tcM->Print(Form("plotsMinBias_correctRate/cells/h_Mip_Phi_%d_%d.root", iL+firstLayer, iR), "root");

      ++savedCout;
    }
  }


  TF1* hfit = new TF1("hfit", "gaus", 0., 2.);
  TCanvas* tcC_Phi = new TCanvas();
  //  h_MPV_values->Rebin(2);
  hfit->SetParameters(h_MPV_values->GetEntries()/2., 1, 0.02);
  h_MPV_values->GetXaxis()->SetTitle("MPV values");
  h_MPV_values->Draw();
  h_MPV_values->Fit("hfit", "R");
  tcC_Phi->Print("plotsMinBias_correctRate/cellsPhi/h_MPV_values.png", "png");

  gPad->SetLogy();
  tcC_Phi->cd();
  h_Count_Phi->GetXaxis()->SetTitle("n counts in rings over #phi");
  h_Count_Phi->Draw();
  tcC_Phi->Print("plotsMinBias_correctRate/cellsPhi/h_Count_Phi.png", "png");
  tcC_Phi->Print("plotsMinBias_correctRate/cellsPhi/h_Count_Phi.root", "root");

  h_minNcounts->GetXaxis()->SetTitle("minimal N. counts");
  h_minNcounts->Draw();
  tcC_Phi->Print("plotsMinBias_correctRate/cellsPhi/h_minNcounts.png", "png");

  gStyle->SetOptStat(0);
  TCanvas* tcC2_Phi = new TCanvas();
  h2_iRvsLayer_MIP->GetXaxis()->SetTitle("layer");
  h2_iRvsLayer_MIP->GetYaxis()->SetTitle("iR");
  h2_iRvsLayer_MIP->GetZaxis()->SetRangeUser(0., 3.);
  h2_iRvsLayer_MIP->Draw("colz");
  tcC2_Phi->Print("plotsMinBias_correctRate/cellsPhi/h2_iRvsLayer_MIP.png", "png");
  tcC2_Phi->Print("plotsMinBias_correctRate/cellsPhi/h2_iRvsLayer_MIP.root", "root");



  std::cout << " all done now print results " << std::endl;


  std::cout << " all done - ciao " << std::endl;

  return 0;



}
