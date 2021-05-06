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
  int iColors[7] = {kOrange-3, kRed, kMagenta, kBlue, kCyan, kGreen+1, kGray+2};
  int iStyle[4] = {20, 21, 22, 23}; //kGray+1};
  //  int iColors[6] = {kGreen+1, kBlue, kRed};
  
  //  TFile* inF = TFile::Open("../test/singleMuon_newGun.root");
  //  TFile* inF = TFile::Open("../test/testMinBias.root");

  // TFile* inF = TFile::Open("../test/MinBias140PU_0.root");

  //TFile* inF = TFile::Open("/tmp/amartell/MinBias_PU140_withNewMB_pt4Cut_vArea.root");
  //  TFile* inF = TFile::Open("/tmp/amartell/MinBias_PU140_withNewMB_pt4Cut_vGoodM.root");
  TFile* inF = TFile::Open("/tmp/amartell/MinBias_PU140_withNewMB_pt4Cut_vGlbM.root");

  //config
  int firstLayer = 37;
  float Nevt = 2.09206e+06;
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
  std::vector<float> *crossZ = 0;
  std::vector<float> *crossL = 0;
  std::vector<float> *crossM = 0;
  std::vector<float> *muonChi = 0;
  std::vector<short int> *muonTrkQ = 0;
  std::vector<float> *recHitX = 0;
  std::vector<float> *recHitY = 0;
  std::vector<float> *recHitZ = 0;
  std::vector<float> *recHitSize = 0;
  std::vector<int> *recHitiPhi = 0;
  std::vector<int> *recHitiR = 0;
  std::vector<int> *recHitL = 0;
  std::vector<float> *recHitEne = 0;
  std::vector<float> *recHitMip = 0;
  std::vector<float> *recHitNoise = 0;


  newT->SetBranchAddress("muonP", &muonP);
  newT->SetBranchAddress("muonPt", &muonPt);
  newT->SetBranchAddress("muonEta", &muonEta);
  newT->SetBranchAddress("muonPhi", &muonPhi);
  newT->SetBranchAddress("crossX", &crossX);
  newT->SetBranchAddress("crossY", &crossY);
  newT->SetBranchAddress("crossZ", &crossZ);
  newT->SetBranchAddress("crossL", &crossL);
  newT->SetBranchAddress("crossM", &crossM);
  newT->SetBranchAddress("muonChi", &muonChi);
  newT->SetBranchAddress("muonTrkQ", &muonTrkQ);
  newT->SetBranchAddress("recHitX", &recHitX);
  newT->SetBranchAddress("recHitY", &recHitY);
  newT->SetBranchAddress("recHitZ", &recHitZ);
  newT->SetBranchAddress("recHitSize", &recHitSize);
  newT->SetBranchAddress("recHitiR", &recHitiR);
  newT->SetBranchAddress("recHitiPhi", &recHitiPhi);
  newT->SetBranchAddress("recHitL", &recHitL);
  newT->SetBranchAddress("recHitEne", &recHitEne);
  newT->SetBranchAddress("recHitMip", &recHitMip);
  newT->SetBranchAddress("recHitNoise", &recHitNoise);

  std::map<std::pair<int, std::pair<int, int> >, float > okChannels;
  std::map<std::pair<int, std::pair<int, int> >, std::vector<float> > enePerRh;
  std::map<std::pair<int, std::pair<int, int> >, std::vector<float> > mipPerRh;
  // std::map<std::pair<int, std::pair<int, int> >, std::vector<float> > sonPerRh;
  // std::map<std::pair<int, std::pair<int, int> >, std::vector<float> > momPerRh;
 
  std::map<std::pair<int, int>, float > okChannelsPhi;
  std::map<std::pair<int, int>, std::vector<float> > enePerPhi;
  std::map<std::pair<int, int>, std::vector<float> > mipPerPhi;
  // std::map<std::pair<int, int>, std::vector<float> > sonPerPhi;
  // std::map<std::pair<int, int>, std::vector<float> > momPerPhi;
  

  TH1F* mipAll = new TH1F("mipAll", "", 80, -1., 19.);
  TH2F* h2_YvsX[14];
  TH2F* h2_iRvsiPhi[14];
  TH1F* h_dX[14];
  TH1F* h_dY[14];
  TH1F* h_dR[14];

  TH2F* h2_dR_vsChi[14];
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
  TH1F* dX_wrt1st = new TH1F("dX_wrt1st", "", 1000, -20, 20);
  TH1F* dY_wrt1st = new TH1F("dY_wrt1st", "", 1000, -20, 20);
  TH1F* dX_wrtPrev = new TH1F("dX_wrtPrev", "", 1000, -20, 20);
  TH1F* dY_wrtPrev = new TH1F("dY_wrtPrev", "", 1000, -20, 20);

   for(int ij=0; ij<14; ++ij){
     h2_YvsX[ij] = new TH2F(Form("h2_YvsX_L%d", ij), "", 600, -300., 300., 600, -300., 300.);
     h2_iRvsiPhi[ij] = new TH2F(Form("h2_iRvsiPhi_L%d", ij), "", 100, -50., 50., 300, 0., 300.);
     h_dX[ij] = new TH1F(Form("h_dX_L%d", ij), "", 100, -10., 10.);
     h_dY[ij] = new TH1F(Form("h_dY_L%d", ij), "", 100, -10., 10.);
     //     h_dR[ij] = new TH1F(Form("h_dR_L%d", ij), "", 100, 0., 10.);
     h_dR[ij] = new TH1F(Form("h_dR_L%d", ij), "", 100, 0., 20.);
     h2_dR_vsChi[ij] = new TH2F(Form("h2_dR_vsChi_L%d", ij), "", 100, 0., 10., 100, 0., 10.);
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
   }

   int nEvents = newT->GetEntries();
   std::cout << " nEvents = " << nEvents << std::endl;
   for(int ij=0; ij<nEvents; ++ij){
   //   for(int ij=3; ij<4; ++ij){
     newT->GetEntry(ij);

     //     std::cout << " evento = " << ij << std::endl;
     std::map<int, std::vector<float>> rhX;
     std::map<int, std::vector<float>> rhY;
     std::map<int, std::vector<float>> rhZ;
     std::map<int, std::vector<float>> rhA;
     std::map<int, std::vector<float>> rhiPhi;
     std::map<int, std::vector<float>> rhiR;
     std::map<int, std::vector<float>> rhMip;
     std::map<int, std::vector<float>> rhEne;
     std::map<int, std::vector<int>> trkN;
     std::map<int, std::vector<float>> trkX;
     std::map<int, std::vector<float>> trkY;
     std::map<int, std::vector<float>> trkEta;
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
       trkEta[trkL].push_back(muonEta->at(iM));
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

       //       std::cout << " trk size = " << trkN[iL].size() << std::endl;
       for(auto iT : trkX[iL]){
	 ++iTc;
	 if(trkChi[iL][iTc] > 1.5) continue;

	 int trackIdx = trkN[iL][iTc];
	 if(foundFirstHit.find(trackIdx) == foundFirstHit.end()) foundFirstHit[trackIdx] = false;
	 //	 std::cout << "  trk " <<  " trkX = " <<  iT << " iL = " << iL << " eta = " << trkEta[iL][iTc] << " ij = " << ij << std::endl;

	 //	 continue;
	 float trackY = trkY[iL][iTc];
	 float trkR = sqrt(iT*iT + trackY*trackY);                                                                                                                        
	 float trkEtaOnLayer = asinh(zVal[iL]/trkR);                                                                                                                      
	 nMuons_vsEta_All[iL]->Fill(std::abs(trkEtaOnLayer));  
	 //std::cout << " iTc = " << iTc << std::endl;
	 // continue;

	 //std::cout << " now recHits" << std::endl;
	 //	 bool alreadyFilled = false;
	 //std::cout << " rhX[iL].size() = " << rhX[iL].size() << std::endl;
	 int iRc = -1;
	 for(auto iR : rhX[iL]){
	   ++iRc;
	   //	   std::cout << " iRc = " << iRc << std::endl; 

	   if(trkEta[iL][iTc] * rhZ[iL][iRc] < 0.) { 
	     //  std::cout << " opposite side " << std::endl; 
	     continue;
	   }

	   float dX = iR - iT;
	   float recY = rhY[iL][iRc];
	   float dY = recY - trackY;
	   float dR2 = dX*dX + dY*dY;
	   float dR = sqrt(dR2);
	   float reciR = rhiR[iL][iRc];

	   h_dX[iL]->Fill(dX);
	   h_dY[iL]->Fill(dY);
	   if(iL == 0)	   h_dR[iL]->Fill(dR);
	   h_dRvsP[iL]->Fill(dR, trkP[iL][iTc]);
	   // std::cout << " iL = " << iL << " recX = " << iR<< " iT = " << iT << " recY = " << recY << " trkY[recL][iTc] = " << trkY[iL][iTc] 
	   // 	     << " muonP = " << trkP[iL][iTc] << std::endl;
	   // std::cout << " dR = " << dR  << std::endl;

	   h2_dR_vsChi[iL]->Fill(dR, trkChi[iL][iTc]);
	   h2_dR_vsTrkQ[iL]->Fill(dR, trkQ[iL][iTc]);
	   // std::cout << " >>> trkChi[recL][iTc] = " << trkChi[recL][iTc] 
	   // 	     << " trkQ[recL][iTc] = " << trkQ[recL][iTc] << std::endl;



	   if(trkChi[iL][iTc] > 1.5){ // || trkQ[iL][iTc] != 5.){
	     //std::cout << " bad chi " << std::endl;
	     //mipPerPhi[channelPhi].push_back(-1);
	     continue;
	   }
	   if(foundFirstHit[iTc] == false && dR2 > rhA[iL][iRc] * 4./pigreco){
	     //std::cout << " bad dR " << std::endl;
	     //	     mipPerPhi[channelPhi].push_back(-2);
	     continue;
	   }
	  
	   /*
	   if(foundFirstHit[trackIdx] == false)	   
	     std::cout << " event = " << ij << " hitX = " << iR << " hitY = " << recY 
	      	       << " trkIdx = " << trackIdx << " iL = " << iL 
	     	       << " foundFirstHit[trackIdx] = " << foundFirstHit[trackIdx] << std::endl;
	   */
	   //	   return 10;

	   float dX_Prev = iR - xPrevHit[trackIdx];
	   float dY_Prev = recY - yPrevHit[trackIdx];
	   float dR2_Prev = dX_Prev * dX_Prev + dY_Prev * dY_Prev;
	   float dX_1st = iR - x1stHit[trackIdx];
	   float dY_1st = recY - y1stHit[trackIdx];

	   // std::cout  << " dR2_Prev = " << dR2_Prev << " rhA[iL][iRc] = " << rhA[iL][iRc] << " foundFirstHit[trackIdx] = " 
	   // 	      << foundFirstHit[trackIdx] << " iL = " << iL << std::endl;
	   if(foundFirstHit[trackIdx] == true){	   
	     dX_wrtPrev->Fill(iR - xPrevHit[trackIdx]);
	     dY_wrtPrev->Fill(recY - yPrevHit[trackIdx]);
	     dX_wrt1st->Fill(iR - x1stHit[trackIdx]);
	     dY_wrt1st->Fill(recY - y1stHit[trackIdx]);
	     
	     if(dR2_Prev >  rhA[iL][iRc] * 4./pigreco){
	       /*
	       std::cout << " other " << " hitX = " << iR << " hitY = " << recY 
			 << " trkIdx = " << trackIdx << " iL = " << iL
			 << " prevX = " << xPrevHit[trackIdx] << " prevY = " << yPrevHit[trackIdx] 
			 << " dxPrev = " << iR - xPrevHit[trackIdx] << " dyPrev = " << recY - yPrevHit[trackIdx] 
			 << " dR2_Prev = " << dR2_Prev << " rhA[iL][iRc] = " << rhA[iL][iRc] << std::endl;
	       */
	       continue;
	     }
	     xPrevHit[trackIdx] = iR;
	     yPrevHit[trackIdx] = recY;	
	     matchedTrack[trackIdx] += 1;       

	     h_TrkChi_passed->Fill(trkChi[iL][iTc]);
	     h_TrkQ_passed->Fill(trkQ[iL][iTc]);

	     h2_YvsX[iL]->Fill(iR, recY);
	     float reciPhi = rhiPhi[iL][iRc];
	     h2_iRvsiPhi[iL]->Fill(reciR, reciPhi);
	     
	     float recR = sqrt(iR*iR + recY*recY);
	     float recEta = asinh(rhZ[iL][iRc]/recR);
	     h2_iR_vsEta_size[iL]->Fill(std::abs(recEta), std::abs(reciR), rhA[iL][iRc]);
	     h2_iR_vsEta[iL]->Fill(std::abs(recEta), std::abs(reciR));
	     h2_iR_vsR[iL]->Fill(recR, std::abs(reciR));
	     mipAll->Fill(rhMip[iL][iRc]);	     
	     /*
	     std::cout << " >>>> OK " << " hitX = " << iR << " hitY = " << recY 
		       << " trkIdx = " << trackIdx << " iL = " << iL
		       << " prevX = " << xPrevHit[trackIdx] << " prevY = " << yPrevHit[trackIdx] 
		       << " dxPrev = " << iR - xPrevHit[trackIdx] << " dyPrev = " << recY - yPrevHit[trackIdx] 
		       << " dR2_Prev = " << dR2_Prev << " rhA[iL][iRc] = " << rhA[iL][iRc] << std::endl;
	     */
	   }
	   if(foundFirstHit[trackIdx] == false) {
	     x1stHit[trackIdx] = iR;
	     y1stHit[trackIdx] = recY;
	     xPrevHit[trackIdx] = iR;
	     yPrevHit[trackIdx] = recY;
	     foundFirstHit[trackIdx] = true; 
	     matchedTrack[trackIdx] = 1;
	     startLMatchedTrack[trackIdx] = iL;
	     //std::cout << " should be true " << std::endl;
	   }
	 }
	//  std::cout << " check :::::: trackIdx = " << trackIdx << " matchedTrack[trackIdx] = " << matchedTrack[trackIdx] 
       // 		   << " foundFirstHit[trackIdx] = " << foundFirstHit[trackIdx] << std::endl;
       }
     }//layer

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
	 //std::cout << " trkIdx = " << trackIdx << " matchedTrack[trackIdx] = " << matchedTrack[trackIdx] << " iL = " << iL << std::endl;
         if(matchedTrack[trackIdx] + startLMatchedTrack[trackIdx] < 5){
	   // std::cout << " missing hits " << std::endl; 
	   continue;
	 }
         if(foundFirstHit.find(trackIdx) == foundFirstHit.end()) foundFirstHit[trackIdx] = false;

	 float trackY = trkY[iL][iTc];
         float trkR = sqrt(iT*iT + trackY*trackY);

         float trkEtaOnLayer = asinh(zVal[iL]/trkR);

         bool alreadyFilled = false;
         int iRc = -1;
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

	   float reciR = rhiR[iL][iRc];
	   std::pair<int, int> channelPhi = std::pair<int, int>(iL, std::abs(reciR));
	   if(okChannelsPhi.find(channelPhi) != okChannelsPhi.end()) { /*std::cout << " channel again " << std::endl;*/ okChannelsPhi[channelPhi] += 1;}
	   else {  okChannelsPhi[channelPhi] = 1;}

	   if(foundFirstHit[trackIdx] == false && dR2 > rhA[iL][iRc] *  4./pigreco){
	     mipPerPhi[channelPhi].push_back(-2);
	     //std::cout << " 1st hit dR "<< std::endl;
	     continue;
	   }

	   float dX_Prev = iR - xPrevHit[trackIdx];
	   float dY_Prev = recY - yPrevHit[trackIdx];
	   float dR2_Prev = dX_Prev * dX_Prev + dY_Prev*dY_Prev;
	   float dX_1st = iR - x1stHit[trackIdx];
	   float dY_1st = recY - y1stHit[trackIdx];

	   if(foundFirstHit[trackIdx] == true){	   
	     
	     if(dR2_Prev > rhA[iL][iRc] * 4./pigreco){
	       mipPerPhi[channelPhi].push_back(-1);
	       //std::cout << " prev hit dR "<< std::endl;
	       continue;
	     }
	     xPrevHit[trackIdx] = iR;
	     yPrevHit[trackIdx] = recY;	
	   }
	   if(foundFirstHit[trackIdx] == false) {
	     x1stHit[trackIdx] = iR;
	     y1stHit[trackIdx] = recY;
	     xPrevHit[trackIdx] = iR;
	     yPrevHit[trackIdx] = recY;
	     foundFirstHit[trackIdx] = true; 
	   }

	   //here rate for trk
	   if(!alreadyFilled){
	     nMuons_vsEta[iL]->Fill(std::abs(trkEtaOnLayer));
	     alreadyFilled = true;
	     //std::cout << " trkEtaOnLayer = " << trkEtaOnLayer << " iL = " << iL << " ij = " << ij << std::endl;
	   }

	   //	   std::cout << " FILLING  channelPhi layer = " << channelPhi.first << " and iR = " << channelPhi.second << std::endl;
	   mipPerPhi[channelPhi].push_back(rhMip[iL][iRc]);
	 }//recHit loop
       }// trk loop
     } // layer loop
   }// events;

   //   return 5;
   std::cout << " okChannels.size = " << okChannels.size() << " okChannelPhi.size() = " << okChannelsPhi.size() << std::endl;

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

     h_dRvsP[ij]->GetXaxis()->SetTitle(Form("dR layer %d", ij));
     h_dRvsP[ij]->GetYaxis()->SetTitle("muon P");
     h_dRvsP[ij]->Draw("colz");
     tc->Print(Form("plotsMinBias_correctRate/h_dRvsdP_L%d.png", ij), "png");

     h2_dR_vsChi[ij]->GetXaxis()->SetTitle(Form("dR layer %d", ij));
     h2_dR_vsChi[ij]->GetYaxis()->SetTitle("muon trk chi2");
     h2_dR_vsChi[ij]->Draw("colz");
     tc->Print(Form("plotsMinBias_correctRate/h2_dR_vsChi_L%d.png", ij), "png");

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
     // tc->Print(Form("plotsMinBias_correctRate/h2_iR_vsR_L%d.root", ij+firstLayer), "root");

     //    tc->Print(Form("plotsMinBias_correctRate/h2_iR_vsEta_L%d.root", ij), "root"); 
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
      tL.DrawLatex(0.55, 0.60, Form("by 1st dR = %d", rejectedBy1stdR));
      tL.DrawLatex(0.55, 0.50, Form("by prev dR = %d", rejectedByPrevdR));
      tL.DrawLatex(0.55, 0.40, Form("found = %d", countFilled));

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
