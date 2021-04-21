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



void lookAtMIP_pt4Cut(){
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
  
  int firstLayer = 37;
  //  TFile* inF = TFile::Open("../test/singleMuon_newGun.root");
  //  TFile* inF = TFile::Open("../test/testMinBias.root");

  // TFile* inF = TFile::Open("../test/MinBias140PU_0.root");

  //  TFile* inF = TFile::Open("/tmp/amartell/MinBias_PU140_withNewMB_pt4Cut_ab.root");
  TFile* inF = TFile::Open("/tmp/amartell/MinBias_PU140_withNewMB_pt4Cut_vChi.root");


  bool doRates = false;
  if(doRates){
  // float xSec = 1947. * 1.e6 * 1.e-24 * 1.e-12; //cm2
  // float Lumi = 5. * 1.e34; //cm-2 s-1
  float Nevt = 2.09206e+06;
  // float timeEqui = Nevt / (Lumi*xSec); //s

  float timeEqui = Nevt / 40.e6; //s

  float weights[51] = {0, 
		       8.89454, 10.9379, 10.9379, 10.9379, 10.9379, 10.9379, 10.9379, 10.9379, 10.9379, 10.9329, 
		       10.9329, 10.9379, 10.9379, 10.9382, 10.9382, 10.9382, 10.9382, 10.9382, 10.9382, 10.9382, 
		       10.9382, 10.9382, 10.9382, 10.9382, 10.9382, 10.9382, 10.9382, 32.3321, 51.5743, 51.4442, 
		       51.4442, 51.4442, 51.4442, 51.4442, 51.4442, 51.4442, 51.4442, 51.4442, 51.4442, 69.5131, 
 /*da 41 a 50*/	       87.582, 87.582, 87.582, 87.582, 87.582, 87.2146, 86.8883, 86.9295, 86.9295, 86.9295}; 



 

  float zVal[14] = {411.29, 416.739, 422.187, 427.636, 436.172, 444.722, 453.263, 461.817, 470.371, 478.925, 487.47, 496.024, 504.577, 513.127};
  float pigreco = 3.1415926535;
  //////// analyze rate

  TH1F* nMuons_vsEta[14];
  TH1F* nMuonsPercm2_vsEta[14];
  TH1F* ratePercm2_vsEta[14];
  TH1F* cm2_vsEta[14];
  for(int ij=0; ij<14; ++ij){
    nMuons_vsEta[ij] = (TH1F*)inF->Get(Form("ana/nMuons_vsEta_L%d", ij+firstLayer));
    //    nMuons_vsEta[ij]->Rebin(4);
    nMuonsPercm2_vsEta[ij] = new TH1F(Form("nMuonsPercm2_vsEta_L%d", ij+firstLayer), "", 400, 0., 4.);
    ratePercm2_vsEta[ij] = new TH1F(Form("ratePercm2_vsEta_L%d", ij+firstLayer), "", 400, 0., 4.);
    cm2_vsEta[ij] = new TH1F(Form("cm2_vsEta_L%d", ij+firstLayer), "", 400, 0., 4.);

    float binW = (4. - 0. ) / 400.;

    for(int iB=2; iB<nMuonsPercm2_vsEta[ij]->GetNbinsX()-1; ++iB){
      float rMax = zVal[ij]/sinh(0. + binW * iB);     
      float rMin = zVal[ij]/sinh(0. + binW * (iB+1));     
      float area = 2. * (rMax*rMax - rMin*rMin) * pigreco;
      std::cout << " rMin = " << rMin << " rMax = " << rMax << " etaL = " << 1 + binW * (iB) << " etaH  = " << (1 + binW * (iB+1)) << " area = " << area << std::endl;
      nMuonsPercm2_vsEta[ij]->SetBinContent(iB, nMuons_vsEta[ij]->GetBinContent(iB)/area);
      ratePercm2_vsEta[ij]->SetBinContent(iB, nMuons_vsEta[ij]->GetBinContent(iB)/area/timeEqui);
      cm2_vsEta[ij]->SetBinContent(iB, area);
    }

    TCanvas* tmuo = new TCanvas();
    tmuo->cd();
    nMuonsPercm2_vsEta[ij]->GetXaxis()->SetTitle(Form("#eta layer %d", ij+firstLayer));
    nMuonsPercm2_vsEta[ij]->GetYaxis()->SetTitle("n Muons / cm2");
    nMuonsPercm2_vsEta[ij]->Draw();
    tmuo->Print(Form("plotsMuon_MinBias_pT4Cut/nMuonsPercm2_vsEta_L%d.png", ij+firstLayer), "png");
    //    tmuo->Print(Form("plotsMuon_MinBias_pT4Cut/nMuonsPercm2_vsEta_L%d.root", ij+firstLayer), "root");

    nMuons_vsEta[ij]->GetXaxis()->SetTitle(Form("#eta layer %d", ij+firstLayer));
    nMuons_vsEta[ij]->Draw();
    tmuo->Print(Form("plotsMuon_MinBias_pT4Cut/nMuons_vsEta_L%d.png", ij+firstLayer), "png");

    gPad->SetLogy();
    cm2_vsEta[ij]->GetYaxis()->SetRangeUser(0.5, 5.e4);
    cm2_vsEta[ij]->GetXaxis()->SetTitle(Form("#eta layer %d", ij+firstLayer));
    cm2_vsEta[ij]->GetYaxis()->SetTitle("cm2");
    cm2_vsEta[ij]->Draw();
    tmuo->Print(Form("plotsMuon_MinBias_pT4Cut/cm2_vsEta_L%d.png", ij+firstLayer), "png");

    ratePercm2_vsEta[ij]->GetYaxis()->SetRangeUser(0.5, 1.e4);
    ratePercm2_vsEta[ij]->GetXaxis()->SetTitle(Form("#eta layer %d", ij+firstLayer));
    ratePercm2_vsEta[ij]->GetYaxis()->SetTitle("n Muons / cm2 / s");
    ratePercm2_vsEta[ij]->Draw();
    tmuo->Print(Form("plotsMuon_MinBias_pT4Cut/ratePercm2_vsEta_L%d.png", ij+firstLayer), "png");
  }
  //  return;
  }


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
  newT->SetBranchAddress("recHitiR", &recHitiR);
  newT->SetBranchAddress("recHitiPhi", &recHitiPhi);
  newT->SetBranchAddress("recHitL", &recHitL);
  newT->SetBranchAddress("recHitEne", &recHitEne);
  newT->SetBranchAddress("recHitMip", &recHitMip);
  newT->SetBranchAddress("recHitNoise", &recHitNoise);

  std::map<std::pair<int, std::pair<int, int> >, float > okChannels;
  std::map<std::pair<int, std::pair<int, int> >, std::vector<float> > enePerRh;
  std::map<std::pair<int, std::pair<int, int> >, std::vector<float> > mipPerRh;
  std::map<std::pair<int, std::pair<int, int> >, std::vector<float> > sonPerRh;
  std::map<std::pair<int, std::pair<int, int> >, std::vector<float> > momPerRh;
 
  std::map<std::pair<int, int>, float > okChannelsPhi;
  std::map<std::pair<int, int>, std::vector<float> > enePerPhi;
  std::map<std::pair<int, int>, std::vector<float> > mipPerPhi;
  std::map<std::pair<int, int>, std::vector<float> > sonPerPhi;
  std::map<std::pair<int, int>, std::vector<float> > momPerPhi;


  TH2F* h2_YvsX[14];
  TH2F* h2_iRvsiPhi[14];
  TH1F* h_dX[14];
  TH1F* h_dY[14];
  TH1F* h_dR[14];

  TH2F* h2_dR_vsChi[14];
  TH2F* h2_dR_vsTrkQ[14];

  TH2F* h_dRvsP[14];
  
  TH2F* h2_iR_vsEta[14];
  TH2F* h2_iR_vsR[14];

  TH1F* h_etaMuon = new TH1F("h_etaMuon", "", 300, 1., 4.);
  TH1F* h_ptMuon = new TH1F("h_ptMuon", "", 100, 0., 100.);

  TH1F* h_TrkChi_all = new TH1F("h_TrkChi_all", "", 1000, 0., 20.);
  TH1F* h_TrkQ_all = new TH1F("h_TrkQ_all", "", 10, 0., 10.);
  TH1F* h_TrkChi_passed = new TH1F("h_TrkChi_passed", "", 1000, 0., 20.);
  TH1F* h_TrkQ_passed = new TH1F("h_TrkQ_passed", "", 10, 0., 10.);

   for(int ij=0; ij<14; ++ij){
     h2_YvsX[ij] = new TH2F(Form("h2_YvsX_L%d", ij), "", 600, -300., 300., 600, -300., 300.);
     h2_iRvsiPhi[ij] = new TH2F(Form("h2_iRvsiPhi_L%d", ij), "", 100, -50., 50., 300, 0., 300.);
     h_dX[ij] = new TH1F(Form("h_dX_L%d", ij), "", 100, -10., 10.);
     h_dY[ij] = new TH1F(Form("h_dY_L%d", ij), "", 100, -10., 10.);
     h_dR[ij] = new TH1F(Form("h_dR_L%d", ij), "", 100, 0., 10.);
     h2_dR_vsChi[ij] = new TH2F(Form("h2_dR_vsChi_L%d", ij), "", 100, 0., 10., 100, 0., 10.);
     h2_dR_vsTrkQ[ij] = new TH2F(Form("h2_dR_vsTrkQ_L%d", ij), "", 100, 0., 10., 10, 0., 10.);
     h_dRvsP[ij] = new TH2F(Form("h_dRvsP_L%d", ij), "", 100, 0., 10., 100, 0., 10.);
     h2_iR_vsEta[ij] = new TH2F(Form("h2_iR_vsEta_L%d", ij), "", 300, 1., 4., 50, 0., 50.);
     h2_iR_vsR[ij] = new TH2F(Form("h2_iR_vsR_L%d", ij), "", 500, 100., 300., 50, 0., 50.);
   }

   int nEvents = newT->GetEntries();
   std::cout << " nEvents = " << nEvents << std::endl;
   for(int ij=0; ij<nEvents; ++ij){
   //   for(int ij=0; ij<20; ++ij){
     newT->GetEntry(ij);

     std::cout << " evento = " << ij << std::endl;
     std::map<int, std::vector<float>> rhX;
     std::map<int, std::vector<float>> rhY;
     std::map<int, std::vector<float>> rhMip;
     std::map<int, std::vector<float>> trkX;
     std::map<int, std::vector<float>> trkY;
     std::map<int, std::vector<float>> trkEta;
     std::map<int, std::vector<float>> trkPt;
     std::map<int, std::vector<float>> trkP;
     std::map<int, std::vector<float>> trkChi;
     std::map<int, std::vector<short int>> trkQ;

     //    std::cout << " muonP->size() = " << muonP->size() << std::endl;
     for(int iM=0; iM<muonP->size(); ++iM){
       trkX[crossL->at(iM)].push_back(crossX->at(iM));
       trkY[crossL->at(iM)].push_back(crossY->at(iM));
       trkEta[crossL->at(iM)].push_back(muonEta->at(iM));
       trkPt[crossL->at(iM)].push_back(muonPt->at(iM));
       trkP[crossL->at(iM)].push_back(muonP->at(iM));
       trkChi[crossL->at(iM)].push_back(muonChi->at(iM));
       trkQ[crossL->at(iM)].push_back(muonTrkQ->at(iM));
       if(crossL->at(iM) == 37){
	 h_etaMuon->Fill(std::abs(muonEta->at(iM)));
	 h_ptMuon->Fill(std::abs(muonPt->at(iM)));

	 h_TrkChi_all->Fill(muonChi->at(iM));
	 h_TrkQ_all->Fill(muonTrkQ->at(iM));
       }
     }



     //    float mindR = 10.;

     //    std::cout << " recHitX->size() = " << recHitX->size() << std::endl;
     for(int iR=0; iR<recHitX->size(); ++iR){
       float recX = recHitX->at(iR);
       float recY = recHitY->at(iR);
       float recZ = recHitZ->at(iR);
       int reciPhi = recHitiPhi->at(iR);
       int reciR = recHitiR->at(iR);
       int recL = recHitL->at(iR);
       //      std::cout << " weights layer = " << recL << " w = " << weights[recL] << std::endl;
       float recEne = recHitEne->at(iR);
       float recMip = recHitMip->at(iR); // * 1.e3 / weights[recL] / 0.9;

       /*
       if(trkX[recL].size() == 0) std::cout << " size trk on layer = 0 " << std::endl;
       else{
	 std::cout << " size trk on layer =  " << trkX[recL].size() << std::endl;
	 std::cout << " recHit >>> recL - firstLayer = " << recL - firstLayer  << " trkP = " << trkP[recL][0] 
		   << " trkPt = " << trkPt[recL][0] << " trkEta = " << trkEta[recL][0] << std::endl;
       }
       */

       int iTc = -1;
       for(auto iT : trkX[recL]){
	   ++iTc;
	   if(trkEta[recL][iTc] * recZ < 0.) continue;
	   //	  std::cout << " recL = " << recL << " recX = " << recX << " iT = " << iT << " recY = " << recY << " trkY[recL][iTc] = " << trkY[recL][iTc] << " muonP = " << trkP[recL][iTc] << std::endl;
	   float dX = recX - iT;
	   float dY = recY - trkY[recL][iTc];
	   float dR = sqrt(dX*dX + dY*dY);
	   h_dX[recL - firstLayer]->Fill(dX);
	   h_dY[recL - firstLayer]->Fill(dY);
	   h_dR[recL - firstLayer]->Fill(dR);
	   h_dRvsP[recL - firstLayer]->Fill(dR, trkP[recL][iTc]);
	   //	  if(recMip < 40.) std::cout << " recMIP = " << recMip << " dR = " << dR << std::endl;

	   //	  std::cout << " dR = " << dR  << std::endl;

	   h2_dR_vsChi[recL - firstLayer]->Fill(dR, trkChi[recL][iTc]);
	   h2_dR_vsTrkQ[recL - firstLayer]->Fill(dR, trkQ[recL][iTc]);
	   // std::cout << " >>> trkChi[recL][iTc] = " << trkChi[recL][iTc] 
	   // 	     << " trkQ[recL][iTc] = " << trkQ[recL][iTc] << std::endl;


	   if(dR > 4.) continue;

	   if(recL == 9){
	     h_TrkChi_passed->Fill(trkChi[recL][iTc]);
	     h_TrkQ_passed->Fill(trkQ[recL][iTc]);
	   }

	   //if(recL - firstLayer == 0){
	   //	   std::cout << " layer " << recL - firstLayer << " mip = " << recMip << " dR = " << dR << std::endl;

	   h2_YvsX[recL - firstLayer]->Fill(recX, recY);
	   h2_iRvsiPhi[recL - firstLayer]->Fill(reciR, reciPhi);

	   float recR = sqrt(recX*recX + recY*recY);
	   float recEta = asinh(recZ/recR);
	   h2_iR_vsEta[recL - firstLayer]->Fill(std::abs(recEta), std::abs(reciR));
	   h2_iR_vsR[recL - firstLayer]->Fill(recR, std::abs(reciR));

	   //	  if(reciR < 0) std::cout << "recEta = " << recEta << " recZ = " << recZ << std::endl;

	   // std::pair<float, float> coordPair = std::pair<float, float>(recX, recY);
	   // std::pair<int, std::pair<float, float> > channel = std::pair<int, std::pair<float, float> >(recL, coordPair);
	   std::pair<int, int> coordPair = std::pair<int, int>(std::abs(reciR), reciPhi);                                                                                      
	   std::pair<int, std::pair<int, int> > channel = std::pair<int, std::pair<int, int> >(recL, coordPair); 	  


	   if(okChannels.find(channel) != okChannels.end()) { /*std::cout << " channel again " << std::endl;*/ okChannels[channel] += 1;}
	   else {  okChannels[channel] = 1;}

	   //	  std::cout << " recMip = " << recMip << " recEne = " << recEne << " SoN = " << 1.*recMip/recHitNoise->at(iR) << std::endl;
	   enePerRh[channel].push_back(recEne);
	   mipPerRh[channel].push_back(recMip);
	   sonPerRh[channel].push_back(recHitNoise->at(iR));
	   momPerRh[channel].push_back(trkP[recL][iTc]);

	   std::pair<int, int> channelPhi = std::pair<int, int>(recL, std::abs(reciR));
	   if(okChannelsPhi.find(channelPhi) != okChannelsPhi.end()) { /*std::cout << " channel again " << std::endl;*/ okChannelsPhi[channelPhi] += 1;}
	   else {  okChannelsPhi[channelPhi] = 1;}
	   enePerPhi[channelPhi].push_back(recEne);
	   mipPerPhi[channelPhi].push_back(recMip);
	   sonPerPhi[channelPhi].push_back(recHitNoise->at(iR));
	   momPerPhi[channelPhi].push_back(trkP[recL][iTc]);
       }
     }//recHits

   }// events;

   std::cout << " okChannels.size = " << okChannels.size() << " okChannelPhi.size() = " << okChannelsPhi.size() << std::endl;

   TCanvas* tChi = new TCanvas();
   tChi->cd();
   h_TrkChi_all->GetXaxis()->SetTitle("trk chi all"); 
   h_TrkChi_all->Draw();
   tChi->Print("plotsMuon_MinBias_pT4Cut/h_TrkChi_all.png", "png");
   tChi->Print("plotsMuon_MinBias_pT4Cut/h_TrkChi_all.root", "root");

   h_TrkQ_all->GetXaxis()->SetTitle("trk Q all"); 
   h_TrkQ_all->Draw();
   tChi->Print("plotsMuon_MinBias_pT4Cut/h_TrkQ_all.png", "png");
   tChi->Print("plotsMuon_MinBias_pT4Cut/h_TrkQ_all.root", "root");

   h_TrkChi_passed->GetXaxis()->SetTitle("trk chi all"); 
   h_TrkChi_passed->Draw();
   tChi->Print("plotsMuon_MinBias_pT4Cut/h_TrkChi_passed.png", "png");
   tChi->Print("plotsMuon_MinBias_pT4Cut/h_TrkChi_passed.root", "root");

   h_TrkQ_passed->GetXaxis()->SetTitle("trk Q all"); 
   h_TrkQ_passed->Draw();
   tChi->Print("plotsMuon_MinBias_pT4Cut/h_TrkQ_passed.png", "png");
   tChi->Print("plotsMuon_MinBias_pT4Cut/h_TrkQ_passed.root", "root");


   return;
   if(1 == 2){
   TCanvas* tMuo = new TCanvas();
   tMuo->cd();
   h_etaMuon->GetXaxis()->SetTitle("muon #eta");
   h_etaMuon->Draw();
   tMuo->Print("plotsMuon_MinBias_pT4Cut/etaMuon.png", "png");
   tMuo->Print("plotsMuon_MinBias_pT4Cut/etaMuon.root", "root");

   gPad->SetLogy();
   h_ptMuon->GetXaxis()->SetTitle("muon pT");
   h_ptMuon->Draw();
   tMuo->Print("plotsMuon_MinBias_pT4Cut/ptMuon.png", "png");
   tMuo->Print("plotsMuon_MinBias_pT4Cut/ptMuon.root", "root");
   }

   //   if(1 == 2){
   //return;
   for(int ij=0; ij<14; ++ij){
     TCanvas* tc = new TCanvas();
     tc->cd();
     h2_YvsX[ij]->GetXaxis()->SetTitle(Form("layer %d", ij));
     h2_YvsX[ij]->Draw("colz");
     tc->Print(Form("plotsMinBias_pT4Cut/h2_YvsX_L%d.png", ij), "png");

     h2_iRvsiPhi[ij]->GetXaxis()->SetTitle(Form("|iR| layer %d", ij));
     h2_iRvsiPhi[ij]->GetYaxis()->SetTitle("iPhi");
     h2_iRvsiPhi[ij]->Draw("colz");
     tc->Print(Form("plotsMinBias_pT4Cut/h2_iRvsiPhi_L%d.png", ij), "png");
     //   continue;
     h_dX[ij]->GetXaxis()->SetTitle(Form("dX layer %d", ij));
     h_dX[ij]->Draw();
     tc->Print(Form("plotsMinBias_pT4Cut/h_dX_L%d.png", ij), "png");

     h_dY[ij]->GetXaxis()->SetTitle(Form("dY layer %d", ij));
     h_dY[ij]->Draw();
     tc->Print(Form("plotsMinBias_pT4Cut/h_dY_L%d.png", ij), "png");

     h_dR[ij]->GetXaxis()->SetTitle(Form("dR layer %d", ij+firstLayer));
     h_dR[ij]->Draw();
     tc->Print(Form("plotsMinBias_pT4Cut/h_dR_L%d.png", ij+firstLayer), "png");

     h_dRvsP[ij]->GetXaxis()->SetTitle(Form("dR layer %d", ij));
     h_dRvsP[ij]->GetYaxis()->SetTitle("muon P");
     h_dRvsP[ij]->Draw("colz");
     tc->Print(Form("plotsMinBias_pT4Cut/h_dRvsdP_L%d.png", ij), "png");

     h2_dR_vsChi[ij]->GetXaxis()->SetTitle(Form("dR layer %d", ij));
     h2_dR_vsChi[ij]->GetYaxis()->SetTitle("muon trk chi2");
     h2_dR_vsChi[ij]->Draw("colz");
     tc->Print(Form("plotsMinBias_pT4Cut/h2_dR_vsChi_L%d.png", ij), "png");

     TCanvas* tc2 = new TCanvas();
     tc2->cd();
     h2_dR_vsTrkQ[ij]->GetXaxis()->SetTitle(Form("dR layer %d", ij));
     h2_dR_vsTrkQ[ij]->GetYaxis()->SetTitle("muon trk quality");
     h2_dR_vsTrkQ[ij]->Draw("colz");
     tc2->Print(Form("plotsMinBias_pT4Cut/h2_dR_vsTrkQ_L%d.png", ij), "png");

     //     if(ij == 1)     return;
     h2_iR_vsEta[ij]->GetXaxis()->SetTitle(Form("eta layer %d", ij+firstLayer));
     h2_iR_vsEta[ij]->GetYaxis()->SetTitle(Form("|iR| layer %d", ij+firstLayer));
     h2_iR_vsEta[ij]->Draw(""); //"colz");
     h2_iR_vsEta[ij]->SetMarkerStyle(7);
     tc->Print(Form("plotsMinBias_pT4Cut/h2_iR_vsEta_L%d.png", ij+firstLayer), "png");

     h2_iR_vsR[ij]->GetXaxis()->SetTitle(Form("R layer %d", ij+firstLayer));
     h2_iR_vsR[ij]->GetYaxis()->SetTitle(Form("|iR| layer %d", ij+firstLayer));
     h2_iR_vsR[ij]->Draw(""); //"colz");
     h2_iR_vsR[ij]->SetMarkerStyle(7);
     tc->Print(Form("plotsMinBias_pT4Cut/h2_iR_vsR_L%d.png", ij+firstLayer), "png");
     // tc->Print(Form("plotsMinBias_pT4Cut/h2_iR_vsR_L%d.root", ij+firstLayer), "root");

     //    tc->Print(Form("plotsMinBias_pT4Cut/h2_iR_vsEta_L%d.root", ij), "root");
   }
   //}//1==2
   return;


  TH1F* h_Ene[300];
  TH1F* h_Mip[300];
  TH1F* h_SoN[300];
  TH2F* h_MipvsP[300];
  TH1F* h_Count = new TH1F("h_Count", "", 20, 0., 20.);

  int iC = -1;
  if(1 == 2){
  for(auto ic : okChannels){
    ++iC;
    
    int iL = ic.first.first;
    int iR = ic.first.second.first;
    int iPhi = ic.first.second.second;

    int iVal = ic.second;
    // if(iVal > 1) std::cout << " iL = " << iL << " iR = " << iR << " iPhi = " << iPhi << " iVal = " << iVal << std::endl;
    h_Count->Fill(iVal);

    //    continue;
    if(iVal < 7) continue;

    if(iL == 37 && iR == 18){
      h_Ene[iC] = new TH1F(Form("h_Ene_%d_%d_%d", iL, iR, iPhi), "", 1000, 0., 1.);
      h_Mip[iC] = new TH1F(Form("h_Mip_%d_%d_%d", iL, iR, iPhi), "", 110, -2., 20);
      h_SoN[iC] = new TH1F(Form("h_SoN_%d_%d_%d", iL, iR, iPhi), "", 100, 0., 1.e3);
      h_MipvsP[iC] = new TH2F(Form("h_MipvsP_%d_%d_%d", iL, iR, iPhi), "", 100, 0., 1.e3, 100, 0., 100.);

      for(int ij=0; ij<iVal; ++ij){
	h_Ene[iC]->Fill(enePerRh[ic.first][ij]);
	//std::cout << " energy = " << enePerRh[ic.first][ij]  << std::endl;
	h_Mip[iC]->Fill(mipPerRh[ic.first][ij]);
	h_SoN[iC]->Fill(sonPerRh[ic.first][ij]);
	h_MipvsP[iC]->Fill(mipPerRh[ic.first][ij], momPerRh[ic.first][ij]);
      }
  
    TCanvas* tcM = new TCanvas();
    tcM->cd();
    // h_Ene[iC]->GetXaxis()->SetTitle("recHit energy (GeV)");
    // h_Ene[iC]->Draw();
    // tcM->Print(Form("plotsMinBias_pT4Cut/cells/h_Ene_%d_%d_%d.png", iL, iR, iPhi), "png");
    // return;
    h_Mip[iC]->GetXaxis()->SetTitle("MIP");
    h_Mip[iC]->Draw();
    tcM->Print(Form("plotsMinBias_pT4Cut/cells/h_Mip_%d_%d_%d.png", iL, iR, iPhi), "png");

    // h_SoN[iC]->GetXaxis()->SetTitle("SoN");
    // h_SoN[iC]->Draw();
    // tcM->Print(Form("plotsMinBias_pT4Cut/cells/h_SoN_%d_%d_%d.png", iL, iR, iPhi), "png");

    // h_MipvsP[iC]->GetXaxis()->SetTitle("MIP");
    // h_MipvsP[iC]->GetYaxis()->SetTitle("muon momentum (GeV)");
    // h_MipvsP[iC]->Draw("colz");
    // tcM->Print(Form("plotsMinBias_pT4Cut/cells/h_MipvsP_%d_%d_%d.png", iL, iR, iPhi), "png");
    }//if(iL == 37 && iR == 18){
  }

  TCanvas* tcC = new TCanvas();
  gPad->SetLogy();
  tcC->cd();
  h_Count->GetXaxis()->SetTitle("n counts in single cells");
  h_Count->Draw();
  tcC->Print("plotsMinBias_pT4Cut/cells/h_Count.png", "png");
  } //1==2

  std::cout << " now plotting " << std::endl;
  //  return;
  // now integrated vs Phi
  TH1F* h_Ene_Phi[200];
  TH1F* h_Mip_Phi[200];
  TH1F* h_SoN_Phi[200];
  TH2F* h_MipvsP_Phi[200];
  TH1F* h_Count_Phi = new TH1F("h_Count_Phi", "", 500, 0., 2000.);

  TH2F* h2_iRvsLayer_MIP = new TH2F("h2_iRvsLayer_MIP", "", 15, 36, 51, 50, 0., 50.);
  TH2F* h2_iRvsLayer_MIPerr = new TH2F("h2_iRvsLayer_MIPerr", "", 15, 36, 51, 50, 0., 50.);
  TH1F* h_minNcounts = new TH1F("h_minNcounts", "", 400, 0., 800.);
  TH1F* h_MPV_values = new TH1F("h_MPV_values", "", 500., 0., 5.);
  int savedCout = 0;
  iC = -1;


  TF1* hfithisto = new TF1("hfithisto", "landaun",-1, 10);
  hfithisto->SetLineColor(kBlue);
  for(auto ic : okChannelsPhi){
    ++iC;
    
    int iL = ic.first.first;
    int iR = ic.first.second;

    int iVal = ic.second;
    //    if(iVal > 500) std::cout << " iL = " << iL << " iR = " << iR << " iPhi = " << iPhi << " iVal = " << iVal << std::endl;
    h_Count_Phi->Fill(iVal);

    TH1F* dummyMIP = new TH1F("dummyMIP", "", 110, -2., 20.);
    for(int ij=0; ij<iVal; ++ij){
      if(mipPerPhi[ic.first].size() != iVal) std::cout << " PROBLEM!!!  " << std::endl;
      dummyMIP->Fill(mipPerPhi[ic.first][ij]);
    }


    float yMax;
    getXmax(dummyMIP, yMax);
    hfithisto->SetParameters(yMax, 1., 0.1);
    hfithisto->SetParLimits(0, yMax/2., yMax * 1.2);
    hfithisto->SetParLimits(1, 0., 10.);
    hfithisto->SetParLimits(2, 0., 0.3);

    TFitResultPtr r = dummyMIP->Fit("hfithisto", "RBS");
    if(r != 0) continue;
    float meanV = hfithisto->GetParameter(1);
    delete dummyMIP;
    if(meanV < 0) continue;
    h_minNcounts->Fill(iVal);
    h2_iRvsLayer_MIP->Fill(iL, iR, hfithisto->GetParameter(1));
    h2_iRvsLayer_MIPerr->Fill(iL, iR, hfithisto->GetParError(1));
    h_MPV_values->Fill(hfithisto->GetParameter(1));

    //    continue;
    //    if(iVal < 500) continue;

    if(savedCout < 200){
      h_Ene_Phi[savedCout] = new TH1F(Form("h_Ene_Phi_%d_%d", iL, iR), "", 1000, 0., 1.);
      h_Mip_Phi[savedCout] = new TH1F(Form("h_Mip_Phi_%d_%d", iL, iR), "", 110, -2., 20.);
      h_SoN_Phi[savedCout] = new TH1F(Form("h_SoN_Phi_%d_%d", iL, iR), "", 100, 0., 1.e3);
      h_MipvsP_Phi[savedCout] = new TH2F(Form("h_MipvsP_Phi_%d_%d", iL, iR), "", 100, 0., 1.e3, 100, 0., 100.);

      for(int ij=0; ij<iVal; ++ij){
	h_Ene_Phi[savedCout]->Fill(enePerPhi[ic.first][ij]);
	//std::cout << " energy = " << enePerRh[ic.first][ij]  << std::endl;


	h_Mip_Phi[savedCout]->Fill(mipPerPhi[ic.first][ij]);
	h_SoN_Phi[savedCout]->Fill(sonPerPhi[ic.first][ij]);
	h_MipvsP_Phi[savedCout]->Fill(mipPerPhi[ic.first][ij], momPerPhi[ic.first][ij]);
      }

      TCanvas* tcM = new TCanvas();
      tcM->cd();
      /*
      h_Ene_Phi[savedCout]->GetXaxis()->SetTitle("recHit energy (GeV)");
      h_Ene_Phi[savedCout]->Draw();
      tcM->Print(Form("plotsMinBias_pT4Cut/cells/h_Ene_Phi_%d_%d.png", iL, iR), "png");
      */

      h_Mip_Phi[savedCout]->GetXaxis()->SetTitle("MIP");
      h_Mip_Phi[savedCout]->Draw();
      hfithisto->SetParameters(yMax, 1., 0.1);
      hfithisto->SetParLimits(0, yMax/2., yMax * 1.2);
      hfithisto->SetParLimits(1, 0., 10.);
      hfithisto->SetParLimits(2, 0., 0.3);
      h_Mip_Phi[savedCout]->Fit("hfithisto", "RB");
      //hfithisto->Draw("same");
      // hfithisto->SetParameter(0, 1.);
      // hfithisto->SetParameter(1, 0.5 );
      // h_Mip_Phi[savedCout]->Fit("hfithisto", "RQ");
      // h2_iRvsLayer_MIP->Fill(iL, iR, hfithisto->GetParameter(0));
      // h2_iRvsLayer_MIPerr->Fill(iL, iR, hfithisto->GetParError(0));
      tcM->Print(Form("plotsMinBias_pT4Cut/cellsPhi/h_Mip_Phi_%d_%d.png", iL, iR), "png");
      //      tcM->Print(Form("plotsMinBias_pT4Cut/cells/h_Mip_Phi_%d_%d_thinnerBins.root", iL, iR), "root");

      /*
      h_SoN_Phi[savedCout]->GetXaxis()->SetTitle("SoN");
      h_SoN_Phi[savedCout]->Draw();
      tcM->Print(Form("plotsQCD_all/cells/h_SoN_Phi_%d_%d.png", iL, iR), "png");      

      h_MipvsP_Phi[savedCout]->GetXaxis()->SetTitle("MIP");
      h_MipvsP_Phi[savedCout]->GetYaxis()->SetTitle("muon momentum (GeV)");
      h_MipvsP_Phi[savedCout]->Draw("colz");
      tcM->Print(Form("plotsQCD_all/cells/h_MipvsP_Phi_%d_%d.png", iL, iR), "png");
      */
      ++savedCout;
    }
  }

  TCanvas* tcC_Phi = new TCanvas();
  h2_iRvsLayer_MIP->GetXaxis()->SetTitle("layer");
  h2_iRvsLayer_MIP->GetYaxis()->SetTitle("iR");
  h2_iRvsLayer_MIP->GetZaxis()->SetRangeUser(0., 5.);
  h2_iRvsLayer_MIP->Draw("colz");
  tcC_Phi->Print("plotsMinBias_pT4Cut/cellsPhi/h2_iRvsLayer_MIP.png", "png");
  tcC_Phi->Print("plotsMinBias_pT4Cut/cellsPhi/h2_iRvsLayer_MIP.root", "root");

  TF1* hfit = new TF1("hfit", "gaus", 0., 2.);
  //  h_MPV_values->Rebin(2);
  hfit->SetParameters(h_MPV_values->GetEntries()/2., 1, 0.02);
  h_MPV_values->GetXaxis()->SetTitle("MPV values");
  h_MPV_values->Draw();
  h_MPV_values->Fit("hfit", "R");
  tcC_Phi->Print("plotsMinBias_pT4Cut/cellsPhi/h_MPV_values.png", "png");

  gPad->SetLogy();
  tcC_Phi->cd();
  h_Count_Phi->GetXaxis()->SetTitle("n counts in rings over #phi");
  h_Count_Phi->Draw();
  tcC_Phi->Print("plotsMinBias_pT4Cut/cellsPhi/h_Count_Phi.png", "png");

  h_minNcounts->GetXaxis()->SetTitle("minimal N. counts");
  h_minNcounts->Draw();
  tcC_Phi->Print("plotsMinBias_pT4Cut/cellsPhi/h_minNcounts.png", "png");


  std::cout << " all done now print results " << std::endl;


  std::cout << " all done - ciao " << std::endl;

  return;



}
