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



void lookAtMIP(){
  gROOT->Macro("/afs/cern.ch/user/a/amartell/public/setStyle.C");

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  
  gROOT->Reset();
  gROOT->Macro("~/public/setStyle.C");
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1);

  std::cout << " inizio ci sono " << std::endl; 


  bool doAllTheFits = true;
  
  //  int iColors[16] = {kRed, kOrange+4, kOrange-3, kOrange-2, kBlue, kBlue-9, kAzure-9, kAzure+10, kCyan, kGreen+1, kCyan-2, kYellow+2}; //kGray+1};
  int iColors[7] = {kOrange-3, kRed, kMagenta, kBlue, kCyan, kGreen+1, kGray+2};
  int iStyle[4] = {20, 21, 22, 23}; //kGray+1};
  //  int iColors[6] = {kGreen+1, kBlue, kRed};
  
  int nBinsRad = 2;
  
  int numberOfBins = 2;
  std::vector<std::string> binNameOK;
  binNameOK.push_back("1.65-1.85");
  binNameOK.push_back("2.60-2.80");
  std::vector<std::string> binName;
  binName.push_back("1.60-2.20");
  binName.push_back("2.20-2.80");
  std::vector<float> binValue;
  binValue.push_back(1.75);
  binValue.push_back(2.70);

 

  TFile* inF = TFile::Open("../test/testMuonNewG.root");
  //  TFile* inF = TFile::Open("../test/testMinBias.root");
 
  TTree* newT = (TTree*)inF->Get("ana/newT");
  std::vector<float> *muonP = 0;
  std::vector<float> *muonEta = 0;
  std::vector<float> *muonPhi = 0;
  std::vector<float> *crossX = 0;
  std::vector<float> *crossY = 0;
  std::vector<float> *crossZ = 0;
  std::vector<float> *crossL = 0;
  std::vector<float> *crossM = 0;
  std::vector<float> *recHitX = 0;
  std::vector<float> *recHitY = 0;
  std::vector<float> *recHitZ = 0;
  std::vector<int> *recHitiPhi = 0;
  std::vector<int> *recHitiR = 0;
  std::vector<int> *recHitL = 0;
  std::vector<float> *recHitEne = 0;
  std::vector<int> *recHitMip = 0;
  std::vector<float> *recHitNoise = 0;


  newT->SetBranchAddress("muonP", &muonP);
  newT->SetBranchAddress("muonEta", &muonEta);
  newT->SetBranchAddress("muonPhi", &muonPhi);
  newT->SetBranchAddress("crossX", &crossX);
  newT->SetBranchAddress("crossY", &crossY);
  newT->SetBranchAddress("crossZ", &crossZ);
  newT->SetBranchAddress("crossL", &crossL);
  newT->SetBranchAddress("crossM", &crossM);
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
  std::map<std::pair<int, std::pair<int, int> >, std::vector<int> > enePerRh;
  std::map<std::pair<int, std::pair<int, int> >, std::vector<int> > mipPerRh;
  std::map<std::pair<int, std::pair<int, int> >, std::vector<float> > sonPerRh;
  std::map<std::pair<int, std::pair<int, int> >, std::vector<float> > momPerRh;
 
  TH2F* h2_YvsX[14];
  TH2F* h2_iRvsiPhi[14];
  TH1F* h_dX[14];
  TH1F* h_dY[14];
  TH1F* h_dR[14];
  TH2F* h_dRvsP[14];
  for(int ij=0; ij<14; ++ij){
    h2_YvsX[ij] = new TH2F(Form("h2_YvsX_L%d", ij), "", 600, -300., 300., 600, -300., 300.);
    h2_iRvsiPhi[ij] = new TH2F(Form("h2_iRvsiPhi_L%d", ij), "", 100, -50., 50., 300, 0., 300.);
    h_dX[ij] = new TH1F(Form("h_dX_L%d", ij), "", 100, -10., 10.);
    h_dY[ij] = new TH1F(Form("h_dY_L%d", ij), "", 100, -10., 10.);
    h_dR[ij] = new TH1F(Form("h_dR_L%d", ij), "", 100, 0., 10.);
    h_dRvsP[ij] = new TH2F(Form("h_dRvsP_L%d", ij), "", 100, 0., 10., 100, 0., 10.);
  }

  int firstLayer = 37;

  int nEvents = newT->GetEntries();
  std::cout << " nEvents = " << nEvents << std::endl;
  for(int ij=0; ij<nEvents; ++ij){
    newT->GetEntry(ij);

    std::map<int, std::vector<float>> rhX;
    std::map<int, std::vector<float>> rhY;
    std::map<int, std::vector<int>> rhMip;
    std::map<int, std::vector<float>> trkX;
    std::map<int, std::vector<float>> trkY;
    std::map<int, std::vector<float>> trkEta;
    std::map<int, std::vector<float>> trkP;

    //    std::cout << " muonP->size() = " << muonP->size() << std::endl;
    for(int iM=0; iM<muonP->size(); ++iM){
      trkX[crossL->at(iM)].push_back(crossX->at(iM));
      trkY[crossL->at(iM)].push_back(crossY->at(iM));
      trkEta[crossL->at(iM)].push_back(muonEta->at(iM));
      trkP[crossL->at(iM)].push_back(muonP->at(iM));
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
      float recMip = recHitMip->at(iR);
      float recEne = recHitEne->at(iR);
      
      //      std::cout << " recL - firstLayer = " << recL - firstLayer  << std::endl;

      if(trkX[recL].size() == 0) continue;
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
	  if(dR > 4.) continue;

	  h2_YvsX[recL - firstLayer]->Fill(recX, recY);
	  h2_iRvsiPhi[recL - firstLayer]->Fill(reciR, reciPhi);

	  // std::pair<float, float> coordPair = std::pair<float, float>(recX, recY);
	  // std::pair<int, std::pair<float, float> > channel = std::pair<int, std::pair<float, float> >(recL, coordPair);
	  std::pair<int, int> coordPair = std::pair<int, int>(reciR, reciPhi);                                                                                         
	  std::pair<int, std::pair<int, int> > channel = std::pair<int, std::pair<int, int> >(recL, coordPair); 	  


	  if(okChannels.find(channel) != okChannels.end()) { std::cout << " channel again " << std::endl; okChannels[channel] += 1;}
	  else {  okChannels[channel] = 1;}
	  
	  std::cout << " recMip = " << recMip << " recEne = " << recEne << " SoN = " << 1.*recMip/recHitNoise->at(iR) << std::endl;
	  enePerRh[channel].push_back(recEne);
	  mipPerRh[channel].push_back(recMip);
	  sonPerRh[channel].push_back(1.*recMip/recHitNoise->at(iR));
	  momPerRh[channel].push_back(trkP[recL][iTc]);
      }
    }//recHits

  }// events;

  std::cout << " okChannels.size = " << okChannels.size() << std::endl;

  //  return;
  for(int ij=0; ij<14; ++ij){
    TCanvas* tc = new TCanvas();
    tc->cd();
    h2_YvsX[ij]->GetXaxis()->SetTitle(Form("layer %d", ij));
    h2_YvsX[ij]->Draw("colz");
    tc->Print(Form("plots/h2_YvsX_L%d.png", ij), "png");

    h2_iRvsiPhi[ij]->GetXaxis()->SetTitle(Form("iR layer %d", ij));
    h2_iRvsiPhi[ij]->GetYaxis()->SetTitle("iPhi");
    h2_iRvsiPhi[ij]->Draw("colz");
    tc->Print(Form("plots/h2_iRvsiPhi_L%d.png", ij), "png");
    //   continue;
    h_dX[ij]->GetXaxis()->SetTitle(Form("dX layer %d", ij));
    h_dX[ij]->Draw();
    tc->Print(Form("plots/h_dX_L%d.png", ij), "png");

    h_dY[ij]->GetXaxis()->SetTitle(Form("dY layer %d", ij));
    h_dY[ij]->Draw();
    tc->Print(Form("plots/h_dY_L%d.png", ij), "png");

    h_dR[ij]->GetXaxis()->SetTitle(Form("dR layer %d", ij));
    h_dR[ij]->Draw();
    tc->Print(Form("plots/h_dR_L%d.png", ij), "png");

    h_dRvsP[ij]->GetXaxis()->SetTitle(Form("dR layer %d", ij));
    h_dRvsP[ij]->GetYaxis()->SetTitle("muon P");
    h_dRvsP[ij]->Draw("colz");
    tc->Print(Form("plots/h_dRvsdP_L%d.png", ij), "png");
  }

  //  return;

  TH1F* h_Ene[100];
  TH1F* h_Mip[100];
  TH1F* h_SoN[100];
  TH2F* h_MipvsP[100];

  int iC = -1;
  for(auto ic : okChannels){
    ++iC;
    if(iC == 100) break;
    int iL = ic.first.first;
    int iR = ic.first.second.first;
    int iPhi = ic.first.second.second;
    h_Ene[iC] = new TH1F(Form("h_Ene_%d_%d_%d", iL, iR, iPhi), "", 100, 0., 10.);
    h_Mip[iC] = new TH1F(Form("h_Mip_%d_%d_%d", iL, iR, iPhi), "", 100, 0., 1.e3);
    h_SoN[iC] = new TH1F(Form("h_SoN_%d_%d_%d", iL, iR, iPhi), "", 100, 0., 1.e3);
    h_MipvsP[iC] = new TH2F(Form("h_MipvsP_%d_%d_%d", iL, iR, iPhi), "", 100, 0., 1.e3, 100, 0., 100.);

    int iVal = ic.second;
    if(iVal > 1) std::cout << " iL = " << iL << " iR = " << iR << " iPhi = " << iPhi << std::endl;
    for(int ij=0; ij<iVal; ++ij){
      h_Ene[iC]->Fill(enePerRh[ic.first][ij]);
      h_Mip[iC]->Fill(mipPerRh[ic.first][ij]);
      h_SoN[iC]->Fill(sonPerRh[ic.first][ij]);
      h_MipvsP[iC]->Fill(mipPerRh[ic.first][ij], momPerRh[ic.first][ij]);
    }

    TCanvas* tcM = new TCanvas();
    tcM->cd();
    h_Ene[iC]->GetXaxis()->SetTitle("recHit energy (GeV)");
    h_Ene[iC]->Draw();
    tcM->Print(Form("plots/h_Ene_%d_%d_%d.png", iL, iR, iPhi), "png");

    h_Mip[iC]->GetXaxis()->SetTitle("MIP");
    h_Mip[iC]->Draw();
    tcM->Print(Form("plots/h_Mip_%d_%d_%d.png", iL, iR, iPhi), "png");

    h_SoN[iC]->GetXaxis()->SetTitle("SoN");
    h_SoN[iC]->Draw();
    tcM->Print(Form("plots/h_SoN_%d_%d_%d.png", iL, iR, iPhi), "png");

    h_MipvsP[iC]->GetXaxis()->SetTitle("MIP");
    h_MipvsP[iC]->GetYaxis()->SetTitle("muon momentum (GeV)");
    h_MipvsP[iC]->Draw("colz");
    tcM->Print(Form("plots/h_MipvsP_%d_%d_%d.png", iL, iR, iPhi), "png");
  }

  std::cout << " all done now print results " << std::endl;


  std::cout << " all done - ciao " << std::endl;

  return;



}
