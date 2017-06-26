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



void drawSoN(){
  gROOT->Macro("/afs/cern.ch/user/a/amartell/public/setStyle.C");

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  
  gROOT->Reset();
  gROOT->Macro("~/public/setStyle.C");
  //  gROOT->LoadMacro("~/public/myStyle.C");
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1);

  std::cout << " inizio ci sono " << std::endl; 


  bool doAllTheFits = true;
  
  //  int iColors[16] = {kRed, kOrange+4, kOrange-3, kOrange-2, kBlue, kBlue-9, kAzure-9, kAzure+10, kCyan, kGreen+1, kCyan-2, kYellow+2}; //kGray+1};
  int iColors[7] = {kOrange-3, kRed, kMagenta, kBlue, kCyan, kGreen+1, kGray+2};
  int iStyle[4] = {20, 21, 22, 23}; //kGray+1};
  //  int iColors[6] = {kGreen+1, kBlue, kRed};
  
  int nBinsEta = 6;
  float binWidth = 0.2;
  int nBinsRad = 3; //refer to 3 thickness
  //int nBinsEta = 3;
  //float binWidth = 0.5;
  float binStart = 1.65;

  int thick = 300;

  //std::string pdgID = "130";
  std::string pdgID = "22";
  //std::string pdgID = "211";


  int nFiles = 7;
  std::vector<std::string> nameFiles;
  std::vector<float> ptValues;
  if(pdgID == "22"){
    nFiles = 3;
    nameFiles.push_back("PDG_"+pdgID+"_Pt2");
    nameFiles.push_back("PDG_"+pdgID+"_Pt5");
    nameFiles.push_back("PDG_"+pdgID+"_Pt60");
    ptValues.push_back(2);
    ptValues.push_back(5);
    ptValues.push_back(60);
  }
  else{
    nameFiles.push_back("PDG_"+pdgID+"_Pt07");
    nameFiles.push_back("PDG_"+pdgID+"_Pt1");
    nameFiles.push_back("PDG_"+pdgID+"_Pt2");
    nameFiles.push_back("PDG_"+pdgID+"_Pt5");
    nameFiles.push_back("PDG_"+pdgID+"_Pt10");
    nameFiles.push_back("PDG_"+pdgID+"_Pt30");
    nameFiles.push_back("PDG_"+pdgID+"_Pt100");

    ptValues.push_back(0.7);
    ptValues.push_back(1);
    ptValues.push_back(2);
    ptValues.push_back(5);
    ptValues.push_back(10);
    ptValues.push_back(30);
    ptValues.push_back(100);
  }

  int nOptions = 2;
  std::vector<std::string> nameOptions;
  nameOptions.push_back("CSF20LBA20");
  nameOptions.push_back("CSF20LEA20");


  TH1F* hDummy[7];
  TH1F* hDummySt[4];

  for(int iT=0; iT<nOptions; ++iT){
    hDummySt[iT] = new TH1F(Form("hDummySt%d", iT), "", 1, 0., 1);
    hDummySt[iT]->SetMarkerStyle(iStyle[iT]);
    hDummySt[iT]->SetMarkerSize(1.5);
  }
  for(int iT=0; iT<nFiles; ++iT){
    hDummy[iT] = new TH1F(Form("hDummy%d", iT), "", 1, 0., 1);
    hDummy[iT]->SetLineColor(iColors[iT]);
  }

  
  for(int ij=0; ij<nameFiles.size(); ++ij) std::cout << " name = " << nameFiles.at(ij) << std::endl;

  std::cout << " >>> ora prendo i file " << std::endl;

  //std::string optionType = "CSF20LBA20"; 
  //std::string optionType = "CSF20LEA20"; 

  
  std::string folder = "SoN";       
  
  
  /*
  std::string folderName = "plotsAllTimeTot_CFD_"+pdgID+"_3hits";
  std::string folder = "plotsTimeTot_CFD_"+pdgID+"_3hits";
  */


  TFile* inF[8][7];
  for(int iR=0; iR<nOptions; ++iR){
    for(int ij=0; ij<nFiles; ++ij){
      std::string inFileName =  "../../test/timingStudies/"+nameOptions.at(iR)+"/JOB_"+nameFiles.at(ij)+"_3fC/OutTimeHGC_RecHits_"+nameFiles.at(ij)+"_3fC.root";
      inF[iR][ij] = TFile::Open(inFileName.c_str());
    }
  }

  std::cout << " >>> fatto = presi " << std::endl;


  TH1F* SoNhisto[2][7];
  for(int iR=0; iR<nOptions; ++iR){
    for(int iT=0; iT<nFiles; ++iT){

      SoNhisto[iR][iT] = (TH1F*)(inF[iR][iT]->Get(Form("ana/SoN_for%d", thick)));

      //      SoNhisto[iR][iT]->Rebin(2);

      SoNhisto[iR][iT]->SetLineColor(iColors[iT]);   
      SoNhisto[iR][iT]->SetMarkerColor(iColors[iT]); 
      SoNhisto[iR][iT]->SetLineWidth(2);
      SoNhisto[iR][iT]->SetMarkerStyle(iStyle[iR]);
      if(iR == 1) SoNhisto[iR][iT]->SetLineStyle(2);
    }
  }


  std::cout << " ci sono ora stampo " << std::endl;


  TLegend *legTGMOpt = new TLegend(0.65,0.75,0.80,0.95,NULL,"brNDC");
  legTGMOpt->SetTextFont(42);
  legTGMOpt->SetTextSize(0.02);
  legTGMOpt->SetFillColor(kWhite);
  legTGMOpt->SetLineColor(kWhite);
  legTGMOpt->SetShadowColor(kWhite);
  for(int iF=0; iF<nOptions; ++iF){
    legTGMOpt->AddEntry(hDummySt[iF], (nameOptions.at(iF)).c_str(), "p");
  }


  TLegend *legTGM = new TLegend(0.80,0.75,0.90,0.95,NULL,"brNDC");
  legTGM->SetTextFont(42);
  legTGM->SetTextSize(0.02);
  legTGM->SetFillColor(kWhite);
  legTGM->SetLineColor(kWhite);
  legTGM->SetShadowColor(kWhite);
  for(int iF=0; iF<nFiles; ++iF){
    legTGM->AddEntry(hDummySt[iF], (nameFiles.at(iF)).c_str(), "l");
  }
  std::cout << " ci sono " << std::endl;

  /*
  TLegend *legTGM2 = new TLegend(0.65,0.75,0.80,0.95,NULL,"brNDC");
  legTGM2->SetTextFont(42);
  legTGM2->SetTextSize(0.05);
  legTGM2->SetFillColor(kWhite);
  legTGM2->SetLineColor(kWhite);
  legTGM2->SetShadowColor(kWhite);
  legTGM2->AddEntry(hDummySt[0], " 100um", "l");
  legTGM2->AddEntry(hDummySt[1], " 200um", "l");
  legTGM2->AddEntry(hDummySt[2], " 300um", "l");
  */

  std::cout << " legends ok  " << std::endl;

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TLatex tL;
  tL.SetNDC();
  tL.SetTextSize(0.05);
  tL.SetTextFont(132);

  TLatex tT;
  tT.SetNDC();
  tT.SetTextSize(0.05);
  tT.SetTextFont(132);

  TCanvas* ch_SoN[7];
  for(int iF=0; iF<nFiles; ++iF){
    ch_SoN[iF] = new TCanvas();
    ch_SoN[iF]->cd();
    gPad->SetLogy();
    SoNhisto[0][iF]->GetXaxis()->SetTitle("S/N");
    SoNhisto[0][iF]->GetYaxis()->SetTitle("events");
    SoNhisto[0][iF]->GetXaxis()->SetRangeUser(0., 100.);
    if(thick == 300)     SoNhisto[0][iF]->GetXaxis()->SetRangeUser(0., 1000.);
    if(thick == 200)     SoNhisto[0][iF]->GetXaxis()->SetRangeUser(0., 500.);
    SoNhisto[0][iF]->GetYaxis()->SetRangeUser(1., 1.e9);
    SoNhisto[0][iF]->Draw("eh");
    for(int iR=1; iR<nOptions; ++iR){
      SoNhisto[iR][iF]->Draw("eh, same");
    }
    tL.DrawLatex(0.7,0.6,nameFiles.at(iF).c_str());
    tT.DrawLatex(0.7,0.5,Form("%d um", thick));
    //    legTGM->Draw("same");
    legTGMOpt->Draw("same");
    ch_SoN[iF]->Print((folder+"/h_SoN_file"+nameFiles.at(iF)+Form("_thick%d.png", thick)).c_str(), "png");
    ch_SoN[iF]->Print((folder+"/h_SoN_file"+nameFiles.at(iF)+Form("_thick%d.pdf", thick)).c_str(), "pdf");
    ch_SoN[iF]->Print((folder+"/h_SoN_file"+nameFiles.at(iF)+Form("_thick%d.root", thick)).c_str(), "root");
  }


  return;
  
}
