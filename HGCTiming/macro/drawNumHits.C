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



void drawNumHits(){
  gROOT->Macro("/afs/cern.ch/user/a/amartell/public/setStyle.C");

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  //  gROOT->Reset();
  //  gROOT->Macro("~/public/setStyle.C");
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
  int nBinsRad = 4;
  //int nBinsEta = 3;
  //float binWidth = 0.5;
  float binStart = 1.65;
  
  //  std::string pdgID = "130";
  std::string pdgID = "22";
  //std::string pdgID = "211";


  int nFiles = 7;
  std::vector<std::string> nameFiles;
  std::vector<int> ptValues;
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
  TH1F* hDummy[7];
  TH1F* hDummySt[4];

  
  for(int ij=0; ij<nameFiles.size(); ++ij) std::cout << " name = " << nameFiles.at(ij) << std::endl;

  std::cout << " >>> ora prendo i file " << std::endl;


  //  std::string optionType = "CSF20LBA50";             
  std::string optionType = "CSF20LBA50keep68";             
  //std::string optionType = "CSF20LEA50";             

  
  std::string folderName = optionType+"/plotsAllTimeTot_CFD_"+pdgID+"_3hits";  
  std::string folder = optionType+"/numHits/";       
  
  
  /*
  std::string folderName = "plotsAllTimeTot_CFD_"+pdgID+"_3hits";
  std::string folder = "plotsTimeTot_CFD_"+pdgID+"_3hits";
  */

  TFile* inF[7];
  for(int ij=0; ij<nFiles; ++ij){
    inF[ij] = TFile::Open(("../test/"+optionType+"/"+"JOB_"+nameFiles.at(ij)+"_3fC/"+"/OutTimeHGC_RecHits_"+nameFiles.at(ij)+"_3fC.root").c_str());
  }
  
  std::cout << " >>> fatto = presi " << std::endl;
  TGraphErrors* nHitsWithTime[7][2];
  TH1F* histoPluto[2];
  for(int iF=0; iF<nFiles; ++iF){
    for(int iRad=0; iRad<nBinsRad; ++iRad){
      if(iRad < 2){
	if(iF == 0){
	  histoPluto[iRad] = new TH1F(Form("histoPluto%d", iRad), "", 1, 0., 1);
	  histoPluto[iRad]->SetMarkerStyle(iStyle[iRad]);
	}
	nHitsWithTime[iF][iRad] = new TGraphErrors();
	nHitsWithTime[iF][iRad]->SetPoint(0, 0., 0.);

	for(int ieta=0; ieta<nBinsEta; ++ieta){
	  //	  std::cout << " >>> name = " << Form("ana/hNumberHitsWithTime_Eta%.2f-%.2f_dRadius%d", (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad) << std::endl;
	  TH1F* dummy = (TH1F*)(inF[iF]->Get(Form("ana/hNumberHitsWithTime_Eta_dRadius_Eta%.2f-%.2f_dRadius%d", (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad)));
	  nHitsWithTime[iF][iRad]->SetPoint(ieta+1, (binStart+ieta*binWidth) + binWidth*0.5, dummy->GetMean());
	  //	  nHitsWithTime[iF][iRad]->SetPointError(ieta+1, 0., dummy->GetRMS());
	}
	nHitsWithTime[iF][iRad]->SetPoint(nBinsEta+1, 3.5, 1.e4);
	nHitsWithTime[iF][iRad]->SetMarkerColor(iColors[iF]);
	nHitsWithTime[iF][iRad]->SetLineColor(iColors[iF]);
	nHitsWithTime[iF][iRad]->SetLineWidth(2);
	nHitsWithTime[iF][iRad]->SetMarkerStyle(iStyle[iRad]);
      }
    }
  }

  TGraphErrors* nHitsWithTimeEtaL[2];
  TH1F* time[7][4];
  TH1F* timeCut[7][4];
  TH1F* energyRadEta_bin1[7][4];
  TH1F* energyRadEta_bin5[7][4];
  TH1F* numHitsWithTime[7][4];
  TH1F* numHitsWithTimeEff[7][4];
  TH1F* numHitsWithTimeDistr[7][4];
  TH1F* testHisto;
  for(int iF=0; iF<nFiles; ++iF){
    hDummy[iF] = new TH1F(Form("hDummy%d", iF), "", 1, 0., 1);
    hDummy[iF]->SetLineColor(iColors[iF]);

    for(int iR=0; iR<nBinsRad; ++iR){
      if(iF == 0){
	hDummySt[iR] = new TH1F(Form("hDummySt%d", iR), "", 1, 0., 1);
	hDummySt[iR]->SetMarkerStyle(iStyle[iR]);
	hDummySt[iR]->SetLineColor(iColors[iR+2]);
      }
      time[iF][iR] = (TH1F*)(inF[iF]->Get(Form("ana/hTime_Eta_dRadius_Eta1.65-1.85_dRadius%d", iR)));
      timeCut[iF][iR] = (TH1F*)(inF[iF]->Get(Form("ana/hTimeCut_Eta_dRadius_Eta1.65-1.85_dRadius%d", iR)));

      energyRadEta_bin1[iF][iR] = (TH1F*)(inF[iF]->Get(Form("ana/hEnergy_Eta_dRadius_Eta1.65-1.85_dRadius%d",iR)));
      energyRadEta_bin5[iF][iR] = (TH1F*)(inF[iF]->Get(Form("ana/hEnergy_Eta_dRadius_Eta2.65-2.85_dRadius%d",iR)));

      testHisto = (TH1F*)(inF[iF]->Get(Form("ana/hNumberHitsWithTime_Eta_dRadius_Eta1.65-1.85_dRadius%d", iR)));
      numHitsWithTimeDistr[iF][iR] = (TH1F*)testHisto->Clone(Form("numHitsWithTimeDistr_iF%d_dRadius%d", iF, iR));
      numHitsWithTime[iF][iR] = (TH1F*)testHisto->Clone(Form("numHitsWithTime_iF%d_dRadius%d", iF, iR));
      numHitsWithTime[iF][iR]->Reset();
      numHitsWithTimeEff[iF][iR] = (TH1F*)testHisto->Clone(Form("numHitsWithTimeEff_iF%d_dRadius%d", iF, iR));
      numHitsWithTimeEff[iF][iR]->Reset();


      if(iF == 0){
	nHitsWithTimeEtaL[iR] = new TGraphErrors();
	nHitsWithTimeEtaL[iR]->SetPoint(0, 0., 0.);
      }
      nHitsWithTimeEtaL[iR]->SetPoint(iF+1, ptValues.at(iF), testHisto->GetMean());
      //      nHitsWithTimeEtaL[iR]->SetPointError(iF+1, 0., testHisto->GetRMS());
      if(iR == 0){
	nHitsWithTimeEtaL[iR]->SetMarkerColor(kBlue);
	nHitsWithTimeEtaL[iR]->SetLineColor(kBlue);
	nHitsWithTimeEtaL[iR]->SetLineWidth(2);
	nHitsWithTimeEtaL[iR]->SetMarkerStyle(20);
      }
      if(iR == 1){
	nHitsWithTimeEtaL[iR]->SetMarkerColor(kRed);
	nHitsWithTimeEtaL[iR]->SetLineColor(kRed);
	nHitsWithTimeEtaL[iR]->SetLineWidth(2);
	nHitsWithTimeEtaL[iR]->SetMarkerStyle(21);
      }


      int totEvt = 0;
      int numEvt = 0;      
      for(int ij=testHisto->GetNbinsX(); ij>0; ij--){
	numEvt += testHisto->GetBinContent(ij);
	numHitsWithTime[iF][iR]->SetBinContent(ij, numEvt);
      }
      int locEvt = 0;
      int locEvt0 = numHitsWithTime[iF][iR]->GetBinContent(1);
      for(int ij=testHisto->GetNbinsX(); ij>0; ij--){
        locEvt = numHitsWithTime[iF][iR]->GetBinContent(ij);
        numHitsWithTimeEff[iF][iR]->SetBinContent(ij, 1.*locEvt/locEvt0);
      }
      //tg[iT][iR]->SetLineColor(iColors[3*iT+iR]);
      //tg[iT][iR]->SetMarkerColor(iColors[3*iT+iR]);

      //      numHitsWithTime[iF][iR]->Rebin(2);

      numHitsWithTimeDistr[iF][iR]->SetLineColor(iColors[iR+2]);   
      numHitsWithTimeDistr[iF][iR]->SetMarkerColor(iColors[iR+2]); 
      numHitsWithTimeDistr[iF][iR]->SetLineWidth(2);
      numHitsWithTimeDistr[iF][iR]->SetMarkerStyle(iStyle[iR]);

      numHitsWithTime[iF][iR]->SetLineColor(iColors[iR+2]);   
      numHitsWithTime[iF][iR]->SetMarkerColor(iColors[iR+2]); 
      numHitsWithTime[iF][iR]->SetLineWidth(2);
      numHitsWithTime[iF][iR]->SetMarkerStyle(iStyle[iR]);

      numHitsWithTimeEff[iF][iR]->SetLineColor(iColors[iR+2]);   
      numHitsWithTimeEff[iF][iR]->SetMarkerColor(iColors[iR+2]); 
      numHitsWithTimeEff[iF][iR]->SetLineWidth(2);
      numHitsWithTimeEff[iF][iR]->SetMarkerStyle(iStyle[iR]);
      
      time[iF][iR]->SetLineColor(iColors[iR+2]);
      time[iF][iR]->SetMarkerColor(iColors[iR+2]);
      time[iF][iR]->SetLineWidth(2);
      time[iF][iR]->SetMarkerStyle(iStyle[iR]);


      timeCut[iF][iR]->SetLineColor(iColors[iR+2]);
      timeCut[iF][iR]->SetMarkerColor(iColors[iR+2]);
      timeCut[iF][iR]->SetLineWidth(2);
      timeCut[iF][iR]->SetLineStyle(2);
      timeCut[iF][iR]->SetMarkerStyle(iStyle[iR]);

      energyRadEta_bin1[iF][iR]->SetLineColor(iColors[iR+2]);
      energyRadEta_bin1[iF][iR]->SetMarkerColor(iColors[iR+2]);
      energyRadEta_bin1[iF][iR]->SetLineWidth(2);
      energyRadEta_bin1[iF][iR]->SetMarkerStyle(iStyle[iR]);


      energyRadEta_bin5[iF][iR]->SetLineColor(iColors[iR+2]);
      energyRadEta_bin5[iF][iR]->SetMarkerColor(iColors[iR+2]);
      energyRadEta_bin5[iF][iR]->SetLineWidth(2);
      energyRadEta_bin5[iF][iR]->SetMarkerStyle(iStyle[iR]);

    }
  }
  nHitsWithTimeEtaL[0]->SetPoint(nFiles+2, 500, 1000);
  nHitsWithTimeEtaL[1]->SetPoint(nFiles+2, 500, 1000);


  std::cout << " ci sono ora stampo " << std::endl;


  TLegend *legTGM = new TLegend(0.80,0.70,0.90,0.90,NULL,"brNDC");
  legTGM->SetTextFont(42);
  legTGM->SetTextSize(0.02);
  legTGM->SetFillColor(kWhite);
  legTGM->SetLineColor(kWhite);
  legTGM->SetShadowColor(kWhite);
  for(int iF=0; iF<nFiles; ++iF){
    legTGM->AddEntry(hDummySt[iF], (nameFiles.at(iF)).c_str(), "l");
  }
  std::cout << " ci sono " << std::endl;

  TLegend *legTGM2 = new TLegend(0.75,0.70,0.90,0.90,NULL,"brNDC");
  legTGM2->SetTextFont(42);
  legTGM2->SetTextSize(0.04);
  legTGM2->SetFillColor(kWhite);
  legTGM2->SetLineColor(kWhite);
  legTGM2->SetShadowColor(kWhite);
  legTGM2->AddEntry(hDummySt[0], " #rho < 2cm", "l");
  legTGM2->AddEntry(hDummySt[1], " #rho < 5cm", "l");
  legTGM2->AddEntry(hDummySt[2], " #rho < 10cm", "l");
  legTGM2->AddEntry(hDummySt[3], " all", "l");


  TLegend *legTGMS = new TLegend(0.80,0.75,0.90,0.95,NULL,"brNDC");
  legTGMS->SetTextFont(42);
  legTGMS->SetTextSize(0.03);
  legTGMS->SetFillColor(kWhite);
  legTGMS->SetLineColor(kWhite);
  legTGMS->SetShadowColor(kWhite);
  for(int iF=0; iF<nFiles; ++iF){
    legTGMS->AddEntry(nHitsWithTime[iF][0], (nameFiles.at(iF)).c_str(), "l");
  }

  TLegend *legTGM2S = new TLegend(0.45,0.80,0.55,0.95,NULL,"brNDC");
  legTGM2S->SetTextFont(42);
  legTGM2S->SetTextSize(0.03);
  legTGM2S->SetFillColor(kWhite);
  legTGM2S->SetLineColor(kWhite);
  legTGM2S->SetShadowColor(kWhite);
  legTGM2S->AddEntry(histoPluto[0], " #rho < 2cm", "p");
  legTGM2S->AddEntry(histoPluto[1], " #rho < 5cm", "p");


  TLegend *legTGMSvsE = new TLegend(0.50,0.80,0.80,0.90,NULL,"brNDC");
  legTGMSvsE->SetTextFont(42);
  legTGMSvsE->SetTextSize(0.04);
  legTGMSvsE->SetFillColor(kWhite);
  legTGMSvsE->SetLineColor(kWhite);
  legTGMSvsE->SetShadowColor(kWhite);
  legTGMSvsE->AddEntry(nHitsWithTimeEtaL[0], "#rho < 2cm", "p");
  legTGMSvsE->AddEntry(nHitsWithTimeEtaL[1], "#rho < 5cm", "p");
  

  std::cout << " legends ok  " << std::endl;

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TLatex tL;
  tL.SetNDC();
  tL.SetTextSize(0.05);
  tL.SetTextFont(132);


  TCanvas* ch_NumberHitsVsPt;
  ch_NumberHitsVsPt = new TCanvas();
  ch_NumberHitsVsPt->cd();
  gPad->SetLogx();
  gPad->SetLogy();
  nHitsWithTimeEtaL[0]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  nHitsWithTimeEtaL[0]->GetYaxis()->SetTitle("<n. hits>");
  //  nHitsWithTimeEtaL[0]->GetXaxis()->SetRangeUser(1.5, 3.);
  nHitsWithTimeEtaL[0]->GetXaxis()->SetRangeUser(0.1, 110.);
  nHitsWithTimeEtaL[0]->GetYaxis()->SetRangeUser(1, 1000.);
  nHitsWithTimeEtaL[0]->Draw("ap");
  nHitsWithTimeEtaL[1]->Draw("same, p");

  if(pdgID == "130") tL.DrawLatex(0.20,0.86, "|#eta| #approx 1.75   K^{0}_{L}");
  if(pdgID == "22")  tL.DrawLatex(0.20,0.86, "|#eta| #approx 1.75   #gamma");

  legTGMSvsE->Draw("same");
  ch_NumberHitsVsPt->Print((folder+"/h_numHitsVsPt_file"+pdgID+".png").c_str(), "png");
  ch_NumberHitsVsPt->Print((folder+"/h_numHitsVsPt_file"+pdgID+".pdf").c_str(), "pdf");
  ch_NumberHitsVsPt->Print((folder+"/h_numHitsVsPt_file"+pdgID+".root").c_str(), "root");


  TCanvas* ch_NumberHitsVseta;
  ch_NumberHitsVseta = new TCanvas();
  ch_NumberHitsVseta->cd();
  gPad->SetLogy();
  nHitsWithTime[0][0]->GetXaxis()->SetTitle("|#eta|");
  nHitsWithTime[0][0]->GetYaxis()->SetTitle("<n. hits>");
  nHitsWithTime[0][0]->GetXaxis()->SetRangeUser(1.5, 3.);
  nHitsWithTime[0][0]->GetYaxis()->SetRangeUser(1., 5000.);
  nHitsWithTime[0][0]->Draw("ap");
  nHitsWithTime[0][1]->Draw("same, p");
  for(int iF=1; iF<nFiles; ++iF){
    for(int iR=0; iR<2; ++iR){
      nHitsWithTime[iF][iR]->Draw("same, p");
    }
  }

  //0.4,0.86                                                                                                                                                                                                                                            
  if(pdgID == "130") tL.DrawLatex(0.20,0.86, "K^{0}_{L}" );
  if(pdgID == "22")  tL.DrawLatex(0.20,0.86, "#gamma   ");

  legTGMS->Draw("same");
  legTGM2S->Draw("same");
  ch_NumberHitsVseta->Print((folder+"/h_numHitsVseta_file"+pdgID+".png").c_str(), "png");
  ch_NumberHitsVseta->Print((folder+"/h_numHitsVseta_file"+pdgID+".pdf").c_str(), "pdf");
  ch_NumberHitsVseta->Print((folder+"/h_numHitsVseta_file"+pdgID+".root").c_str(), "root");
  ch_NumberHitsVseta->SaveAs((folder+"/h_numHitsVseta_file"+pdgID+".C").c_str(), "C");

  bool doAll = false;
  if(doAll){
  TCanvas* ch_NumberHits[7];
  for(int iF=0; iF<nFiles; ++iF){
    ch_NumberHits[iF] = new TCanvas();
    ch_NumberHits[iF]->cd();
    numHitsWithTime[iF][0]->GetXaxis()->SetTitle("#geq n. hits with time (1.65 < |#eta| < 1.85)");
    numHitsWithTime[iF][0]->GetYaxis()->SetTitle("events");
    float Ymax = 0.;
    getXmax(numHitsWithTime[iF][3], Ymax);
    numHitsWithTime[iF][0]->GetXaxis()->SetRangeUser(0., 50.);
    numHitsWithTime[iF][0]->GetYaxis()->SetRangeUser(0., Ymax * 4.);
    numHitsWithTime[iF][0]->Draw("");

    for(int iR=1; iR<nBinsRad; ++iR){
      numHitsWithTime[iF][iR]->Draw("same");
    }

    tL.DrawLatex(0.7,0.6,nameFiles.at(iF).c_str());
    //    legTGM->Draw("same");
    legTGM2->Draw("same");
    ch_NumberHits[iF]->Print((folder+"/h_numHits_file"+nameFiles.at(iF)+".png").c_str(), "png");
    ch_NumberHits[iF]->Print((folder+"/h_numHits_file"+nameFiles.at(iF)+".pdf").c_str(), "pdf");
    ch_NumberHits[iF]->Print((folder+"/h_numHits_file"+nameFiles.at(iF)+".root").c_str(), "root");
  }

  TCanvas* ch_NumberHitsEff[7];
  for(int iF=0; iF<nFiles; ++iF){
    ch_NumberHitsEff[iF] = new TCanvas();
    ch_NumberHitsEff[iF]->cd();
    numHitsWithTimeEff[iF][0]->GetXaxis()->SetTitle("#geq n. hits with time (1.65 < |#eta| < 1.85)");
    numHitsWithTimeEff[iF][0]->GetYaxis()->SetTitle("fraction events");
    numHitsWithTimeEff[iF][0]->GetXaxis()->SetRangeUser(0., 50.);
    numHitsWithTimeEff[iF][0]->GetYaxis()->SetRangeUser(0., 1.01);
    numHitsWithTimeEff[iF][0]->Draw("");

    for(int iR=1; iR<nBinsRad; ++iR){
      numHitsWithTimeEff[iF][iR]->Draw("same");
    }
    tL.DrawLatex(0.7,0.6,nameFiles.at(iF).c_str());
    //    legTGM->Draw("same");
    legTGM2->Draw("same");
    ch_NumberHitsEff[iF]->Print((folder+"/h_numHitsEff_file"+nameFiles.at(iF)+".png").c_str(), "png");
    ch_NumberHitsEff[iF]->Print((folder+"/h_numHitsEff_file"+nameFiles.at(iF)+".pdf").c_str(), "pdf");
    ch_NumberHitsEff[iF]->Print((folder+"/h_numHitsEff_file"+nameFiles.at(iF)+".root").c_str(), "root");
  }


  TCanvas* ch_NumberHitsDistr[7];
  for(int iF=0; iF<nFiles; ++iF){
    ch_NumberHitsDistr[iF] = new TCanvas();
    ch_NumberHitsDistr[iF]->cd();
    numHitsWithTimeDistr[iF][0]->GetXaxis()->SetTitle("#n. hits with time (1.65 < |#eta| < 1.85)");
    numHitsWithTimeDistr[iF][0]->GetYaxis()->SetTitle(" events ");
    numHitsWithTimeDistr[iF][0]->GetXaxis()->SetRangeUser(0., 200.);
    if(iF == nFiles-1) numHitsWithTimeDistr[iF][0]->GetXaxis()->SetRangeUser(0., 1000.);
    numHitsWithTimeDistr[iF][0]->Draw("");

    for(int iR=1; iR<nBinsRad; ++iR){
      numHitsWithTimeDistr[iF][iR]->Draw("same");
    }

    tL.DrawLatex(0.7,0.6,nameFiles.at(iF).c_str());
    //    legTGM->Draw("same");
    legTGM2->Draw("same");
    ch_NumberHitsDistr[iF]->Print((folder+"/h_numHitsDistr_file"+nameFiles.at(iF)+".png").c_str(), "png");
    ch_NumberHitsDistr[iF]->Print((folder+"/h_numHitsDistr_file"+nameFiles.at(iF)+".pdf").c_str(), "pdf");
    ch_NumberHitsDistr[iF]->Print((folder+"/h_numHitsDistr_file"+nameFiles.at(iF)+".root").c_str(), "root");
  }


  TCanvas* ch_Time[7];
  for(int iF=0; iF<nFiles; ++iF){
    ch_Time[iF] = new TCanvas();
    ch_Time[iF]->cd();
    if(iF != 2)    gPad->SetLogy();
    if(iF == 2) time[iF][0]->GetXaxis()->SetTitle("time of the hits (ns)");
    else time[iF][0]->GetXaxis()->SetTitle("recHit time (1.65 < |#eta| < 1.85)");
    time[iF][0]->GetXaxis()->SetRangeUser(-0.5, 2.);
    if(iF == 2) time[iF][0]->GetYaxis()->SetRangeUser(0., 4000.);
    //    time[iF][0]->GetXaxis()->SetRangeUser(-1.5, 4.5);
    time[iF][0]->GetYaxis()->SetTitle(" events");
    //    numHitsWithTimeDistr[iF][0]->GetXaxis()->SetRangeUser(0., 20.);
    //    numHitsWithTimeDistr[iF][0]->GetYaxis()->SetRangeUser(0., 1.01);
    time[iF][0]->Draw("");
    timeCut[iF][0]->Draw("same");

    for(int iR=1; iR<nBinsRad; ++iR){
      time[iF][iR]->Draw("same");
      timeCut[iF][iR]->Draw("same");
    }
    if(iF == 2 && pdgID == "130") tL.DrawLatex(0.4,0.86, "|#eta| #approx 1.75   K^{0}_{L}  p_{T} = 2 GeV/c");
    else if(iF == 2 && pdgID == "22") tL.DrawLatex(0.6,0.86, "|#eta| #approx 1.75   #gamma   p_{T} = 2 GeV/c");
    else tL.DrawLatex(0.7,0.6,nameFiles.at(iF).c_str());
    //    legTGM->Draw("same");
    legTGM2->Draw("same");
    ch_Time[iF]->Print((folder+"/h_recTime_file"+nameFiles.at(iF)+".png").c_str(), "png");
    ch_Time[iF]->Print((folder+"/h_recTime_file"+nameFiles.at(iF)+".pdf").c_str(), "pdf");
    ch_Time[iF]->Print((folder+"/h_recTime_file"+nameFiles.at(iF)+".root").c_str(), "root");
  }
  }

  /////////
  TCanvas* ch_TimeDistr[7];
  for(int iF=0; iF<nFiles; ++iF){
    ch_TimeDistr[iF] = new TCanvas();
    ch_TimeDistr[iF]->cd();
    time[iF][0]->GetXaxis()->SetTitle("time of the hits (ns)");
    time[iF][0]->GetXaxis()->SetRangeUser(-0.5, 2.);

    // For the axis titles: 
    gStyle->SetTitleColor(1, "XYZ");
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetTitleSize(0.06, "XYZ");
    gStyle->SetTitleXOffset(0.9);
    // gStyle->SetTitleYOffset(1.5);
    // gStyle->SetTitleYOffset(1.5); 
    gStyle->SetTitleYOffset(1.4);
    gStyle->SetTitleYOffset(1.4);
    //  gStyle->SetTitleOffset(1.3,"Y"); 
    
    // For the axis labels:   
    gStyle->SetLabelColor(1, "XYZ");
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetLabelOffset(0.007, "XYZ");
    gStyle->SetLabelSize(0.05, "XYZ");

    float Ymax = 0.;
    getXmax(time[iF][0], Ymax);
    time[iF][0]->GetYaxis()->SetRangeUser(0., Ymax*2.5);
    if(iF == 2) time[iF][0]->GetYaxis()->SetRangeUser(0., 4000.);
    else time[iF][0]->GetYaxis()->SetRangeUser(0., 30000.);
    time[iF][0]->GetYaxis()->SetTitle("");

    time[iF][0]->Draw("");
    for(int iR=1; iR<nBinsRad; ++iR){
      time[iF][iR]->Draw("same");
    }
    //0.4,0.86
    if(pdgID == "130") tL.DrawLatex(0.25,0.86, Form("|#eta| #approx 1.75   K^{0}_{L}   p_{T} = %d GeV/c", ptValues.at(iF)) );
    if(pdgID == "22")  tL.DrawLatex(0.25,0.86, Form("|#eta| #approx 1.75   #gamma      p_{T} = %d GeV/c", ptValues.at(iF)));


    legTGM2->Draw("same");
    ch_TimeDistr[iF]->Print((folder+"/h_recTimeDistr_file"+nameFiles.at(iF)+".png").c_str(), "png");
    ch_TimeDistr[iF]->Print((folder+"/h_recTimeDistr_file"+nameFiles.at(iF)+".pdf").c_str(), "pdf");
    ch_TimeDistr[iF]->Print((folder+"/h_recTimeDistr_file"+nameFiles.at(iF)+".root").c_str(), "root");

    //    time[iF][0]->GetYaxis()->SetRangeUser(0.1, Ymax*1.2);
    // gPad->SetLogy();
    // ch_TimeDistr[iF]->Print((folder+"/h_recTimeDistr_file"+nameFiles.at(iF)+"_log.png").c_str(), "png");
    // ch_TimeDistr[iF]->Print((folder+"/h_recTimeDistr_file"+nameFiles.at(iF)+"_log.pdf").c_str(), "pdf");
    // ch_TimeDistr[iF]->Print((folder+"/h_recTimeDistr_file"+nameFiles.at(iF)+"_log.root").c_str(), "root");
  }

  return;

  TCanvas* ch_Energy_bin1[7];
  for(int iF=0; iF<nFiles; ++iF){
    ch_Energy_bin1[iF] = new TCanvas();
    ch_Energy_bin1[iF]->cd();
    //    gPad->SetLogy();

    energyRadEta_bin1[iF][0]->GetXaxis()->SetTitle("recHits Sum(pT) (1.65 < |#eta| < 1.85)");
    energyRadEta_bin1[iF][0]->GetXaxis()->SetRangeUser(0., ptValues.at(iF)*2.);
    energyRadEta_bin1[iF][0]->GetYaxis()->SetTitle(" events");
    energyRadEta_bin1[iF][0]->Draw("");

    for(int iR=1; iR<nBinsRad; ++iR){
      energyRadEta_bin1[iF][iR]->Draw("same");
      energyRadEta_bin1[iF][iR]->Draw("same");
    }

    tL.DrawLatex(0.7,0.6,nameFiles.at(iF).c_str());
    legTGM2->Draw("same");
    ch_Energy_bin1[iF]->Print((folder+"/h_energy_bin1_file"+nameFiles.at(iF)+".png").c_str(), "png");
    ch_Energy_bin1[iF]->Print((folder+"/h_energy_bin1_file"+nameFiles.at(iF)+".pdf").c_str(), "pdf");
    ch_Energy_bin1[iF]->Print((folder+"/h_energy_bin1_file"+nameFiles.at(iF)+".root").c_str(), "root");
  }

  TCanvas* ch_Energy_bin5[7];
  for(int iF=0; iF<nFiles; ++iF){
    ch_Energy_bin5[iF] = new TCanvas();
    ch_Energy_bin5[iF]->cd();
    //    gPad->SetLogy();

    energyRadEta_bin5[iF][0]->GetXaxis()->SetTitle("recHits Sum(pT) (2.65 < |#eta| < 2.85)");
    energyRadEta_bin5[iF][0]->GetXaxis()->SetRangeUser(0., ptValues.at(iF)*2.);
    energyRadEta_bin5[iF][0]->GetYaxis()->SetTitle(" events");
    energyRadEta_bin5[iF][0]->Draw("");

    for(int iR=1; iR<nBinsRad; ++iR){
      energyRadEta_bin5[iF][iR]->Draw("same");
      energyRadEta_bin5[iF][iR]->Draw("same");
    }

    tL.DrawLatex(0.7,0.6,nameFiles.at(iF).c_str());
    legTGM2->Draw("same");
    ch_Energy_bin5[iF]->Print((folder+"/h_energy_bin5_file"+nameFiles.at(iF)+".png").c_str(), "png");
    ch_Energy_bin5[iF]->Print((folder+"/h_energy_bin5_file"+nameFiles.at(iF)+".pdf").c_str(), "pdf");
    ch_Energy_bin5[iF]->Print((folder+"/h_energy_bin5_file"+nameFiles.at(iF)+".root").c_str(), "root");
  }


  return;

}
