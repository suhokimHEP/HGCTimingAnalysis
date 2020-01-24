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



void doPlots_fitTiming(std::string pdgID, std::string optionType = "in11X_200PU_allEta"){
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

  //  std::string pdgID = "130";
  //std::string pdgID = "22";
  

  int nFiles = 3;
  std::vector<std::string> nameFiles;
  std::vector<float> ptValues;
  nameFiles.push_back("PDG_"+pdgID+"_pt5");
  nameFiles.push_back("PDG_"+pdgID+"_pt5");
  nameFiles.push_back("PDG_"+pdgID+"_pt5");
  ptValues.push_back(5);
  ptValues.push_back(5);
  ptValues.push_back(5);
  TH1F* hDummy[7];
  TH1F* hDummySt[4];

  
  for(int ij=0; ij<nameFiles.size(); ++ij) std::cout << " name = " << nameFiles.at(ij) << std::endl;

  TGraphErrors* tg[7][4];
  TGraphErrors* tgM[7][4];
  TGraphErrors* tgFrEvt[7][4];
  for(int iT=0; iT<nFiles; ++iT){
    hDummy[iT] = new TH1F(Form("hDummy%d", iT), "", 1, 0., 1);
    hDummy[iT]->SetLineColor(iColors[iT]);

    for(int iR=0; iR<nBinsRad; ++iR){
      if(iT == 0){
	hDummySt[iR] = new TH1F(Form("hDummySt%d", iR), "", 1, 0., 1);
	hDummySt[iR]->SetMarkerStyle(iStyle[iR]);
      }

      tg[iT][iR] = new TGraphErrors();
      tg[iT][iR]->SetName(Form("resolution_%d_rad%d", iT, iR));
      tg[iT][iR]->SetPoint(0, -1, -1);
      
      tg[iT][iR]->SetLineColor(iColors[iT]);   
      tg[iT][iR]->SetMarkerColor(iColors[iT]); 
      tg[iT][iR]->SetMarkerSize(1.2); 
      tg[iT][iR]->SetLineWidth(2);
      tg[iT][iR]->SetMarkerStyle(iStyle[iR]);
      ////////      
      tgM[iT][iR] = new TGraphErrors();
      tgM[iT][iR]->SetName(Form("mean_%d_rad%d", iT, iR));
      tgM[iT][iR]->SetPoint(0, -1, -1);
      
      tgM[iT][iR]->SetLineColor(iColors[iT]);
      tgM[iT][iR]->SetLineWidth(2);
      tgM[iT][iR]->SetMarkerSize(1.2); 
      tgM[iT][iR]->SetMarkerColor(iColors[iT]);
      tgM[iT][iR]->SetMarkerStyle(iStyle[iR]);

      //////////
           
      tgFrEvt[iT][iR] = new TGraphErrors();
      tgFrEvt[iT][iR]->SetName(Form("fractionEvt_%d_rad%d", iT, iR));
      tgFrEvt[iT][iR]->SetPoint(0, -1, -1);
      
      tgFrEvt[iT][iR]->SetLineColor(iColors[iT]);
      tgFrEvt[iT][iR]->SetLineWidth(2);
      tgFrEvt[iT][iR]->SetMarkerSize(1.2); 
      tgFrEvt[iT][iR]->SetMarkerColor(iColors[iT]);
      tgFrEvt[iT][iR]->SetMarkerStyle(iStyle[iR]);
  
    }
  }  

  std::cout << " >>> now load files " << std::endl;
  
  //  std::string optionType = "TDRset_0PU_allEta";

  std::string singleFitsFolderName = optionType+"/plotsAllHistos_"+pdgID;  
  std::string plotFolder = optionType+"/plotsTime_"+pdgID;       
  

  TFile* inF[3];
  for(int ij=0; ij<nFiles; ++ij){
    inF[ij] = TFile::Open(("../test/"+optionType+"/"+"/OutTimeHGC_RecHits_"+nameFiles.at(ij)+"_200PU_AE.root").c_str());
  }
  
  std::cout << " >>> files loaded " << std::endl;

  std::ofstream outFileLong[7];
  for(int ij=0; ij<nFiles; ++ij){
    outFileLong[ij].open((optionType+"/timingResults_"+nameFiles.at(ij)+".txt").c_str(), std::ios::out);
    std::cout << (optionType+"/timingResults_"+nameFiles.at(ij)+".txt").c_str() << std::endl; 
  }

  TH1F* hFractionEvents_HitsWithTime_Eta_dRadius;
  TH1F* hAverageTime_Eta_dRadius;
  TF1* hfithisto = new TF1("hfithisto", "gaus", -0.5, 1.);
  TF1* hfithisto2 = new TF1("hfithisto2", "gaus", -0.5, 1.);
  hfithisto2->SetLineColor(kBlue);
  TCanvas* c1;

  for(int iF=0; iF<nFiles; ++iF){
    std::cout << " >>> iF = " << iF << std::endl;

    outFileLong[iF] << " eta " << " \t" << " rad " << " \t" << " res " <<" \t" << " eff " << std::endl;   

    for(int ieta=0; ieta<numberOfBins; ++ieta){
      for(int iRad=0; iRad<nBinsRad; ++iRad){
	
	std::cout << " iF = " << iF  << " ieta = " << ieta << " irad = " << iRad << " etaBin = " << binValue.at(ieta) << std::endl;
	if(iRad < nBinsRad){
	  hFractionEvents_HitsWithTime_Eta_dRadius = (TH1F*)(inF[iF]->Get(Form(("ana/hFractionEvents_Hits_Eta"+binName.at(ieta)+"_dRadius%d").c_str(), iRad)));
	  float efficiency = hFractionEvents_HitsWithTime_Eta_dRadius->GetMean();
	  tgFrEvt[iF][iRad]->SetPoint(ieta, binValue.at(ieta),  efficiency);
	  hFractionEvents_HitsWithTime_Eta_dRadius->Delete();

	  
	  hAverageTime_Eta_dRadius = (TH1F*)(inF[iF]->Get(Form(("ana/hAverageTime_Eta"+binName.at(ieta)+"_dRadius%d_AvgInt").c_str(), iRad)));
	  hAverageTime_Eta_dRadius->Rebin(10);

	  hfithisto->SetRange(hAverageTime_Eta_dRadius->GetMean() - 2.*hAverageTime_Eta_dRadius->GetRMS(), 
			      hAverageTime_Eta_dRadius->GetMean() + 2.*hAverageTime_Eta_dRadius->GetRMS());
	  hfithisto->SetParameter(1, hAverageTime_Eta_dRadius->GetMean());
	  hfithisto->SetParameter(2, hAverageTime_Eta_dRadius->GetRMS());

	  hAverageTime_Eta_dRadius->Fit("hfithisto", "RQ");
	  hfithisto2->SetRange(hfithisto->GetParameter(1) - 2.*hfithisto->GetParameter(2), hfithisto->GetParameter(1) + 2.*hfithisto->GetParameter(2));
	  hfithisto2->SetParameters(hfithisto->GetParameter(0), hfithisto->GetParameter(1), hfithisto->GetParameter(2));
	  
	  hAverageTime_Eta_dRadius->Fit("hfithisto2", "RQ");  
	  tg[iF][iRad]->SetPoint(ieta, binValue.at(ieta),  hfithisto2->GetParameter(2));                        
	  tg[iF][iRad]->SetPointError(ieta, 0, hfithisto2->GetParError(2));     
	    
	  outFileLong[iF] << binValue.at(ieta) << " \t " << iRad << " \t " << hfithisto2->GetParameter(2) <<"+/-"<< hfithisto2->GetParError(2)
			  << " \t " << efficiency<<"+/-" << std::endl;

	  tgM[iF][iRad]->SetPoint(ieta, binValue.at(ieta),  hfithisto2->GetParameter(1));
	  tgM[iF][iRad]->SetPointError(ieta, 0, hfithisto2->GetParError(1));
	  
	  c1 = new TCanvas();
	  c1->cd();
	  hAverageTime_Eta_dRadius->GetXaxis()->SetRangeUser(-0.5, 1.);
	  hAverageTime_Eta_dRadius->Draw();
	  c1->Print(Form((singleFitsFolderName+"/avgTime_Eta"+binNameOK.at(ieta)+"_dRadius%d_"+nameFiles.at(iF)+".png").c_str(), iRad), "png");
	  c1->Delete();
	  hAverageTime_Eta_dRadius->Delete();
	}//rad < 2

      }//loop over radius
    }//loop over eta
  }//loop over files

  
  for(int iF=0; iF<nFiles; ++iF)
    outFileLong[iF].close();


  std::cout << " all done now print results " << std::endl;


  TLegend *legTGM = new TLegend(0.65,0.70,0.85,0.85,NULL,"brNDC");
  legTGM->SetTextFont(42);
  legTGM->SetTextSize(0.05);
  legTGM->SetFillColor(kWhite);
  legTGM->SetLineColor(kWhite);
  legTGM->SetShadowColor(kWhite);
  for(int iF=0; iF<nFiles; ++iF){
    legTGM->AddEntry(hDummy[iF], Form("p_{T} = %.1f GeV",ptValues.at(iF)), "l");
  }

  TLegend *legTGM2s = new TLegend(0.40,0.70,0.60,0.85,NULL,"brNDC");
  legTGM2s->SetTextFont(42);
  legTGM2s->SetTextSize(0.05);
  legTGM2s->SetFillColor(kWhite);
  legTGM2s->SetLineColor(kWhite);
  legTGM2s->SetShadowColor(kWhite);
  legTGM2s->AddEntry(hDummySt[0], " #rho < 2cm", "p");
  legTGM2s->AddEntry(hDummySt[1], " #rho < 5cm", "p");
  if(nBinsRad > 2){
    legTGM2s->AddEntry(hDummySt[2], " #rho < 10cm", "p");
    legTGM2s->AddEntry(hDummySt[3], " all", "p");
  }

  ////////// down
  TLegend *legTGMd = new TLegend(0.65,0.15,0.85,0.30,NULL,"brNDC");
  legTGMd->SetTextFont(42);
  legTGMd->SetTextSize(0.05);
  legTGMd->SetFillColor(kWhite);
  legTGMd->SetLineColor(kWhite);
  legTGMd->SetShadowColor(kWhite);
  for(int iF=0; iF<nFiles; ++iF){
    legTGMd->AddEntry(hDummy[iF], Form("p_{T} = %.1f GeV",ptValues.at(iF)), "l");
  }
  TLegend *legTGM2d = new TLegend(0.40,0.15,0.60,0.30,NULL,"brNDC");
  legTGM2d->SetTextFont(42);
  legTGM2d->SetTextSize(0.05);
  legTGM2d->SetFillColor(kWhite);
  legTGM2d->SetLineColor(kWhite);
  legTGM2d->SetShadowColor(kWhite);
  legTGM2d->AddEntry(hDummySt[0], " #rho < 2cm", "p");
  legTGM2d->AddEntry(hDummySt[1], " #rho < 5cm", "p");
  if(nBinsRad > 2){
    legTGM2d->AddEntry(hDummySt[2], " #rho < 10cm", "p");
    legTGM2d->AddEntry(hDummySt[3], " all", "p");
  }
  std::cout << " legends ok  " << std::endl;


  TLatex tL;
  tL.SetNDC();
  tL.SetTextSize(0.05);
  tL.SetTextFont(132);

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  
  TCanvas* ch_R = new TCanvas();
  ch_R->cd();
  tg[0][0]->GetXaxis()->SetTitle("|#eta|");
  tg[0][0]->GetYaxis()->SetTitle("#sigma_{t} ns");
  if(pdgID == "22")tg[0][0]->GetYaxis()->SetRangeUser(0., 0.05);
  if(pdgID == "130")tg[0][0]->GetYaxis()->SetRangeUser(0., 0.1);
  tg[0][0]->Draw("apl");
  tg[0][1]->Draw("pl, same");
  if(nBinsRad > 2){
    tg[0][2]->Draw("pl, same");
    tg[0][3]->Draw("pl, same");
  }
  for(int iF=1; iF<nFiles; ++iF){
    tg[iF][0]->Draw("pl, same");
    tg[iF][1]->Draw("pl, same");
    if(nBinsRad > 2){
      tg[iF][2]->Draw("pl, same");
      tg[iF][3]->Draw("pl, same");
    }
  }
  if(pdgID == "130")  tL.DrawLatex(0.25, 0.78, "K^{0}_{L}");
  if(pdgID == "22")   tL.DrawLatex(0.25, 0.78, "#gamma");
  legTGM->Draw("same");
  legTGM2s->Draw("same");
  ch_R->Print((plotFolder+"/h_timeResoVsEta.png").c_str(), "png");
  ch_R->Print((plotFolder+"/h_timeResoVsEta.pdf").c_str(), "pdf");  
  ch_R->Print((plotFolder+"/h_timeResoVsEta.root").c_str(), "root");
  tg[0][0]->GetYaxis()->SetRangeUser(0., 0.2);
  ch_R->Print((plotFolder+"/h_timeResoVsEta_zoomIn.png").c_str(), "png");
  ch_R->Print((plotFolder+"/h_timeResoVsEta_zoomIn.pdf").c_str(), "pdf");
  ch_R->Print((plotFolder+"/h_timeResoVsEta_zoomIn.root").c_str(), "root");


  TCanvas* ch_M = new TCanvas();
  ch_M->cd();
  tgM[0][0]->GetXaxis()->SetTitle("|#eta|");
  tgM[0][0]->GetYaxis()->SetTitle("mean(t) ns");
  tgM[0][0]->GetYaxis()->SetRangeUser(-0.02, 0.6);
  tgM[0][0]->Draw("apl");
  tgM[0][1]->Draw("pl, same");
  if(nBinsRad > 2){
     tgM[0][2]->Draw("pl, same");
     tgM[0][3]->Draw("pl, same");
  }
  for(int iF=1; iF<nFiles; ++iF){
    tgM[iF][0]->Draw("pl, same");
    tgM[iF][1]->Draw("pl, same");
    if(nBinsRad > 2){
      tgM[iF][2]->Draw("pl, same");
      tgM[iF][3]->Draw("pl, same");
    }
  }
  if(pdgID == "130")  tL.DrawLatex(0.25, 0.78, "K^{0}_{L}");
  if(pdgID == "22")   tL.DrawLatex(0.25, 0.78, "#gamma");
  legTGM->Draw("same");
  legTGM2s->Draw("same");
  ch_M->Print((plotFolder+"/h_timeMeanVsEta.png").c_str(), "png");
  ch_M->Print((plotFolder+"/h_timeMeanVsEta.pdf").c_str(), "pdf");  
  ch_M->Print((plotFolder+"/h_timeMeanVsEta.root").c_str(), "root");


  std::cout <<  " now efficiency " << std::endl;
  TCanvas* ch_FrEvt = new TCanvas();
  ch_FrEvt->cd();
  tgFrEvt[0][0]->GetXaxis()->SetTitle("|#eta|");
  tgFrEvt[0][0]->GetYaxis()->SetTitle("efficiency"); //#geq 3 hits with time / all events>");
  tgFrEvt[0][0]->GetYaxis()->SetRangeUser(0., 1.1);
  tgFrEvt[0][0]->Draw("apl");
  tgFrEvt[0][1]->Draw("pl, same");
  if(nBinsRad > 2){
    tgFrEvt[0][2]->Draw("pl, same");
    tgFrEvt[0][3]->Draw("pl, same");
  }
  for(int iF=1; iF<nFiles; ++iF){
    tgFrEvt[iF][0]->Draw("pl, same");
    tgFrEvt[iF][1]->Draw("pl, same");
    if(nBinsRad > 2){
      tgFrEvt[iF][2]->Draw("pl, same");
      tgFrEvt[iF][3]->Draw("pl, same");
    }
  }
  if(pdgID == "130")  tL.DrawLatex(0.25, 0.23, "K^{0}_{L}");
  if(pdgID == "22")   tL.DrawLatex(0.25, 0.23, "#gamma");
  legTGMd->Draw("same");
  legTGM2d->Draw("same");
  ch_FrEvt->Print((plotFolder+"/h_fractionEvtWithTime.png").c_str(), "png");
  ch_FrEvt->Print((plotFolder+"/h_fractionEvtWithTime.pdf").c_str(), "pdf");  
  ch_FrEvt->Print((plotFolder+"/h_fractionEvtWithTime.root").c_str(), "root");
 

  std::cout << " all done - ciao " << std::endl;

  return;



}
