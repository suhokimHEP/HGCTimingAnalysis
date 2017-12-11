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

  //  std::cout << histo->GetName() << " " << histo->GetBinCenter(xBin) << std::endl; 
  return histo->GetBinCenter(xBin);
}



void comparePlots_Timing_PU(){
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
  int iColors[7] = {kBlue, kRed, kBlack};
  int iStyle[7] = {20, 21, 22, 24, 25, 26};
  //  int iColors[6] = {kGren+1, kBlue, kRed};
  
  int nBinsEta = 6;
  float binWidth = 0.2;
  int nBinsRad = 1;
  //int nBinsEta = 3;
  //float binWidth = 0.5;
  float binStart = 1.65;
  
  std::string pdgID = "130";
  //std::string pdgID = "211";
  //std::string pdgID = "22";

  int rad1 = 0;

  std::vector<int> PUvalue;
  PUvalue.push_back(200);
  PUvalue.push_back(140);
  PUvalue.push_back(0);

  std::vector<std::string> PUname;
  PUname.push_back("200");
  PUname.push_back("140");
  PUname.push_back("0");



  int nFiles = 3;
  std::vector<std::string> nameFiles;
  std::vector<float> ptValues;
  nFiles = 3;
  for(int ij=0; ij<nFiles; ++ij){
    nameFiles.push_back("PDG_"+pdgID+"_pt5_"+PUname.at(ij)+"PU");
    ptValues.push_back(5);
  }


  int numberOfBins = 2;
  std::vector<std::string> binName;
  binName.push_back("1.65-1.85");
  binName.push_back("2.65-2.85");
  std::vector<float> binValue;
  binValue.push_back(1.65);
  binValue.push_back(2.65);


  int nOptions = 2;
  std::vector<std::string> nameOptions;
  std::vector<std::string> nameOptionsLeg;
  //  nameOptions.push_back("210Fix_noAged");
  //  nameOptions.push_back("210Fix_Aged");
  nameOptions.push_back("210Fix_noAged_VtxZPUfixed");
  nameOptions.push_back("210Fix_Aged_VtxZPUfixed");
  nameOptionsLeg.push_back("0 fb^{-1}");
  nameOptionsLeg.push_back("3000 fb^{-1}");


  TH1F* hDummy[3];
  TH1F* hDummySt[2];

  for(int iT=0; iT<nOptions; ++iT){
    hDummySt[iT] = new TH1F(Form("hDummySt_ageing%d", iT), "", 1, 0., 1);
    hDummySt[iT]->SetMarkerStyle(iStyle[iT+nFiles]);    
    if(iT == 0)     hDummySt[iT]->SetMarkerStyle(20);    
    if(iT == 1)     hDummySt[iT]->SetMarkerStyle(24);    
    hDummySt[iT]->SetMarkerColor(kBlack);    
    hDummySt[iT]->SetMarkerSize(1.5);    
  }
  for(int iT=0; iT<nFiles; ++iT){
    hDummy[iT] = new TH1F(Form("hDummy_PU%d", iT), "", 1, 0., 1);
    hDummy[iT]->SetLineColor(iColors[iT]);
    hDummy[iT]->SetMarkerColor(iColors[iT]);
    hDummy[iT]->SetMarkerStyle(iStyle[iT]);
    hDummy[iT]->SetMarkerSize(1.5);
  }



  for(int ij=0; ij<nameFiles.size(); ++ij) std::cout << " name = " << nameFiles.at(ij) << std::endl;

  std::cout << " >>> ora riempio TGraph " << std::endl;
  
  TGraphErrors* tgM[8][3];
  TGraphErrors* tg[8][3];
  TGraphErrors* tgFrEvt[8][3];
  for(int iR=0; iR<nOptions; ++iR){
    for(int iT=0; iT<nFiles; ++iT){
      tgM[iR][iT] = new TGraphErrors();
      tgM[iR][iT]->SetName(Form("mean_%d_opt%d", iT, iR));
      tgM[iR][iT]->SetPoint(0, -1, -1);

      tgM[iR][iT]->SetLineColor(iColors[iT]);
      tgM[iR][iT]->SetLineWidth(2);
      tgM[iR][iT]->SetMarkerColor(iColors[iT]);
      tgM[iR][iT]->SetMarkerStyle(iStyle[(3 * iR)+iT]);
      tgM[iR][iT]->SetMarkerSize(1.5);
      ///////////
      tg[iR][iT] = new TGraphErrors();
      tg[iR][iT]->SetName(Form("reso_%d_opt%d", iT, iR));
      tg[iR][iT]->SetPoint(0, -1, -1);

      tg[iR][iT]->SetLineColor(iColors[iT]);
      tg[iR][iT]->SetLineWidth(2);
      tg[iR][iT]->SetMarkerColor(iColors[iT]);
      tg[iR][iT]->SetMarkerStyle(iStyle[(3 * iR)+iT]);
      tg[iR][iT]->SetMarkerSize(1.5);
      /////////
      tgFrEvt[iR][iT] = new TGraphErrors();
      tgFrEvt[iR][iT]->SetName(Form("fractionEvt_%d_opt%d", iT, iR));
      tgFrEvt[iR][iT]->SetPoint(0, -1, -1);
      
      tgFrEvt[iR][iT]->SetLineColor(iColors[iT]);
      tgFrEvt[iR][iT]->SetLineWidth(2);
      tgFrEvt[iR][iT]->SetMarkerColor(iColors[iT]);
      tgFrEvt[iR][iT]->SetMarkerStyle(iStyle[(3 * iR)+iT]);
      tgFrEvt[iR][iT]->SetMarkerSize(1.5);
    }
  }

  std::string folder = "combinationOption_PU_agedANDnon";  

  std::cout << " >>> ora prendo i file " << std::endl;   
  TFile* inF[8][3];
  for(int iR=0; iR<nOptions; ++iR){
    for(int ij=0; ij<nFiles; ++ij){
      std::string inFileName;
      if(iR == 0) inFileName = "../../test/testPU/"+nameOptions.at(iR)+"/OutTimeHGC_RecHits_"+nameFiles.at(ij)+"_allEta.root";
      if(iR == 1) inFileName = "../../test/testPU/"+nameOptions.at(iR)+"/OutTimeHGC_RecHits_"+nameFiles.at(ij)+"_allEta_aged.root";
      inF[iR][ij] = TFile::Open(inFileName.c_str());
    }
  }
  

  std::string folderName = "combinationOption_PU_agedANDnon_singleFit_"+pdgID;
  std::cout << " >>> fatto = presi " << std::endl;

  // std::ofstream outFileLong[7];
  // for(int ij=0; ij<nFiles; ++ij)
  //   outFileLong[ij].open(Form("myTiming_CFD_%d.txt", ij), std::ios::out);
  
  //  std::ofstream outFile("myTiming_CFD.txt", std::ios::app);


  TH1F* hFractionEvents_HitsWithTime_Eta_dRadius;
  TH1F* hAverageTime_Eta_dRadius;
  TF1* hfithisto = new TF1("hfithisto", "gaus", -5., 5.);
  TF1* hfithisto2 = new TF1("hfithisto2", "gaus", -5., 5.);
  hfithisto2->SetLineColor(kBlue);
  TF1* hradiusFit = new TF1("hradiusFit", "landau", 0., 5.);
  TCanvas* c1;

  for(int iT=0; iT<nOptions; ++iT){
    for(int iF=0; iF<nFiles; ++iF){
      std::cout << " >>> iF = " << iF << " iT = " << iT  << " name file = " << inF[iT][iF]->GetName() << std::endl;

      for(int ieta=0; ieta<numberOfBins; ++ieta){
	for(int iRad=0; iRad<nBinsRad; ++iRad){	
	  //	  std::cout << " iF = " << iF  << " ieta = " << ieta << " irad = " << iRad << " name file = " << inF[iT][iF]->GetName() << std::endl;
	  if(iRad < 2){
	  if(doAllTheFits == true){
	    hAverageTime_Eta_dRadius = (TH1F*)(inF[iT][iF]->Get(Form(("ana/hAverageTime_Eta"+binName.at(ieta)+"_dRadius%d_AvgInt").c_str(), iRad)) );
	    hAverageTime_Eta_dRadius->Rebin(5);


	    hfithisto->SetRange(hAverageTime_Eta_dRadius->GetMean() - 2.*hAverageTime_Eta_dRadius->GetRMS(), hAverageTime_Eta_dRadius->GetMean() + 2.*hAverageTime_Eta_dRadius->GetRMS());
	    hfithisto->SetParameter(1, hAverageTime_Eta_dRadius->GetMean());
	    hfithisto->SetParameter(2, hAverageTime_Eta_dRadius->GetRMS());
	      	      
	    hAverageTime_Eta_dRadius->Fit("hfithisto", "RQ");
	    hfithisto2->SetRange(hfithisto->GetParameter(1) - 2.*hfithisto->GetParameter(2), hfithisto->GetParameter(1) + 2.*hfithisto->GetParameter(2));
	    hfithisto2->SetParameters(hfithisto->GetParameter(0), hfithisto->GetParameter(1), hfithisto->GetParameter(2));

	    hAverageTime_Eta_dRadius->Fit("hfithisto2", "RQ");  
	    if(iT == 1) tg[iT][iF]->SetPoint(ieta, binValue.at(ieta)+0.05,  hfithisto2->GetParameter(2));                        
	    else tg[iT][iF]->SetPoint(ieta, binValue.at(ieta),  hfithisto2->GetParameter(2));                        
	    tg[iT][iF]->SetPointError(ieta, 0, hfithisto2->GetParError(2));     
	    if(iT == 1) tgM[iT][iF]->SetPoint(ieta, binValue.at(ieta)+0.05,  hfithisto2->GetParameter(1));
	    else tgM[iT][iF]->SetPoint(ieta, binValue.at(ieta),  hfithisto2->GetParameter(1));
	    tgM[iT][iF]->SetPointError(ieta, 0, hfithisto2->GetParError(1));
	  }
	  /*
	  c1 = new TCanvas();
	  c1->cd();
	  hAverageTime_Eta_dRadius->GetXaxis()->SetRangeUser(-1., 1.);
	  hAverageTime_Eta_dRadius->Draw();
	  //	hfithisto2->Draw("same");
	  c1->Print(Form((folderName+"/avgTime_Eta"+binName.at(ieta)+"_dRadius%d_"+nameFiles.at(iF)+".png").c_str(), iRad), "png");
	  c1->Delete();
	  hAverageTime_Eta_dRadius->Delete();
	  */
	  }	  
	  hFractionEvents_HitsWithTime_Eta_dRadius = (TH1F*)(inF[iT][iF]->Get(Form(("ana/hFractionEvents_Hits_Eta"+binName.at(ieta)+"_dRadius%d").c_str(), iRad)));
	  float efficiency = hFractionEvents_HitsWithTime_Eta_dRadius->GetMean();
	  if(iT == 1)tgFrEvt[iT][iF]->SetPoint(ieta, binValue.at(ieta)+0.05,  efficiency);
	  else tgFrEvt[iT][iF]->SetPoint(ieta, binValue.at(ieta),  efficiency);
	  std::cout << " eta = " << binValue.at(ieta) << " efficiency = " << efficiency << std::endl;
	  hFractionEvents_HitsWithTime_Eta_dRadius->Delete();
	}//rad

      }//eta
    }//file
  }//type


  std::cout << " ci sono " << std::endl;

  // gStyle->SetStatStyle(0);
  // gStyle->SetTitleStyle(0);
  // gROOT->ForceStyle();

  TLegend *legTGM = new TLegend(0.77,0.7,0.87,0.93,NULL,"brNDC");
  legTGM->SetTextFont(42);
  legTGM->SetTextSize(0.04);
  legTGM->SetFillStyle(0); 
  //  legTGM->SetFillColor(kWhite);
  legTGM->SetLineColor(kWhite);
  legTGM->SetShadowColor(kWhite);
  for(int iF=0; iF<nFiles; ++iF){
    legTGM->AddEntry(hDummy[iF], Form(" %d PU",PUvalue.at(iF)), "p");
  }

  TLegend *legTGMD = new TLegend(0.77,0.30,0.87,0.53,NULL,"brNDC");
  legTGMD->SetTextFont(42);
  legTGMD->SetTextSize(0.04);
  legTGMD->SetFillStyle(0); 
  //  legTGMD->SetFillColor(kWhite);
  legTGMD->SetLineColor(kWhite);
  legTGMD->SetShadowColor(kWhite);
  for(int iF=0; iF<nFiles; ++iF){
    legTGMD->AddEntry(hDummy[iF], Form(" %d PU",PUvalue.at(iF)), "p");
  }



  TLegend *legTGMOpt = new TLegend(0.57,0.70,0.72,0.93,NULL,"brNDC");
  legTGMOpt->SetTextFont(42);
  legTGMOpt->SetTextSize(0.04);
  legTGMOpt->SetFillColor(-1);
  legTGMOpt->SetFillStyle(0); 
  legTGMOpt->SetLineColor(kWhite);
  legTGMOpt->SetShadowColor(kWhite);
  for(int iF=0; iF<nOptions; ++iF){
    legTGMOpt->AddEntry(hDummySt[iF], (nameOptionsLeg.at(iF)).c_str(), "p");
  }


  TLegend *legTGMOptD = new TLegend(0.57,0.30,0.72,0.53,NULL,"brNDC");
  legTGMOptD->SetTextFont(42);
  legTGMOptD->SetTextSize(0.04);
  legTGMOptD->SetFillColor(-1);
  legTGMOptD->SetFillStyle(0); 
  legTGMOptD->SetLineColor(kWhite);
  legTGMOptD->SetShadowColor(kWhite);
  for(int iF=0; iF<nOptions; ++iF){
    legTGMOptD->AddEntry(hDummySt[iF], (nameOptionsLeg.at(iF)).c_str(), "p");
  }
  
  TLatex tCMS;
  tCMS.SetNDC();
  tCMS.SetTextSize(0.05);
  tCMS.SetTextFont(132);

  TLatex latexLabel;
  latexLabel.SetNDC();
  latexLabel.SetTextSize(0.045);
  //  latexLabel.DrawLatex(0.18, 0.96, "#font[132]{CMS Phase-2 Simulation}");
  //latexLabel.DrawLatex(0.77, 0.955, "#font[132]{#sqrt{s} = 7 TeV}");
  //latexLabel.DrawLatex(0.2, 0.85, "#font[132]{#intL dt= 36.1 pb^{-1}}");


  TLatex tL;
  tL.SetNDC();
  tL.SetTextSize(0.05);
  tL.SetTextFont(132);
  TLatex tL2;
  tL2.SetNDC();
  tL2.SetTextSize(0.05);
  tL2.SetTextFont(132);

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  if(doAllTheFits == true){
  TCanvas* ch_M = new TCanvas();
  ch_M->cd();
  tgM[0][0]->GetXaxis()->SetTitle("|#eta|");
  tgM[0][0]->GetYaxis()->SetTitle("<time> (ns)");
  tgM[0][0]->GetXaxis()->SetRangeUser(1.5, 3.);
  tgM[0][0]->GetYaxis()->SetRangeUser(0., 0.6);
  if(pdgID == "22") tgM[0][0]->GetYaxis()->SetRangeUser(0., 0.05);
  tgM[0][0]->Draw("ap");
  for(int iF=1; iF<nFiles; ++iF){
    tgM[0][iF]->Draw("p, same");
  }

  for(int iT=1; iT<nOptions; ++iT){
    tgM[iT][0]->Draw("p, same");
    for(int iF=1; iF<nFiles; ++iF){
      tgM[iT][iF]->Draw("p, same");
    } 
  }
  legTGM->Draw("same");
  legTGMOpt->Draw("same");

  latexLabel.DrawLatex(0.20, 0.96, "#font[62]{CMS} #font[52]{Phase-2 Simulation}");  

  if(rad1 == 0 && pdgID == "130") {
    tL.DrawLatex(0.23, 0.85, "K^{0}_{L}    p_{T} = 5 GeV");
    tL2.DrawLatex(0.23, 0.75, "#rho #leq 2cm");
  }
  if(rad1 == 0 && pdgID == "22") {
    tL.DrawLatex(0.23, 0.85, "#gamma    p_{T} = 5 GeV");
    tL2.DrawLatex(0.23, 0.75, "#rho #leq 2cm");
  }

  if(rad1 == 1 && pdgID == "130")  tL.DrawLatex(0.25, 0.85, "#rho #leq 5cm  K^{0}_{L}");
  if(rad1 == 1 && pdgID == "22")  tL.DrawLatex(0.25, 0.85, "#rho #leq 5cm  #gamma");

  if(rad1 == 1){
  ch_M->Print((folder+"/h_timeMeanVsEta_pdg"+pdgID+"_rad1.png").c_str(), "png");
  ch_M->Print((folder+"/h_timeMeanVsEta_pdg"+pdgID+"_rad1.pdf").c_str(), "pdf");  
  ch_M->Print((folder+"/h_timeMeanVsEta_pdg"+pdgID+"_rad1.root").c_str(), "root");
  }
  else{
  ch_M->Print((folder+"/h_timeMeanVsEta_pdg"+pdgID+".png").c_str(), "png");
  ch_M->Print((folder+"/h_timeMeanVsEta_pdg"+pdgID+".pdf").c_str(), "pdf");  
  ch_M->Print((folder+"/h_timeMeanVsEta_pdg"+pdgID+".root").c_str(), "root");

  }

  TCanvas* ch_R = new TCanvas();
  ch_R->cd();
  tg[0][0]->GetXaxis()->SetTitle("|#eta|");
  tg[0][0]->GetYaxis()->SetTitle("#sigma_{time} (ns)");
  tg[0][0]->GetXaxis()->SetRangeUser(1.5, 3.);
  tg[0][0]->GetYaxis()->SetRangeUser(0., 0.1);
  tg[0][0]->Draw("ap");

  tg[0][1]->Draw("p, same");
  tg[1][1]->Draw("p, same");
  tg[0][0]->Draw("p, same");
  tg[1][0]->Draw("p, same");
  tg[0][2]->Draw("p, same");
  tg[1][2]->Draw("p, same");

  /*
  for(int iF=1; iF<nFiles; ++iF){
    tg[0][iF]->Draw("p, same");
  }

  for(int iT=1; iT<nOptions; ++iT){
    tg[iT][0]->Draw("p, same");
    for(int iF=1; iF<nFiles; ++iF){
      tg[iT][iF]->Draw("p, same");
    } 
  }
  */

  legTGM->Draw("same");
  legTGMOpt->Draw("same");

  latexLabel.DrawLatex(0.20, 0.96, "#font[62]{CMS} #font[52]{Phase-2 Simulation}");    

  if(rad1 == 0 && pdgID == "130") {
    tL.DrawLatex(0.23, 0.85, "K^{0}_{L}    p_{T} = 5 GeV");
    tL2.DrawLatex(0.23, 0.75, "#rho #leq 2cm");
  }
  if(rad1 == 0 && pdgID == "22") {
    tL.DrawLatex(0.23, 0.85, "#gamma    p_{T} = 5 GeV");
    tL2.DrawLatex(0.23, 0.75, "#rho #leq 2cm");
  }
  if(rad1 == 1 && pdgID == "130")  tL.DrawLatex(0.25, 0.85, "#rho #leq 5cm  K^{0}_{L}");
  if(rad1 == 1 && pdgID == "22")  tL.DrawLatex(0.25, 0.85, "#rho #leq 5cm  #gamma");

  if(rad1 == 1){
  ch_R->Print((folder+"/h_timeResoVsEta_pdg"+pdgID+"_rad1.png").c_str(), "png");
  ch_R->Print((folder+"/h_timeResoVsEta_pdg"+pdgID+"_rad1.pdf").c_str(), "pdf");  
  ch_R->Print((folder+"/h_timeResoVsEta_pdg"+pdgID+"_rad1.root").c_str(), "root");
  }
  else{
  ch_R->Print((folder+"/h_timeResoVsEta_pdg"+pdgID+".png").c_str(), "png");
  ch_R->Print((folder+"/h_timeResoVsEta_pdg"+pdgID+".pdf").c_str(), "pdf");  
  ch_R->Print((folder+"/h_timeResoVsEta_pdg"+pdgID+".root").c_str(), "root");
  if(pdgID == "22") tg[0][0]->GetYaxis()->SetRangeUser(0., 0.05);
  ch_R->Print((folder+"/h_timeResoVsEta_pdg"+pdgID+"_zoomIn.png").c_str(), "png");
  ch_R->Print((folder+"/h_timeResoVsEta_pdg"+pdgID+"_zoomIn.pdf").c_str(), "pdf");  
  ch_R->Print((folder+"/h_timeResoVsEta_pdg"+pdgID+"_zoomIn.root").c_str(), "root");
  }

  TCanvas* ch_FrEvt = new TCanvas();
  ch_FrEvt->cd();
  tgFrEvt[0][0]->GetXaxis()->SetTitle("|#eta|");
  tgFrEvt[0][0]->GetYaxis()->SetTitle("efficiency"); //#geq 3 hits with time / all events>");
  tgFrEvt[0][0]->GetXaxis()->SetRangeUser(1.5, 3.);
  tgFrEvt[0][0]->GetYaxis()->SetRangeUser(0., 1.1);
  //tgFrEvt[0][0]->GetYaxis()->SetRangeUser(0., 1.6);
  tgFrEvt[0][0]->Draw("ap");

  tgFrEvt[0][1]->Draw("p, same");
  tgFrEvt[1][1]->Draw("p, same");
  tgFrEvt[0][0]->Draw("p, same");
  tgFrEvt[1][0]->Draw("p, same");
  tgFrEvt[0][2]->Draw("p, same");
  tgFrEvt[1][2]->Draw("p, same");


  /*
  for(int iF=1; iF<nFiles; ++iF){
    tgFrEvt[0][iF]->Draw("p, same");
  }
  for(int iT=1; iT<nOptions; ++iT){
    tgFrEvt[iT][0]->Draw("p, same");
    for(int iF=1; iF<nFiles; ++iF){
      tgFrEvt[iT][iF]->Draw("p, same");
    }
  }
  */

  legTGMD->Draw("same");
  //legTGM->Draw("same");
  legTGMOptD->Draw("same");

  latexLabel.DrawLatex(0.20, 0.96, "#font[62]{CMS} #font[52]{Phase-2 Simulation}");  

  if(rad1 == 0 && pdgID == "130") {
    tL.DrawLatex(0.23, 0.45, "K^{0}_{L}    p_{T} = 5 GeV");
    tL2.DrawLatex(0.23, 0.35, "#rho #leq 2cm");
  }
  if(rad1 == 0 && pdgID == "22") {
    tL.DrawLatex(0.23, 0.45, "#gamma    p_{T} = 5 GeV");
    tL2.DrawLatex(0.23, 0.35, "#rho #leq 2cm");
  }
  if(rad1 == 1 && pdgID == "130")  tL.DrawLatex(0.25, 0.25, "#rho #leq 5cm  K^{0}_{L}");
  if(rad1 == 1 && pdgID == "22")  tL.DrawLatex(0.25, 0.25, "#rho #leq 5cm  #gamma");


  if(rad1 == 1){
  ch_FrEvt->Print((folder+"/h_fractionEvtWithTime"+pdgID+"_rad1.png").c_str(), "png");
  ch_FrEvt->Print((folder+"/h_fractionEvtWithTime"+pdgID+"_rad1.pdf").c_str(), "pdf");  
  ch_FrEvt->Print((folder+"/h_fractionEvtWithTime"+pdgID+"_rad1.root").c_str(), "root");
  }
  else{
  ch_FrEvt->Print((folder+"/h_fractionEvtWithTime"+pdgID+".png").c_str(), "png");
  ch_FrEvt->Print((folder+"/h_fractionEvtWithTime"+pdgID+".pdf").c_str(), "pdf");  
  ch_FrEvt->Print((folder+"/h_fractionEvtWithTime"+pdgID+".root").c_str(), "root");
  }
  
  }  

  return;

}
