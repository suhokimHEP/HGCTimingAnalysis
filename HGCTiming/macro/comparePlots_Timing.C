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



void comparePlots_Timing(){
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
  int iStyle[8] = {20, 24, 21, 25, 22, 26, 23, 27}; //kGray+1};
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

  int nFiles = 3;
  std::vector<std::string> nameFiles;
  std::vector<float> ptValues;
  nFiles = 3;
  nameFiles.push_back("PDG_"+pdgID+"_pt2");
  nameFiles.push_back("PDG_"+pdgID+"_pt5");
  nameFiles.push_back("PDG_"+pdgID+"_pt10");
  ptValues.push_back(2);
  ptValues.push_back(5);
  ptValues.push_back(10);

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
  nameOptions.push_back("ACfix_bkp210FixTime_double_0PU_allEta");
  nameOptions.push_back("ACfix_bkp210FixTime_double_0PU_allEta_aged");
  nameOptionsLeg.push_back("0 fb^{-1}");
  nameOptionsLeg.push_back("3000 fb^{-1}");

   // nameOptions.push_back("CSF20LBA40");
   // nameOptions.push_back("CSF20LEA40");
   // nameOptions.push_back("CSF30LBA40");
   // nameOptions.push_back("CSF30LEA40");

  // nameOptions.push_back("CSF30LBA20");   
  // nameOptions.push_back("CSF30LEA20");   


  TH1F* hDummy[3];
  TH1F* hDummySt[8];

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
      tgM[iR][iT]->SetMarkerStyle(iStyle[iR]);
      //      tgM[iR][iT]->SetMarkerSize(1.5);
      ///////////
      tg[iR][iT] = new TGraphErrors();
      tg[iR][iT]->SetName(Form("reso_%d_opt%d", iT, iR));
      tg[iR][iT]->SetPoint(0, -1, -1);

      tg[iR][iT]->SetLineColor(iColors[iT]);
      tg[iR][iT]->SetLineWidth(2);
      tg[iR][iT]->SetMarkerColor(iColors[iT]);
      tg[iR][iT]->SetMarkerStyle(iStyle[iR]);
      //      tg[iR][iT]->SetMarkerSize(1.5);
      /////////
      tgFrEvt[iR][iT] = new TGraphErrors();
      tgFrEvt[iR][iT]->SetName(Form("fractionEvt_%d_opt%d", iT, iR));
      tgFrEvt[iR][iT]->SetPoint(0, -1, -1);
      
      tgFrEvt[iR][iT]->SetLineColor(iColors[iT]);
      tgFrEvt[iR][iT]->SetLineWidth(2);
      tgFrEvt[iR][iT]->SetMarkerColor(iColors[iT]);
      tgFrEvt[iR][iT]->SetMarkerStyle(iStyle[iR]);
      if(iR%2 == 1) tgFrEvt[iR][iT]->SetMarkerSize(1.5);
    }
  }

  std::string folder = "combinationOption_"+nameOptions.at(0);  
  //  std::cout << "../test/timingStudies/ " << nameOptions.at(0) << " /JOB_ " << nameFiles.at(0) << " _3fC/OutTimeHGC_RecHits_ " << nameFiles.at(0) << " _3fC.root " << std::endl;

  std::cout << " >>> ora prendo i file " << std::endl;   
  TFile* inF[8][3];
  for(int iR=0; iR<nOptions; ++iR){
    for(int ij=0; ij<nFiles; ++ij){
      std::string inFileName;
      if(iR == 0) inFileName = "../test/testPU/"+nameOptions.at(iR)+"/OutTimeHGC_RecHits_"+nameFiles.at(ij)+"_0PU_allEta.root";
      if(iR == 1) inFileName = "../test/testPU/"+nameOptions.at(iR)+"/OutTimeHGC_RecHits_"+nameFiles.at(ij)+"_0PU_allEta_aged.root";

      // if(nameOptions.at(iR) == "CSF30LE" && ij == 6 && pdgID == "130") 
      // 	inFileName =  "../test/timingStudies/"+nameOptions.at(iR)+"/JOB_"+nameFiles.at(ij)+"_3fC/OutTimeHGC_RecHits_"+nameFiles.at(ij)+"_3fC.root_bkp";

      // else if(nameOptions.at(iR) == "CHF20LE" && ij == 6) 
      // 	inFileName = "../test/timingStudies/CHF30LE/JOB_"+nameFiles.at(ij)+"_3fC/OutTimeHGC_RecHits_"+nameFiles.at(ij)+"_3fC.root";

      // else      if(nameOptions.at(iR) == "CHF20LE" && ij == 1) 
      // 	inFileName = "../test/timingStudies/CHF30LE/JOB_"+nameFiles.at(ij)+"_3fC/OutTimeHGC_RecHits_"+nameFiles.at(ij)+"_3fC.root";

      // else
      // 	inFileName =  "../test/timingStudies/"+nameOptions.at(iR)+"/JOB_"+nameFiles.at(ij)+"_3fC/OutTimeHGC_RecHits_"+nameFiles.at(ij)+"_3fC.root";


      // std::cout << inFileName << std::endl;
      inF[iR][ij] = TFile::Open(inFileName.c_str());
    }
  }
  

  std::string folderName = "compareAgeing_singleFit_"+pdgID;
  std::cout << " >>> fatto = presi " << std::endl;

  // std::ofstream outFileLong[7];
  // for(int ij=0; ij<nFiles; ++ij)
  //   outFileLong[ij].open(Form("myTiming_CFD_%d.txt", ij), std::ios::out);
  
  //  std::ofstream outFile("myTiming_CFD.txt", std::ios::app);


  TH1F* hFractionEvents_HitsWithTime_Eta_dRadius;
  TH1F* hAverageTime_Eta_dRadius;
  TF1* hfithisto = new TF1("hfithisto", "gaus", -0.5, 1.);
  TF1* hfithisto2 = new TF1("hfithisto2", "gaus", -0.5, 1.);
  hfithisto2->SetLineColor(kBlue);
  TF1* hradiusFit = new TF1("hradiusFit", "landau", 0., 5.);
  TCanvas* c1;

  for(int iT=0; iT<nOptions; ++iT){
    for(int iF=0; iF<nFiles; ++iF){
      std::cout << " >>> iF = " << iF << " iT = " << iT<< std::endl;

      for(int ieta=0; ieta<numberOfBins; ++ieta){
	for(int iRad=0; iRad<nBinsRad; ++iRad){	
	  std::cout << " iF = " << iF  << " ieta = " << ieta << " irad = " << iRad << std::endl;
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
	    tg[iT][iF]->SetPoint(ieta, binValue.at(ieta),  hfithisto2->GetParameter(2));                        
	    tg[iT][iF]->SetPointError(ieta, 0, hfithisto2->GetParError(2));     
	    
	    tgM[iT][iF]->SetPoint(ieta, binValue.at(ieta),  hfithisto2->GetParameter(1));
	    tgM[iT][iF]->SetPointError(ieta, 0, hfithisto2->GetParError(1));
	  }

	  c1 = new TCanvas();
	  c1->cd();
	  hAverageTime_Eta_dRadius->GetXaxis()->SetRangeUser(-0.05, 1.);
	  hAverageTime_Eta_dRadius->Draw();
	  //	hfithisto2->Draw("same");
	  c1->Print(Form((folderName+"/avgTime_Eta"+binName.at(ieta)+"_dRadius%d_"+nameFiles.at(iF)+".png").c_str(), iRad), "png");
	  c1->Delete();
	  hAverageTime_Eta_dRadius->Delete();

	  }	  
	  hFractionEvents_HitsWithTime_Eta_dRadius = (TH1F*)(inF[iT][iF]->Get(Form(("ana/hFractionEvents_Hits_Eta"+binName.at(ieta)+"_dRadius%d").c_str(), iRad)));
	  float efficiency = hFractionEvents_HitsWithTime_Eta_dRadius->GetMean();
	  tgFrEvt[iT][iF]->SetPoint(ieta, binValue.at(ieta),  efficiency);
	  hFractionEvents_HitsWithTime_Eta_dRadius->Delete();
	}//rad

      }//eta
    }//file
  }//type


  std::cout << " ci sono " << std::endl;

  TLegend *legTGM = new TLegend(0.7,0.7,0.90,0.93,NULL,"brNDC");
  legTGM->SetTextFont(42);
  legTGM->SetTextSize(0.03);
  legTGM->SetFillColor(kWhite);
  legTGM->SetLineColor(kWhite);
  legTGM->SetShadowColor(kWhite);
  for(int iF=0; iF<nFiles; ++iF){
    //legTGM->AddEntry(hDummy[iF], (nameFiles.at(iF)).c_str(), "l");
    legTGM->AddEntry(hDummy[iF], Form("p_{T} = %.1f GeV",ptValues.at(iF)), "l");
  }

  TLegend *legTGMShort = new TLegend(0.7,0.7,0.90,0.93,NULL,"brNDC");
  legTGMShort->SetTextFont(42);
  legTGMShort->SetTextSize(0.03);
  legTGMShort->SetFillColor(kWhite);
  legTGMShort->SetLineColor(kWhite);
  legTGMShort->SetShadowColor(kWhite);
  for(int iF=2; iF<nFiles; ++iF){
    //legTGM->AddEntry(hDummy[iF], (nameFiles.at(iF)).c_str(), "l");
    legTGMShort->AddEntry(hDummy[iF], Form("p_{T} = %.1f GeV",ptValues.at(iF)), "l");
  }

  TLegend *legTGMD = new TLegend(0.7,0.15,0.90,0.35,NULL,"brNDC");
  legTGMD->SetTextFont(42);
  legTGMD->SetTextSize(0.05);
  legTGMD->SetFillColor(kWhite);
  legTGMD->SetLineColor(kWhite);
  legTGMD->SetShadowColor(kWhite);
  for(int iF=0; iF<nFiles; ++iF){
    //legTGM->AddEntry(hDummy[iF], (nameFiles.at(iF)).c_str(), "l");
    legTGMD->AddEntry(hDummy[iF], Form("p_{T} = %.1f GeV",ptValues.at(iF)), "l");
  }


  TLegend *legTGMShortD = new TLegend(0.7,0.13,0.90,0.33,NULL,"brNDC");
  legTGMShortD->SetTextFont(42);
  legTGMShortD->SetTextSize(0.05);
  legTGMShortD->SetFillColor(kWhite);
  legTGMShortD->SetLineColor(kWhite);
  legTGMShortD->SetShadowColor(kWhite);
  for(int iF=2; iF<nFiles; ++iF){
    //legTGM->AddEntry(hDummy[iF], (nameFiles.at(iF)).c_str(), "l");
    legTGMShortD->AddEntry(hDummy[iF], Form("p_{T} = %.1f GeV",ptValues.at(iF)), "l");
  }


  TLegend *legTGMOpt = new TLegend(0.5,0.75,0.65,0.93,NULL,"brNDC");
  legTGMOpt->SetTextFont(42);
  legTGMOpt->SetTextSize(0.05);
  legTGMOpt->SetFillColor(kWhite);
  legTGMOpt->SetLineColor(kWhite);
  legTGMOpt->SetShadowColor(kWhite);
  for(int iF=0; iF<nOptions; ++iF){
    legTGMOpt->AddEntry(hDummySt[iF], (nameOptionsLeg.at(iF)).c_str(), "p");
  }


  TLegend *legTGMOptD = new TLegend(0.5,0.15,0.65,0.33,NULL,"brNDC");
  legTGMOptD->SetTextFont(42);
  legTGMOptD->SetTextSize(0.05);
  legTGMOptD->SetFillColor(kWhite);
  legTGMOptD->SetLineColor(kWhite);
  legTGMOptD->SetShadowColor(kWhite);
  for(int iF=0; iF<nOptions; ++iF){
    legTGMOptD->AddEntry(hDummySt[iF], (nameOptionsLeg.at(iF)).c_str(), "p");
  }
  
  TLatex tL;
  tL.SetNDC();
  tL.SetTextSize(0.05);
  tL.SetTextFont(132);

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  if(doAllTheFits == true){
  TCanvas* ch_M = new TCanvas();
  ch_M->cd();
  tgM[0][0]->GetXaxis()->SetTitle("|#eta|");
  tgM[0][0]->GetYaxis()->SetTitle("<time> (ns)");
  tgM[0][0]->GetYaxis()->SetRangeUser(0., 0.6);
  if(pdgID == "22") tgM[0][0]->GetYaxis()->SetRangeUser(0., 0.05);
  tgM[0][0]->Draw("apl");
  for(int iF=1; iF<nFiles; ++iF){
    tgM[0][iF]->Draw("pl, same");
  }


  for(int iT=1; iT<nOptions; ++iT){
    tgM[iT][0]->Draw("pl, same");
    for(int iF=1; iF<nFiles; ++iF){
      tgM[iT][iF]->Draw("pl, same");
    } 
  }
  legTGM->Draw("same");
  legTGMOpt->Draw("same");

  if(rad1 == 0 && pdgID == "130")  tL.DrawLatex(0.25, 0.85, "#rho #leq 2cm  K^{0}_{L}");
  if(rad1 == 1 && pdgID == "130")  tL.DrawLatex(0.25, 0.85, "#rho #leq 5cm  K^{0}_{L}");
  if(rad1 == 0 && pdgID == "22")  tL.DrawLatex(0.25, 0.85, "#rho #leq 2cm  #gamma");
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
  tg[0][0]->GetYaxis()->SetRangeUser(0., 0.2);
  if(pdgID == "22") tg[0][0]->GetYaxis()->SetRangeUser(0., 0.05);
  tg[0][0]->Draw("apl");
  for(int iF=1; iF<nFiles; ++iF){
    tg[0][iF]->Draw("pl, same");
  }

  for(int iT=1; iT<nOptions; ++iT){
    tg[iT][0]->Draw("pl, same");
    for(int iF=1; iF<nFiles; ++iF){
      tg[iT][iF]->Draw("pl, same");
    } 
  }
  legTGM->Draw("same");
  legTGMOpt->Draw("same");

  if(rad1 == 0 && pdgID == "130")  tL.DrawLatex(0.25, 0.85, "#rho #leq 2cm  K^{0}_{L}");
  if(rad1 == 1 && pdgID == "130")  tL.DrawLatex(0.25, 0.85, "#rho #leq 5cm  K^{0}_{L}");
  if(rad1 == 0 && pdgID == "22")  tL.DrawLatex(0.25, 0.85, "#rho #leq 2cm  #gamma");
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
  }

  TCanvas* ch_Rcut = new TCanvas();
  ch_Rcut->cd();
  tg[0][2]->GetXaxis()->SetTitle("|#eta|");
  tg[0][2]->GetYaxis()->SetTitle("#sigma_{time} (ns)");
  tg[0][2]->GetYaxis()->SetRangeUser(0., 0.2);
  if(pdgID == "22") tg[0][2]->GetYaxis()->SetRangeUser(0., 0.05);
  tg[0][2]->Draw("apl");
  for(int iF=3; iF<nFiles; ++iF){
    tg[0][iF]->Draw("pl, same");
  } 
  for(int iT=1; iT<nOptions; ++iT){
    tg[iT][2]->Draw("pl, same");
    for(int iF=3; iF<nFiles; ++iF){
      tg[iT][iF]->Draw("pl, same");
    } 
  }
  legTGMShort->Draw("same");
  //legTGM->Draw("same");
  legTGMOpt->Draw("same");

  if(rad1 == 0 && pdgID == "130")  tL.DrawLatex(0.25, 0.85, "#rho #leq 2cm  K^{0}_{L}");
  if(rad1 == 1 && pdgID == "130")  tL.DrawLatex(0.25, 0.85, "#rho #leq 5cm  K^{0}_{L}");
  if(rad1 == 0 && pdgID == "22")  tL.DrawLatex(0.25, 0.85, "#rho #leq 2cm  #gamma");
  if(rad1 == 1 && pdgID == "22")  tL.DrawLatex(0.25, 0.85, "#rho #leq 5cm  #gamma");


  if(rad1 == 1){
  ch_Rcut->Print((folder+"/h_timeResoVsEtaCut_pdg"+pdgID+"_rad1.png").c_str(), "png");
  ch_Rcut->Print((folder+"/h_timeResoVsEtaCut_pdg"+pdgID+"_rad1.pdf").c_str(), "pdf");  
  ch_Rcut->Print((folder+"/h_timeResoVsEtaCut_pdg"+pdgID+"_rad1.root").c_str(), "root");
  }
  else{
  ch_Rcut->Print((folder+"/h_timeResoVsEtaCut_pdg"+pdgID+".png").c_str(), "png");
  ch_Rcut->Print((folder+"/h_timeResoVsEtaCut_pdg"+pdgID+".pdf").c_str(), "pdf");  
  ch_Rcut->Print((folder+"/h_timeResoVsEtaCut_pdg"+pdgID+".root").c_str(), "root");
  }

  TCanvas* ch_Rcut2 = new TCanvas();
  ch_Rcut2->cd();
  tg[0][2]->GetXaxis()->SetTitle("|#eta|");
  tg[0][2]->GetYaxis()->SetTitle("#sigma_{time} (ns)");
  tg[0][2]->GetYaxis()->SetRangeUser(0., 0.1);
  if(pdgID == "22") tg[0][2]->GetYaxis()->SetRangeUser(0., 0.05);
  tg[0][2]->Draw("apl");
  for(int iF=3; iF<nFiles; ++iF){
    tg[0][iF]->Draw("pl, same");
  } 
  for(int iT=1; iT<nOptions; ++iT){
    tg[iT][2]->Draw("pl, same");
    for(int iF=3; iF<nFiles; ++iF){
      tg[iT][iF]->Draw("pl, same");
    } 
  }
  legTGMShort->Draw("same");
  //legTGM->Draw("same");
  legTGMOpt->Draw("same");

  if(rad1 == 0 && pdgID == "130")  tL.DrawLatex(0.25, 0.85, "#rho #leq 2cm  K^{0}_{L}");
  if(rad1 == 1 && pdgID == "130")  tL.DrawLatex(0.25, 0.85, "#rho #leq 5cm  K^{0}_{L}");
  if(rad1 == 0 && pdgID == "22")  tL.DrawLatex(0.25, 0.85, "#rho #leq 2cm  #gamma");
  if(rad1 == 1 && pdgID == "22")  tL.DrawLatex(0.25, 0.85, "#rho #leq 5cm  #gamma");


  if(rad1 == 1){
  ch_Rcut2->Print((folder+"/h_timeResoVsEtaCut2_pdg"+pdgID+"_rad1.png").c_str(), "png");
  ch_Rcut2->Print((folder+"/h_timeResoVsEtaCut2_pdg"+pdgID+"_rad1.pdf").c_str(), "pdf");  
  ch_Rcut2->Print((folder+"/h_timeResoVsEtaCut2_pdg"+pdgID+"_rad1.root").c_str(), "root");
  }
  else{
  ch_Rcut2->Print((folder+"/h_timeResoVsEtaCut2_pdg"+pdgID+".png").c_str(), "png");
  ch_Rcut2->Print((folder+"/h_timeResoVsEtaCut2_pdg"+pdgID+".pdf").c_str(), "pdf");  
  ch_Rcut2->Print((folder+"/h_timeResoVsEtaCut2_pdg"+pdgID+".root").c_str(), "root");
  }

  TCanvas* ch_FrEvt = new TCanvas();
  ch_FrEvt->cd();
  tgFrEvt[0][0]->GetXaxis()->SetTitle("|#eta|");
  //  tgFrEvt[0][0]->GetYaxis()->SetTitle("<n. events with time/all>"); //#geq 3 hits with time / all events>");
  tgFrEvt[0][0]->GetYaxis()->SetTitle("efficiency"); //#geq 3 hits with time / all events>");
  tgFrEvt[0][0]->GetYaxis()->SetRangeUser(0., 1.1);
  //tgFrEvt[0][0]->GetYaxis()->SetRangeUser(0., 1.6);
  tgFrEvt[0][0]->Draw("apl");
  for(int iF=1; iF<nFiles; ++iF){
    tgFrEvt[0][iF]->Draw("pl, same");
  }
  for(int iT=1; iT<nOptions; ++iT){
    tgFrEvt[iT][0]->Draw("pl, same");
    for(int iF=1; iF<nFiles; ++iF){
      tgFrEvt[iT][iF]->Draw("pl, same");
      //      std::cout << " file = " << nameFiles.at(iF) << " options = " << nameOptions.at(iT) << std::endl; 
    }
  }
  legTGMD->Draw("same");
  //legTGM->Draw("same");
  legTGMOptD->Draw("same");
  
  if(rad1 == 0 && pdgID == "130")  tL.DrawLatex(0.25, 0.25, "#rho #leq 2cm  K^{0}_{L}");
  if(rad1 == 1 && pdgID == "130")  tL.DrawLatex(0.25, 0.25, "#rho #leq 5cm  K^{0}_{L}");
  if(rad1 == 0 && pdgID == "22")  tL.DrawLatex(0.25, 0.25, "#rho #leq 2cm  #gamma");
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
  
  TCanvas* ch_FrEvtCut = new TCanvas();
  ch_FrEvtCut->cd();
  tgFrEvt[0][2]->GetXaxis()->SetTitle("|#eta|");
  tgFrEvt[0][2]->GetYaxis()->SetTitle("efficiency"); //#geq 3 hits with time / all events>");
  tgFrEvt[0][2]->GetYaxis()->SetRangeUser(0., 1.1);
  //tgFrEvt[0][0]->GetYaxis()->SetRangeUser(0., 1.6);
  tgFrEvt[0][2]->Draw("apl");
  for(int iF=3; iF<nFiles; ++iF){
    tgFrEvt[0][iF]->Draw("pl, same");
  }
  for(int iT=1; iT<nOptions; ++iT){
    tgFrEvt[iT][2]->Draw("pl, same");
    for(int iF=3; iF<nFiles; ++iF){
      tgFrEvt[iT][iF]->Draw("pl, same");
      //      std::cout << " file = " << nameFiles.at(iF) << " options = " << nameOptions.at(iT) << std::endl; 
    }
  }
  legTGMShortD->Draw("same");
  //legTGM->Draw("same");
  legTGMOptD->Draw("same");

  if(rad1 == 0 && pdgID == "130")  tL.DrawLatex(0.25, 0.25, "#rho #leq 2cm  K^{0}_{L}");
  if(rad1 == 1 && pdgID == "130")  tL.DrawLatex(0.25, 0.25, "#rho #leq 5cm  K^{0}_{L}");
  if(rad1 == 0 && pdgID == "22")  tL.DrawLatex(0.25, 0.25, "#rho #leq 2cm  #gamma");
  if(rad1 == 1 && pdgID == "22")  tL.DrawLatex(0.25, 0.25, "#rho #leq 5cm  #gamma");


  if(rad1 == 1){
  ch_FrEvtCut->Print((folder+"/h_fractionEvtWithTimeCut"+pdgID+"_rad1.png").c_str(), "png");
  ch_FrEvtCut->Print((folder+"/h_fractionEvtWithTimeCut"+pdgID+"_rad1.pdf").c_str(), "pdf");  
  ch_FrEvtCut->Print((folder+"/h_fractionEvtWithTimeCut"+pdgID+"_rad1.root").c_str(), "root");
  }
  else{
  ch_FrEvtCut->Print((folder+"/h_fractionEvtWithTimeCut"+pdgID+".png").c_str(), "png");
  ch_FrEvtCut->Print((folder+"/h_fractionEvtWithTimeCut"+pdgID+".pdf").c_str(), "pdf");  
  ch_FrEvtCut->Print((folder+"/h_fractionEvtWithTimeCut"+pdgID+".root").c_str(), "root");
  }
  }  

  return;

}
