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
  int nBinsRad = 4;
  //int nBinsEta = 3;
  //float binWidth = 0.5;
  float binStart = 1.65;

  //std::string pdgID = "130";
  //std::string pdgID = "211";
  std::string pdgID = "22";

  int rad1 = 1;

  int nFiles = 7;
  std::vector<std::string> nameFiles;
  if(pdgID == "22"){
    nFiles = 3;
    nameFiles.push_back("PDG_"+pdgID+"_Pt2");
    nameFiles.push_back("PDG_"+pdgID+"_Pt5");
    nameFiles.push_back("PDG_"+pdgID+"_Pt60");
  }
  else{
    nameFiles.push_back("PDG_"+pdgID+"_Pt07");
    nameFiles.push_back("PDG_"+pdgID+"_Pt1");
    nameFiles.push_back("PDG_"+pdgID+"_Pt2");
    nameFiles.push_back("PDG_"+pdgID+"_Pt5");
    nameFiles.push_back("PDG_"+pdgID+"_Pt10");
    nameFiles.push_back("PDG_"+pdgID+"_Pt30");
    nameFiles.push_back("PDG_"+pdgID+"_Pt100");
  }


  int nOptions = 2;
  std::vector<std::string> nameOptions;
  nameOptions.push_back("CSF20LBA20");
  nameOptions.push_back("CSF20LEA20");


  TH1F* hDummy[7];
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
  
  TGraphErrors* tgM[8][7];
  TGraphErrors* tg[8][7];
  TGraphErrors* tgFrEvt[8][7];
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

  std::string folder = "combinationOption";  
  //  std::cout << "../test/timingStudies/ " << nameOptions.at(0) << " /JOB_ " << nameFiles.at(0) << " _3fC/OutTimeHGC_RecHits_ " << nameFiles.at(0) << " _3fC.root " << std::endl;

  std::cout << " >>> ora prendo i file " << std::endl;   
  TFile* inF[8][7];
  for(int iR=0; iR<nOptions; ++iR){
    for(int ij=0; ij<nFiles; ++ij){
      std::string inFileName =  "../../test/timingStudies/"+nameOptions.at(iR)+"/JOB_"+nameFiles.at(ij)+"_3fC/OutTimeHGC_RecHits_"+nameFiles.at(ij)+"_3fC.root";
      

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

      for(int ieta=0; ieta<nBinsEta; ++ieta){
	for(int iRad=0; iRad<nBinsRad; ++iRad){	
	  std::cout << " iF = " << iF  << " ieta = " << ieta << " irad = " << iRad << std::endl;
	  if((rad1 == 0  && iRad < 1) || (rad1 == 1 && iRad == 1)){
	if(doAllTheFits == true){
	  //hAverageTime_Eta_dRadius = (TH1F*)(inF[iT][iF]->Get(Form("ana/hAverageTime_Eta%.1f-%.1f_dRadius%d_Avg68", (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad)));
	  //hAverageTime_Eta_dRadius = (TH1F*)(inF[iT][iF]->Get(Form("ana/hAverageTime_Eta%.1f-%.1f_dRadius%d_AvgArm", (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad)));
	  //hAverageTime_Eta_dRadius = (TH1F*)(inF[iT][iF]->Get(Form("ana/hAverageTime_Eta%.2f-%.2f_dRadius%d", (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad)));
	  hAverageTime_Eta_dRadius = (TH1F*)(inF[iT][iF]->Get(Form("ana/hAverageTime_Eta%.2f-%.2f_dRadius%d_ResoWe", (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad)));
	  if(hAverageTime_Eta_dRadius->GetEntries() < 100) 	  hAverageTime_Eta_dRadius->Rebin(10);
	  else	  hAverageTime_Eta_dRadius->Rebin(5);
	  //	  tgRMS[iF][iRad]->SetPoint(ieta, (binStart+(2*ieta+1)*binWidth/2.),  hAverageTime_Eta_dRadius->GetRMS());
	  //tgMean[iF][iRad]->SetPoint(ieta, (binStart+(2*ieta+1)*binWidth/2.), hAverageTime_Eta_dRadius->GetMean());

	  if(hAverageTime_Eta_dRadius->GetRMS() < 0.01){
	    tg[iT][iF]->SetPoint(ieta, (binStart+(2*ieta+1)*binWidth/2.),  hAverageTime_Eta_dRadius->GetRMS());
	    tgM[iT][iF]->SetPoint(ieta, (binStart+(2*ieta+1)*binWidth/2.), hAverageTime_Eta_dRadius->GetMean());
	  }
	  else{	  
	    if(pdgID == "130" && iF == nFiles-1 && ieta == 2 && iRad == 0){
	      hfithisto2->SetRange(-0.05, 0.05);
	      hfithisto2->SetParameters(44, 0.03, 0.01);
	    }
	    if(pdgID == "130" && iF == nFiles-1 && ieta == 2 && iRad == 1){
	      hfithisto2->SetRange(-0.05, 0.05);
	      hfithisto2->SetParameters(44, 0.03, 0.01);
	    }
	    // if(pdgID == "130" && iF == nFiles-1 && ieta == nBinsEta-1 && iRad == 0){
	    //   hfithisto2->SetRange(-0.2, 0.2);
	    //   hfithisto2->SetParameters(32, 0.02, 0.01);
	    // }
	    // if(pdgID == "130" && iF == nFiles-1 && ieta == nBinsEta-1 && iRad == 0){
	    //   hfithisto2->SetRange(-0.2, 0.2);
	    //   hfithisto2->SetParameters(32, 0.02, 0.01);
	    // }    
	    else{
	      float YMax = 0.;
	      float xMaxVal = getXmax(hAverageTime_Eta_dRadius, YMax);
	      hfithisto->SetRange(xMaxVal - 2.*hAverageTime_Eta_dRadius->GetRMS(), xMaxVal + 2.*hAverageTime_Eta_dRadius->GetRMS());
	      hfithisto->SetParameter(1, xMaxVal);
	      hfithisto->SetParameter(2, hAverageTime_Eta_dRadius->GetRMS());
	      
	      hAverageTime_Eta_dRadius->Fit("hfithisto", "RQ");
	      hfithisto2->SetRange(xMaxVal - 2.*hfithisto->GetParameter(2), xMaxVal + 2.*hfithisto->GetParameter(2));
	      hfithisto2->SetParameters(YMax, hfithisto->GetParameter(1), hfithisto->GetParameter(2));
	    }

	    hAverageTime_Eta_dRadius->Fit("hfithisto2", "RQ");  
	    tg[iT][iF]->SetPoint(ieta, (binStart+(2*ieta+1)*binWidth/2.),  hfithisto2->GetParameter(2));                        
	    tg[iT][iF]->SetPointError(ieta, 0, hfithisto2->GetParError(2));     
	    
	    tgM[iT][iF]->SetPoint(ieta, (binStart+(2*ieta+1)*binWidth/2.),  hfithisto2->GetParameter(1));
	    tgM[iT][iF]->SetPointError(ieta, 0, hfithisto2->GetParError(1));
	  }
	  /*
	  c1 = new TCanvas();
	  c1->cd();
	  hAverageTime_Eta_dRadius->GetXaxis()->SetRangeUser(-0.05, 1.);
	  hAverageTime_Eta_dRadius->Draw();
	  //	hfithisto2->Draw("same");
	  c1->Print(Form((folderName+"/avgTime_Eta%.1f-%.1f_dRadius%d_"+nameFiles.at(iF)+".png").c_str(), (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "png");
	  c1->Delete();
	  hAverageTime_Eta_dRadius->Delete();
	  */
	}

	hFractionEvents_HitsWithTime_Eta_dRadius = (TH1F*)(inF[iT][iF]->Get(Form("ana/hFractionEvents_Hits_Eta%.2f-%.2f_dRadius%d", (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad)));
	tgFrEvt[iT][iF]->SetPoint(ieta, (binStart+(2*ieta+1)*binWidth/2.),  hFractionEvents_HitsWithTime_Eta_dRadius->GetMean());
	hFractionEvents_HitsWithTime_Eta_dRadius->Delete();

	}
	}//rad
      }//eta
    }//file
  }//type


  std::cout << " ci sono " << std::endl;

  TLegend *legTGM = new TLegend(0.80,0.75,0.90,0.95,NULL,"brNDC");
  legTGM->SetTextFont(42);
  legTGM->SetTextSize(0.02);
  legTGM->SetFillColor(kWhite);
  legTGM->SetLineColor(kWhite);
  legTGM->SetShadowColor(kWhite);
  for(int iF=0; iF<nFiles; ++iF){
    legTGM->AddEntry(hDummy[iF], (nameFiles.at(iF)).c_str(), "l");
  }
  TLegend *legTGMOpt = new TLegend(0.65,0.75,0.80,0.95,NULL,"brNDC");
  legTGMOpt->SetTextFont(42);
  legTGMOpt->SetTextSize(0.02);
  legTGMOpt->SetFillColor(kWhite);
  legTGMOpt->SetLineColor(kWhite);
  legTGMOpt->SetShadowColor(kWhite);
  for(int iF=0; iF<nOptions; ++iF){
    legTGMOpt->AddEntry(hDummySt[iF], (nameOptions.at(iF)).c_str(), "p");
  }
  
  

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  if(doAllTheFits == true){
  TCanvas* ch_M = new TCanvas();
  ch_M->cd();
  tgM[0][0]->GetXaxis()->SetTitle("#eta");
  tgM[0][0]->GetYaxis()->SetTitle("mean(t) ns");
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
  tg[0][0]->GetXaxis()->SetTitle("#eta");
  tg[0][0]->GetYaxis()->SetTitle("#sigma ns");
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
  tg[0][2]->GetXaxis()->SetTitle("#eta");
  tg[0][2]->GetYaxis()->SetTitle("#sigma ns");
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
  legTGM->Draw("same");
  legTGMOpt->Draw("same");
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

  TCanvas* ch_FrEvt = new TCanvas();
  ch_FrEvt->cd();
  tgFrEvt[0][0]->GetXaxis()->SetTitle("#eta");
  //  tgFrEvt[0][0]->GetYaxis()->SetTitle("<n. events with time/all>"); //#geq 3 hits with time / all events>");
  tgFrEvt[0][0]->GetYaxis()->SetTitle("fraction events (#geq 3 hits with time)"); //#geq 3 hits with time / all events>");
  tgFrEvt[0][0]->GetYaxis()->SetRangeUser(0., 1.1);
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
  //  legTGMd->Draw("same");
  //legTGM2o->Draw("same");
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
