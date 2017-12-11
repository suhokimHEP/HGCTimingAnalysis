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



void doPlots_fitTiming(){
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
  int nBinsRad = 4;
  //int nBinsEta = 3;
  //float binWidth = 0.5;
  float binStart = 1.65;
  
  std::string pdgID = "130";
  //  std::string pdgID = "22";
  //std::string pdgID = "211";


  int nFiles = 3;
  std::vector<std::string> nameFiles;
  std::vector<float> ptValues;
  nameFiles.push_back("PDG_"+pdgID+"_pt2");
  nameFiles.push_back("PDG_"+pdgID+"_pt5");
  nameFiles.push_back("PDG_"+pdgID+"_pt10");
  ptValues.push_back(2);
  ptValues.push_back(5);
  ptValues.push_back(10);
  TH1F* hDummy[7];
  TH1F* hDummySt[4];

  
  for(int ij=0; ij<nameFiles.size(); ++ij) std::cout << " name = " << nameFiles.at(ij) << std::endl;

  TGraphErrors* tgRMS[7][4];
  TGraphErrors* tgMean[7][4];
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
      
      //tg[iT][iR]->SetLineColor(iColors[3*iT+iR]);
      //tg[iT][iR]->SetMarkerColor(iColors[3*iT+iR]);

      tg[iT][iR]->SetLineColor(iColors[iT]);   
      tg[iT][iR]->SetMarkerColor(iColors[iT]); 
      tg[iT][iR]->SetLineWidth(2);
      tg[iT][iR]->SetMarkerStyle(iStyle[iR]);
      
      tgM[iT][iR] = new TGraphErrors();
      tgM[iT][iR]->SetName(Form("mean_%d_rad%d", iT, iR));
      tgM[iT][iR]->SetPoint(0, -1, -1);
      
      tgM[iT][iR]->SetLineColor(iColors[iT]);
      tgM[iT][iR]->SetLineWidth(2);
      tgM[iT][iR]->SetMarkerColor(iColors[iT]);
      tgM[iT][iR]->SetMarkerStyle(iStyle[iR]);

      /////////
      tgRMS[iT][iR] = new TGraphErrors();
      tgRMS[iT][iR]->SetName(Form("resolutionRMS_%d_rad%d", iT, iR));
      tgRMS[iT][iR]->SetPoint(0, -1, -1);
      
      tgRMS[iT][iR]->SetLineColor(iColors[iT]);   
      tgRMS[iT][iR]->SetMarkerColor(iColors[iT]); 
      tgRMS[iT][iR]->SetLineWidth(2);
      tgRMS[iT][iR]->SetMarkerStyle(iStyle[iR]);
      
      tgMean[iT][iR] = new TGraphErrors();
      tgMean[iT][iR]->SetName(Form("meanMEAN_%d_rad%d", iT, iR));
      tgMean[iT][iR]->SetPoint(0, -1, -1);
      
      tgMean[iT][iR]->SetLineColor(iColors[iT]);
      tgMean[iT][iR]->SetLineWidth(2);
      tgMean[iT][iR]->SetMarkerColor(iColors[iT]);
      tgMean[iT][iR]->SetMarkerStyle(iStyle[iR]);

      //////////
           
      tgFrEvt[iT][iR] = new TGraphErrors();
      tgFrEvt[iT][iR]->SetName(Form("fractionEvt_%d_rad%d", iT, iR));
      tgFrEvt[iT][iR]->SetPoint(0, -1, -1);
      
      tgFrEvt[iT][iR]->SetLineColor(iColors[iT]);
      tgFrEvt[iT][iR]->SetLineWidth(2);
      tgFrEvt[iT][iR]->SetMarkerColor(iColors[iT]);
      tgFrEvt[iT][iR]->SetMarkerStyle(iStyle[iR]);
  
    }
  }  

  std::cout << " >>> ora prendo i file " << std::endl;
  
  //  std::string optionType = "bkp310FixTime_double_0PU_thr12fC";
  std::string optionType = "ACfix_bkp310FixTime_double_0PU_thr12fC";

  std::string folderName = optionType+"/plotsAllTimeTot_CFD_"+pdgID;  
  std::string folder = optionType+"/plotsTimeTot_CFD_"+pdgID;       
  

  TFile* inF[3];
  for(int ij=0; ij<nFiles; ++ij){
    //inF[ij] = TFile::Open(("../test/testPU/"+optionType+"/"+nameFiles.at(ij)+"_0PU_lowEta/"+"JOB_"+nameFiles.at(ij)+"_0PU_lowEta/"+"/OutTimeHGC_RecHits_"+nameFiles.at(ij)+"_0PU_lowEta.root").c_str());
    inF[ij] = TFile::Open(("../test/testPU/"+optionType+"/"+"/OutTimeHGC_RecHits_"+nameFiles.at(ij)+"_0PU_allEta.root").c_str());
  }
  
  std::cout << " >>> fatto = presi " << std::endl;

  std::ofstream outFileLong[7];
  for(int ij=0; ij<nFiles; ++ij){
    outFileLong[ij].open((optionType+"/myTiming_CFD_"+nameFiles.at(ij)+".txt").c_str(), std::ios::out);
    std::cout << (optionType+"/myTiming_CFD_"+nameFiles.at(ij)+".txt").c_str() << std::endl; 
  }
  //  std::ofstream outFile("myTiming_CFD.txt", std::ios::app);



  TH1F* hFractionEvents_HitsWithTime_Eta_dRadius;
  TH1F* hAverageTime_Eta_dRadius;
  TF1* hfithisto = new TF1("hfithisto", "gaus", -0.5, 1.);
  TF1* hfithisto2 = new TF1("hfithisto2", "gaus", -0.5, 1.);
  hfithisto2->SetLineColor(kBlue);
  TCanvas* c1;

  for(int iF=0; iF<nFiles; ++iF){
    std::cout << " >>> iF = " << iF << std::endl;

    outFileLong[iF] << " eta " << " \t" << " rad " << " \t" << " res " <<" \t" << " eff " << std::endl;   

    for(int ieta=0; ieta<nBinsEta; ++ieta){
      if(ieta > 0 && ieta < nBinsEta-1) continue;
      for(int iRad=0; iRad<nBinsRad; ++iRad){
	
	std::cout << " iF = " << iF  << " ieta = " << ieta << " irad = " << iRad << std::endl;
	if(iRad < 2){
	hFractionEvents_HitsWithTime_Eta_dRadius = (TH1F*)(inF[iF]->Get(Form("ana/hFractionEvents_Hits_Eta%.2f-%.2f_dRadius%d", (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad)));
	float efficiency = hFractionEvents_HitsWithTime_Eta_dRadius->GetMean();
	tgFrEvt[iF][iRad]->SetPoint(ieta, (binStart+(2*ieta+1)*binWidth/2.),  efficiency);
	//	std::cout << " file  = " << iF << " eta BIn = " << (binStart+(2*ieta+1)*binWidth/2.) << " value = " << hFractionEvents_HitsWithTime_Eta_dRadius->GetMean() << std::endl;
	hFractionEvents_HitsWithTime_Eta_dRadius->Delete();


	if(doAllTheFits == true && iRad < 2){
	  //	  hAverageTime_Eta_dRadius = (TH1F*)(inF[iF]->Get(Form("ana/hAverageTime_Eta%.2f-%.2f_dRadius%d_Avg68NoRW", (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad)));
	  hAverageTime_Eta_dRadius = (TH1F*)(inF[iF]->Get(Form("ana/hAverageTime_Eta%.2f-%.2f_dRadius%d_AvgInt", (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad)));
	  //	  hAverageTime_Eta_dRadius = (TH1F*)(inF[iF]->Get(Form("ana/hAverageTime_Eta%.2f-%.2f_dRadius%d_Avg68Corr", (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad)));
	  //hAverageTime_Eta_dRadius = (TH1F*)(inF[iF]->Get(Form("ana/hAverageTime_Eta%.2f-%.2f_dRadius%d_AvgCutH", (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad)));
	  //hAverageTime_Eta_dRadius = (TH1F*)(inF[iF]->Get(Form("ana/hAverageTime_Eta%.2f-%.2f_dRadius%d_AvgCutHCorr", (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad)));
	  // if(hAverageTime_Eta_dRadius->GetEntries() < 100) 	  hAverageTime_Eta_dRadius->Rebin(10);
	  // else	  
	  hAverageTime_Eta_dRadius->Rebin(10);

	  tgFrEvt[iF][iRad]->SetPointError(ieta, 0.,  1./sqrt(hAverageTime_Eta_dRadius->GetEntries()) );

	  tgRMS[iF][iRad]->SetPoint(ieta, (binStart+(2*ieta+1)*binWidth/2.),  hAverageTime_Eta_dRadius->GetRMS());
	  tgMean[iF][iRad]->SetPoint(ieta, (binStart+(2*ieta+1)*binWidth/2.), hAverageTime_Eta_dRadius->GetMean());

	  /*
	  if(hAverageTime_Eta_dRadius->GetRMS() < 0.01){
	    tg[iF][iRad]->SetPoint(ieta, (binStart+(2*ieta+1)*binWidth/2.),  hAverageTime_Eta_dRadius->GetRMS());
	    tgM[iF][iRad]->SetPoint(ieta, (binStart+(2*ieta+1)*binWidth/2.), hAverageTime_Eta_dRadius->GetMean());
	    outFileLong[iF] << (binStart+(2*ieta+1)*binWidth/2.) << " \t " << iRad << " \t " << hAverageTime_Eta_dRadius->GetRMS() << " \t " << efficiency << std::endl;
	  }
	  else{	  
	    
	    if(pdgID == "130" && iF == nFiles-1 && ieta == 4 && iRad == 0){
	      hfithisto2->SetRange(-0.02, 0.02);
	      hfithisto2->SetParameters(80, 0.001, 0.003);
	    }
	    else if(pdgID == "130" && iF == nFiles-1 && ieta == 4 && iRad == 1){
              hfithisto2->SetRange(-0.05, 0.05);
	      hfithisto2->SetParameters(28, 0.01, 0.007);
	    }
	    else if(pdgID == "130" && iF == nFiles-1 && ieta == 5 && iRad == 0){
              hfithisto2->SetRange(-0.02, 0.02);
              hfithisto2->SetParameters(80, 0.001, 0.003);
            }
	    else if(pdgID == "130" && iF == nFiles-1 && ieta == 5 && iRad == 1){
              hfithisto2->SetRange(-0.05, 0.05);
	      hfithisto2->SetParameters(28, 0.01, 0.007);
	    }
	    else{
	  */
	      /*
	      float YMax = 0.;
	      float xMaxVal = getXmax(hAverageTime_Eta_dRadius, YMax);
	      //    hfithisto->SetRange(xMaxVal - 2.*hAverageTime_Eta_dRadius->GetRMS(), xMaxVal + 2.*hAverageTime_Eta_dRadius->GetRMS());
	      hfithisto->SetRange(hAverageTime_Eta_dRadius->GetMean() - 2.*hAverageTime_Eta_dRadius->GetRMS(), hAverageTime_Eta_dRadius->GetMean() + 2.*hAverageTime_Eta_dRadius->GetRMS());
	      hfithisto->SetParameter(1, xMaxVal);
	      hfithisto->SetParameter(2, hAverageTime_Eta_dRadius->GetRMS());
	      */

	      hfithisto->SetRange(hAverageTime_Eta_dRadius->GetMean() - 2.*hAverageTime_Eta_dRadius->GetRMS(), hAverageTime_Eta_dRadius->GetMean() + 2.*hAverageTime_Eta_dRadius->GetRMS());
	      hfithisto->SetParameter(1, hAverageTime_Eta_dRadius->GetMean());
	      hfithisto->SetParameter(2, hAverageTime_Eta_dRadius->GetRMS());

	      hAverageTime_Eta_dRadius->Fit("hfithisto", "RQ");
	      hfithisto2->SetRange(hfithisto->GetParameter(1) - 2.*hfithisto->GetParameter(2), hfithisto->GetParameter(1) + 2.*hfithisto->GetParameter(2));
	      hfithisto2->SetParameters(hfithisto->GetParameter(0), hfithisto->GetParameter(1), hfithisto->GetParameter(2));
	      // }

	    hAverageTime_Eta_dRadius->Fit("hfithisto2", "RQ");  
	    tg[iF][iRad]->SetPoint(ieta, (binStart+(2*ieta+1)*binWidth/2.),  hfithisto2->GetParameter(2));                        
	    tg[iF][iRad]->SetPointError(ieta, 0, hfithisto2->GetParError(2));     
	    
	    outFileLong[iF] << (binStart+(2*ieta+1)*binWidth/2.) << " \t " << iRad << " \t " << hfithisto2->GetParameter(2) << " \t " << efficiency << std::endl;

	    tgM[iF][iRad]->SetPoint(ieta, (binStart+(2*ieta+1)*binWidth/2.),  hfithisto2->GetParameter(1));
	    tgM[iF][iRad]->SetPointError(ieta, 0, hfithisto2->GetParError(1));
	    //	  }

	  c1 = new TCanvas();
	  c1->cd();
	  hAverageTime_Eta_dRadius->GetXaxis()->SetRangeUser(-0.5, 1.);
	  hAverageTime_Eta_dRadius->Draw();
	  //	hfithisto2->Draw("same");
	  c1->Print(Form((folderName+"/avgTime_Eta%.2f-%.2f_dRadius%d_"+nameFiles.at(iF)+".png").c_str(), (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "png");
	  c1->Delete();
	  hAverageTime_Eta_dRadius->Delete();
	}
	
	}

      }
    }
  }

  std::cout << " ci sono ora stampo " << std::endl;

  for(int iF=0; iF<nFiles; ++iF)
    outFileLong[iF].close();
  //  outFile.close();


  std::cout << " ci sono " << std::endl;

  TLegend *legTGM = new TLegend(0.65,0.70,0.85,0.85,NULL,"brNDC");
  legTGM->SetTextFont(42);
  legTGM->SetTextSize(0.04);
  legTGM->SetFillColor(kWhite);
  legTGM->SetLineColor(kWhite);
  legTGM->SetShadowColor(kWhite);
  for(int iF=0; iF<nFiles; ++iF){
    //legTGM->AddEntry(hDummy[iF], (nameFiles.at(iF)).c_str(), "l");
    legTGM->AddEntry(hDummy[iF], Form("p_{T} = %.1f GeV/c",ptValues.at(iF)), "l");
  }
  TLegend *legTGM2 = new TLegend(0.45,0.75,0.70,0.90,NULL,"brNDC");
  legTGM2->SetTextFont(42);
  legTGM2->SetTextSize(0.04);
  legTGM2->SetFillColor(kWhite);
  legTGM2->SetLineColor(kWhite);
  legTGM2->SetShadowColor(kWhite);
  legTGM2->AddEntry(hDummySt[0], " #rho < 2cm", "p");
  legTGM2->AddEntry(hDummySt[1], " #rho < 5cm", "p");
  legTGM2->AddEntry(hDummySt[2], " #rho < 10cm", "p");
  legTGM2->AddEntry(hDummySt[3], " all", "p");
  
  TLegend *legTGM2s = new TLegend(0.40,0.70,0.60,0.85,NULL,"brNDC");
  legTGM2s->SetTextFont(42);
  legTGM2s->SetTextSize(0.04);
  legTGM2s->SetFillColor(kWhite);
  legTGM2s->SetLineColor(kWhite);
  legTGM2s->SetShadowColor(kWhite);
  legTGM2s->AddEntry(hDummySt[0], " #rho < 2cm", "p");
  legTGM2s->AddEntry(hDummySt[1], " #rho < 5cm", "p");


  TLegend *legTGM2o = new TLegend(0.40,0.70,0.60,0.85,NULL,"brNDC");
  legTGM2o->SetTextFont(42);
  legTGM2o->SetTextSize(0.04);
  legTGM2o->SetFillColor(kWhite);
  legTGM2o->SetLineColor(kWhite);
  legTGM2o->SetShadowColor(kWhite);
  legTGM2o->AddEntry(hDummySt[0], " #rho < 2cm", "p");

  ////////// down
  TLegend *legTGMd = new TLegend(0.65,0.15,0.85,0.30,NULL,"brNDC");
  legTGMd->SetTextFont(42);
  legTGMd->SetTextSize(0.04);
  legTGMd->SetFillColor(kWhite);
  legTGMd->SetLineColor(kWhite);
  legTGMd->SetShadowColor(kWhite);
  for(int iF=0; iF<nFiles; ++iF){
    //legTGMd->AddEntry(hDummy[iF], (nameFiles.at(iF)).c_str(), "l");
    legTGMd->AddEntry(hDummy[iF], Form("p_{T} = %.1f GeV/c",ptValues.at(iF)), "l");
  }
  TLegend *legTGM2d = new TLegend(0.40,0.15,0.60,0.30,NULL,"brNDC");
  //TLegend *legTGM2d = new TLegend(0.65,0.15,0.80,0.35,NULL,"brNDC");
  legTGM2d->SetTextFont(42);
  legTGM2d->SetTextSize(0.04);
  legTGM2d->SetFillColor(kWhite);
  legTGM2d->SetLineColor(kWhite);
  legTGM2d->SetShadowColor(kWhite);
  legTGM2d->AddEntry(hDummySt[0], " #rho < 2cm", "p");
  legTGM2d->AddEntry(hDummySt[1], " #rho < 5cm", "p");
  //  legTGM2d->AddEntry(hDummySt[2], " R < 10cm", "p");
  //  legTGM2d->AddEntry(hDummySt[3], " all", "p");
  

  std::cout << " legends ok  " << std::endl;


  TLatex tL;
  tL.SetNDC();
  tL.SetTextSize(0.05);
  tL.SetTextFont(132);

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  if(doAllTheFits == true){
  TCanvas* ch_R = new TCanvas();
  ch_R->cd();
  tg[0][0]->GetXaxis()->SetTitle("|#eta|");
  tg[0][0]->GetYaxis()->SetTitle("#sigma_{t} ns");
  //  tg[0][0]->GetYaxis()->SetTitle("RMS ns");
  tg[0][0]->GetYaxis()->SetRangeUser(0., 0.1);
  tg[0][0]->Draw("apl");
  tg[0][1]->Draw("pl, same");
  //  tg[0][2]->Draw("pl, same");
  // tg[0][3]->Draw("pl, same");
  for(int iF=1; iF<nFiles; ++iF){
    // if((pdgID == "130" || pdgID == "211") &&
    //    iF >= nFiles - 3) continue;
    tg[iF][0]->Draw("pl, same");
    tg[iF][1]->Draw("pl, same");
    //    tg[iF][2]->Draw("pl, same");
    //    tg[iF][3]->Draw("pl, same");
  }


  if(pdgID == "130")  tL.DrawLatex(0.25, 0.85, "K^{0}_{L}");
  if(pdgID == "22")   tL.DrawLatex(0.25, 0.85, "#gamma");

  legTGM->Draw("same");
  legTGM2s->Draw("same");
  //  legTGM2o->Draw("same");
  ch_R->Print((folder+"/h_timeResoVsEta.png").c_str(), "png");
  ch_R->Print((folder+"/h_timeResoVsEta.pdf").c_str(), "pdf");  
  ch_R->Print((folder+"/h_timeResoVsEta.root").c_str(), "root");
  tg[0][0]->GetYaxis()->SetRangeUser(0., 0.2);
  ch_R->Print((folder+"/h_timeResoVsEta_zoomIn.png").c_str(), "png");
  ch_R->Print((folder+"/h_timeResoVsEta_zoomIn.pdf").c_str(), "pdf");
  ch_R->Print((folder+"/h_timeResoVsEta_zoomIn.root").c_str(), "root");

  TCanvas* ch_M = new TCanvas();
  ch_M->cd();
  tgM[0][0]->GetXaxis()->SetTitle("|#eta|");
  tgM[0][0]->GetYaxis()->SetTitle("mean(t) ns");
  //RAFIXME  tgM[0][0]->GetYaxis()->SetRangeUser(0., 0.6);
  tgM[0][0]->GetYaxis()->SetRangeUser(-0.02, 0.6);
  tgM[0][0]->Draw("apl");
  tgM[0][1]->Draw("pl, same");
  // tgM[0][2]->Draw("pl, same");
  // tgM[0][3]->Draw("pl, same");
  for(int iF=1; iF<nFiles; ++iF){
    tgM[iF][0]->Draw("pl, same");
    tgM[iF][1]->Draw("pl, same");
    // tgM[iF][2]->Draw("pl, same");
    // tgM[iF][3]->Draw("pl, same");
  }


  if(pdgID == "130")  tL.DrawLatex(0.25, 0.85, "K^{0}_{L}");
  if(pdgID == "22")   tL.DrawLatex(0.25, 0.85, "#gamma");


  legTGM->Draw("same");
  //legTGM2o->Draw("same");
  legTGM2s->Draw("same");
  ch_M->Print((folder+"/h_timeMeanVsEta.png").c_str(), "png");
  ch_M->Print((folder+"/h_timeMeanVsEta.pdf").c_str(), "pdf");  
  ch_M->Print((folder+"/h_timeMeanVsEta.root").c_str(), "root");

  /////////////
  TCanvas* ch_Rms = new TCanvas();
  ch_Rms->cd();
  tgRMS[0][0]->GetXaxis()->SetTitle("|#eta|");
  tgRMS[0][0]->GetYaxis()->SetTitle("RMS ns");
  tgRMS[0][0]->GetYaxis()->SetRangeUser(0., 0.1);
  tgRMS[0][0]->Draw("apl");
  for(int iF=1; iF<nFiles; ++iF){
    tg[iF][0]->Draw("pl, same");
    //tg[iF][1]->Draw("pl, same");
    //    tg[iF][2]->Draw("pl, same");
    //    tg[iF][3]->Draw("pl, same");
  }

  if(pdgID == "130")  tL.DrawLatex(0.25, 0.85, "K^{0}_{L}");
  if(pdgID == "22")   tL.DrawLatex(0.25, 0.85, "#gamma");


  legTGM->Draw("same");
  legTGM2o->Draw("same");
  ch_Rms->Print((folder+"/h_timeRMSVsEta.png").c_str(), "png");
  ch_Rms->Print((folder+"/h_timeRMSVsEta.pdf").c_str(), "pdf");  
  ch_Rms->Print((folder+"/h_timeRMSVsEta.root").c_str(), "root");
  tgRMS[0][0]->GetYaxis()->SetRangeUser(0., 0.2);
  ch_Rms->Print((folder+"/h_timeRMSVsEta_zoomIn.png").c_str(), "png");
  ch_Rms->Print((folder+"/h_timeRMSVsEta_zoomIn.pdf").c_str(), "pdf");
  ch_Rms->Print((folder+"/h_timeRMSVsEta_zoomIn.root").c_str(), "root");

  TCanvas* ch_Mean = new TCanvas();
  ch_Mean->cd();
  tgMean[0][0]->GetXaxis()->SetTitle("|#eta|");
  tgMean[0][0]->GetYaxis()->SetTitle("mean(t) ns");
  tgMean[0][0]->GetYaxis()->SetRangeUser(-0.02, 0.6);
  tgMean[0][0]->Draw("apl");
  // tgM[0][1]->Draw("pl, same");
  // tgM[0][2]->Draw("pl, same");
  // tgM[0][3]->Draw("pl, same");
  for(int iF=1; iF<nFiles; ++iF){
    tgM[iF][0]->Draw("pl, same");
    // tgM[iF][1]->Draw("pl, same");
    // tgM[iF][2]->Draw("pl, same");
    // tgM[iF][3]->Draw("pl, same");
  }

  if(pdgID == "130")  tL.DrawLatex(0.25, 0.85, "#rho #leq 2cm  K^{0}_{L}");
  if(pdgID == "22")   tL.DrawLatex(0.25, 0.85, "#rho #leq 2cm  #gamma");


  legTGM->Draw("same");
  //  legTGM2o->Draw("same");
  ch_Mean->Print((folder+"/h_timeMEANVsEta.png").c_str(), "png");
  ch_Mean->Print((folder+"/h_timeMEANVsEta.pdf").c_str(), "pdf");  
  ch_Mean->Print((folder+"/h_timeMEANVsEta.root").c_str(), "root");

  ///solo per rho == 2cm
  TCanvas* ch_R2cm = new TCanvas();
  ch_R2cm->cd();
  tg[0][0]->GetXaxis()->SetTitle("|#eta|");
  tg[0][0]->GetYaxis()->SetTitle("#sigma_{t} ns");
  tg[0][0]->GetYaxis()->SetRangeUser(0., 0.1);
  tg[0][0]->Draw("apl");
  for(int iF=1; iF<nFiles; ++iF){
    tg[iF][0]->Draw("pl, same");
  }

  if(pdgID == "130")  tL.DrawLatex(0.25, 0.85, "#rho #leq 2cm  K^{0}_{L}");
  if(pdgID == "22")   tL.DrawLatex(0.25, 0.85, "#rho #leq 2cm  #gamma");


  legTGM->Draw("same");
  //  legTGM2o->Draw("same");
  ch_R2cm->Print((folder+"/h_timeResoVsEta_2cm.png").c_str(), "png");
  ch_R2cm->Print((folder+"/h_timeResoVsEta_2cm.pdf").c_str(), "pdf");  
  ch_R2cm->Print((folder+"/h_timeResoVsEta_2cm.root").c_str(), "root");
  tg[0][0]->GetYaxis()->SetRangeUser(0., 0.2);
  ch_R2cm->Print((folder+"/h_timeResoVsEta_zoomIn_2cm.png").c_str(), "png");
  ch_R2cm->Print((folder+"/h_timeResoVsEta_zoomIn_2cm.pdf").c_str(), "pdf");
  ch_R2cm->Print((folder+"/h_timeResoVsEta_zoomIn_2cm.root").c_str(), "root");


  TCanvas* ch_M2cm = new TCanvas();
  ch_M2cm->cd();
  tgM[0][0]->GetXaxis()->SetTitle("|#eta|");
  tgM[0][0]->GetYaxis()->SetTitle("mean(t) ns");
  tgM[0][0]->GetYaxis()->SetRangeUser(-0.02, 0.6);
  tgM[0][0]->Draw("apl");
  for(int iF=1; iF<nFiles; ++iF){
    tgM[iF][0]->Draw("pl, same");
  }


  if(pdgID == "130")  tL.DrawLatex(0.25, 0.85, "#rho #leq 2cm  K^{0}_{L}");
  if(pdgID == "22")   tL.DrawLatex(0.25, 0.85, "#rho #leq 2cm  #gamma");


  legTGM->Draw("same");
  //  legTGM2o->Draw("same");
  ch_M2cm->Print((folder+"/h_timeMeanVsEta_2cm.png").c_str(), "png");
  ch_M2cm->Print((folder+"/h_timeMeanVsEta_2cm.pdf").c_str(), "pdf");  
  ch_M2cm->Print((folder+"/h_timeMeanVsEta_2cm.root").c_str(), "root");



  }

  std::cout <<  " ok ora riprendo " << std::endl;
  TCanvas* ch_FrEvt = new TCanvas();
  ch_FrEvt->cd();
  tgFrEvt[0][0]->GetXaxis()->SetTitle("|#eta|");
  tgFrEvt[0][0]->GetYaxis()->SetTitle("efficiency"); //#geq 3 hits with time / all events>");
  tgFrEvt[0][0]->GetYaxis()->SetRangeUser(0., 1.1);
  tgFrEvt[0][0]->Draw("apl");
  tgFrEvt[0][1]->Draw("pl, same");
  tgFrEvt[0][2]->Draw("pl, same");
  tgFrEvt[0][3]->Draw("pl, same");
  for(int iF=1; iF<nFiles; ++iF){
    tgFrEvt[iF][0]->Draw("pl, same");
    tgFrEvt[iF][1]->Draw("pl, same");
    tgFrEvt[iF][2]->Draw("pl, same");
    tgFrEvt[iF][3]->Draw("pl, same");
  }
  std::cout << " fin qui ok " << std::endl;
  //return;

  if(pdgID == "130")  tL.DrawLatex(0.25, 0.25, "K^{0}_{L}");
  if(pdgID == "22")   tL.DrawLatex(0.25, 0.25, "#gamma");


  legTGMd->Draw("same");
  legTGM2d->Draw("same");

  ch_FrEvt->Print((folder+"/h_fractionEvtWithTime.png").c_str(), "png");
  ch_FrEvt->Print((folder+"/h_fractionEvtWithTime.pdf").c_str(), "pdf");  
  ch_FrEvt->Print((folder+"/h_fractionEvtWithTime.root").c_str(), "root");
  // ch_FrEvt2->Print((folder+"/h_Pippo.png").c_str(), "png");        
  // ch_FrEvt2->Print((folder+"/h_Pippo.pdf").c_str(), "pdf");        
  // ch_FrEvt2->Print((folder+"/h_Pippo.root").c_str(), "root");      


  TCanvas* ch_FrEvt2cm = new TCanvas();
  ch_FrEvt2cm->cd();
  tgFrEvt[0][0]->GetXaxis()->SetTitle("|#eta|");
  tgFrEvt[0][0]->GetYaxis()->SetTitle("efficiency"); //#geq 3 hits with time / all events>");
  tgFrEvt[0][0]->GetYaxis()->SetRangeUser(0., 1.1);
  tgFrEvt[0][0]->Draw("apl");
  for(int iF=1; iF<nFiles; ++iF){
    tgFrEvt[iF][0]->Draw("pl, same");
  }
  std::cout << " fin qui ok " << std::endl;
  //return;

  if(pdgID == "130")  tL.DrawLatex(0.25, 0.25, "#rho #leq 2cm  K^{0}_{L}");
  if(pdgID == "22")   tL.DrawLatex(0.25, 0.25, "#rho #leq 2cm  #gamma");


  legTGMd->Draw("same");
  //legTGM2d->Draw("same");

  ch_FrEvt2cm->Print((folder+"/h_fractionEvtWithTime_2cm.png").c_str(), "png");
  ch_FrEvt2cm->Print((folder+"/h_fractionEvtWithTime_2cm.pdf").c_str(), "pdf");  
  ch_FrEvt2cm->Print((folder+"/h_fractionEvtWithTime_2cm.root").c_str(), "root");
  // ch_FrEvt2->Print((folder+"/h_Pippo.png").c_str(), "png");        
  // ch_FrEvt2->Print((folder+"/h_Pippo.pdf").c_str(), "pdf");        
  // ch_FrEvt2->Print((folder+"/h_Pippo.root").c_str(), "root");      


 
  std::cout << " qui finito " << std::endl;

  return;



}
