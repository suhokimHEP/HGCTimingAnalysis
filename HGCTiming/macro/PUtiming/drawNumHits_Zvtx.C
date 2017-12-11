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



void drawNumHits_Zvtx(){
  gROOT->Macro("/afs/cern.ch/user/a/amartell/public/setStyle.C");
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
  int iStyle[4] = {20, 21, 22, 23}; //kGray+1};
  //  int iColors[6] = {kGreen+1, kBlue, kRed};
  
  int nBinsEta = 6;
  float binWidth = 0.2;
  int nBinsRad = 4;
  //int nBinsEta = 3;
  //float binWidth = 0.5;
  float binStart = 1.65;
  
  //std::string pdgID = "130";
  std::string pdgID = "22";
  //std::string pdgID = "211";

    std::vector<int> PUvalue;
    PUvalue.push_back(-10);
    PUvalue.push_back(+10);
    PUvalue.push_back(0);

    std::vector<std::string> PUname;
    PUname.push_back("Min10");
    PUname.push_back("Plu10");
    PUname.push_back("At0");


  std::string radius = "0";
  std::vector<int> radiusR;
  radiusR.push_back(2);
  radiusR.push_back(5);

  std::vector<std::string> etaBin;
  std::vector<float> etaValue;
  etaBin.push_back("1.65-1.85");
  etaBin.push_back("2.65-2.85");
  etaValue.push_back(1.75);
  etaValue.push_back(2.75);


  TFile* inF[3];
  inF[0] = TFile::Open(("../../test/testPU/bkp310FixTime_double_0PU_Vtx_thr12fC/PDG_"+pdgID+"_Pt10_0PU_thr12_allEta_ZvtsMin10/JOB_PDG_"+pdgID+"_Pt10_0PU_thr12_allEta_ZvtsMin10/OutTimeHGC_RecHits_PDG_"+pdgID+"_Pt10_0PU_thr12_allEta_ZvtsMin10.root").c_str());
  inF[1] = TFile::Open(("../../test/testPU/bkp310FixTime_double_0PU_Vtx_thr12fC/PDG_"+pdgID+"_Pt10_0PU_thr12_allEta_ZvtsPlu10/JOB_PDG_"+pdgID+"_Pt10_0PU_thr12_allEta_ZvtsPlu10/OutTimeHGC_RecHits_PDG_"+pdgID+"_Pt10_0PU_thr12_allEta_ZvtsPlu10.root").c_str());
  //inF[1] = TFile::Open(("../../test/testPU/bkp310FixTime_double_0PU_Vtx_thr12fC/PDG_"+pdgID+"_Pt10_0PU_thr12_allEta_ZvtsMin10/JOB_PDG_"+pdgID+"_Pt10_0PU_thr12_allEta_ZvtsMin10/OutTimeHGC_RecHits_PDG_"+pdgID+"_Pt10_0PU_thr12_allEta_ZvtsMin10.root").c_str());
  inF[2] = TFile::Open(("../../test/testPU/bkp310FixTime_double_0PU_thr12fC/PDG_"+pdgID+"_Pt5_0PU_thr12_allEta/JOB_PDG_"+pdgID+"_Pt5_0PU_thr12_allEta/OutTimeHGC_RecHits_PDG_"+pdgID+"_Pt5_0PU_thr12_allEta.root").c_str());


  for(int ij=0; ij<3; ++ij) std::cout << " file name = " << inF[ij]->GetName() << std::endl;

  std::cout << " >>> fatto = presi " << std::endl;


  TH1F* hDummy[7];
  TH1F* hDummySt[4];

  
  //  for(int ij=0; ij<nameFiles.size(); ++ij) std::cout << " name = " << nameFiles.at(ij) << std::endl;

  std::cout << " >>> ora prendo i file " << std::endl;


  std::string folder = "plotsPU_VtxPos_"+pdgID;


  TH1F* hTime_Eta_dRadius[2][3];
  TH1F* hTimeCP_Eta_dRadius[2][3];
  TH1F* hTimeCut_Eta_dRadius[2][3];
  TH1F* hTimeIntCut_Eta_dRadius[2][3];

  TH1F* hAverageTime_Eta_dRadius[2][3];
  TH1F* hAverageTimeCP_Eta_dRadius[2][3];
  TH1F* hAverageTime_Eta_dRadius_Avg68[2][3];
  TH1F* hAverageTime_Eta_dRadius_AvgInt[2][3];

  TGraphErrors* tgFrEvt[3];
  TGraphErrors* tgReso[3];

  for(int iF=0; iF<3; ++iF){
    tgFrEvt[iF] = new TGraphErrors();
    tgFrEvt[iF]->SetName(Form("tgFrEvt_PU%d", PUvalue.at(iF)) );
    tgFrEvt[iF]->SetPoint(0, -1., -1.);
    tgFrEvt[iF]->SetMarkerColor(iColors[iF]);
    tgFrEvt[iF]->SetMarkerStyle(iStyle[iF]);
    tgFrEvt[iF]->SetMarkerSize(1.1);

    tgReso[iF] = new TGraphErrors();
    tgReso[iF]->SetName(Form("tgReso_PU%d", PUvalue.at(iF)) );
    tgReso[iF]->SetPoint(0, -1., -1.);
    tgReso[iF]->SetMarkerColor(iColors[iF]);
    tgReso[iF]->SetMarkerStyle(iStyle[iF]);

    for(int ij=0; ij<2; ++ij){
      TH1F* test  = (TH1F*)inF[iF]->Get(("ana/hFractionEvents_Hits_Eta"+etaBin.at(ij)+"_dRadius"+radius).c_str());

      std::cout << "prendo da  file name = " << inF[iF]->GetName() 
		<< " test name = " << test->GetName() << " test mean = " << test->GetMean() << std::endl;

      if(ij == 0) tgFrEvt[iF]->SetPoint(ij+1, 1.75, test->GetMean());
      if(ij == 1) tgFrEvt[iF]->SetPoint(ij+1, 2.75, test->GetMean());

      //      test2->Draw();

      test->Delete();
      
      hTime_Eta_dRadius[ij][iF] = (TH1F*)inF[iF]->Get(("ana/hTime_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius).c_str());
      hTimeCP_Eta_dRadius[ij][iF] = (TH1F*)inF[iF]->Get(("ana/hTimeCP_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius).c_str());
      hTimeCut_Eta_dRadius[ij][iF] = (TH1F*)inF[iF]->Get(("ana/hTimeCut_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius).c_str());
      hTimeIntCut_Eta_dRadius[ij][iF] = (TH1F*)inF[iF]->Get(("ana/hTimeIntCut_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius).c_str());

      hAverageTime_Eta_dRadius[ij][iF] = (TH1F*)(inF[iF]->Get(("ana/hAverageTime_Eta"+etaBin.at(ij)+"_dRadius"+radius+"").c_str()) );
      hAverageTimeCP_Eta_dRadius[ij][iF] = (TH1F*)(inF[iF]->Get(("ana/hAverageTimeCP_Eta"+etaBin.at(ij)+"_dRadius"+radius+"").c_str()) );
      hAverageTime_Eta_dRadius_Avg68[ij][iF] = (TH1F*)(inF[iF]->Get(("ana/hAverageTime_Eta"+etaBin.at(ij)+"_dRadius"+radius+"_Avg68NoRW").c_str()) );
      hAverageTime_Eta_dRadius_AvgInt[ij][iF] = (TH1F*)(inF[iF]->Get(("ana/hAverageTime_Eta"+etaBin.at(ij)+"_dRadius"+radius+"_AvgInt").c_str()) );

      hTime_Eta_dRadius[ij][iF]->SetLineColor(kBlack);
      hTimeCP_Eta_dRadius[ij][iF]->SetLineColor(kBlue);
      hTimeCut_Eta_dRadius[ij][iF]->SetLineColor(kGreen+2);
      hTimeIntCut_Eta_dRadius[ij][iF]->SetLineColor(kViolet+1);

      hAverageTime_Eta_dRadius[ij][iF]->SetLineColor(kBlack);
      hAverageTimeCP_Eta_dRadius[ij][iF]->SetLineColor(kBlue);
      hAverageTime_Eta_dRadius_Avg68[ij][iF]->SetLineColor(kGreen+2);
      hAverageTime_Eta_dRadius_AvgInt[ij][iF]->SetLineColor(kViolet+1);

      hTime_Eta_dRadius[ij][iF]->SetLineWidth(2);
      hTimeCP_Eta_dRadius[ij][iF]->SetLineWidth(2);
      hTimeCut_Eta_dRadius[ij][iF]->SetLineWidth(2);
      hTimeIntCut_Eta_dRadius[ij][iF]->SetLineWidth(2);

      hAverageTime_Eta_dRadius[ij][iF]->SetLineWidth(2);
      hAverageTimeCP_Eta_dRadius[ij][iF]->SetLineWidth(2);
      hAverageTime_Eta_dRadius_Avg68[ij][iF]->SetLineWidth(2);
      hAverageTime_Eta_dRadius_AvgInt[ij][iF]->SetLineWidth(2);
    }
  }

  for(int iF=0; iF<3; ++iF){
    tgFrEvt[iF]->SetPoint(3, 10, 10);
  }

  std::cout << " ci sono ora stampo " << std::endl;


  TLegend *legTGM = new TLegend(0.70,0.70,0.80,0.93,NULL,"brNDC");
  legTGM->SetTextFont(42);
  legTGM->SetTextSize(0.04);
  legTGM->SetFillColor(kWhite);
  legTGM->SetLineColor(kWhite);
  legTGM->SetShadowColor(kWhite);
  //  for(int iF=0; iF<nFiles; ++iF){
  legTGM->AddEntry(hTime_Eta_dRadius[0][0], "all", "l");
  legTGM->AddEntry(hTimeCP_Eta_dRadius[0][0], "all signal", "l");
  legTGM->AddEntry(hTimeCut_Eta_dRadius[0][0], "selected 68%",  "l");
  legTGM->AddEntry(hTimeIntCut_Eta_dRadius[0][0], "selected 300ps", "l");
  //}

  std::cout << " legends ok  " << std::endl;

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TLatex tL;
  tL.SetNDC();
  tL.SetTextSize(0.04);
  tL.SetTextFont(132);


  TLatex tLPU;
  tLPU.SetNDC();
  tLPU.SetTextSize(0.04);
  tLPU.SetTextFont(132);


  TCanvas* ch_hitsTime[2][3];
  for(int iP=0; iP<3; ++iP){
  for(int ij=0; ij<2; ++ij){
    ch_hitsTime[ij][iP] = new TCanvas();
    ch_hitsTime[ij][iP]->cd();
  // gPad->SetLogx();
    gPad->SetLogy();

    hTime_Eta_dRadius[ij][iP]->GetXaxis()->SetTitle("recHits time (ns)");
    hTime_Eta_dRadius[ij][iP]->GetXaxis()->SetRangeUser(-3., 22.);
    hTime_Eta_dRadius[ij][iP]->GetYaxis()->SetRangeUser(0.1, 1.e+5);
    hTime_Eta_dRadius[ij][iP]->Draw("h");

    hTimeCP_Eta_dRadius[ij][iP]->Draw("h, same");
    hTimeCut_Eta_dRadius[ij][iP]->Draw("h, same");
    hTimeIntCut_Eta_dRadius[ij][iP]->Draw("h, same");

    if(pdgID == "130") tL.DrawLatex(0.25,0.90, Form("K^{0}_{L}  #rho #leq%dcm  |#eta| = %.2f", radiusR.at(std::stoi(radius)), etaValue.at(ij) ) );
    if(pdgID == "22")  tL.DrawLatex(0.25,0.90, Form("#gamma     #rho #leq%dcm  |#eta| = %.2f", radiusR.at(std::stoi(radius)), etaValue.at(ij)  ) );

    tLPU.DrawLatex(0.25,0.85, Form("vtx_{Z} = %dcm   p_{T} = 10GeV/c", PUvalue.at(iP)) );

    legTGM->Draw("same");
    ch_hitsTime[ij][iP]->Print((folder+"/h_hitsTime_pdg_"+pdgID+"_rad_"+radius+"_eta_"+etaBin.at(ij)+"Vtx_"+PUname.at(iP)+".png").c_str(), "png");
    ch_hitsTime[ij][iP]->Print((folder+"/h_hitsTime_pdg_"+pdgID+"_rad_"+radius+"_eta_"+etaBin.at(ij)+"Vtx_"+PUname.at(iP)+".pdf").c_str(), "pdf");
    ch_hitsTime[ij][iP]->Print((folder+"/h_hitsTime_pdg_"+pdgID+"_rad_"+radius+"_eta_"+etaBin.at(ij)+"Vtx_"+PUname.at(iP)+".root").c_str(), "root");

    hTime_Eta_dRadius[ij][iP]->GetXaxis()->SetRangeUser(-3., 5.);
    ch_hitsTime[ij][iP]->Print((folder+"/h_hitsTime_pdg_"+pdgID+"_rad_"+radius+"_eta_"+etaBin.at(ij)+"Vtx_"+PUname.at(iP)+"zoom.png").c_str(), "png");
    ch_hitsTime[ij][iP]->Print((folder+"/h_hitsTime_pdg_"+pdgID+"_rad_"+radius+"_eta_"+etaBin.at(ij)+"Vtx_"+PUname.at(iP)+"zoom.pdf").c_str(), "pdf");
    ch_hitsTime[ij][iP]->Print((folder+"/h_hitsTime_pdg_"+pdgID+"_rad_"+radius+"_eta_"+etaBin.at(ij)+"Vtx_"+PUname.at(iP)+"zoom.root").c_str(), "root");

  }
  }


  TLegend *legTGMresoAll = new TLegend(0.70,0.70,0.80,0.93,NULL,"brNDC");
  legTGMresoAll->SetTextFont(42);
  legTGMresoAll->SetTextSize(0.04);
  legTGMresoAll->SetFillColor(kWhite);
  legTGMresoAll->SetLineColor(kWhite);
  legTGMresoAll->SetShadowColor(kWhite);
  //  for(int iF=0; iF<nFiles; ++iF){
  legTGMresoAll->AddEntry(hAverageTime_Eta_dRadius[0][0], "all", "l");
  legTGMresoAll->AddEntry(hAverageTimeCP_Eta_dRadius[0][0], "all signal", "l");
  // legTGMresoAll->AddEntry(hAverageTime_Eta_dRadius_Avg68[0], "selected 68% (140PU)", "l");
  // legTGMresoAll->AddEntry(hAverageTime_Eta_dRadius_AvgInt[0], "selected 300ps (140PU)", "l");
  legTGMresoAll->AddEntry(hAverageTime_Eta_dRadius_Avg68[0][0], "signal 68%", "l");
  legTGMresoAll->AddEntry(hAverageTime_Eta_dRadius_AvgInt[0][0], "signal 300ps", "l");
  //  legTGM->AddEntry(timeCut0PU[0], "signal (0PU) new sim", "l");
  //}

  TLegend *legTGMresoFit = new TLegend(0.70,0.70,0.80,0.93,NULL,"brNDC");
  legTGMresoFit->SetTextFont(42);
  legTGMresoFit->SetTextSize(0.04);
  legTGMresoFit->SetFillColor(kWhite);
  legTGMresoFit->SetLineColor(kWhite);
  legTGMresoFit->SetShadowColor(kWhite);
  //  for(int iF=0; iF<nFiles; ++iF){
  //legTGMresoFit->AddEntry(hAverageTime_Eta_dRadius[0], "all", "l");
  //legTGMresoFit->AddEntry(hAverageTimeCP_Eta_dRadius[0], "all signal (140PU)", "l");
  legTGMresoFit->AddEntry(hAverageTime_Eta_dRadius_Avg68[0][0], "selected 68%",  "l");
  legTGMresoFit->AddEntry(hAverageTime_Eta_dRadius_AvgInt[0][0], "selected 300ps",  "l");
  //legTGMresoFit->AddEntry(hAverageTime_Eta_dRadius_Avg680PU[0], "signal (0PU) old sim", "l");
  //  legTGM->AddEntry(timeCut0PU[0], "signal (0PU) new sim", "l");
  //}



  TCanvas* ch_Reso_all[2][3];
  for(int iP=0; iP<3; ++iP){
  for(int ij=0; ij<2; ++ij){
    ch_Reso_all[ij][iP] = new TCanvas();
    ch_Reso_all[ij][iP]->cd();
  // gPad->SetLogx();
    gPad->SetLogy();

    hAverageTime_Eta_dRadius[ij][iP]->Rebin(5);
    hAverageTimeCP_Eta_dRadius[ij][iP]->Rebin(5);
    hAverageTime_Eta_dRadius_Avg68[ij][iP]->Rebin(5);
    hAverageTime_Eta_dRadius_AvgInt[ij][iP]->Rebin(5);
    
    hAverageTime_Eta_dRadius[ij][iP]->GetXaxis()->SetTitle("average time per shower (ns)");
    hAverageTime_Eta_dRadius[ij][iP]->GetXaxis()->SetRangeUser(-0.5, 0.5);
    hAverageTime_Eta_dRadius[ij][iP]->GetYaxis()->SetRangeUser(0.1, 100.);
    hAverageTime_Eta_dRadius[ij][iP]->Draw("h");

    hAverageTimeCP_Eta_dRadius[ij][iP]->Draw("h, same");
    hAverageTime_Eta_dRadius_Avg68[ij][iP]->Draw("h, same");
    hAverageTime_Eta_dRadius_AvgInt[ij][iP]->Draw("h, same");
    

    if(pdgID == "130") tL.DrawLatex(0.25,0.90, Form("K^{0}_{L}  #rho #leq%dcm  |#eta| = %.2f", radiusR.at(std::stoi(radius)), etaValue.at(ij) ) );
    if(pdgID == "22")  tL.DrawLatex(0.25,0.90, Form("#gamma     #rho #leq%dcm  |#eta| = %.2f", radiusR.at(std::stoi(radius)), etaValue.at(ij) ) );

    tLPU.DrawLatex(0.25,0.85, Form("vtx_{Z} = %dcm   p_{T} = 10GeV/c", PUvalue.at(iP)) );

    legTGMresoAll->Draw("same");
    ch_Reso_all[ij][iP]->Print((folder+"/h_Reso_all_pdg_"+pdgID+"_rad_"+radius+"_eta_"+etaBin.at(ij)+"Vtx_"+PUname.at(iP)+".png").c_str(), "png");
    ch_Reso_all[ij][iP]->Print((folder+"/h_Reso_all_pdg_"+pdgID+"_rad_"+radius+"_eta_"+etaBin.at(ij)+"Vtx_"+PUname.at(iP)+".pdf").c_str(), "pdf");
    ch_Reso_all[ij][iP]->Print((folder+"/h_Reso_all_pdg_"+pdgID+"_rad_"+radius+"_eta_"+etaBin.at(ij)+"Vtx_"+PUname.at(iP)+".root").c_str(), "root");
  }
  }



  TF1* hfithisto68[2][3];
  TF1* hfithisto68_2[2][3];
  TF1* hfithistoInt[2][3];
  TF1* hfithistoInt_2[2][3];
  
  //  gStyle->SetOptStat("rm");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  // gStyle->SetOptFit(0001);
  // gStyle->SetOptFit(0001);

  TCanvas* ch_Reso_fit[2][3];
  for(int iP=0; iP<3; ++iP){
  for(int ij=0; ij<2; ++ij){
    ch_Reso_fit[ij][iP] = new TCanvas();
    ch_Reso_fit[ij][iP]->cd();
  // gPad->SetLogx();
    
    //hAverageTime_Eta_dRadius[ij][iP]->Rebin(5);
    //hAverageTimeCP_Eta_dRadius[ij][iP]->Rebin(5);
    // hAverageTime_Eta_dRadius_Avg68[ij][iP]->Rebin(5);
    // hAverageTime_Eta_dRadius_AvgInt[ij][iP]->Rebin(5);
    //hAverageTime_Eta_dRadius_Avg680PU[ij][iP]->Rebin(5);
    

    hAverageTime_Eta_dRadius_Avg68[ij][iP]->GetXaxis()->SetTitle("average time per shower (ns)");
    hAverageTime_Eta_dRadius_Avg68[ij][iP]->GetXaxis()->SetRangeUser(-0.5, 0.5);
    hAverageTime_Eta_dRadius_Avg68[ij][iP]->GetYaxis()->SetRangeUser(0.1, 100.);
    hAverageTime_Eta_dRadius_Avg68[ij][iP]->Draw("h");

    hAverageTime_Eta_dRadius_AvgInt[ij][iP]->Draw("h, sames");


    if(pdgID == "130") tL.DrawLatex(0.25,0.90, Form("K^{0}_{L}  #rho #leq%dcm   |#eta| = %.2f", radiusR.at(std::stoi(radius)), etaValue.at(ij) ) );
    if(pdgID == "22")  tL.DrawLatex(0.25,0.90, Form("#gamma     #rho #leq%dcm   |#eta| = %.2f", radiusR.at(std::stoi(radius)), etaValue.at(ij) ) );

    tLPU.DrawLatex(0.25,0.85, Form("vtx_{Z} = %dcm   p_{T} = 10GeV/c", PUvalue.at(iP)) );

    /////////////

    hfithisto68[ij][iP] = new TF1(Form("hfithisto68_%d_PU%d", ij, PUvalue.at(iP)), "gaus", -0.5, 0.5);
    hfithisto68[ij][iP]->SetLineColor(kGreen+2);
    hfithisto68_2[ij][iP] = new TF1(Form("hfithisto68_2_%d_PU%d", ij, PUvalue.at(iP)), "gaus", -0.5, 0.5);
    hfithisto68_2[ij][iP]->SetLineColor(kGreen+2);
    hfithisto68_2[ij][iP]->SetLineWidth(2);

    hfithisto68[ij][iP]->SetRange(hAverageTime_Eta_dRadius_Avg68[ij][iP]->GetMean() - 2.*hAverageTime_Eta_dRadius_Avg68[ij][iP]->GetRMS(),
				  hAverageTime_Eta_dRadius_Avg68[ij][iP]->GetMean() + 2.*hAverageTime_Eta_dRadius_Avg68[ij][iP]->GetRMS());
    hfithisto68[ij][iP]->SetParameter(1, hAverageTime_Eta_dRadius_Avg68[ij][iP]->GetMean());
    hfithisto68[ij][iP]->SetParameter(2, hAverageTime_Eta_dRadius_Avg68[ij][iP]->GetRMS());
    hAverageTime_Eta_dRadius_Avg68[ij][iP]->Fit(Form("hfithisto68_%d_PU%d", ij, PUvalue.at(iP)), "RQ");
    hfithisto68_2[ij][iP]->SetRange(hfithisto68[ij][iP]->GetParameter(1) - 2.*hfithisto68[ij][iP]->GetParameter(2), 
				    hfithisto68[ij][iP]->GetParameter(1) + 2.*hfithisto68[ij][iP]->GetParameter(2));
    hfithisto68_2[ij][iP]->SetParameters(hfithisto68[ij][iP]->GetParameter(0), hfithisto68[ij][iP]->GetParameter(1), hfithisto68[ij][iP]->GetParameter(2));
    hAverageTime_Eta_dRadius_Avg68[ij][iP]->Fit(Form("hfithisto68_2_%d_PU%d", ij, PUvalue.at(iP)), "R");
    hfithisto68_2[ij][iP]->Draw("same");


    hfithistoInt[ij][iP] = new TF1(Form("hfithistoInt_%d_PU%d", ij, PUvalue.at(iP)), "gaus", -0.5, 0.5);
    hfithistoInt[ij][iP]->SetLineColor(kViolet+1);
    hfithistoInt_2[ij][iP] = new TF1(Form("hfithistoInt_2_%d_PU%d", ij, PUvalue.at(iP)), "gaus", -0.5, 0.5);
    hfithistoInt_2[ij][iP]->SetLineColor(kViolet+1);
    hfithistoInt_2[ij][iP]->SetLineWidth(2);

    hfithistoInt[ij][iP]->SetRange(hAverageTime_Eta_dRadius_AvgInt[ij][iP]->GetMean() - 2.*hAverageTime_Eta_dRadius_AvgInt[ij][iP]->GetRMS(),
				   hAverageTime_Eta_dRadius_AvgInt[ij][iP]->GetMean() + 2.*hAverageTime_Eta_dRadius_AvgInt[ij][iP]->GetRMS());
    hfithistoInt[ij][iP]->SetParameter(1, hAverageTime_Eta_dRadius_AvgInt[ij][iP]->GetMean());
    hfithistoInt[ij][iP]->SetParameter(2, hAverageTime_Eta_dRadius_AvgInt[ij][iP]->GetRMS());
    hAverageTime_Eta_dRadius_AvgInt[ij][iP]->Fit(Form("hfithistoInt_%d_PU%d", ij, PUvalue.at(iP)), "RQ");
    hfithistoInt_2[ij][iP]->SetRange(hfithistoInt[ij][iP]->GetParameter(1) - 2.*hfithistoInt[ij][iP]->GetParameter(2), 
				     hfithistoInt[ij][iP]->GetParameter(1) + 2.*hfithistoInt[ij][iP]->GetParameter(2));
    hfithistoInt_2[ij][iP]->SetParameters(hfithistoInt[ij][iP]->GetParameter(0), hfithistoInt[ij][iP]->GetParameter(1), hfithistoInt[ij][iP]->GetParameter(2));
    hAverageTime_Eta_dRadius_AvgInt[ij][iP]->Fit(Form("hfithistoInt_2_%d_PU%d", ij, PUvalue.at(iP)), "R");
    hfithistoInt_2[ij][iP]->Draw("same");


    ch_Reso_all[ij][iP]->Update();


    TLatex tPS68_a;
    tPS68_a.SetNDC();
    tPS68_a.SetTextSize(0.04);
    tPS68_a.SetTextFont(132);
    tPS68_a.SetTextColor(kGreen+2);
    tPS68_a.DrawLatex(0.70, 0.65, Form("mean = %.3f #pm %.3f",  hfithisto68_2[ij][iP]->GetParameter(1),  hfithisto68_2[ij][iP]->GetParError(1)) );

    TLatex tPS68_b;
    tPS68_b.SetNDC();
    tPS68_b.SetTextSize(0.04);
    tPS68_b.SetTextFont(132);
    tPS68_b.SetTextColor(kGreen+2);
    tPS68_b.DrawLatex(0.70, 0.60, Form("#sigma = %.3f #pm %.3f", hfithisto68_2[ij][iP]->GetParameter(2),  hfithisto68_2[ij][iP]->GetParError(2)) );


    TLatex tPSInt_a;
    tPSInt_a.SetNDC();
    tPSInt_a.SetTextSize(0.04);
    tPSInt_a.SetTextFont(132);
    tPSInt_a.SetTextColor(kViolet+1);
    tPSInt_a.DrawLatex(0.70, 0.55, Form("mean = %.3f #pm %.3f",  hfithistoInt_2[ij][iP]->GetParameter(1),  hfithistoInt_2[ij][iP]->GetParError(1)) );

    TLatex tPSInt_b;
    tPSInt_b.SetNDC();
    tPSInt_b.SetTextSize(0.04);
    tPSInt_b.SetTextFont(132);
    tPSInt_b.SetTextColor(kViolet+1);
    tPSInt_b.DrawLatex(0.70, 0.50, Form("#sigma = %.3f #pm %.3f", hfithistoInt_2[ij][iP]->GetParameter(2),  hfithistoInt_2[ij][iP]->GetParError(2)) );


    if(ij == 0){
      tgReso[iP]->SetPoint(iP+1, 1.75, hfithistoInt_2[ij][iP]->GetParameter(1));
      tgReso[iP]->SetPointError(iP+1, 0, sqrt(pow(hfithistoInt_2[ij][iP]->GetParError(1), 2) + pow(hfithistoInt_2[ij][iP]->GetParameter(2), 2) ));
    }

    if(ij == 1){
      tgReso[iP]->SetPoint(iP+2, 2.75, hfithistoInt_2[ij][iP]->GetParameter(1));
      tgReso[iP]->SetPointError(iP+2, 0, sqrt( pow(hfithistoInt_2[ij][iP]->GetParError(1), 2) + pow(hfithistoInt_2[ij][iP]->GetParameter(2), 2) ) );
    }

    std::cout << " ij = " << ij << " iP = " << iP << " reso = " << hfithistoInt_2[ij][iP]->GetParameter(2) << std::endl;

    legTGMresoFit->Draw("same");
    ch_Reso_fit[ij][iP]->Print((folder+"/h_Reso_fit_pdg_"+pdgID+"_rad_"+radius+"_eta_"+etaBin.at(ij)+"Vtx_"+PUname.at(iP)+".png").c_str(), "png");
    ch_Reso_fit[ij][iP]->Print((folder+"/h_Reso_fit_pdg_"+pdgID+"_rad_"+radius+"_eta_"+etaBin.at(ij)+"Vtx_"+PUname.at(iP)+".pdf").c_str(), "pdf");
    ch_Reso_fit[ij][iP]->Print((folder+"/h_Reso_fit_pdg_"+pdgID+"_rad_"+radius+"_eta_"+etaBin.at(ij)+"Vtx_"+PUname.at(iP)+".root").c_str(), "root");

    hAverageTime_Eta_dRadius_Avg68[ij][iP]->GetYaxis()->SetRangeUser(0.1, 500.);
    gPad->SetLogy();
    ch_Reso_fit[ij][iP]->Print((folder+"/h_Reso_fit_pdg_"+pdgID+"_rad_"+radius+"_eta_"+etaBin.at(ij)+"Vtx_"+PUname.at(iP)+"logY.png").c_str(), "png");
    ch_Reso_fit[ij][iP]->Print((folder+"/h_Reso_fit_pdg_"+pdgID+"_rad_"+radius+"_eta_"+etaBin.at(ij)+"Vtx_"+PUname.at(iP)+"logY.pdf").c_str(), "pdf");
    ch_Reso_fit[ij][iP]->Print((folder+"/h_Reso_fit_pdg_"+pdgID+"_rad_"+radius+"_eta_"+etaBin.at(ij)+"Vtx_"+PUname.at(iP)+"logY.root").c_str(), "root");
  }
  }

  TLegend *legReso = new TLegend(0.70,0.30,0.90,0.53,NULL,"brNDC");
  legReso->SetTextFont(42);
  legReso->SetTextSize(0.04);
  legReso->SetFillColor(kWhite);
  legReso->SetLineColor(kWhite);
  legReso->SetShadowColor(kWhite);
  for(int iF=0; iF<3; ++iF)
    legReso->AddEntry(tgFrEvt[iF], Form("vtx_{Z} = %dcm",PUvalue.at(iF)), "p");

  TLegend *legResoD = new TLegend(0.70,0.70,0.90,0.93,NULL,"brNDC");
  legResoD->SetTextFont(42);
  legResoD->SetTextSize(0.04);
  legResoD->SetFillColor(kWhite);
  legResoD->SetLineColor(kWhite);
  legResoD->SetShadowColor(kWhite);
  for(int iF=0; iF<3; ++iF)
    legResoD->AddEntry(tgFrEvt[iF], Form("vtx_{Z} = %dcm",PUvalue.at(iF)), "p");



  TCanvas* ch_Eff = new TCanvas();
  ch_Eff->cd();
  for(int iF=0; iF<3; ++iF){
  tgFrEvt[iF]->GetXaxis()->SetTitle("|#eta|");
  tgFrEvt[iF]->GetXaxis()->SetRangeUser(1.5, 3.);
  tgFrEvt[iF]->GetYaxis()->SetTitle("efficiency");
  tgFrEvt[iF]->GetYaxis()->SetRangeUser(0., 1.2);
  if(iF == 0)  tgFrEvt[iF]->Draw("ap");
  else tgFrEvt[iF]->Draw("p, same");
  }

  if(pdgID == "130") tL.DrawLatex(0.25,0.90, Form("K^{0}_{L}  #rho #leq%dcm ", radiusR.at(std::stoi(radius)) ) );
  if(pdgID == "22")  tL.DrawLatex(0.25,0.90, Form("#gamma     #rho #leq%dcm ", radiusR.at(std::stoi(radius)) ) );

  tLPU.DrawLatex(0.25,0.85, "p_{T} = 10GeV/c" );

  legReso->Draw("same");
  ch_Eff->Print((folder+"/h_Eff"+pdgID+"_rad_"+radius+".png").c_str(), "png");
  ch_Eff->Print((folder+"/h_Eff"+pdgID+"_rad_"+radius+".pdf").c_str(), "pdf");
  ch_Eff->Print((folder+"/h_Eff"+pdgID+"_rad_"+radius+".root").c_str(), "root");

  //reso
  TCanvas* ch_Reso = new TCanvas();
  ch_Reso->cd();
  for(int iF=0; iF<3; ++iF){
  tgReso[iF]->GetXaxis()->SetTitle("|#eta|");
  tgReso[iF]->GetXaxis()->SetRangeUser(1.5, 3.);
  //  tgReso[iF]->GetYaxis()->SetTitle("#sigma (ns)");
  tgReso[iF]->GetYaxis()->SetTitle("<time> (ns)");
  tgReso[iF]->GetYaxis()->SetRangeUser(-0.5, 1.);
  if(iF == 0)  tgReso[iF]->Draw("ap");
  else tgReso[iF]->Draw("p, same");
  }

  if(pdgID == "130") tL.DrawLatex(0.25,0.90, Form("K^{0}_{L}  #rho #leq%dcm ", radiusR.at(std::stoi(radius)) ) );
  if(pdgID == "22")  tL.DrawLatex(0.25,0.90, Form("#gamma     #rho #leq%dcm ", radiusR.at(std::stoi(radius)) ) );

  tLPU.DrawLatex(0.25,0.85, "p_{T} = 10GeV/c" );

  legResoD->Draw("same");
  ch_Reso->Print((folder+"/h_Resp"+pdgID+"_rad_"+radius+".png").c_str(), "png");
  ch_Reso->Print((folder+"/h_Reso"+pdgID+"_rad_"+radius+".pdf").c_str(), "pdf");
  ch_Reso->Print((folder+"/h_Reso"+pdgID+"_rad_"+radius+".root").c_str(), "root");


  //////////
  TCanvas* ch_TimeAll[2];
  for(int ij=0; ij<2; ++ij){
    ch_TimeAll[ij] = new TCanvas();
    ch_TimeAll[ij]->cd();
    gPad->SetLogy();
    hTime_Eta_dRadius[ij][0]->GetXaxis()->SetTitle(("hit time (ns) eta bin "+etaBin.at(ij)).c_str());
    hTime_Eta_dRadius[ij][0]->GetXaxis()->SetRangeUser(-3., 5.);
    hTime_Eta_dRadius[ij][0]->GetYaxis()->SetRangeUser(0.1, 5.e+3);

    for(int iF=0; iF<3; ++iF) {
      hTime_Eta_dRadius[ij][iF]->SetLineColor(iColors[iF]);
      hTime_Eta_dRadius[ij][iF]->Scale(1.*hTime_Eta_dRadius[ij][2]->Integral()/hTime_Eta_dRadius[ij][iF]->Integral());
    }

    hTime_Eta_dRadius[ij][0]->Draw("h");
    hTime_Eta_dRadius[ij][1]->Draw("h, same");
    hTime_Eta_dRadius[ij][2]->Draw("h, same");
    
    // hAverageTime_Eta_dRadius_Avg68[ij][iP]->Rebin(5);                                                                                                                                                            
    // hAverageTime_Eta_dRadius_AvgInt[ij][iP]->Rebin(5);

    
    if(pdgID == "130") tL.DrawLatex(0.25,0.90, Form("K^{0}_{L}  #rho #leq%dcm ", radiusR.at(std::stoi(radius)) ) );
    if(pdgID == "22")  tL.DrawLatex(0.25,0.90, Form("#gamma     #rho #leq%dcm ", radiusR.at(std::stoi(radius)) ) );

    tLPU.DrawLatex(0.25,0.85, "p_{T} = 10GeV/c" );

    legResoD->SetHeader(" all hits ");
    legResoD->Draw("same");
    ch_TimeAll[ij]->Print((folder+"/h_TimeAll"+pdgID+"_rad_"+radius+"_eta_"+etaBin.at(ij)+".png").c_str(), "png");
    ch_TimeAll[ij]->Print((folder+"/h_TimeAll"+pdgID+"_rad_"+radius+"_eta_"+etaBin.at(ij)+".pdf").c_str(), "pdf");
    ch_TimeAll[ij]->Print((folder+"/h_TimeAll"+pdgID+"_rad_"+radius+"_eta_"+etaBin.at(ij)+".root").c_str(), "root");
  }

  TCanvas* ch_TimeInt[2];
  for(int ij=0; ij<2; ++ij){
    ch_TimeInt[ij] = new TCanvas();
    ch_TimeInt[ij]->cd();
    gPad->SetLogy();
    hTimeIntCut_Eta_dRadius[ij][0]->GetXaxis()->SetTitle(("hit time (ns) eta eta "+etaBin.at(ij)).c_str());
    hTimeIntCut_Eta_dRadius[ij][0]->GetXaxis()->SetRangeUser(-3., 5.);
    hTimeIntCut_Eta_dRadius[ij][0]->GetYaxis()->SetRangeUser(0.1, 5.e+3);

    for(int iF=0; iF<3; ++iF){
      hTimeIntCut_Eta_dRadius[ij][iF]->SetLineColor(iColors[iF]);
      hTimeIntCut_Eta_dRadius[ij][iF]->Scale(1.*hTimeIntCut_Eta_dRadius[ij][2]->Integral()/hTimeIntCut_Eta_dRadius[ij][iF]->Integral());
    }


    hTimeIntCut_Eta_dRadius[ij][0]->Draw("h");
    hTimeIntCut_Eta_dRadius[ij][1]->Draw("h, same");
    hTimeIntCut_Eta_dRadius[ij][2]->Draw("h, same");
    
    // hAverageTime_Eta_dRadius_Avg68[ij][iP]->Rebin(5);                                                                                                                                                            
    // hAverageTime_Eta_dRadius_AvgInt[ij][iP]->Rebin(5);

    
    if(pdgID == "130") tL.DrawLatex(0.25,0.90, Form("K^{0}_{L}  #rho #leq%dcm ", radiusR.at(std::stoi(radius)) ) );
    if(pdgID == "22")  tL.DrawLatex(0.25,0.90, Form("#gamma     #rho #leq%dcm ", radiusR.at(std::stoi(radius)) ) );

    tLPU.DrawLatex(0.25,0.85, "p_{T} = 10GeV/c" );
    
    legResoD->SetHeader(" selected hits ");
    legResoD->Draw("same");
    ch_TimeInt[ij]->Print((folder+"/h_TimeInt"+pdgID+"_rad_"+radius+"_eta_"+etaBin.at(ij)+".png").c_str(), "png");
    ch_TimeInt[ij]->Print((folder+"/h_TimeInt"+pdgID+"_rad_"+radius+"_eta_"+etaBin.at(ij)+".pdf").c_str(), "pdf");
    ch_TimeInt[ij]->Print((folder+"/h_TimeInt"+pdgID+"_rad_"+radius+"_eta_"+etaBin.at(ij)+".root").c_str(), "root");
  }


  return;

}
