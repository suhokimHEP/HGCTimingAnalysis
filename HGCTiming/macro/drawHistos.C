#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TFile.h>
#include <TCanvas.h>
#include "TRandom3.h"
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





void drawHistos(){

  int type = 2;
  TFile *_file0;
  if(type == 0) _file0 =  TFile::Open("../test/CSF20LBA50/JOB_PDG_22_Pt2_3fC/OutTimeHGC_RecHits_PDG_22_Pt2_3fC.root");  
  if(type == 1) _file0 =  TFile::Open("../test/CSF20LBA50/JOB_PDG_22_Pt60_3fC/OutTimeHGC_RecHits_PDG_22_Pt60_3fC.root");
    if(type == 2) _file0 = TFile::Open("../test/CSF20LBA50/JOB_PDG_130_Pt2_3fC/OutTimeHGC_RecHits_PDG_130_Pt2_3fC.root");
  //  if(type == 2) _file0 = TFile::Open("../test/CSF20LEA50/JOB_PDG_130_Pt2_3fC/OutTimeHGC_RecHits_PDG_130_Pt2_3fC.root");
  if(type == 3) _file0 = TFile::Open("../test/CSF20LBA50/JOB_PDG_130_Pt100_3fC/OutTimeHGC_RecHits_PDG_130_Pt100_3fC.root");
  if(type == 4) _file0 = TFile::Open("../test/CSF20LBA50/JOB_PDG_130_Pt1_3fC/OutTimeHGC_RecHits_PDG_130_Pt1_3fC.root");
  if(type == 5) _file0 = TFile::Open("../test/CSF20LBA50/JOB_PDG_130_Pt07_3fC/OutTimeHGC_RecHits_PDG_130_Pt07_3fC.root");

  std::vector<std::string> names;
  names.push_back("");
  names.push_back("_ResoWe");
  //names.push_back("_AvgCutH");
  names.push_back("_Avg68");

  int iColors[4] = {kBlack, kRed, kBlue, kGreen+2};

  TH1F* hN[4];
  TLatex tL[4];
  for(int ij=0; ij<3; ++ij){
    hN[ij] = (TH1F*)_file0->Get(("ana/hAverageTime_Eta1.65-1.85_dRadius0"+names.at(ij)).c_str());
    //  hN[ij] = (TH1F*)_file0->Get(("ana/hAverageTime_Eta2.45-2.65_dRadius0"+names.at(ij)).c_str());
    std::cout <<("ana/hAverageTime_Eta1.65-1.85_dRadius0"+names.at(ij)).c_str() << std::endl;
    // if(hN[ij]->GetEntries() < 100) hN[ij]->Rebin(10);
    if(type == 2 || type > 3) hN[ij]->Rebin(5);
    hN[ij]->Rebin(5);
    hN[ij]->SetLineColor(iColors[ij]);
    hN[ij]->SetLineWidth(2);


    //    if(type == 2)  hN[ij]->Rebin(4);
    tL[ij].SetNDC();
    tL[ij].SetTextSize(0.05);
    tL[ij].SetTextFont(132);
    tL[ij].SetTextColor(iColors[ij]);
  }

  TLegend *legTGM = new TLegend(0.5,0.65,0.70,0.85,NULL,"brNDC");
  legTGM->SetTextFont(42);
  legTGM->SetTextSize(0.04);
  legTGM->SetFillColor(kWhite);
  legTGM->SetLineColor(kWhite);
  legTGM->SetShadowColor(kWhite);
  legTGM->AddEntry(hN[0], "all hits", "l");
  legTGM->AddEntry(hN[1], "res weig", "l");
  //  legTGM->AddEntry(hN[2], "trunc 30\%", "l");
  legTGM->AddEntry(hN[2], "avg 60\% res weig", "l");

  gStyle->SetOptStat(0);

  TCanvas* c1 = new TCanvas();
  c1->cd();
  hN[1]->GetXaxis()->SetTitle("average time (ns)");
  hN[1]->GetYaxis()->SetTitle("events");
  hN[1]->GetXaxis()->SetRangeUser(-0.2, 0.5); 
  if(type == 2 || type > 3)  hN[1]->GetXaxis()->SetRangeUser(-0.2, 1.5);
  //  hN[3]->GetYaxis()->SetRangeUser(0., 400.);
  hN[1]->Draw();
  hN[0]->Draw("same");
  tL[0].DrawLatex(0.5,0.6,Form("m = %.2e rms = %.2e", hN[0]->GetMean(), hN[0]->GetRMS()) );
  for(int ij=1; ij<3; ++ij){
    hN[ij]->Draw("same");
    tL[ij].DrawLatex(0.5,0.6-(0.05*ij),Form("m = %.2e rms = %.2e", hN[ij]->GetMean(), hN[ij]->GetRMS()) );
  }
  legTGM->Draw("same");
  if(type == 0){
    c1->Print("averageHistos/averageHistos_22_Pt2_CSF20LBA50.png", "png");
    c1->Print("averageHistos/averageHistos_22_Pt2_CSF20LBA50.pdf", "pdf");
    c1->Print("averageHistos/averageHistos_22_Pt2_CSF20LBA50.root", "root");
  }
  if(type == 1){
    c1->Print("averageHistos/averageHistos_22_Pt60_CSF20LBA50.png", "png");
    c1->Print("averageHistos/averageHistos_22_Pt60_CSF20LBA50.pdf", "pdf");
    c1->Print("averageHistos/averageHistos_22_Pt60_CSF20LBA50.root", "root");
  }
  if(type == 2){
    c1->Print("averageHistos/averageHistos_130_Pt2_CSF20LBA50.png", "png");
    c1->Print("averageHistos/averageHistos_130_Pt2_CSF20LBA50.pdf", "pdf");
    c1->Print("averageHistos/averageHistos_130_Pt2_CSF20LBA50.root", "root");
  }
  if(type == 3){
    c1->Print("averageHistos/averageHistos_130_Pt100_CSF20LBA50.png", "png");
    c1->Print("averageHistos/averageHistos_130_Pt100_CSF20LBA50.pdf", "pdf");
    c1->Print("averageHistos/averageHistos_130_Pt100_CSF20LBA50.root", "root");
  }
  if(type == 4){
    c1->Print("averageHistos/averageHistos_130_Pt1_CSF20LBA50.png", "png");
    c1->Print("averageHistos/averageHistos_130_Pt1_CSF20LBA50.pdf", "pdf");
    c1->Print("averageHistos/averageHistos_130_Pt1_CSF20LBA50.root", "root");
  }
  if(type == 5){
    c1->Print("averageHistos/averageHistos_130_Pt07_CSF20LBA50.png", "png");
    c1->Print("averageHistos/averageHistos_130_Pt07_CSF20LBA50.pdf", "pdf");
    c1->Print("averageHistos/averageHistos_130_Pt07_CSF20LBA50.root", "root");
  }



}
