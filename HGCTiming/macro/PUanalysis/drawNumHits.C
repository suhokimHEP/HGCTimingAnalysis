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
  
  gROOT->Reset();
  gROOT->Macro("~/public/setStyle.C");
  //  gROOT->LoadMacro("~/public/myStyle.C");
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1);

  std::cout << " inizio ci sono " << std::endl; 


  bool doAllTheFits = true;
  
  //  int iColors[16] = {kRed, kOrange+4, kOrange-3, kOrange-2, kBlue, kBlue-9, kAzure-9, kAzure+10, kCyan, kGreen+1, kCyan-2, kYellow+2}; //kGray+1};
  //                 all     PU     ME   shared    sharedME   sharedPU
  int iColors[6] = {kBlack, kBlue, kRed, kGreen+2, kMagenta+2, kBlue-7};
  int iStyle[4] = {20, 21, 22, 23}; //kGray+1};
  //  int iColors[6] = {kGreen+1, kBlue, kRed};
  
  int nBinsEta = 6;
  float binWidth = 0.2;
  int nBinsRad = 4;
  //int nBinsEta = 3;
  //float binWidth = 0.5;
  float binStart = 1.65;
  
  std::string pdgID = "130";
  //std::string pdgID = "22";
  //std::string pdgID = "211";

  std::string radius = "0";
  std::vector<int> radiusR;
  radiusR.push_back(2);
  radiusR.push_back(5);

  std::vector<std::string> etaBin;
  etaBin.push_back("1.65-1.85");
  etaBin.push_back("2.65-2.85");

  std::cout << " >>> ora prendo i file " << std::endl;


  
  /*
  std::string folderName = "plotsAllTimeTot_CFD_"+pdgID+"_3hits";
  */

  std::string folder = "plotsPU_"+pdgID;


  //  return;
  TFile* inF[2];
  //inF[0] = TFile::Open(("../../test/testPU/ROOTFILES/OutTimeHGC_RecHits_PDG_"+pdgID+"_Pt5_200PU_lowEta.root").c_str());
  //inF[1] = TFile::Open(("../../test/testPU/ROOTFILES/OutTimeHGC_RecHits_PDG_"+pdgID+"_Pt5_200PU_highEta.root").c_str());
  inF[0] = TFile::Open(("../../test/testPU/OutTimeHGC_RecHits_PDG_"+pdgID+"_Pt5_200PU_lowEta.root").c_str());
  inF[1] = TFile::Open(("../../test/testPU/OutTimeHGC_RecHits_PDG_"+pdgID+"_Pt5_200PU_highEta.root").c_str());
  TFile* inFanyPU[2];
  // inFanyPU[0] =  TFile::Open(("../../test/testPU/OutTimeHGC_RecHits_PDG_"+pdgID+"_Pt5_AllPU_lowEta.root").c_str());
  // inFanyPU[1] =  TFile::Open(("../../test/testPU/OutTimeHGC_RecHits_PDG_"+pdgID+"_Pt5_AllPU_highEta.root").c_str());
  inFanyPU[0] =  TFile::Open("../../test/testPU/OutTimeHGC_RecHits_PDG_allPDG_Pt5_AllPU_lowEta.root");
  inFanyPU[1] =  TFile::Open("../../test/testPU/OutTimeHGC_RecHits_PDG_allPDG_Pt5_AllPU_highEta.root");
  TFile* inFanyME[2];  
  inFanyME[0] =  TFile::Open(("../../test/testPU/OutTimeHGC_RecHits_PDG_"+pdgID+"_Pt5_0PU.root").c_str());
  inFanyME[1] =  TFile::Open(("../../test/testPU/OutTimeHGC_RecHits_PDG_"+pdgID+"_Pt5_0PU.root").c_str());

  std::cout << " >>> fatto = presi " << std::endl;
  TH1F* hFractionEvents_PU_Eta_dRadius_3MIP[2];
  TH1F* hFractionEvents_PUany_Eta_dRadius_3MIP[2][2];
  TH1F* hFractionEvents_MainEvt_Eta_dRadius_3MIP[2];
  TH1F* hFractionEvents_MainEvtany_Eta_dRadius_3MIP[2][2];
  TH1F* hFractionEvents_Shared_Eta_dRadius_3MIP[2];
  TH1F* hFractionEvents_SharedME_Eta_dRadius_3MIP[2];
  TH1F* hFractionEvents_SharedPU_Eta_dRadius_3MIP[2];
  
  TH1F* h_NumberHits_Eta_dRadius_3MIP[2];
  TH1F* h_NumberHits_PU_Eta_dRadius_3MIP[2];
  TH1F* h_NumberHits_PUany_Eta_dRadius_3MIP[2][2];
  TH1F* h_NumberHits_MainEvt_Eta_dRadius_3MIP[2];
  TH1F* h_NumberHits_MainEvtany_Eta_dRadius_3MIP[2][2];
  TH1F* h_NumberHits_Shared_Eta_dRadius_3MIP[2];
  TH1F* h_NumberHits_SharedME_Eta_dRadius_3MIP[2];
  TH1F* h_NumberHits_SharedPU_Eta_dRadius_3MIP[2];

  TH1F* h_FractionHits_PU_Eta_dRadius_3MIP[2];
  TH1F* h_FractionHits_PUany_Eta_dRadius_3MIP[2];
  TH1F* h_FractionHits_MainEvt_Eta_dRadius_3MIP[2];
  TH1F* h_FractionHits_MainEvtany_Eta_dRadius_3MIP[2];
  TH1F* h_FractionHits_Shared_Eta_dRadius_3MIP[2];
  TH1F* h_FractionHits_SharedME_Eta_dRadius_3MIP[2];
  TH1F* h_FractionHits_SharedPU_Eta_dRadius_3MIP[2];
  TH1F* h_Energy_Eta_dRadius_3MIP[2];
  TH1F* h_Energy_PU_Eta_dRadius_3MIP[2];
  TH1F* h_Energy_PUany_Eta_dRadius_3MIP[2][2];
  TH1F* h_Energy_MainEvt_Eta_dRadius_3MIP[2];
  TH1F* h_Energy_MainEvtany_Eta_dRadius_3MIP[2][2];
  TH1F* h_Energy_Shared_Eta_dRadius_3MIP[2];
  TH1F* h_Energy_SharedME_Eta_dRadius_3MIP[2];
  TH1F* h_Energy_SharedPU_Eta_dRadius_3MIP[2];
  TH1F* h_FractionEnergy_PU_Eta_dRadius_3MIP[2];
  TH1F* h_FractionEnergy_PUany_Eta_dRadius_3MIP[2];
  TH1F* h_FractionEnergy_MainEvt_Eta_dRadius_3MIP[2];
  TH1F* h_FractionEnergy_MainEvtany_Eta_dRadius_3MIP[2];
  TH1F* h_FractionEnergy_Shared_Eta_dRadius_3MIP[2];
  TH1F* h_FractionEnergy_SharedME_Eta_dRadius_3MIP[2];
  TH1F* h_FractionEnergy_SharedPU_Eta_dRadius_3MIP[2];



  for(int ij=0; ij<2; ++ij){
    /*
    hFractionEvents_PU_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/hFractionEvents_PU_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"_3MIP").c_str()));
    hFractionEvents_PUany_Eta_dRadius_3MIP[ij][0] = (TH1F*)(inF[ij]->Get(("ana/hFractionEvents_PUany_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"_3MIP").c_str()));
    hFractionEvents_PUany_Eta_dRadius_3MIP[ij][1] = (TH1F*)(inFanyPU[ij]->Get(("ana/hFractionEvents_PUany_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"_3MIP").c_str()));
    hFractionEvents_MainEvt_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/hFractionEvents_MainEvt_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"_3MIP").c_str()));
    hFractionEvents_MainEvtany_Eta_dRadius_3MIP[ij][0] = (TH1F*)(inF[ij]->Get(("ana/hFractionEvents_MainEvtany_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"_3MIP").c_str()));
    hFractionEvents_MainEvtany_Eta_dRadius_3MIP[ij][1] = (TH1F*)(inFanyME[ij]->Get(("ana/hFractionEvents_MainEvtany_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"_3MIP").c_str()));
    hFractionEvents_Shared_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/hFractionEvents_Shared_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"_3MIP").c_str()));
    hFractionEvents_SharedME_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/hFractionEvents_SharedME_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"_3MIP").c_str()));
    hFractionEvents_SharedPU_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/hFractionEvents_SharedPU_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"_3MIP").c_str()));
    */
    hFractionEvents_PU_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/hFractionEvents_PU_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"").c_str()));
    hFractionEvents_PUany_Eta_dRadius_3MIP[ij][0] = (TH1F*)(inF[ij]->Get(("ana/hFractionEvents_PUany_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"").c_str()));
    hFractionEvents_PUany_Eta_dRadius_3MIP[ij][1] = (TH1F*)(inFanyPU[ij]->Get(("ana/hFractionEvents_PUany_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"").c_str()));
    hFractionEvents_MainEvt_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/hFractionEvents_MainEvt_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"").c_str()));
    hFractionEvents_MainEvtany_Eta_dRadius_3MIP[ij][0] = (TH1F*)(inF[ij]->Get(("ana/hFractionEvents_MainEvtany_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"").c_str()));
    hFractionEvents_MainEvtany_Eta_dRadius_3MIP[ij][1] = (TH1F*)(inFanyME[ij]->Get(("ana/hFractionEvents_MainEvtany_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"").c_str()));
    hFractionEvents_Shared_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/hFractionEvents_Shared_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"").c_str()));
    hFractionEvents_SharedME_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/hFractionEvents_SharedME_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"").c_str()));
    hFractionEvents_SharedPU_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/hFractionEvents_SharedPU_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"").c_str()));


    std::cout << " bin " << etaBin.at(ij) << " fraction of events (counting hits > 3 MIP): " << std::endl;
    std::cout << " from PU =                      " << hFractionEvents_PU_Eta_dRadius_3MIP[ij]->GetMean() 
	      << " pm " << hFractionEvents_PU_Eta_dRadius_3MIP[ij]->GetRMS() / sqrt(hFractionEvents_PU_Eta_dRadius_3MIP[ij]->GetEntries()) << std::endl;
    // std::cout << " from PU any (200PU +ME)  =     " << hFractionEvents_PUany_Eta_dRadius_3MIP[ij][0]->GetMean() << std::endl;
    // std::cout << " from PU any (200PU)  =         " << hFractionEvents_PUany_Eta_dRadius_3MIP[ij][1]->GetMean() << std::endl;
    std::cout << " from mainEvt =                 " << hFractionEvents_MainEvt_Eta_dRadius_3MIP[ij]->GetMean() 
	      << " pm " << hFractionEvents_MainEvt_Eta_dRadius_3MIP[ij]->GetRMS() / sqrt(hFractionEvents_MainEvt_Eta_dRadius_3MIP[ij]->GetEntries()) << std::endl;
    // std::cout << " from mainEvt any (200PU +ME) = " << hFractionEvents_MainEvtany_Eta_dRadius_3MIP[ij][0]->GetMean() << std::endl;
    // std::cout << " from mainEvt any (200PU) =     " << hFractionEvents_MainEvtany_Eta_dRadius_3MIP[ij][1]->GetMean() << std::endl;
    std::cout << " shared  =                      " << hFractionEvents_Shared_Eta_dRadius_3MIP[ij]->GetMean() 
	      << " pm " << hFractionEvents_Shared_Eta_dRadius_3MIP[ij]->GetRMS() / sqrt(hFractionEvents_Shared_Eta_dRadius_3MIP[ij]->GetEntries()) << std::endl;
    std::cout << " shared ME>3Mip  =              " << hFractionEvents_SharedME_Eta_dRadius_3MIP[ij]->GetMean() 
	      << " pm " << hFractionEvents_SharedME_Eta_dRadius_3MIP[ij]->GetRMS() / sqrt(hFractionEvents_SharedME_Eta_dRadius_3MIP[ij]->GetEntries()) << std::endl;
    std::cout << " shared PU>3Mip  =              " << hFractionEvents_SharedPU_Eta_dRadius_3MIP[ij]->GetMean() 
	      << " pm " << hFractionEvents_SharedPU_Eta_dRadius_3MIP[ij]->GetRMS() / sqrt(hFractionEvents_SharedPU_Eta_dRadius_3MIP[ij]->GetEntries()) << std::endl;

    /*
    h_NumberHits_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/h_NumberHits_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"_3MIP").c_str()));
    h_NumberHits_PU_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/h_NumberHits_PU_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"_3MIP").c_str()));
    h_NumberHits_PUany_Eta_dRadius_3MIP[ij][0] = (TH1F*)(inF[ij]->Get(("ana/h_NumberHits_PUany_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"_3MIP").c_str()));
    h_NumberHits_PUany_Eta_dRadius_3MIP[ij][1] = (TH1F*)(inFanyPU[ij]->Get(("ana/h_NumberHits_PUany_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"_3MIP").c_str()));
    h_NumberHits_MainEvt_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/h_NumberHits_MainEvt_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"_3MIP").c_str()));
    h_NumberHits_MainEvtany_Eta_dRadius_3MIP[ij][0] = (TH1F*)(inF[ij]->Get(("ana/h_NumberHits_MainEvtany_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"_3MIP").c_str()));
    h_NumberHits_MainEvtany_Eta_dRadius_3MIP[ij][1] = (TH1F*)(inFanyME[ij]->Get(("ana/h_NumberHits_MainEvtany_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"_3MIP").c_str()));
    */
    h_NumberHits_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/h_NumberHits_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"").c_str()));
    h_NumberHits_PU_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/h_NumberHits_PU_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"").c_str()));
    h_NumberHits_PUany_Eta_dRadius_3MIP[ij][0] = (TH1F*)(inF[ij]->Get(("ana/h_NumberHits_PUany_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"").c_str()));
    h_NumberHits_PUany_Eta_dRadius_3MIP[ij][1] = (TH1F*)(inFanyPU[ij]->Get(("ana/h_NumberHits_PUany_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"").c_str()));
    h_NumberHits_MainEvt_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/h_NumberHits_MainEvt_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"").c_str()));
    h_NumberHits_MainEvtany_Eta_dRadius_3MIP[ij][0] = (TH1F*)(inF[ij]->Get(("ana/h_NumberHits_MainEvtany_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"").c_str()));
    h_NumberHits_MainEvtany_Eta_dRadius_3MIP[ij][1] = (TH1F*)(inFanyME[ij]->Get(("ana/h_NumberHits_MainEvtany_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"").c_str()));


    std::cout << " bin " << etaBin.at(ij) << " number of hits (counting hits > 3 MIP): " << std::endl;
    std::cout << " from PU any 200 + ME           " << h_NumberHits_PUany_Eta_dRadius_3MIP[ij][0]->GetMean() 
	      << " pm " << h_NumberHits_PUany_Eta_dRadius_3MIP[ij][0]->GetRMS() / sqrt(h_NumberHits_PUany_Eta_dRadius_3MIP[ij][0]->GetEntries()) << std::endl;
    std::cout << " from PU any 200                " << h_NumberHits_PUany_Eta_dRadius_3MIP[ij][1]->GetMean() 
	      << " pm " << h_NumberHits_PUany_Eta_dRadius_3MIP[ij][1]->GetRMS() / sqrt(h_NumberHits_PUany_Eta_dRadius_3MIP[ij][1]->GetEntries()) << std::endl;

    std::cout << " mainEvt any 200 + ME =         " << h_NumberHits_MainEvtany_Eta_dRadius_3MIP[ij][0]->GetMean()
	      << " pm " << h_NumberHits_MainEvtany_Eta_dRadius_3MIP[ij][0]->GetRMS() / sqrt(h_NumberHits_MainEvtany_Eta_dRadius_3MIP[ij][0]->GetEntries()) << std::endl;
    std::cout << " mainEvt any 200                " << h_NumberHits_MainEvtany_Eta_dRadius_3MIP[ij][1]->GetMean()
	      << " pm " << h_NumberHits_MainEvtany_Eta_dRadius_3MIP[ij][1]->GetRMS() / sqrt(h_NumberHits_MainEvtany_Eta_dRadius_3MIP[ij][1]->GetEntries()) << std::endl;


    /*
    h_NumberHits_Shared_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/h_NumberHits_Shared_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"_3MIP").c_str()));
    h_NumberHits_SharedME_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/h_NumberHits_SharedME_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"_3MIP").c_str()));
    h_NumberHits_SharedPU_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/h_NumberHits_SharedPU_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"_3MIP").c_str()));
    h_FractionHits_PU_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/h_FractionHits_PU_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"_3MIP").c_str()));
    h_FractionHits_PUany_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/h_FractionHits_PUany_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"_3MIP").c_str()));
    h_FractionHits_MainEvt_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/h_FractionHits_MainEvt_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"_3MIP").c_str()));
    h_FractionHits_MainEvtany_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/h_FractionHits_MainEvtany_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"_3MIP").c_str()));
    h_FractionHits_Shared_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/h_FractionHits_Shared_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"_3MIP").c_str()));
    h_FractionHits_Shared_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/h_FractionHits_Shared_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"_3MIP").c_str()));
    h_FractionHits_SharedME_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/h_FractionHits_SharedME_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"_3MIP").c_str()));
    h_FractionHits_SharedPU_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/h_FractionHits_SharedPU_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"_3MIP").c_str()));

    h_Energy_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/h_Energy_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"_3MIP").c_str()));
    h_Energy_PU_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/h_Energy_PU_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"_3MIP").c_str()));
    h_Energy_PUany_Eta_dRadius_3MIP[ij][0] = (TH1F*)(inF[ij]->Get(("ana/h_Energy_PUany_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"_3MIP").c_str()));
    h_Energy_PUany_Eta_dRadius_3MIP[ij][1] = (TH1F*)(inFanyPU[ij]->Get(("ana/h_Energy_PUany_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"_3MIP").c_str()));
    h_Energy_MainEvt_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/h_Energy_MainEvt_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"_3MIP").c_str()));
    h_Energy_MainEvtany_Eta_dRadius_3MIP[ij][0] = (TH1F*)(inF[ij]->Get(("ana/h_Energy_MainEvtany_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"_3MIP").c_str()));
    h_Energy_MainEvtany_Eta_dRadius_3MIP[ij][1] = (TH1F*)(inFanyME[ij]->Get(("ana/h_Energy_MainEvtany_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"_3MIP").c_str()));
    h_Energy_Shared_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/h_Energy_Shared_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"_3MIP").c_str()));
    h_Energy_SharedME_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/h_Energy_SharedME_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"_3MIP").c_str()));
    h_Energy_SharedPU_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/h_Energy_SharedPU_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"_3MIP").c_str()));
    h_FractionEnergy_PU_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/h_FractionEnergy_PU_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"_3MIP").c_str()));
    h_FractionEnergy_PUany_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/h_FractionEnergy_PUany_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"_3MIP").c_str()));
    h_FractionEnergy_MainEvt_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/h_FractionEnergy_MainEvt_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"_3MIP").c_str()));
    h_FractionEnergy_MainEvtany_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/h_FractionEnergy_MainEvtany_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"_3MIP").c_str()));
    h_FractionEnergy_Shared_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/h_FractionEnergy_Shared_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"_3MIP").c_str()));
    h_FractionEnergy_SharedME_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/h_FractionEnergy_SharedME_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"_3MIP").c_str()));
    h_FractionEnergy_SharedPU_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/h_FractionEnergy_SharedPU_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"_3MIP").c_str()));
    */
    h_NumberHits_Shared_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/h_NumberHits_Shared_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"").c_str()));
    h_NumberHits_SharedME_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/h_NumberHits_SharedME_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"").c_str()));
    h_NumberHits_SharedPU_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/h_NumberHits_SharedPU_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"").c_str()));
    h_FractionHits_PU_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/h_FractionHits_PU_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"").c_str()));
    h_FractionHits_PUany_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/h_FractionHits_PUany_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"").c_str()));
    h_FractionHits_MainEvt_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/h_FractionHits_MainEvt_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"").c_str()));
    h_FractionHits_MainEvtany_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/h_FractionHits_MainEvtany_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"").c_str()));
    h_FractionHits_Shared_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/h_FractionHits_Shared_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"").c_str()));
    h_FractionHits_Shared_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/h_FractionHits_Shared_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"").c_str()));
    h_FractionHits_SharedME_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/h_FractionHits_SharedME_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"").c_str()));
    h_FractionHits_SharedPU_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/h_FractionHits_SharedPU_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"").c_str()));

    h_Energy_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/h_Energy_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"").c_str()));
    h_Energy_PU_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/h_Energy_PU_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"").c_str()));
    h_Energy_PUany_Eta_dRadius_3MIP[ij][0] = (TH1F*)(inF[ij]->Get(("ana/h_Energy_PUany_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"").c_str()));
    h_Energy_PUany_Eta_dRadius_3MIP[ij][1] = (TH1F*)(inFanyPU[ij]->Get(("ana/h_Energy_PUany_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"").c_str()));
    h_Energy_MainEvt_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/h_Energy_MainEvt_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"").c_str()));
    h_Energy_MainEvtany_Eta_dRadius_3MIP[ij][0] = (TH1F*)(inF[ij]->Get(("ana/h_Energy_MainEvtany_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"").c_str()));
    h_Energy_MainEvtany_Eta_dRadius_3MIP[ij][1] = (TH1F*)(inFanyME[ij]->Get(("ana/h_Energy_MainEvtany_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"").c_str()));
    h_Energy_Shared_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/h_Energy_Shared_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"").c_str()));
    h_Energy_SharedME_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/h_Energy_SharedME_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"").c_str()));
    h_Energy_SharedPU_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/h_Energy_SharedPU_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"").c_str()));
    h_FractionEnergy_PU_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/h_FractionEnergy_PU_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"").c_str()));
    h_FractionEnergy_PUany_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/h_FractionEnergy_PUany_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"").c_str()));
    h_FractionEnergy_MainEvt_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/h_FractionEnergy_MainEvt_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"").c_str()));
    h_FractionEnergy_MainEvtany_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/h_FractionEnergy_MainEvtany_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"").c_str()));
    h_FractionEnergy_Shared_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/h_FractionEnergy_Shared_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"").c_str()));
    h_FractionEnergy_SharedME_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/h_FractionEnergy_SharedME_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"").c_str()));
    h_FractionEnergy_SharedPU_Eta_dRadius_3MIP[ij] = (TH1F*)(inF[ij]->Get(("ana/h_FractionEnergy_SharedPU_Eta_dRadius_Eta"+etaBin.at(ij)+"_dRadius"+radius+"").c_str()));



    h_NumberHits_Eta_dRadius_3MIP[ij]->SetLineColor(iColors[0]);
    h_NumberHits_PU_Eta_dRadius_3MIP[ij]->SetLineColor(iColors[1]);
    h_NumberHits_PUany_Eta_dRadius_3MIP[ij][1]->SetLineColor(iColors[1]);
    h_NumberHits_MainEvt_Eta_dRadius_3MIP[ij]->SetLineColor(iColors[2]);
    h_NumberHits_MainEvtany_Eta_dRadius_3MIP[ij][1]->SetLineColor(iColors[2]);
    h_NumberHits_Shared_Eta_dRadius_3MIP[ij]->SetLineColor(iColors[3]);
    h_NumberHits_SharedME_Eta_dRadius_3MIP[ij]->SetLineColor(iColors[4]);
    h_NumberHits_SharedPU_Eta_dRadius_3MIP[ij]->SetLineColor(iColors[5]);

    h_FractionHits_PU_Eta_dRadius_3MIP[ij]->SetLineColor(iColors[1]);
    h_FractionHits_PUany_Eta_dRadius_3MIP[ij]->SetLineColor(iColors[1]);
    h_FractionHits_MainEvt_Eta_dRadius_3MIP[ij]->SetLineColor(iColors[2]);
    h_FractionHits_MainEvtany_Eta_dRadius_3MIP[ij]->SetLineColor(iColors[2]);
    h_FractionHits_Shared_Eta_dRadius_3MIP[ij]->SetLineColor(iColors[3]);
    h_FractionHits_SharedME_Eta_dRadius_3MIP[ij]->SetLineColor(iColors[4]);
    h_FractionHits_SharedPU_Eta_dRadius_3MIP[ij]->SetLineColor(iColors[5]);

    h_Energy_Eta_dRadius_3MIP[ij]->SetLineColor(iColors[0]);
    h_Energy_PU_Eta_dRadius_3MIP[ij]->SetLineColor(iColors[1]);
    h_Energy_PUany_Eta_dRadius_3MIP[ij][1]->SetLineColor(iColors[1]);
    h_Energy_MainEvt_Eta_dRadius_3MIP[ij]->SetLineColor(iColors[2]);
    h_Energy_MainEvtany_Eta_dRadius_3MIP[ij][1]->SetLineColor(iColors[2]);
    h_Energy_Shared_Eta_dRadius_3MIP[ij]->SetLineColor(iColors[3]);
    h_Energy_SharedME_Eta_dRadius_3MIP[ij]->SetLineColor(iColors[4]);
    h_Energy_SharedPU_Eta_dRadius_3MIP[ij]->SetLineColor(iColors[5]);
    h_FractionEnergy_PU_Eta_dRadius_3MIP[ij]->SetLineColor(iColors[1]);
    h_FractionEnergy_PUany_Eta_dRadius_3MIP[ij]->SetLineColor(iColors[1]);
    h_FractionEnergy_MainEvt_Eta_dRadius_3MIP[ij]->SetLineColor(iColors[2]);
    h_FractionEnergy_MainEvtany_Eta_dRadius_3MIP[ij]->SetLineColor(iColors[2]);
    h_FractionEnergy_Shared_Eta_dRadius_3MIP[ij]->SetLineColor(iColors[3]);
    h_FractionEnergy_SharedME_Eta_dRadius_3MIP[ij]->SetLineColor(iColors[4]);
    h_FractionEnergy_SharedPU_Eta_dRadius_3MIP[ij]->SetLineColor(iColors[5]);


    h_NumberHits_Eta_dRadius_3MIP[ij]->SetLineWidth(2);
    h_NumberHits_PU_Eta_dRadius_3MIP[ij]->SetLineWidth(2);
    h_NumberHits_PUany_Eta_dRadius_3MIP[ij][1]->SetLineWidth(2);
    h_NumberHits_MainEvt_Eta_dRadius_3MIP[ij]->SetLineWidth(2);
    h_NumberHits_MainEvtany_Eta_dRadius_3MIP[ij][1]->SetLineWidth(2);
    h_NumberHits_Shared_Eta_dRadius_3MIP[ij]->SetLineWidth(2);
    h_NumberHits_SharedME_Eta_dRadius_3MIP[ij]->SetLineWidth(2);
    h_NumberHits_SharedPU_Eta_dRadius_3MIP[ij]->SetLineWidth(2);
    h_FractionHits_PU_Eta_dRadius_3MIP[ij]->SetLineWidth(2);
    h_FractionHits_PUany_Eta_dRadius_3MIP[ij]->SetLineWidth(2);
    h_FractionHits_MainEvt_Eta_dRadius_3MIP[ij]->SetLineWidth(2);
    h_FractionHits_MainEvtany_Eta_dRadius_3MIP[ij]->SetLineWidth(2);
    h_FractionHits_Shared_Eta_dRadius_3MIP[ij]->SetLineWidth(2);
    h_FractionHits_SharedME_Eta_dRadius_3MIP[ij]->SetLineWidth(2);
    h_FractionHits_SharedPU_Eta_dRadius_3MIP[ij]->SetLineWidth(2);
    h_Energy_Eta_dRadius_3MIP[ij]->SetLineWidth(2);
    h_Energy_PU_Eta_dRadius_3MIP[ij]->SetLineWidth(2);
    h_Energy_PUany_Eta_dRadius_3MIP[ij][1]->SetLineWidth(2);
    h_Energy_MainEvt_Eta_dRadius_3MIP[ij]->SetLineWidth(2);
    h_Energy_MainEvtany_Eta_dRadius_3MIP[ij][1]->SetLineWidth(2);
    h_Energy_Shared_Eta_dRadius_3MIP[ij]->SetLineWidth(2);
    h_Energy_SharedME_Eta_dRadius_3MIP[ij]->SetLineWidth(2);
    h_Energy_SharedPU_Eta_dRadius_3MIP[ij]->SetLineWidth(2);
    h_FractionEnergy_PU_Eta_dRadius_3MIP[ij]->SetLineWidth(2);
    h_FractionEnergy_PUany_Eta_dRadius_3MIP[ij]->SetLineWidth(2);
    h_FractionEnergy_MainEvt_Eta_dRadius_3MIP[ij]->SetLineWidth(2);
    h_FractionEnergy_MainEvtany_Eta_dRadius_3MIP[ij]->SetLineWidth(2);
    h_FractionEnergy_Shared_Eta_dRadius_3MIP[ij]->SetLineWidth(2);
    h_FractionEnergy_SharedME_Eta_dRadius_3MIP[ij]->SetLineWidth(2);
    h_FractionEnergy_SharedPU_Eta_dRadius_3MIP[ij]->SetLineWidth(2);

    h_FractionHits_PU_Eta_dRadius_3MIP[ij]->Rebin(10);
    h_FractionHits_PUany_Eta_dRadius_3MIP[ij]->Rebin(10);
    h_FractionHits_MainEvt_Eta_dRadius_3MIP[ij]->Rebin(10);
    h_FractionHits_MainEvtany_Eta_dRadius_3MIP[ij]->Rebin(10);
    h_FractionHits_Shared_Eta_dRadius_3MIP[ij]->Rebin(10);
    h_FractionHits_SharedME_Eta_dRadius_3MIP[ij]->Rebin(10);
    h_FractionHits_SharedPU_Eta_dRadius_3MIP[ij]->Rebin(10);

    h_FractionEnergy_PU_Eta_dRadius_3MIP[ij]->Rebin(10);
    h_FractionEnergy_PUany_Eta_dRadius_3MIP[ij]->Rebin(10);
    h_FractionEnergy_MainEvt_Eta_dRadius_3MIP[ij]->Rebin(10);
    h_FractionEnergy_MainEvtany_Eta_dRadius_3MIP[ij]->Rebin(10);
    h_FractionEnergy_Shared_Eta_dRadius_3MIP[ij]->Rebin(10);
    h_FractionEnergy_SharedME_Eta_dRadius_3MIP[ij]->Rebin(15);
    h_FractionEnergy_SharedPU_Eta_dRadius_3MIP[ij]->Rebin(15);


    h_Energy_PUany_Eta_dRadius_3MIP[ij][1]->Rebin(2);
    h_Energy_MainEvtany_Eta_dRadius_3MIP[ij][1]->Rebin(2);

    
    h_NumberHits_PUany_Eta_dRadius_3MIP[ij][1]->Rebin(8);
    h_NumberHits_MainEvtany_Eta_dRadius_3MIP[ij][1]->Rebin(8);

    h_NumberHits_Eta_dRadius_3MIP[ij]->Rebin(5);
    h_NumberHits_PU_Eta_dRadius_3MIP[ij]->Rebin(5);
    h_NumberHits_MainEvt_Eta_dRadius_3MIP[ij]->Rebin(5);
    h_NumberHits_Shared_Eta_dRadius_3MIP[ij]->Rebin(5);
    h_NumberHits_SharedME_Eta_dRadius_3MIP[ij]->Rebin(5);
    h_NumberHits_SharedPU_Eta_dRadius_3MIP[ij]->Rebin(5);
    


  }

  std::cout << " ci sono ora stampo " << std::endl;


  TLegend *legTGMany = new TLegend(0.70,0.78,0.90,0.95,NULL,"brNDC");
  legTGMany->SetTextFont(42);
  legTGMany->SetTextSize(0.03);
  legTGMany->SetFillColor(kWhite);
  legTGMany->SetLineColor(kWhite);
  legTGMany->SetShadowColor(kWhite);
  legTGMany->AddEntry(h_NumberHits_PUany_Eta_dRadius_3MIP[0][1], "PU", "l");
  legTGMany->AddEntry(h_NumberHits_MainEvtany_Eta_dRadius_3MIP[0][1], "ME", "l");


  TLegend *legTGM = new TLegend(0.70,0.78,0.90,0.95,NULL,"brNDC");
  legTGM->SetTextFont(42);
  legTGM->SetTextSize(0.03);
  legTGM->SetFillColor(kWhite);
  legTGM->SetLineColor(kWhite);
  legTGM->SetShadowColor(kWhite);
  legTGM->AddEntry(h_NumberHits_Eta_dRadius_3MIP[0], "total hits.", "l");
  legTGM->AddEntry(h_NumberHits_PU_Eta_dRadius_3MIP[0], "PU only", "l");
  legTGM->AddEntry(h_NumberHits_MainEvt_Eta_dRadius_3MIP[0], "ME only", "l");
  legTGM->AddEntry(h_NumberHits_Shared_Eta_dRadius_3MIP[0], "shared", "l");
  legTGM->AddEntry(h_NumberHits_SharedME_Eta_dRadius_3MIP[0], "shared ME >3MIP", "l");
  legTGM->AddEntry(h_NumberHits_SharedPU_Eta_dRadius_3MIP[0], "shared PU >3MIP", "l");


  TLegend *legTGM_fr = new TLegend(0.70,0.78,0.90,0.95,NULL,"brNDC");
  legTGM_fr->SetTextFont(42);
  legTGM_fr->SetTextSize(0.03);
  legTGM_fr->SetFillColor(kWhite);
  legTGM_fr->SetLineColor(kWhite);
  legTGM_fr->SetShadowColor(kWhite);
  legTGM_fr->AddEntry(h_NumberHits_PU_Eta_dRadius_3MIP[0], "PU only", "l");
  legTGM_fr->AddEntry(h_NumberHits_MainEvt_Eta_dRadius_3MIP[0], "ME only", "l");
  legTGM_fr->AddEntry(h_NumberHits_Shared_Eta_dRadius_3MIP[0], "shared", "l");



  TLegend *legTGM_fr2 = new TLegend(0.60,0.78,0.85,0.90,NULL,"brNDC");
  legTGM_fr2->SetTextFont(42);
  legTGM_fr2->SetTextSize(0.04);
  legTGM_fr2->SetFillColor(kWhite);
  legTGM_fr2->SetLineColor(kWhite);
  legTGM_fr2->SetShadowColor(kWhite);
  legTGM_fr2->AddEntry(h_NumberHits_SharedME_Eta_dRadius_3MIP[0], "shared ME >3MIP", "l");
  legTGM_fr2->AddEntry(h_NumberHits_SharedPU_Eta_dRadius_3MIP[0], "shared PU >3MIP", "l");

  std::cout << " legends ok  " << std::endl;

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TLatex tL;
  tL.SetNDC();
  tL.SetTextSize(0.04);
  tL.SetTextFont(132);

  TLatex tFr;
  tFr.SetNDC();
  tFr.SetTextSize(0.04);
  tFr.SetTextColor(iColors[0]);
  tFr.SetTextFont(132);

  TLatex tFrPU;
  tFrPU.SetNDC();
  tFrPU.SetTextSize(0.04);
  tFrPU.SetTextColor(iColors[1]);
  tFrPU.SetTextFont(132);

  TLatex tFrME;
  tFrME.SetNDC();
  tFrME.SetTextSize(0.04);
  tFrME.SetTextColor(iColors[2]);
  tFrME.SetTextFont(132);

  TLatex tFrS;
  tFrS.SetNDC();
  tFrS.SetTextSize(0.04);
  tFrS.SetTextColor(iColors[3]);
  tFrS.SetTextFont(132);

  TLatex tFrSM;
  tFrSM.SetNDC();
  tFrSM.SetTextSize(0.04);
  tFrSM.SetTextColor(iColors[4]);
  tFrSM.SetTextFont(132);

  TLatex tFrSP;
  tFrSP.SetNDC();
  tFrSP.SetTextSize(0.04);
  tFrSP.SetTextColor(iColors[5]);
  tFrSP.SetTextFont(132);


  TCanvas* ch_NumberHitsAny[2];
  for(int iF=0; iF<2; ++iF){
    ch_NumberHitsAny[iF] = new TCanvas();
    ch_NumberHitsAny[iF]->cd();
    h_NumberHits_PUany_Eta_dRadius_3MIP[iF][1]->GetXaxis()->SetTitle(("n. hits per event etaBin "+etaBin.at(iF)).c_str() );
    if(radius == "0")h_NumberHits_PUany_Eta_dRadius_3MIP[iF][1]->GetXaxis()->SetRangeUser(0., 300.);    
    else if(radius == "1")h_NumberHits_PUany_Eta_dRadius_3MIP[iF][1]->GetXaxis()->SetRangeUser(0., 600.);    
    float yMax = 0;
    getXmax(h_NumberHits_MainEvtany_Eta_dRadius_3MIP[iF][1], yMax);
    h_NumberHits_PUany_Eta_dRadius_3MIP[iF][1]->GetYaxis()->SetRangeUser(0., yMax*2.);

    h_NumberHits_PUany_Eta_dRadius_3MIP[iF][1]->Draw();
    h_NumberHits_MainEvtany_Eta_dRadius_3MIP[iF][1]->Draw("same");


    tL.DrawLatex(0.2,0.8, Form("hits > 3 MIP and  #rho #leq%dcm", radiusR.at(std::stoi(radius))) );

    tFrPU.DrawLatex(0.65,0.7, Form("mean = %.2f +/- %.2f", h_NumberHits_PUany_Eta_dRadius_3MIP[iF][1]->GetMean(), 
				   h_NumberHits_PUany_Eta_dRadius_3MIP[iF][1]->GetRMS() / sqrt(h_NumberHits_PUany_Eta_dRadius_3MIP[iF][1]->GetEntries())));
    tFrME.DrawLatex(0.65,0.65, Form("mean = %.2f +/- %.2f", h_NumberHits_MainEvtany_Eta_dRadius_3MIP[iF][1]->GetMean(),
				    h_NumberHits_MainEvtany_Eta_dRadius_3MIP[iF][1]->GetRMS()/sqrt(h_NumberHits_MainEvtany_Eta_dRadius_3MIP[iF][1]->GetEntries())));

    legTGMany->Draw("same");
    ch_NumberHitsAny[iF]->Print((folder+"/h_numHitsAny_pdg_"+pdgID+"_rad_"+radius+"_eta_"+etaBin.at(iF)+".png").c_str(), "png");
    ch_NumberHitsAny[iF]->Print((folder+"/h_numHitsAny_pdg_"+pdgID+"_rad_"+radius+"_eta_"+etaBin.at(iF)+".pdf").c_str(), "pdf");
    ch_NumberHitsAny[iF]->Print((folder+"/h_numHitsAny_pdg_"+pdgID+"_rad_"+radius+"_eta_"+etaBin.at(iF)+".root").c_str(), "root");
  }




  TCanvas* ch_NumberHits[2];
  for(int iF=0; iF<2; ++iF){
    ch_NumberHits[iF] = new TCanvas();
    ch_NumberHits[iF]->cd();
    h_NumberHits_Eta_dRadius_3MIP[iF]->GetXaxis()->SetTitle(("n. hits per event etaBin "+etaBin.at(iF)).c_str() );
    if(radius == "0")h_NumberHits_Eta_dRadius_3MIP[iF]->GetXaxis()->SetRangeUser(0., 300.);    
    else if(radius == "1")h_NumberHits_Eta_dRadius_3MIP[iF]->GetXaxis()->SetRangeUser(0., 600.);    
    float yMax = 0;
    getXmax(h_NumberHits_SharedPU_Eta_dRadius_3MIP[iF], yMax);
    h_NumberHits_Eta_dRadius_3MIP[iF]->GetYaxis()->SetRangeUser(0., yMax*2.);

    h_NumberHits_Eta_dRadius_3MIP[iF]->Draw();
    h_NumberHits_PU_Eta_dRadius_3MIP[iF]->Draw("same");
    h_NumberHits_MainEvt_Eta_dRadius_3MIP[iF]->Draw("same");
    h_NumberHits_Shared_Eta_dRadius_3MIP[iF]->Draw("same");
    h_NumberHits_SharedME_Eta_dRadius_3MIP[iF]->Draw("same");
    h_NumberHits_SharedPU_Eta_dRadius_3MIP[iF]->Draw("same");


    tL.DrawLatex(0.2,0.8, Form("hits > 3 MIP and  #rho #leq%dcm", radiusR.at(std::stoi(radius))) );

    tFr.DrawLatex(0.65,0.75, Form("mean = %.2f +/- %.2f", h_NumberHits_Eta_dRadius_3MIP[iF]->GetMean(), 
				  h_NumberHits_Eta_dRadius_3MIP[iF]->GetRMS()/sqrt(h_NumberHits_Eta_dRadius_3MIP[iF]->GetEntries())));
    tFrPU.DrawLatex(0.65,0.7, Form("mean = %.2f +/- %.2f", h_NumberHits_PU_Eta_dRadius_3MIP[iF]->GetMean(), 
				   h_NumberHits_PU_Eta_dRadius_3MIP[iF]->GetRMS()/sqrt(h_NumberHits_PU_Eta_dRadius_3MIP[iF]->GetEntries())));
    tFrME.DrawLatex(0.65,0.65, Form("mean = %.2f +/- %.2f", h_NumberHits_MainEvt_Eta_dRadius_3MIP[iF]->GetMean(), 
				    h_NumberHits_MainEvt_Eta_dRadius_3MIP[iF]->GetRMS()/sqrt(h_NumberHits_MainEvt_Eta_dRadius_3MIP[iF]->GetEntries())));
    tFrS.DrawLatex(0.65,0.60, Form("mean = %.2f +/- %.2f", h_NumberHits_Shared_Eta_dRadius_3MIP[iF]->GetMean(), 
				   h_NumberHits_Shared_Eta_dRadius_3MIP[iF]->GetRMS()/sqrt(h_NumberHits_Shared_Eta_dRadius_3MIP[iF]->GetEntries())));
    tFrSM.DrawLatex(0.65,0.55, Form("mean = %.2f +/- %.2f", h_NumberHits_SharedME_Eta_dRadius_3MIP[iF]->GetMean(), 
				    h_NumberHits_SharedME_Eta_dRadius_3MIP[iF]->GetRMS()/sqrt(h_NumberHits_SharedME_Eta_dRadius_3MIP[iF]->GetEntries())));
    tFrSP.DrawLatex(0.65,0.50, Form("mean = %.2f +/- %.2f", h_NumberHits_SharedPU_Eta_dRadius_3MIP[iF]->GetMean(), 
				    h_NumberHits_SharedPU_Eta_dRadius_3MIP[iF]->GetRMS()/sqrt(h_NumberHits_SharedPU_Eta_dRadius_3MIP[iF]->GetEntries())));

    legTGM->Draw("same");
    ch_NumberHits[iF]->Print((folder+"/h_numHits_pdg_"+pdgID+"_rad_"+radius+"_eta_"+etaBin.at(iF)+".png").c_str(), "png");
    ch_NumberHits[iF]->Print((folder+"/h_numHits_pdg_"+pdgID+"_rad_"+radius+"_eta_"+etaBin.at(iF)+".pdf").c_str(), "pdf");
    ch_NumberHits[iF]->Print((folder+"/h_numHits_pdg_"+pdgID+"_rad_"+radius+"_eta_"+etaBin.at(iF)+".root").c_str(), "root");
  }


  TCanvas* ch_NumberHitsCut[2];
  for(int iF=0; iF<2; ++iF){
    ch_NumberHitsCut[iF] = new TCanvas();
    ch_NumberHitsCut[iF]->cd();
    h_NumberHits_PU_Eta_dRadius_3MIP[iF]->GetXaxis()->SetTitle(("n. hits per event etaBin "+etaBin.at(iF)).c_str() );
    if(radius == "0")h_NumberHits_PU_Eta_dRadius_3MIP[iF]->GetXaxis()->SetRangeUser(0., 300.);    
    else if(radius == "1")h_NumberHits_PU_Eta_dRadius_3MIP[iF]->GetXaxis()->SetRangeUser(0., 600.);    
    float yMax = 0;
    getXmax(h_NumberHits_Shared_Eta_dRadius_3MIP[iF], yMax);
    h_NumberHits_PU_Eta_dRadius_3MIP[iF]->GetYaxis()->SetRangeUser(0., yMax*2.);

    h_NumberHits_PU_Eta_dRadius_3MIP[iF]->Draw();
    h_NumberHits_MainEvt_Eta_dRadius_3MIP[iF]->Draw("same");
    h_NumberHits_Shared_Eta_dRadius_3MIP[iF]->Draw("same");


    tL.DrawLatex(0.2,0.8, Form("hits > 3 MIP and  #rho #leq%dcm", radiusR.at(std::stoi(radius))) );

    tFrPU.DrawLatex(0.65,0.7, Form("mean = %.2f +/- %.2f", h_NumberHits_PU_Eta_dRadius_3MIP[iF]->GetMean(), 
				   h_NumberHits_PU_Eta_dRadius_3MIP[iF]->GetRMS()/sqrt(h_NumberHits_PU_Eta_dRadius_3MIP[iF]->GetEntries())));
    tFrME.DrawLatex(0.65,0.65, Form("mean = %.2f +/- %.2f", h_NumberHits_MainEvt_Eta_dRadius_3MIP[iF]->GetMean(), 
				    h_NumberHits_MainEvt_Eta_dRadius_3MIP[iF]->GetRMS()/sqrt(h_NumberHits_MainEvt_Eta_dRadius_3MIP[iF]->GetEntries())));
    tFrS.DrawLatex(0.65,0.60, Form("mean = %.2f +/- %.2f", h_NumberHits_Shared_Eta_dRadius_3MIP[iF]->GetMean(), 
				   h_NumberHits_Shared_Eta_dRadius_3MIP[iF]->GetRMS()/sqrt(h_NumberHits_Shared_Eta_dRadius_3MIP[iF]->GetEntries())));

    legTGM_fr->Draw("same");
    ch_NumberHitsCut[iF]->Print((folder+"/h_numHitsCut_pdg_"+pdgID+"_rad_"+radius+"_eta_"+etaBin.at(iF)+".png").c_str(), "png");
    ch_NumberHitsCut[iF]->Print((folder+"/h_numHitsCut_pdg_"+pdgID+"_rad_"+radius+"_eta_"+etaBin.at(iF)+".pdf").c_str(), "pdf");
    ch_NumberHitsCut[iF]->Print((folder+"/h_numHitsCut_pdg_"+pdgID+"_rad_"+radius+"_eta_"+etaBin.at(iF)+".root").c_str(), "root");

    h_NumberHits_PU_Eta_dRadius_3MIP[iF]->GetYaxis()->SetRangeUser(0.5, yMax*200.);
    gPad->SetLogy();
    ch_NumberHitsCut[iF]->Print((folder+"/h_numHitsCut_pdg_"+pdgID+"_rad_"+radius+"_eta_"+etaBin.at(iF)+"_log.png").c_str(), "png");
    ch_NumberHitsCut[iF]->Print((folder+"/h_numHitsCut_pdg_"+pdgID+"_rad_"+radius+"_eta_"+etaBin.at(iF)+"_log.pdf").c_str(), "pdf");
    ch_NumberHitsCut[iF]->Print((folder+"/h_numHitsCut_pdg_"+pdgID+"_rad_"+radius+"_eta_"+etaBin.at(iF)+"_log.root").c_str(), "root");
  }



  //fraction
  TCanvas* ch_FractionHits[2];
  for(int iF=0; iF<2; ++iF){
    ch_FractionHits[iF] = new TCanvas();
    ch_FractionHits[iF]->cd();
    h_FractionHits_PU_Eta_dRadius_3MIP[iF]->GetXaxis()->SetTitle(("fraction hits per event etaBin "+etaBin.at(iF)).c_str() );
    //h_FractionHits_Eta_dRadius_3MIP[iF]->GetXaxis()->SetRangeUser(0., 500.);    
    float yMax = 0;
    getXmax(h_FractionHits_Shared_Eta_dRadius_3MIP[iF], yMax);
    h_FractionHits_PU_Eta_dRadius_3MIP[iF]->GetYaxis()->SetRangeUser(0., yMax*2.);

    h_FractionHits_PU_Eta_dRadius_3MIP[iF]->Draw();
    h_FractionHits_MainEvt_Eta_dRadius_3MIP[iF]->Draw("same");
    h_FractionHits_Shared_Eta_dRadius_3MIP[iF]->Draw("same");


    tL.DrawLatex(0.2,0.8, Form("hits > 3 MIP and  #rho #leq%dcm", radiusR.at(std::stoi(radius))) );

    tFrPU.DrawLatex(0.65,0.7, Form("mean = %.2f +/- %.2f", h_FractionHits_PU_Eta_dRadius_3MIP[iF]->GetMean(), 
				   h_FractionHits_PU_Eta_dRadius_3MIP[iF]->GetRMS()/sqrt(h_FractionHits_PU_Eta_dRadius_3MIP[iF]->GetEntries())));
    tFrME.DrawLatex(0.65,0.65, Form("mean = %.2f +/- %.2f", h_FractionHits_MainEvt_Eta_dRadius_3MIP[iF]->GetMean(), 
				    h_FractionHits_MainEvt_Eta_dRadius_3MIP[iF]->GetRMS()/sqrt(h_FractionHits_MainEvt_Eta_dRadius_3MIP[iF]->GetEntries())));
    tFrS.DrawLatex(0.65,0.60, Form("mean = %.2f +/- %.2f", h_FractionHits_Shared_Eta_dRadius_3MIP[iF]->GetMean(), 
				   h_FractionHits_Shared_Eta_dRadius_3MIP[iF]->GetRMS()/sqrt(h_FractionHits_Shared_Eta_dRadius_3MIP[iF]->GetEntries())));

    legTGM_fr->Draw("same");
    ch_FractionHits[iF]->Print((folder+"/h_fractionHits_pdg_"+pdgID+"_rad_"+radius+"_eta_"+etaBin.at(iF)+".png").c_str(), "png");
    ch_FractionHits[iF]->Print((folder+"/h_fractionHits_pdg_"+pdgID+"_rad_"+radius+"_eta_"+etaBin.at(iF)+".pdf").c_str(), "pdf");
    ch_FractionHits[iF]->Print((folder+"/h_fractionHits_pdg_"+pdgID+"_rad_"+radius+"_eta_"+etaBin.at(iF)+".root").c_str(), "root");
  }


  TCanvas* ch_FractionHitsAny[2];
  for(int iF=0; iF<2; ++iF){
    ch_FractionHitsAny[iF] = new TCanvas();
    ch_FractionHitsAny[iF]->cd();
    h_FractionHits_PUany_Eta_dRadius_3MIP[iF]->GetXaxis()->SetTitle(("fraction hits per event etaBin "+etaBin.at(iF)).c_str() );
    //h_FractionHits_Eta_dRadius_3MIP[iF]->GetXaxis()->SetRangeUser(0., 500.);    
    float yMax = 0;
    getXmax(h_FractionHits_PUany_Eta_dRadius_3MIP[iF], yMax);
    h_FractionHits_PUany_Eta_dRadius_3MIP[iF]->GetYaxis()->SetRangeUser(0., yMax*2.);

    h_FractionHits_PUany_Eta_dRadius_3MIP[iF]->Draw();
    h_FractionHits_MainEvtany_Eta_dRadius_3MIP[iF]->Draw("same");


    tL.DrawLatex(0.2,0.8, Form("hits > 3 MIP and  #rho #leq%dcm", radiusR.at(std::stoi(radius))) );

    tFrPU.DrawLatex(0.65,0.7, Form("mean = %.2f +/- %.2f", h_FractionHits_PUany_Eta_dRadius_3MIP[iF]->GetMean(), 
				   h_FractionHits_PUany_Eta_dRadius_3MIP[iF]->GetRMS()/sqrt(h_FractionHits_PUany_Eta_dRadius_3MIP[iF]->GetEntries())));
    tFrME.DrawLatex(0.65,0.65, Form("mean = %.2f +/- %.2f", h_FractionHits_MainEvtany_Eta_dRadius_3MIP[iF]->GetMean(), 
				    h_FractionHits_MainEvtany_Eta_dRadius_3MIP[iF]->GetRMS()/sqrt(h_FractionHits_MainEvtany_Eta_dRadius_3MIP[iF]->GetEntries())));

    legTGMany->Draw("same");
    ch_FractionHitsAny[iF]->Print((folder+"/h_fractionHitsAny_pdg_"+pdgID+"_rad_"+radius+"_eta_"+etaBin.at(iF)+".png").c_str(), "png");
    ch_FractionHitsAny[iF]->Print((folder+"/h_fractionHitsAny_pdg_"+pdgID+"_rad_"+radius+"_eta_"+etaBin.at(iF)+".pdf").c_str(), "pdf");
    ch_FractionHitsAny[iF]->Print((folder+"/h_fractionHitsAny_pdg_"+pdgID+"_rad_"+radius+"_eta_"+etaBin.at(iF)+".root").c_str(), "root");
  }

  
  TCanvas* ch_FractionHits_2[2];
  for(int iF=0; iF<2; ++iF){
    ch_FractionHits_2[iF] = new TCanvas();
    ch_FractionHits_2[iF]->cd();
    h_FractionHits_SharedME_Eta_dRadius_3MIP[iF]->GetXaxis()->SetTitle(("fraction of shared hits per event etaBin "+etaBin.at(iF)).c_str() );
    //h_FractionHits_Eta_dRadius_3MIP[iF]->GetXaxis()->SetRangeUser(0., 500.);    
    float yMax = 0;
    getXmax(h_FractionHits_SharedPU_Eta_dRadius_3MIP[iF], yMax);
    h_FractionHits_SharedME_Eta_dRadius_3MIP[iF]->GetYaxis()->SetRangeUser(0., yMax*2.);

    //    h_FractionHits_Shared_Eta_dRadius_3MIP[iF]->Draw();
    h_FractionHits_SharedME_Eta_dRadius_3MIP[iF]->Draw();
    h_FractionHits_SharedPU_Eta_dRadius_3MIP[iF]->Draw("same");


    tL.DrawLatex(0.2,0.8, Form("hits > 3 MIP and  #rho #leq%dcm", radiusR.at(std::stoi(radius))) );

    //    tFrS.DrawLatex(0.65,0.56, Form("mean = %.2f +/- %.2f", h_FractionHits_Shared_Eta_dRadius_3MIP[iF]->GetMean()));
    tFrSM.DrawLatex(0.65,0.55, Form("mean = %.2f +/- %.2f", h_FractionHits_SharedME_Eta_dRadius_3MIP[iF]->GetMean(), 
				    h_FractionHits_SharedME_Eta_dRadius_3MIP[iF]->GetRMS()/sqrt(h_FractionHits_SharedME_Eta_dRadius_3MIP[iF]->GetEntries())));
    tFrSP.DrawLatex(0.65,0.50, Form("mean = %.2f +/- %.2f", h_FractionHits_SharedPU_Eta_dRadius_3MIP[iF]->GetMean(), 
				    h_FractionHits_SharedPU_Eta_dRadius_3MIP[iF]->GetRMS()/sqrt(h_FractionHits_SharedPU_Eta_dRadius_3MIP[iF]->GetEntries())));

    legTGM_fr2->Draw("same");
    ch_FractionHits_2[iF]->Print((folder+"/h_fractionHitsShared_pdg_"+pdgID+"_rad_"+radius+"_eta_"+etaBin.at(iF)+".png").c_str(), "png");
    ch_FractionHits_2[iF]->Print((folder+"/h_fractionHitsShared_pdg_"+pdgID+"_rad_"+radius+"_eta_"+etaBin.at(iF)+".pdf").c_str(), "pdf");
    ch_FractionHits_2[iF]->Print((folder+"/h_fractionHitsShared_pdg_"+pdgID+"_rad_"+radius+"_eta_"+etaBin.at(iF)+".root").c_str(), "root");
  }

  ////////////energy
  TCanvas* ch_EnergyAny[2];
  for(int iF=0; iF<2; ++iF){
    ch_EnergyAny[iF] = new TCanvas();
    ch_EnergyAny[iF]->cd();
    h_Energy_PUany_Eta_dRadius_3MIP[iF][1]->GetXaxis()->SetTitle(("pT sum etaBin "+etaBin.at(iF)).c_str() );
    h_Energy_PUany_Eta_dRadius_3MIP[iF][1]->GetXaxis()->SetRangeUser(0., 20.);
    float yMax = 0;
    getXmax(h_Energy_MainEvtany_Eta_dRadius_3MIP[iF][1], yMax);
    h_Energy_PUany_Eta_dRadius_3MIP[iF][1]->GetYaxis()->SetRangeUser(0., yMax*2.);

    h_Energy_PUany_Eta_dRadius_3MIP[iF][1]->Draw();
    h_Energy_MainEvtany_Eta_dRadius_3MIP[iF][1]->Draw("same");


    tL.DrawLatex(0.2,0.8, Form("hits > 3 MIP and  #rho #leq%dcm", radiusR.at(std::stoi(radius))) );

    tFrPU.DrawLatex(0.65,0.7, Form("mean = %.2f +/- %.2f", h_Energy_PUany_Eta_dRadius_3MIP[iF][1]->GetMean(), 
				   h_Energy_PUany_Eta_dRadius_3MIP[iF][1]->GetRMS()/sqrt(h_Energy_PUany_Eta_dRadius_3MIP[iF][1]->GetEntries())));
    tFrME.DrawLatex(0.65,0.65, Form("mean = %.2f +/- %.2f", h_Energy_MainEvtany_Eta_dRadius_3MIP[iF][1]->GetMean(), 
				    h_Energy_MainEvtany_Eta_dRadius_3MIP[iF][1]->GetRMS()/sqrt(h_Energy_MainEvtany_Eta_dRadius_3MIP[iF][1]->GetEntries())));

    legTGMany->Draw("same");
    ch_EnergyAny[iF]->Print((folder+"/h_energyAny_pdg_"+pdgID+"_rad_"+radius+"_eta_"+etaBin.at(iF)+".png").c_str(), "png");
    ch_EnergyAny[iF]->Print((folder+"/h_energyAny_pdg_"+pdgID+"_rad_"+radius+"_eta_"+etaBin.at(iF)+".pdf").c_str(), "pdf");
    ch_EnergyAny[iF]->Print((folder+"/h_energyAny_pdg_"+pdgID+"_rad_"+radius+"_eta_"+etaBin.at(iF)+".root").c_str(), "root");
  }

  TCanvas* ch_Energy[2];
  for(int iF=0; iF<2; ++iF){
    ch_Energy[iF] = new TCanvas();
    ch_Energy[iF]->cd();
    h_Energy_Eta_dRadius_3MIP[iF]->GetXaxis()->SetTitle(("pT sum  etaBin "+etaBin.at(iF)).c_str() );
    h_Energy_Eta_dRadius_3MIP[iF]->GetXaxis()->SetRangeUser(0., 20.);    
    float yMax = 0;
    getXmax(h_Energy_SharedPU_Eta_dRadius_3MIP[iF], yMax);
    h_Energy_Eta_dRadius_3MIP[iF]->GetYaxis()->SetRangeUser(0., yMax*2.);

    h_Energy_Eta_dRadius_3MIP[iF]->Draw();
    h_Energy_PU_Eta_dRadius_3MIP[iF]->Draw("same");
    h_Energy_MainEvt_Eta_dRadius_3MIP[iF]->Draw("same");
    h_Energy_Shared_Eta_dRadius_3MIP[iF]->Draw("same");
    h_Energy_SharedME_Eta_dRadius_3MIP[iF]->Draw("same");
    h_Energy_SharedPU_Eta_dRadius_3MIP[iF]->Draw("same");


    tL.DrawLatex(0.2,0.8, Form("hits > 3 MIP and  #rho #leq%dcm", radiusR.at(std::stoi(radius))) );

    tFr.DrawLatex(0.65,0.75, Form("mean = %.2f +/- %.2f", h_Energy_Eta_dRadius_3MIP[iF]->GetMean(), 
				  h_Energy_Eta_dRadius_3MIP[iF]->GetRMS()/sqrt(h_Energy_Eta_dRadius_3MIP[iF]->GetEntries())));
    tFrPU.DrawLatex(0.65,0.7, Form("mean = %.2f +/- %.2f", h_Energy_PU_Eta_dRadius_3MIP[iF]->GetMean(), 
				   h_Energy_PU_Eta_dRadius_3MIP[iF]->GetRMS()/sqrt(h_Energy_PU_Eta_dRadius_3MIP[iF]->GetEntries())));
    tFrME.DrawLatex(0.65,0.65, Form("mean = %.2f +/- %.2f", h_Energy_MainEvt_Eta_dRadius_3MIP[iF]->GetMean(), 
				    h_Energy_MainEvt_Eta_dRadius_3MIP[iF]->GetRMS()/sqrt(h_Energy_MainEvt_Eta_dRadius_3MIP[iF]->GetEntries())));
    tFrS.DrawLatex(0.65,0.60, Form("mean = %.2f +/- %.2f", h_Energy_Shared_Eta_dRadius_3MIP[iF]->GetMean(),
				   h_Energy_Shared_Eta_dRadius_3MIP[iF]->GetRMS()/sqrt(h_Energy_Shared_Eta_dRadius_3MIP[iF]->GetEntries())));
    tFrSM.DrawLatex(0.65,0.55, Form("mean = %.2f +/- %.2f", h_Energy_SharedME_Eta_dRadius_3MIP[iF]->GetMean(), 
				    h_Energy_SharedME_Eta_dRadius_3MIP[iF]->GetRMS()/sqrt(h_Energy_SharedME_Eta_dRadius_3MIP[iF]->GetEntries())));
    tFrSP.DrawLatex(0.65,0.50, Form("mean = %.2f +/- %.2f", h_Energy_SharedPU_Eta_dRadius_3MIP[iF]->GetMean(), 
				    h_Energy_SharedPU_Eta_dRadius_3MIP[iF]->GetRMS()/sqrt(h_Energy_SharedPU_Eta_dRadius_3MIP[iF]->GetEntries())));

    legTGM->Draw("same");
    ch_Energy[iF]->Print((folder+"/h_energy_pdg_"+pdgID+"_rad_"+radius+"_eta_"+etaBin.at(iF)+".png").c_str(), "png");
    ch_Energy[iF]->Print((folder+"/h_energy_pdg_"+pdgID+"_rad_"+radius+"_eta_"+etaBin.at(iF)+".pdf").c_str(), "pdf");
    ch_Energy[iF]->Print((folder+"/h_energy_pdg_"+pdgID+"_rad_"+radius+"_eta_"+etaBin.at(iF)+".root").c_str(), "root");
  }

  //fraction
  TCanvas* ch_FractionEnergy[2];
  for(int iF=0; iF<2; ++iF){
    ch_FractionEnergy[iF] = new TCanvas();
    ch_FractionEnergy[iF]->cd();
    h_FractionEnergy_PU_Eta_dRadius_3MIP[iF]->GetXaxis()->SetTitle(("fraction energy per event etaBin "+etaBin.at(iF)).c_str() );
    //h_FractionEnergy_Eta_dRadius_3MIP[iF]->GetXaxis()->SetRangeUser(0., 500.);    
    float yMax = 0;
    getXmax(h_FractionEnergy_Shared_Eta_dRadius_3MIP[iF], yMax);
    h_FractionEnergy_PU_Eta_dRadius_3MIP[iF]->GetYaxis()->SetRangeUser(0., yMax*2.);

    h_FractionEnergy_PU_Eta_dRadius_3MIP[iF]->Draw();
    h_FractionEnergy_MainEvt_Eta_dRadius_3MIP[iF]->Draw("same");
    h_FractionEnergy_Shared_Eta_dRadius_3MIP[iF]->Draw("same");

    tL.DrawLatex(0.2,0.8, Form("hits > 3 MIP and  #rho #leq%dcm", radiusR.at(std::stoi(radius))) );

    tFrPU.DrawLatex(0.65,0.7, Form("mean = %.2f +/- %.2f", h_FractionEnergy_PU_Eta_dRadius_3MIP[iF]->GetMean(), 
				   h_FractionEnergy_PU_Eta_dRadius_3MIP[iF]->GetRMS()/sqrt(h_FractionEnergy_PU_Eta_dRadius_3MIP[iF]->GetEntries())));
    tFrME.DrawLatex(0.65,0.65, Form("mean = %.2f +/- %.2f", h_FractionEnergy_MainEvt_Eta_dRadius_3MIP[iF]->GetMean(), 
				    h_FractionEnergy_MainEvt_Eta_dRadius_3MIP[iF]->GetRMS()/sqrt(h_FractionEnergy_MainEvt_Eta_dRadius_3MIP[iF]->GetEntries())));
    tFrS.DrawLatex(0.65,0.60, Form("mean = %.2f +/- %.2f", h_FractionEnergy_Shared_Eta_dRadius_3MIP[iF]->GetMean(), 
				   h_FractionEnergy_Shared_Eta_dRadius_3MIP[iF]->GetRMS()/sqrt(h_FractionEnergy_Shared_Eta_dRadius_3MIP[iF]->GetEntries())));

    legTGM_fr->Draw("same");
    ch_FractionEnergy[iF]->Print((folder+"/h_fractionEnergy_pdg_"+pdgID+"_rad_"+radius+"_eta_"+etaBin.at(iF)+".png").c_str(), "png");
    ch_FractionEnergy[iF]->Print((folder+"/h_fractionEnergy_pdg_"+pdgID+"_rad_"+radius+"_eta_"+etaBin.at(iF)+".pdf").c_str(), "pdf");
    ch_FractionEnergy[iF]->Print((folder+"/h_fractionEnergy_pdg_"+pdgID+"_rad_"+radius+"_eta_"+etaBin.at(iF)+".root").c_str(), "root");

    h_FractionEnergy_PU_Eta_dRadius_3MIP[iF]->GetYaxis()->SetRangeUser(0.5, yMax*200.);
    gPad->SetLogy();
    ch_FractionEnergy[iF]->Print((folder+"/h_fractionEnergy_pdg_"+pdgID+"_rad_"+radius+"_eta_"+etaBin.at(iF)+"_log.png").c_str(), "png");
    ch_FractionEnergy[iF]->Print((folder+"/h_fractionEnergy_pdg_"+pdgID+"_rad_"+radius+"_eta_"+etaBin.at(iF)+"_log.pdf").c_str(), "pdf");
    ch_FractionEnergy[iF]->Print((folder+"/h_fractionEnergy_pdg_"+pdgID+"_rad_"+radius+"_eta_"+etaBin.at(iF)+"_log.root").c_str(), "root");
  }


  TCanvas* ch_FractionEnergyAny[2];
  for(int iF=0; iF<2; ++iF){
    ch_FractionEnergyAny[iF] = new TCanvas();
    ch_FractionEnergyAny[iF]->cd();
    h_FractionEnergy_PUany_Eta_dRadius_3MIP[iF]->GetXaxis()->SetTitle(("fraction energy per event etaBin "+etaBin.at(iF)).c_str() );
    //h_FractionEnergy_Eta_dRadius_3MIP[iF]->GetXaxis()->SetRangeUser(0., 500.);    
    float yMax = 0;
    getXmax(h_FractionEnergy_PUany_Eta_dRadius_3MIP[iF], yMax);
    h_FractionEnergy_PUany_Eta_dRadius_3MIP[iF]->GetYaxis()->SetRangeUser(0., yMax*2.);

    h_FractionEnergy_PUany_Eta_dRadius_3MIP[iF]->Draw();
    h_FractionEnergy_MainEvtany_Eta_dRadius_3MIP[iF]->Draw("same");


    tL.DrawLatex(0.2,0.8, Form("hits > 3 MIP and  #rho #leq%dcm", radiusR.at(std::stoi(radius))) );

    tFrPU.DrawLatex(0.65,0.7, Form("mean = %.2f +/- %.2f", h_FractionEnergy_PUany_Eta_dRadius_3MIP[iF]->GetMean(), 
				   h_FractionEnergy_PUany_Eta_dRadius_3MIP[iF]->GetRMS()/sqrt(h_FractionEnergy_PUany_Eta_dRadius_3MIP[iF]->GetEntries())));
    tFrME.DrawLatex(0.65,0.65, Form("mean = %.2f +/- %.2f", h_FractionEnergy_MainEvtany_Eta_dRadius_3MIP[iF]->GetMean(), 
				    h_FractionEnergy_MainEvtany_Eta_dRadius_3MIP[iF]->GetRMS()/sqrt(h_FractionEnergy_MainEvtany_Eta_dRadius_3MIP[iF]->GetEntries())));

    legTGMany->Draw("same");
    ch_FractionEnergyAny[iF]->Print((folder+"/h_fractionEnergyAny_pdg_"+pdgID+"_rad_"+radius+"_eta_"+etaBin.at(iF)+".png").c_str(), "png");
    ch_FractionEnergyAny[iF]->Print((folder+"/h_fractionEnergyAny_pdg_"+pdgID+"_rad_"+radius+"_eta_"+etaBin.at(iF)+".pdf").c_str(), "pdf");
    ch_FractionEnergyAny[iF]->Print((folder+"/h_fractionEnergyAny_pdg_"+pdgID+"_rad_"+radius+"_eta_"+etaBin.at(iF)+".root").c_str(), "root");
  }


  TCanvas* ch_FractionEnergy_2[2];
  for(int iF=0; iF<2; ++iF){
    ch_FractionEnergy_2[iF] = new TCanvas();
    ch_FractionEnergy_2[iF]->cd();
    h_FractionEnergy_SharedME_Eta_dRadius_3MIP[iF]->GetXaxis()->SetTitle(("fraction of shared energy per event etaBin "+etaBin.at(iF)).c_str() );
    //h_FractionEnergy_Eta_dRadius_3MIP[iF]->GetXaxis()->SetRangeUser(0., 500.);    
    float yMax = 0;
    getXmax(h_FractionEnergy_SharedPU_Eta_dRadius_3MIP[iF], yMax);
    h_FractionEnergy_SharedME_Eta_dRadius_3MIP[iF]->GetYaxis()->SetRangeUser(0., yMax*2.);

    //    h_FractionEnergy_Shared_Eta_dRadius_3MIP[iF]->Draw();
    h_FractionEnergy_SharedME_Eta_dRadius_3MIP[iF]->Draw();
    h_FractionEnergy_SharedPU_Eta_dRadius_3MIP[iF]->Draw("same");


    tL.DrawLatex(0.2,0.8, Form("hits > 3 MIP and  #rho #leq%dcm", radiusR.at(std::stoi(radius))) );

    //tFrS.DrawLatex(0.65,0.56, Form("mean = %.2f +/- %.2f", h_FractionEnergy_Shared_Eta_dRadius_3MIP[iF]->GetMean()));
    tFrSM.DrawLatex(0.65,0.55, Form("mean = %.2f +/- %.2f", h_FractionEnergy_SharedME_Eta_dRadius_3MIP[iF]->GetMean(), 
				    h_FractionEnergy_SharedME_Eta_dRadius_3MIP[iF]->GetRMS()/sqrt(h_FractionEnergy_SharedME_Eta_dRadius_3MIP[iF]->GetEntries())));
    tFrSP.DrawLatex(0.65,0.50, Form("mean = %.2f +/- %.2f", h_FractionEnergy_SharedPU_Eta_dRadius_3MIP[iF]->GetMean(), 
				    h_FractionEnergy_SharedPU_Eta_dRadius_3MIP[iF]->GetRMS()/sqrt(h_FractionEnergy_SharedPU_Eta_dRadius_3MIP[iF]->GetEntries())));

    legTGM_fr2->Draw("same");
    ch_FractionEnergy_2[iF]->Print((folder+"/h_fractionEnergyShared_pdg_"+pdgID+"_rad_"+radius+"_eta_"+etaBin.at(iF)+".png").c_str(), "png");
    ch_FractionEnergy_2[iF]->Print((folder+"/h_fractionEnergyShared_pdg_"+pdgID+"_rad_"+radius+"_eta_"+etaBin.at(iF)+".pdf").c_str(), "pdf");
    ch_FractionEnergy_2[iF]->Print((folder+"/h_fractionEnergyShared_pdg_"+pdgID+"_rad_"+radius+"_eta_"+etaBin.at(iF)+".root").c_str(), "root");

    h_FractionEnergy_SharedME_Eta_dRadius_3MIP[iF]->GetYaxis()->SetRangeUser(0.5, yMax*200.);
    gPad->SetLogy();
    ch_FractionEnergy_2[iF]->Print((folder+"/h_fractionEnergyShared_pdg_"+pdgID+"_rad_"+radius+"_eta_"+etaBin.at(iF)+"_log.png").c_str(), "png");
    ch_FractionEnergy_2[iF]->Print((folder+"/h_fractionEnergyShared_pdg_"+pdgID+"_rad_"+radius+"_eta_"+etaBin.at(iF)+"_log.pdf").c_str(), "pdf");
    ch_FractionEnergy_2[iF]->Print((folder+"/h_fractionEnergyShared_pdg_"+pdgID+"_rad_"+radius+"_eta_"+etaBin.at(iF)+"_log.root").c_str(), "root");
  }



  return;

}
