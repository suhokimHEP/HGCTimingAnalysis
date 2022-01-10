//g++  -o fitMIP  fitMIP.cpp `root-config --cflags --glibs`

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>

#include "TROOT.h"
#include "TSystem.h"
#include "TKey.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TH2F.h"
#include "TF1.h"
#include "TString.h"
#include "TCut.h"
#include "TMath.h"
#include "TApplication.h"
#include "TError.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TPad.h"
#include "TMultiGraph.h"
#include "TStyle.h"
#include "TChain.h"
#include "TPaveStats.h"

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <iostream>
#include <string>
#include "TStopwatch.h"
#include <vector>
#include <ROOT/RDF/RInterface.hxx>
#include <ROOT/RDF/InterfaceUtils.hxx>

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooCBShape.h"
#include "RooArgusBG.h"
#include "RooFFTConvPdf.h"
#include "RooAddPdf.h"
#include "RooWorkspace.h"
#include "RooFitResult.h"
#include "TText.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"

using namespace RooFit;
using namespace std;
//using namespace ROOT::VecOps;

void fitRegLine(){
 TText* title = new TText(1,1,"") ;
 title->SetTextSize(0.04);
 title->SetTextColor(kBlack);
 title->SetTextAlign(11);
 title->SetTextFont(62);
 
 TText* extra = new TText(1,1,"") ;
 extra->SetTextSize(0.03);
 extra->SetTextColor(kBlack);
 extra->SetTextAlign(11);
 extra->SetTextFont(52);
   
 TText* extra2 = new TText(1,1,"") ;
 extra2->SetTextSize(0.025);
 extra2->SetTextColor(kBlack);
 extra2->SetTextAlign(11);
 extra2->SetTextFont(62);
 TText* Punzi = new TText(1,1,"") ;
 Punzi->SetTextSize(0.03);
 Punzi->SetTextColor(kBlack);
 Punzi->SetTextAlign(31);
 Punzi->SetTextFont(62);
 
 TText* Ent = new TText(1,1,"") ;
 Ent->SetTextSize(0.03);
 Ent->SetTextColor(kBlack);
 Ent->SetTextAlign(31);
 Ent->SetTextFont(62);
 
 TText* SEnt = new TText(1,1,"") ;
 SEnt->SetTextSize(0.03);
 SEnt->SetTextColor(kBlack);
 SEnt->SetTextAlign(31);
 SEnt->SetTextFont(62);
  gROOT->Reset();
  gROOT->Macro("./setStyle.C");

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptTitle(0);
 
  std::string nametag = "startup";
  //std::string nametag = "eol";
  std::string line;
  std::string block;
  std::vector<vector<float>> meanmap;
  std::vector<vector<float>> errmap;

  std::string meanstr = "tex/mlandau_"+nametag+".tex";
  std::ifstream meanfile(meanstr);
  if (!meanfile.is_open()) {
    std::cout << "Unable to open '" << meanstr << "'" << std::endl;
  }
  while (getline(meanfile, line)) {
    std::stringstream linestream(line);
  std::vector<float> meanline;
  while (linestream>>block) {
	meanline.push_back(stof(block)); 
	}
	meanmap.push_back(meanline); 
	}

  std::string errstr = "tex/mlandauerr_"+nametag+".tex";
  std::ifstream errfile(errstr);
  if (!errfile.is_open()) {
    std::cout << "Unable to open '" << errstr << "'" << std::endl;
  }
  while (getline(errfile, line)) {
    std::stringstream linestream(line);
  std::vector<float> errline;
  while (linestream>>block) {
	errline.push_back(stof(block)); 
	}
	errmap.push_back(errline); 
	}
  int start=1;
   TCanvas *c1 = new TCanvas("c1","A Simple Graph Example",700,700);
    Double_t x[100], y[100], e[100], r[100], m[100];  
  for (Int_t i=0;i<meanmap.size();i++) {
   int refnum =meanmap[i][0];
   int graphcount =0;
   for (Int_t j=start;j<meanmap[i].size();j++) {
     y[j-start] = meanmap[i][j];
     e[j-start] = errmap[i][j];
     r[j-start] = e[j-start]/y[j-start];
     m[j-start] = r[j-start]/r[0];
     x[j-start] = refnum+j-start;
	graphcount++;
   	}
   TGraph* gr = new TGraph(graphcount,x,y);
	gr->Fit("pol1");
	gr->GetYaxis()->SetTitleOffset(1.75);
	gr->GetYaxis()->SetTitleSize(0.04);
	gr->GetYaxis()->SetLabelSize(0.04);
	gr->GetXaxis()->SetTitleSize(0.04);
	gr->GetXaxis()->SetLabelSize(0.04);
	gr->GetYaxis()->SetTitle("m_landau(MIP peak)");
	gr->GetXaxis()->SetTitle("i#eta");
	float correlation = gr->GetCorrelationFactor();
	if(correlation>0.3) gr->Draw("AP");
	else gr->Draw("P");
     title->DrawTextNDC(0.2,0.96,"CMS");
     extra->DrawTextNDC(0.28,0.96,"Work In Progress");
     extra2->DrawTextNDC(0.51,0.96,Form("Phase2 m_landau vs Ring %s",nametag.c_str()));
	Ent->DrawTextNDC(0.8,0.81,Form("L=%d, R^{2}:%f",i+37,correlation));	
      c1->Print(Form("Corr/Corr_%s_L%d.png",nametag.c_str(),i+37), ".png"); 
   TGraph* gr2 = new TGraph(graphcount,x,r);
	gr2->Fit("pol1");
	gr2->GetYaxis()->SetTitleOffset(2.25);
	gr2->GetYaxis()->SetTitleSize(0.04);
	gr2->GetYaxis()->SetLabelSize(0.04);
	gr2->GetXaxis()->SetTitleSize(0.04);
	gr2->GetXaxis()->SetLabelSize(0.04);
	gr2->GetYaxis()->SetTitle("#frac{e_landau}{m_landau}(MIP peak)");
	gr2->GetXaxis()->SetTitle("i#eta");
	float ecorrelation = gr2->GetCorrelationFactor();
	if(ecorrelation>0.3) gr2->Draw("AP");
	else gr2->Draw("P");
     title->DrawTextNDC(0.2,0.96,"CMS");
     extra->DrawTextNDC(0.28,0.96,"Work In Progress");
     extra2->DrawTextNDC(0.51,0.96,Form("Phase2 Resolution vs Ring %s",nametag.c_str()));
	Ent->DrawTextNDC(0.8,0.81,Form("L=%d, R^{2}:%f",i+37,ecorrelation));	
      c1->Print(Form("Corr/ErrCorr_%s_L%d.png",nametag.c_str(),i+37), ".png"); 
   TGraph* gr3 = new TGraph(graphcount,x,m);
	gr3->Fit("pol1");
	gr3->GetYaxis()->SetTitleOffset(3.25);
	gr3->GetYaxis()->SetTitleSize(0.02);
	gr3->GetYaxis()->SetLabelSize(0.04);
	gr3->GetXaxis()->SetTitleSize(0.04);
	gr3->GetXaxis()->SetLabelSize(0.04);
	gr3->GetYaxis()->SetTitle("#frac{#frac{e_landau}{m_landau}}{#frac{e_landau}{m_landau}_{InnMost}} (MIP peak)");
	gr3->GetXaxis()->SetTitle("i#eta");
	float mcorrelation = gr3->GetCorrelationFactor();
	if(mcorrelation>0.3) gr3->Draw("AP");
	else gr3->Draw("P");
     title->DrawTextNDC(0.2,0.96,"CMS");
     extra->DrawTextNDC(0.28,0.96,"Work In Progress");
     extra2->DrawTextNDC(0.51,0.96,Form("Phase2 Resolution(unit-ref) vs Ring %s",nametag.c_str()));
	Ent->DrawTextNDC(0.8,0.81,Form("L=%d, R^{2}:%f",i+37,mcorrelation));	
      c1->Print(Form("Corr/RefErrCorr_%s_L%d.png",nametag.c_str(),i+37), ".png"); 
  	// gr->Draw("AC*");
   }



 /*
  std::string errstr = "tex/mlandau_"+nametag+".tex";
  std::ifstream errfile(errstr);
  if (!errfile.is_open()) {
    std::cout << "Unable to open '" << errstr << "'" << std::endl;
  }


 TTree* newT = (TTree*)inF->Get("MIPtree");
 std::map<std::pair<int,int>, int> rentries;
 std::map<std::pair<int,int>, int> MIPentries;
rentries.clear();
	TCanvas*canvas1 = new TCanvas();
 FILE * outfulltable;
      outfulltable = fopen(Form("%s.tex",outtag.c_str()),"w");

  for(int iL=0; iL <= 4; ++iL){
  for(int Rn=16; Rn <= 38; ++Rn){
        int entries = newT->GetEntries(Form("MIP_iR==%d&&MIP_layer==%d",Rn,iL));
        //int Mentries = newT->GetEntries(Form("MIP_iR==%d&&MIP_layer==%d&&MIP_val<0.543",Rn,iL));
	newT->Draw("MIP_val>>temp(150,-3,3)",Form("MIP_iR==%d&&MIP_layer==%d",Rn,iL));
        TH1F* fithist = (TH1F*)newT->GetHistogram();
	fithist->Draw("hist");
	fithist->GetYaxis()->SetTitle("Events/0.04");
	fithist->GetXaxis()->SetTitle("MIP value");
	fithist->Draw("hist");
	if (fithist->Integral(0,-1)==0) {
        fprintf (outfulltable, "0.  ");		
	continue;

	}
	TF1 *func = new TF1("fit","gaus",-3,3);
	func->SetParameters(1,0.);
	func->SetParLimits(1,-.05,.05);
	fithist->Fit("fit");
	double c0 = func->GetParameter(0);	
	double c1 = func->GetParameter(1);
	double c2 = func->GetParameter(2);	
        fprintf (outfulltable, "%f   ", c2);		
	fithist->GetYaxis()->SetRangeUser(0,c0);
	fithist->Draw("same");
	printf("%d:%d:%d\n",iL,Rn,entries);
	printf("%d:%d:%f\n",iL,Rn,c0);
	printf("%d:%d:%f\n",iL,Rn,c1);
	printf("%d:%d:%f\n",iL,Rn,c2);






     char lumistring2 [50] = {0};
     int dummy2; 
     char lumistring3 [50] = {0};
     int dummy3; 
     char lumistring4 [50] = {0};
     int dummy4; 
	dummy2=0;
	dummy3=0;
	dummy4=0;




	dummy2=sprintf (lumistring2, "%3.5f", c0);
	dummy3=sprintf (lumistring3, "%3.5f", c1);
	dummy4=sprintf (lumistring4, "%3.5f", c2);
	Ent->DrawTextNDC(0.8,0.31,"amplitude:"+(TString)lumistring2);	
	Punzi->DrawTextNDC(0.8,0.21,"mean:"+(TString)lumistring3);	
	SEnt->DrawTextNDC(0.8,0.26,"spread:"+(TString)lumistring4);	
     title->DrawTextNDC(0.2,0.96,"CMS");
     extra->DrawTextNDC(0.3,0.96,"Work In Progress");
     extra2->DrawTextNDC(0.51,0.96,"Phase2 Noise fit");

      canvas1->Print(Form("NeutrinoFits/NeutrinoFits_%s_L%d_R%d_nEvts%d.png",nametag.c_str(),iL,Rn,entries), ".png"); 
	}
        fprintf (outfulltable, "\n");		
   }*/
}








	//fithist->Fit("gaus");
	//TF1 *g = (TF1*)fithist->GetListOfFunctions()->FindObject("gaus") ;
        ////std::pair<int, int> channelPhi = std::pair<int, int>(iL,Rn);
	////rentries[channelPhi] = entries;
	////MIPentries[channelPhi] = Mentries;
