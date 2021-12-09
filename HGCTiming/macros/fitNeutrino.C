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

void fitNeutrino(){
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
 
  std::string inputFileList = "outMIP_Neutfit_eol.root";
  std::string outtag = inputFileList.substr(0,inputFileList.find(".root"));
  printf("%s",outtag.c_str());
 TFile* inF = TFile::Open(TString(inputFileList));
 TTree* newT = (TTree*)inF->Get("MIPtree");
 std::map<std::pair<int,int>, int> rentries;
 std::map<std::pair<int,int>, int> MIPentries;
rentries.clear();
	TCanvas*canvas1 = new TCanvas();
 FILE * outfulltable;
      outfulltable = fopen("try.tex","w");

  for(int iL=0; iL <= 14; ++iL){
  for(int Rn=16; Rn <= 38; ++Rn){
        int entries = newT->GetEntries(Form("MIP_iR==%d&&MIP_layer==%d",Rn,iL));
        //int Mentries = newT->GetEntries(Form("MIP_iR==%d&&MIP_layer==%d&&MIP_val<0.543",Rn,iL));
	newT->Draw("MIP_val>>temp(150,-3,3)",Form("MIP_iR==%d&&MIP_layer==%d",Rn,iL));
        TH1F* fithist = (TH1F*)newT->GetHistogram();
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

      canvas1->Print(Form("singleFits/fitNeutrino_mip_L%d_R%d_nEvts%d.png",iL,Rn,entries), ".png"); 
	}
        fprintf (outfulltable, "\n");		
  }
}








	//fithist->Fit("gaus");
	//TF1 *g = (TF1*)fithist->GetListOfFunctions()->FindObject("gaus") ;
        ////std::pair<int, int> channelPhi = std::pair<int, int>(iL,Rn);
	////rentries[channelPhi] = entries;
	////MIPentries[channelPhi] = Mentries;
