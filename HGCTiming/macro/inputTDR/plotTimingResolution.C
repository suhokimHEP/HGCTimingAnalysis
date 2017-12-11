#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TCanvas.h>
#include "TRandom3.h"


void plotTimingResolution(){

  gROOT->Macro("/afs/cern.ch/user/a/amartell/public/setStyle.C");

  int iColors[3] = {kBlack, kRed+2, kBlue+2};
  //  int iColorsHG[3] = {kGray+3, kRed+2, kAzure-2};
  int iColorsHG[3] = {kGreen+2, kRed+2, kAzure+2};
  int iColorsHG2[3] = {kBlack, kRed+2, kBlack};

  //from TB
  float paramA[3];
  float parErrA[3];
  float paramC[3];
  float parErrC[3];
  //100um
  paramA[0] = 1.00;
  parErrA[0] = 0.01;
  paramC[0] = 0.02;
  parErrC[0] = 0.001;
  //200um  
  paramA[1] = 1.06;
  parErrA[1] = 0.02;
  paramC[1] = 0.02;
  parErrC[1] = 0.001;
  //300um 
  paramA[2] = 1.11;
  parErrA[2] = 0.02;
  paramC[2] = 0.02;
  parErrC[2] = 0.001;

  float noise[3];


  float paramA_MIP[3];
  float parErrA_MIP[3];
  float paramC_MIP[3];
  float parErrC_MIP[3];
  paramA_MIP[0] = 0.69;
  parErrA_MIP[0] = 0.01;
  paramC_MIP[0] = 0.020;
  parErrC_MIP[0] = 0.001;

  paramA_MIP[1] = 0.38;
  parErrA_MIP[1] = 0.01;
  paramC_MIP[1] = 0.02;
  parErrC_MIP[1] = 0.001;

  paramA_MIP[2] = 0.34;
  parErrA_MIP[2] = 0.01;
  paramC_MIP[2] = 0.020;
  parErrC_MIP[2] = 0.001;

  float parX[3];
  parX[0] = 5.;
  parX[1] = 5;
  parX[2] = 5.;

  float parQ[3];
  parQ[0] = 1.69;
  parQ[1] = 1.68;
  parQ[2] = 1.28;



  float fromTBtoHGC[3];
  fromTBtoHGC[0] = sqrt(2.5) * 3. * 0.5;
  fromTBtoHGC[1] = sqrt(2.5) * 60./14. * 0.5;
  fromTBtoHGC[2] = sqrt(2.5) * 4. * 0.5;


  float SoverNperMIP[3];
  SoverNperMIP[0] = 1./0.269;
  SoverNperMIP[1] = 1./0.131;
  SoverNperMIP[2] = 1./0.066;


  TF1* TB_MIP_100 = new TF1("TB_MIP_100", "sqrt(pow([0]/x, 2) + pow([1], 2) )", 1., 1000.);
  TF1* TB_MIP_200 = new TF1("TB_MIP_200", "sqrt(pow([0]/x, 2) + pow([1], 2) )", 1., 1000.);
  TF1* TB_MIP_300 = new TF1("TB_MIP_300", "sqrt(pow([0]/x, 2) + pow([1], 2) )", 1., 1000.);

  TF1* TB_SoN_100 = new TF1("TB_SoN_100", "sqrt(pow([0]/x, 2) + pow([1], 2) )", 1., 1000.);
  TF1* TB_SoN_200 = new TF1("TB_SoN_200", "sqrt(pow([0]/x, 2) + pow([1], 2) )", 1., 1000.);
  TF1* TB_SoN_300 = new TF1("TB_SoN_300", "sqrt(pow([0]/x, 2) + pow([1], 2) )", 1., 1000.);


  TF1* HGC_MIP_100 = new TF1("HGC_MIP_100", "sqrt(pow([0]/x, 2) + pow([1], 2) )", 1., 1000.);
  TF1* HGC_MIP_200 = new TF1("HGC_MIP_200", "sqrt(pow([0]/x, 2) + pow([1], 2) )", 1., 1000.);
  TF1* HGC_MIP_300 = new TF1("HGC_MIP_300", "sqrt(pow([0]/x, 2) + pow([1], 2) )", 1., 1000.);

  TF1* HGC_SoN_100 = new TF1("HGC_SoN_100", "sqrt(pow([0]/x, 2) + pow([1], 2) )", 1., 1000.);
  TF1* HGC_SoN_200 = new TF1("HGC_SoN_200", "sqrt(pow([0]/x, 2) + pow([1], 2) )", 1., 1000.);
  TF1* HGC_SoN_300 = new TF1("HGC_SoN_300", "sqrt(pow([0]/x, 2) + pow([1], 2) )", 1., 1000.);


  TF1* X_SoN_100 = new TF1("X_SoN_100", "sqrt(pow([0]/x, 2) + pow([1], 2) )", 1., 1000.);
  TF1* X_SoN_200 = new TF1("X_SoN_200", "sqrt(pow([0]/x, 2) + pow([1], 2) )", 1., 1000.);
  TF1* X_SoN_300 = new TF1("X_SoN_300", "sqrt(pow([0]/x, 2) + pow([1], 2) )", 1., 1000.);
  TF1* X_MIP_100 = new TF1("X_MIP_100", "sqrt(pow([0]/x, 2) + pow([1], 2) )", 1., 1000.);
  TF1* X_MIP_200 = new TF1("X_MIP_200", "sqrt(pow([0]/x, 2) + pow([1], 2) )", 1., 1000.);
  TF1* X_MIP_300 = new TF1("X_MIP_300", "sqrt(pow([0]/x, 2) + pow([1], 2) )", 1., 1000.);
  TF1* X_Q_100 = new TF1("X_Q_100", "sqrt(pow([0]/x, 2) + pow([1], 2) )", 1., 1000.);
  TF1* X_Q_200 = new TF1("X_Q_200", "sqrt(pow([0]/x, 2) + pow([1], 2) )", 1., 1000.);
  TF1* X_Q_300 = new TF1("X_Q_300", "sqrt(pow([0]/x, 2) + pow([1], 2) )", 1., 1000.);



  X_SoN_100->SetParameters(parX[0], paramC_MIP[0]);
  X_SoN_200->SetParameters(parX[1], paramC_MIP[1]);
  X_SoN_300->SetParameters(parX[2], paramC_MIP[2]);
  X_MIP_100->SetParameters(parX[0]/SoverNperMIP[0], paramC_MIP[0]);
  X_MIP_200->SetParameters(parX[1]/SoverNperMIP[1], paramC_MIP[1]);
  X_MIP_300->SetParameters(parX[2]/SoverNperMIP[2], paramC_MIP[2]);
  X_Q_100->SetParameters(parQ[0], paramC_MIP[0]);
  X_Q_200->SetParameters(parQ[1], paramC_MIP[1]);
  X_Q_300->SetParameters(parQ[2], paramC_MIP[2]);


  X_SoN_100->SetLineColor(iColorsHG2[0]);
  X_SoN_200->SetLineColor(iColorsHG2[1]);
  X_SoN_300->SetLineColor(iColorsHG2[2]);
  X_SoN_100->SetLineWidth(2);
  X_SoN_200->SetLineWidth(2);
  X_SoN_300->SetLineWidth(2);

  X_Q_100->SetLineColor(iColorsHG[0]);
  X_Q_200->SetLineColor(iColorsHG[1]);
  X_Q_300->SetLineColor(iColorsHG[2]);
  X_Q_100->SetLineWidth(2);
  X_Q_200->SetLineWidth(2);
  X_Q_300->SetLineWidth(2);

  X_MIP_100->SetLineColor(iColorsHG[0]);
  X_MIP_200->SetLineColor(iColorsHG[1]);
  X_MIP_300->SetLineColor(iColorsHG[2]);
  X_MIP_100->SetLineWidth(2);
  X_MIP_200->SetLineWidth(2);
  X_MIP_300->SetLineWidth(2);


  TB_MIP_100->SetParameters(paramA_MIP[0]/sqrt(2.), paramC_MIP[0]);
  TB_MIP_200->SetParameters(paramA_MIP[1]/sqrt(2.), paramC_MIP[1]);
  TB_MIP_300->SetParameters(paramA_MIP[2]/sqrt(2.), paramC_MIP[2]);
  TB_MIP_100->SetLineColor(iColors[0]);
  TB_MIP_200->SetLineColor(iColors[1]);
  TB_MIP_300->SetLineColor(iColors[2]);


  HGC_MIP_100->SetParameters(paramA_MIP[0]/sqrt(2.)*fromTBtoHGC[0], paramC_MIP[0]);
  HGC_MIP_200->SetParameters(paramA_MIP[1]/sqrt(2.)*fromTBtoHGC[1], paramC_MIP[1]);
  HGC_MIP_300->SetParameters(paramA_MIP[2]/sqrt(2.)*fromTBtoHGC[2], paramC_MIP[2]);
  HGC_MIP_100->SetLineColor(iColorsHG[0]);
  HGC_MIP_200->SetLineColor(iColorsHG[1]);
  HGC_MIP_300->SetLineColor(iColorsHG[2]);


  TB_SoN_100->SetParameters(paramA[0]/sqrt(2.), paramC[0]);
  TB_SoN_200->SetParameters(paramA[1]/sqrt(2.), paramC[1]);
  TB_SoN_300->SetParameters(paramA[2]/sqrt(2.), paramC[2]);
  TB_SoN_100->SetLineColor(iColors[0]);
  TB_SoN_200->SetLineColor(iColors[1]);
  TB_SoN_300->SetLineColor(iColors[2]);


  HGC_SoN_100->SetParameters(paramA_MIP[0]/sqrt(2.)*SoverNperMIP[0]*fromTBtoHGC[0], paramC[0]);
  HGC_SoN_200->SetParameters(paramA_MIP[1]/sqrt(2.)*SoverNperMIP[1]*fromTBtoHGC[1], paramC[1]);
  HGC_SoN_300->SetParameters(paramA_MIP[2]/sqrt(2.)*SoverNperMIP[2]*fromTBtoHGC[2], paramC[2]);
  HGC_SoN_100->SetLineColor(iColorsHG[0]);
  HGC_SoN_200->SetLineColor(iColorsHG[1]);
  HGC_SoN_300->SetLineColor(iColorsHG[2]);


  gStyle->SetOptTitle(0);


  TLatex tBX;
  tBX.SetNDC();
  tBX.SetTextSize(0.05);
  tBX.SetTextFont(132);

  TLatex tB100;
  tB100.SetNDC();
  tB100.SetTextSize(0.05);
  tB100.SetTextFont(132);
  tB100.SetTextColor(iColors[0]);

  TLatex tB200;
  tB200.SetNDC();
  tB200.SetTextSize(0.05);
  tB200.SetTextFont(132);
  tB200.SetTextColor(iColors[1]);

  TLatex tB300;
  tB300.SetNDC();
  tB300.SetTextSize(0.05);
  tB300.SetTextFont(132);
  tB300.SetTextColor(iColors[2]);

  TLatex tH100;
  tH100.SetNDC();
  tH100.SetTextSize(0.05);
  tH100.SetTextFont(132);
  tH100.SetTextColor(iColorsHG[0]);

  TLatex tH200;
  tH200.SetNDC();
  tH200.SetTextSize(0.05);
  tH200.SetTextFont(132);
  tH200.SetTextColor(iColorsHG[1]);

  TLatex tH300;
  tH300.SetNDC();
  tH300.SetTextSize(0.05);
  tH300.SetTextFont(132);
  tH300.SetTextColor(iColorsHG[2]);


  TCanvas* cMIP = new TCanvas(); 
  gPad->SetLogx();
  gPad->SetLogy();
  gPad->SetGrid();
  X_MIP_100->SetRange(1, 1000);
  X_MIP_100->GetXaxis()->SetTitle("signal (MIP)");
  X_MIP_100->GetYaxis()->SetRangeUser(0.01, 10.);
  X_MIP_100->GetYaxis()->SetTitle("#sigma_{t} (ns)");
  cMIP->cd();
  // TB_MIP_100->Draw();
  // TB_MIP_200->Draw("same");
  // TB_MIP_300->Draw("same");
  // HGC_MIP_100->Draw("same");
  // HGC_MIP_200->Draw("same");
  // HGC_MIP_300->Draw("same");

  X_MIP_100->Draw("");
  X_MIP_200->Draw("same");
  X_MIP_300->Draw("same");

  tBX.DrawLatex(0.50,0.85, "#sigma(t) (ns) = A/x #oplus C");
  tH100.DrawLatex(0.4,0.75, Form("Si 100#mum  A = %.2f  C = %.2f", X_MIP_100->GetParameter(0), HGC_MIP_100->GetParameter(1)));
  tH200.DrawLatex(0.4,0.70, Form("Si 200#mum  A = %.2f  C = %.2f", X_MIP_200->GetParameter(0), HGC_MIP_200->GetParameter(1)));
  tH300.DrawLatex(0.4,0.65, Form("Si 300#mum  A = %.2f  C = %.2f", X_MIP_300->GetParameter(0), HGC_MIP_300->GetParameter(1)));

  // tB100.DrawLatex(0.5,0.60, Form("TB 100um A = %.2f C = %.3f", TB_MIP_100->GetParameter(0), TB_MIP_100->GetParameter(1)));
  // tB200.DrawLatex(0.5,0.55, Form("TB 200um A = %.2f C = %.3f", TB_MIP_200->GetParameter(0), TB_MIP_200->GetParameter(1)));
  // tB300.DrawLatex(0.5,0.50, Form("TB 300um A = %.2f C = %.3f", TB_MIP_300->GetParameter(0), TB_MIP_300->GetParameter(1)));
  cMIP->Print("timeResolution_singleSensor_MIP.png", "png");
  cMIP->Print("timeResolution_singleSensor_MIP.pdf", "pdf");
  cMIP->Print("timeResolution_singleSensor_MIP.root", "root");

  std::cout << " 100 X A = " << X_MIP_100->GetParameter(0) << std::endl;
  std::cout << " 200 X A = " << X_MIP_200->GetParameter(0) << std::endl;
  std::cout << " 300 X A = " << X_MIP_300->GetParameter(0) << std::endl;
  ////////

  TCanvas* cQfC = new TCanvas(); 
  gPad->SetLogx();
  gPad->SetLogy();
  gPad->SetGrid();
  X_Q_100->SetRange(1, 1000);
  X_Q_100->GetXaxis()->SetTitle("Q (fC)   ");
  X_Q_100->GetYaxis()->SetRangeUser(0.01, 10.);
  X_Q_100->GetYaxis()->SetTitle("#sigma_{t} (ns)");
  cQfC->cd();
  // TB_MIP_100->Draw();
  // TB_MIP_200->Draw("same");
  // TB_MIP_300->Draw("same");
  // HGC_MIP_100->Draw("same");
  // HGC_MIP_200->Draw("same");
  // HGC_MIP_300->Draw("same");

  X_Q_100->Draw("");
  X_Q_200->Draw("same");
  X_Q_300->Draw("same");

  tBX.DrawLatex(0.50,0.85, "#sigma(t) (ns) = #frac{A}{S/N} #oplus C");
  // tH100.DrawLatex(0.4,0.75, Form("Si 100#mum  A = %.2f  C = %.2f", X_Q_100->GetParameter(0), HGC_MIP_100->GetParameter(1)));
  // tH200.DrawLatex(0.4,0.70, Form("Si 200#mum  A = %.2f  C = %.2f", X_Q_200->GetParameter(0), HGC_MIP_200->GetParameter(1)));
  // tH300.DrawLatex(0.4,0.65, Form("Si 300#mum  A = %.2f  C = %.2f", X_Q_300->GetParameter(0), HGC_MIP_300->GetParameter(1)));

  //  tH100.DrawLatex(0.4,0.75, "Si 100#mum");
  tH200.DrawLatex(0.5,0.70, "Si 100#mum and 200#mum");
  tH300.DrawLatex(0.5,0.65, "Si 300#mum");

  // tB100.DrawLatex(0.5,0.60, Form("TB 100um A = %.2f C = %.3f", TB_MIP_100->GetParameter(0), TB_MIP_100->GetParameter(1)));
  // tB200.DrawLatex(0.5,0.55, Form("TB 200um A = %.2f C = %.3f", TB_MIP_200->GetParameter(0), TB_MIP_200->GetParameter(1)));
  // tB300.DrawLatex(0.5,0.50, Form("TB 300um A = %.2f C = %.3f", TB_MIP_300->GetParameter(0), TB_MIP_300->GetParameter(1)));
  cQfC->Print("timeResolution_singleSensor_QfC.png", "png");
  cQfC->Print("timeResolution_singleSensor_QfC.pdf", "pdf");
  cQfC->Print("timeResolution_singleSensor_QfC.root", "root");

  std::cout << " 100 X A = " << X_MIP_100->GetParameter(0) << std::endl;
  std::cout << " 200 X A = " << X_MIP_200->GetParameter(0) << std::endl;
  std::cout << " 300 X A = " << X_MIP_300->GetParameter(0) << std::endl;


  TCanvas* cSoN = new TCanvas();
  gPad->SetLogx();
  gPad->SetLogy();
  gPad->SetGrid();
  cSoN->cd();
  X_SoN_100->SetRange(1, 1000);
  X_SoN_100->GetXaxis()->SetTitle("S/N   ");
  X_SoN_100->GetYaxis()->SetRangeUser(0.01, 10.);
  X_SoN_100->GetYaxis()->SetTitle("#sigma_{t} (ns)");
  cSoN->cd();
  // TB_SoN_100->Draw();
  // TB_SoN_200->Draw("same");
  // TB_SoN_300->Draw("same");
  // HGC_SoN_100->Draw("");
  // HGC_SoN_200->Draw("same");
  // HGC_SoN_300->Draw("same");
  X_SoN_100->Draw("");
  //  X_SoN_200->Draw("same");
  //  X_SoN_300->Draw("same");

  tBX.DrawLatex(0.50,0.85, "#sigma_{t} (ns) = #frac{A}{S/N} #oplus C");
  tB100.DrawLatex(0.5,0.75, Form("A = %.2f ns,  C = %.2f ns", X_SoN_100->GetParameter(0), X_SoN_100->GetParameter(1)));
  //tB200.DrawLatex(0.5,0.65, Form("A = %.2f C = %.2f", X_SoN_200->GetParameter(0), X_SoN_200->GetParameter(1)));
  // tH300.DrawLatex(0.4,0.70, Form("HGC 300um A = %.2f C = %.3f", X_SoN_300->GetParameter(0), X_SoN_300->GetParameter(1)));

  // tB100.DrawLatex(0.5,0.60, Form("TB 100um A = %.2f C = %.3f", TB_SoN_100->GetParameter(0), TB_SoN_100->GetParameter(1)));
  // tB200.DrawLatex(0.5,0.55, Form("TB 200um A = %.2f C = %.3f", TB_SoN_200->GetParameter(0), TB_SoN_200->GetParameter(1)));
  // tB300.DrawLatex(0.5,0.50, Form("TB 300um A = %.2f C = %.3f", TB_SoN_300->GetParameter(0), TB_SoN_300->GetParameter(1)));
  cSoN->Print("timeResolution_singleSensor_SoN.png", ".png");
  cSoN->Print("timeResolution_singleSensor_SoN.pdf", ".pdf");
  cSoN->Print("timeResolution_singleSensor_SoN.root", ".root");


  /*
  TLegend *legTGM = new TLegend(0.65,0.5,0.80,0.85,NULL,"brNDC");
  legTGM->SetTextFont(42);
  legTGM->SetTextSize(0.05);
  legTGM->SetFillColor(kWhite);
  legTGM->SetLineColor(kWhite);
  legTGM->SetShadowColor(kWhite);
  legTGM->AddEntry(timeResolution_300, "#sigma(t) from TB - cell 1cm^{2}", "l");
  legTGM->AddEntry(timeMIP_300, "j = 1.2ns/fC - 300um", "l");
  legTGM->AddEntry(timeResolution05_300, "#sigma(t) - cell 0.5cm^{2}", "l");
  legTGM->AddEntry(timeMIP05_300, "j = 0.6ns/fC - 300um", "l");
  TCanvas* cTB_100 = new TCanvas();
  gPad->SetLogx();
  gPad->SetLogy();
  cTB_100->cd();
  timeResolution_300->GetYaxis()->SetRangeUser(0.005, 1.);
  timeResolution_300->GetYaxis()->SetTitle("#sigma(t) (ns)");
  timeResolution_300->GetXaxis()->SetTitle("S/N");
  timeResolution_300->Draw();
  timeMIP_300->Draw("same");
  timeResolution05_300->Draw("same");
  timeMIP05_300->Draw("same");
  legTGM->Draw("same");
  cTB_100->Print("timeResolution_TB_300.png", "png");
  */





}
