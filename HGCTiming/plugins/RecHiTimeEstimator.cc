// user include files
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TCanvas.h>
#include "TRandom3.h"
#include <cstdlib> 
#include <ctime>

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "Geometry/HGCalCommonData/interface/HGCalDDDConstants.h"
#include "Geometry/HcalCommonData/interface/HcalDDDRecConstants.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include "HGCTimingAnalysis/HGCTiming/plugins/RecHiTimeEstimator.h"

RecHiTimeEstimator::RecHiTimeEstimator(const edm::ParameterSet& ps){
  //  std::cout << " >>> constructor " << std::endl;
  //section type
  keV2fC[0] =  ps.getParameter<double>("HGCEE_keV2fC");
  keV2fC[1] =  ps.getParameter<double>("HGCHEF_keV2fC");
  keV2MIP = ps.getParameter<double>("HGCHB_keV2MIP");

  //cell type
  noiseEndOfLife[0] = (ps.getParameter<std::vector<double> >("endOfLifeNoises")).at(0);
  noiseEndOfLife[1] = (ps.getParameter<std::vector<double> >("endOfLifeNoises")).at(1);
  noiseEndOfLife[2] = (ps.getParameter<std::vector<double> >("endOfLifeNoises")).at(2);

  noiseBegOfLife[0] = (ps.getParameter<std::vector<double> >("nonAgedNoises")).at(0);
  noiseBegOfLife[1] = (ps.getParameter<std::vector<double> >("nonAgedNoises")).at(1);
  noiseBegOfLife[2] = (ps.getParameter<std::vector<double> >("nonAgedNoises")).at(2);

  noisefC_bkp[0] = (ps.getParameter<std::vector<double> >("HGCEE_noisefC")).at(0);
  noisefC_bkp[1] = (ps.getParameter<std::vector<double> >("HGCEE_noisefC")).at(1);
  noisefC_bkp[2] = (ps.getParameter<std::vector<double> >("HGCEE_noisefC")).at(2);

  noisefC[0] = (ps.getParameter<std::vector<double> >("HGCEE_noisefC")).at(0);
  noisefC[1] = (ps.getParameter<std::vector<double> >("HGCEE_noisefC")).at(1);
  noisefC[2] = (ps.getParameter<std::vector<double> >("HGCEE_noisefC")).at(2);

  noiseMIP = ps.getParameter<double>("HGCBH_noiseMIP");

  fCPerMIP[0] =  (ps.getParameter<std::vector<double> >("HGCEE_fCPerMIP")).at(0);
  fCPerMIP[1] =  (ps.getParameter<std::vector<double> >("HGCEE_fCPerMIP")).at(1);
  fCPerMIP[2] =  (ps.getParameter<std::vector<double> >("HGCEE_fCPerMIP")).at(2);

  const auto& rcorr = ps.getParameter<std::vector<double> >("thicknessCorrection");
  scaleCorrection.clear();
  for( auto corr : rcorr ) {
    scaleCorrection.push_back(1.0/corr);
  }

  const auto& dweights = ps.getParameter<std::vector<double> >("dEdXweights");
  for( auto weight : dweights ) {
    weights.push_back(weight);
  }

  keV2GeV = 1e-6;
  keV2MeV = 1e-3;

  //100um
  paramA[0] = 1.;
  parErrA[0] = 0.01;
  paramC[0] = 0.020;
  parErrC[0] = 0.001;

  //200um
  paramA[1] = 1.;
  parErrA[1] = 0.01;
  paramC[1] = 0.02;
  parErrC[1] = 0.001;

  //300um
  paramA[2] = 1.;
  parErrA[2] = 0.01;
  paramC[2] = 0.020;
  parErrC[2] = 0.001;

  SoverNperMIP[0] = 1./(noisefC[0]/fCPerMIP[0]);
  SoverNperMIP[1] = 1./(noisefC[1]/fCPerMIP[1]);
  SoverNperMIP[2] = 1./(noisefC[2]/fCPerMIP[2]);

  fromTBtoHGC[0] = sqrt(2.5) * 3. * 0.5;
  fromTBtoHGC[1] = sqrt(2.5) * 60./14. * 0.5;
  fromTBtoHGC[2] = sqrt(2.5) * 4. * 0.5;

  floorValue = 0.02;
  absoluteTrend = 1.;

  chargeCollEff[0] = 1.;
  chargeCollEff[1] = 1.;
  chargeCollEff[2] = 1.;

  cellSize[0] = 1.;
  cellSize[1] = 1.;
  cellSize[2] = 1.;

  timeResolution = new TF1("timeSi100", "sqrt(pow([0]/x, 2) + pow([1], 2) )", 1., 1000.);

  SoN_for100 = fs->make<TH1F>("SoN_for100", "", 1000, 0., 1000.);
  SoN_for200 = fs->make<TH1F>("SoN_for200", "", 1000, 0., 1000.);
  SoN_for300 = fs->make<TH1F>("SoN_for300", "", 1000, 0., 1000.);
}


void RecHiTimeEstimator::setEventSetup(const edm::EventSetup& es){
  recHitTools.getEventSetup(es);
}


void RecHiTimeEstimator::setOptions(int cellType, float floor, int lifeAge, float absTrend){
  if(cellType == 1){
    cellSize[0] = 1.;
    cellSize[1] = 0.5;
    cellSize[2] = 0.5;
  }
  if(floor == 0.03) {
    floorValue = floor;
  }
  if(lifeAge == 1){
    chargeCollEff[0] = 0.5;
    chargeCollEff[1] = 0.5;
    chargeCollEff[2] = 0.7;

    noisefC[0] = noisefC_bkp[0] * noiseEndOfLife[0] / noiseBegOfLife[0];
    noisefC[1] = noisefC_bkp[1] * noiseEndOfLife[1] / noiseBegOfLife[1];
    noisefC[2] = noisefC_bkp[2] * noiseEndOfLife[2] / noiseBegOfLife[2];

    SoverNperMIP[0] = 1./(noisefC[0]/fCPerMIP[0]);
    SoverNperMIP[1] = 1./(noisefC[1]/fCPerMIP[1]);
    SoverNperMIP[2] = 1./(noisefC[2]/fCPerMIP[2]);
  }
  if(absTrend == 1.5){
    absoluteTrend = absTrend;
  }
  if(absTrend == 2.){
    absoluteTrend = absTrend;
  }

  paramC[0] = floor;
  paramC[1] = floor;
  paramC[2] = floor;

  paramA[0] = absoluteTrend;
  paramA[1] = absoluteTrend;
  paramA[2] = absoluteTrend;

}


double RecHiTimeEstimator::getExpectedReso(const DetId detid, double energyHit){
  int thick = (detid.det() != DetId::Forward) ? -1 : recHitTools.getSiThickness(detid) / 100. - 1.;
  timeResolution->SetParameters(paramA[thick]*cellSize[thick], paramC[thick]);

  int sectionType = -1;
  if(detid.subdetId() == HGCEE) sectionType = 0;
  else if(detid.subdetId() == HGCHEF) sectionType = 1;
  else if(detid.subdetId() == HGCHEB) sectionType = 2;

  float energy = energyHit*chargeCollEff[thick];
  unsigned int layer = recHitTools.getLayerWithOffset(detid);

  double sigmaNoiseMIP = noiseMIP;
  if(sectionType != 2) sigmaNoiseMIP = noisefC[thick]/fCPerMIP[thick];

  int energyMIP = 0.;
  if(sectionType == 2) energyMIP = energy/keV2GeV * keV2MIP;
  else if(sectionType == 0 || sectionType == 1) energyMIP = energy/scaleCorrection.at(thick)/keV2GeV / (weights.at(layer)/keV2MeV);

  float SoverN = energyMIP / sigmaNoiseMIP;

  timeResolution->SetParameters(paramA[thick]*cellSize[thick], paramC[thick]);
  double sigma = 0.2;
  if(SoverN > 1) sigma = timeResolution->Eval(SoverN);
  if(sigma < floorValue) sigma = floorValue;

  return sigma;
}


double RecHiTimeEstimator::getTimeHit(int thick, double SoverN){
  timeResolution->SetParameters(paramA[thick]*cellSize[thick], paramC[thick]);

  //resolution from TB results with floor of 20ps at high S/N
  double sigma = 0.2;
  if(SoverN > 1) sigma = timeResolution->Eval(SoverN);
  if(sigma < floorValue || SoverN > 1000.) sigma = floorValue;

  TRandom3* rand = new TRandom3();
  rand->SetSeed(0);
  double smearing = rand->Gaus(0., sigma);
  return smearing;
}


double RecHiTimeEstimator::getTimeHitFixThr(){
  //flat resolution at 50ps
  double sigma = 0.05;

  TRandom3* rand = new TRandom3();
  rand->SetSeed(0);
  double smearing = rand->Gaus(0., sigma);
  return smearing;
}


void RecHiTimeEstimator::correctTime(const HGCRecHitCollection& rechits, HGCRecHitCollection* Newrechits){

  for(HGCRecHitCollection::const_iterator it_hit = rechits.begin(); it_hit < rechits.end(); ++it_hit) {
    const DetId detid = it_hit->detid();

    int thick = (detid.det() != DetId::Forward) ? -1 : recHitTools.getSiThickness(detid) / 100. - 1.;

    int sectionType = -1;
    if(detid.subdetId() == HGCEE) sectionType = 0;
    else if(detid.subdetId() == HGCHEF) sectionType = 1;
    else if(detid.subdetId() == HGCHEB) sectionType = 2;

    HGCRecHit myrechit(*it_hit);
    float energy = it_hit->energy();
    float time = it_hit->time();

    if(sectionType == -1 || thick == -1){
      myrechit.setTime(-1.);
      Newrechits->push_back(myrechit);
      continue;
    }

    energy = it_hit->energy()*chargeCollEff[thick];

    unsigned int layer = recHitTools.getLayerWithOffset(detid);

    double sigmaNoiseMIP = noiseMIP;
    if(sectionType != 2) sigmaNoiseMIP = noisefC[thick]/fCPerMIP[thick];

    int energyMIP = 0.;
    if(sectionType == 2) energyMIP = energy/keV2GeV * keV2MIP;
    else if(sectionType == 0 || sectionType == 1) energyMIP = energy/scaleCorrection.at(thick)/keV2GeV / (weights.at(layer)/keV2MeV);

    float SoverN = energyMIP / sigmaNoiseMIP;

    if(thick == 0) SoN_for100->Fill(SoverN);
    if(thick == 1) SoN_for200->Fill(SoverN);
    if(thick == 2) SoN_for300->Fill(SoverN);

    if(energyMIP > 3. && SoverN > 10.){
      double smearedTime = getTimeHit(thick, SoverN);
      myrechit.setTime(time + smearedTime);
    }
    else myrechit.setTime(-1.); 

    Newrechits->push_back(myrechit);
  }
  return;
}


void RecHiTimeEstimator::correctTimeFixThr(const HGCRecHitCollection& rechits, HGCRecHitCollection* Newrechits){

  for(HGCRecHitCollection::const_iterator it_hit = rechits.begin(); it_hit < rechits.end(); ++it_hit) {
    const DetId detid = it_hit->detid();
    int thick = (detid.det() != DetId::Forward) ? -1 : recHitTools.getSiThickness(detid) / 100. - 1.;

    int sectionType = -1;
    if(detid.subdetId() == HGCEE) sectionType = 0;
    else if(detid.subdetId() == HGCHEF) sectionType = 1;
    else if(detid.subdetId() == HGCHEB) sectionType = 2;

    HGCRecHit myrechit(*it_hit);
    float energy = it_hit->energy();
    float time = it_hit->time();

    if(sectionType == -1 || thick == -1){
      myrechit.setTime(-1.);
      Newrechits->push_back(myrechit);

      continue;
    }

    energy = it_hit->energy()*chargeCollEff[thick];

    unsigned int layer = recHitTools.getLayerWithOffset(detid);

    int energyMIP = 0.;
    if(sectionType == 2) energyMIP = energy/keV2GeV * keV2MIP;
    else if(sectionType == 1 || sectionType == 0) energyMIP = energy/scaleCorrection.at(thick)/keV2GeV / (weights.at(layer)/keV2MeV);

    float energyCharge = 0.;
    // from SimCalorimetry/HGCalSimProducers/src/HGCHEbackDigitizer.cc L 58
    if(sectionType == 2) energyCharge = energyMIP * 1.; 
    else if(sectionType == 1 || sectionType == 0) energyCharge = energyMIP * fCPerMIP[thick];


    if(energyCharge > 60){
      double smearedTime = getTimeHitFixThr();
      myrechit.setTime(time + smearedTime);
    }
    else myrechit.setTime(-1.); 

    Newrechits->push_back(myrechit);
  }
  return;
}


