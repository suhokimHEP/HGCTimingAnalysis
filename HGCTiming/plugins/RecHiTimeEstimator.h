#ifndef RecHiTimeEstimator_h
#define RecHiTimeEstimator_h

// user include files
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TCanvas.h>

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
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include "TRandom3.h"
#include "TProfile.h"
using namespace std;




//legend
/*
cellType 0 = cell standard (1cm2)
cellType 1 = cell half (0.5cm2) => effect ~half timing resolution
floor    0.02 = Floor 20ps time at high S/N
floor    0.03 = Floor 30ps time at high S/N
lifeAge   0 = life at beginning for charge collection efficiency
lifeAge   1 = life end for charge collection efficiency (70% 50% 50% for 300, 200, 100)
*/


class RecHiTimeEstimator
{
  
public:
  explicit RecHiTimeEstimator(const edm::ParameterSet& ps);
  ~RecHiTimeEstimator();

  double getTimeHit(int thick, double SoverN);
  double getTimeHitFixThr();
  double getExpectedReso(const DetId detid, double energy);
  void correctTime(const HGCRecHitCollection& rechits, HGCRecHitCollection* Newrechits);
  void correctTimeFixThr(const HGCRecHitCollection& rechits, HGCRecHitCollection* Newrechits);

  void setOptions(int cellType = 0, float floor = 0.02, int lifeAge = 0, float absTrend=1.);

  /* void timeCSF20LB(const HGCRecHitCollection& rechits, HGCRecHitCollection* Newrechits); */
  /* void timeCSF30LB(const HGCRecHitCollection& rechits, HGCRecHitCollection* Newrechits); */
  /* void timeCHF20LB(const HGCRecHitCollection& rechits, HGCRecHitCollection* Newrechits); */
  /* void timeCHF30LB(const HGCRecHitCollection& rechits, HGCRecHitCollection* Newrechits); */
  /* void timeCSF20LE(const HGCRecHitCollection& rechits, HGCRecHitCollection* Newrechits); */
  /* void timeCSF30LE(const HGCRecHitCollection& rechits, HGCRecHitCollection* Newrechits); */
  /* void timeCHF20LE(const HGCRecHitCollection& rechits, HGCRecHitCollection* Newrechits); */
  /* void timeCHF30LE(const HGCRecHitCollection& rechits, HGCRecHitCollection* Newrechits); */

  void setEventSetup(const edm::EventSetup& es);

private:
  hgcal::RecHitTools recHitTools;

  TF1* timeResolution;

  edm::Service<TFileService> fs;

  TH1F* SoN_for100;
  TH1F* SoN_for200;
  TH1F* SoN_for300;

  TH1F* originalTime_100;
  TH1F* originalTime_200;
  TH1F* originalTime_300;
  TH2F* originalTime_vsMip_100;
  TH2F* originalTime_vsMip_200;
  TH2F* originalTime_vsMip_300;

  float fromTBtoHGC[3];
  float SoverNperMIP[3];

  float paramA[3];
  float parErrA[3];
  float paramC[3];
  float parErrC[3];

  float floorValue;
  float absoluteTrend;

  float chargeCollEff[3];

  float cellSize[3];

  std::vector<float> scaleCorrection;
  std::vector<float> weights;

  float keV2GeV;
  float keV2MeV;

  double keV2fC[2];
  double keV2MIP;

  double noiseEndOfLife[3];
  double noiseBegOfLife[3];

  //for cell type
  double noisefC[3];
  double noisefC_bkp[3];
  double noiseMIP;

  double fCPerMIP[3];
};

#endif
