// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"



#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCHEDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"

#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"

#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/TruncatedPyramid.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"


#include "TRandom3.h"
#include "TTree.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TH2F.h"
#include "TH3F.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"


#include "RecoLocalCalo/HGCalRecAlgos/interface/ClusterTools.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/HGCalDepthPreClusterer.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include "HGCTimingAnalysis/HGCTiming/plugins/RecHiTimeEstimator.h"
#include "HGCTimingAnalysis/HGCTiming/interface/UtilClasses.h"



#include <string>
#include <map>

class HGCalTimingAnalyzerPU : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
//
// constructors and destructor
//
  explicit HGCalTimingAnalyzerPU(const edm::ParameterSet&);
  ~HGCalTimingAnalyzerPU();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  static float getXmax(TH1F* histo, float& YMax);

  static bool comparePairs(const std::pair<float, float>& i, const std::pair<float, float>& j);

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;


 // ----------member data ---------------------------

  edm::EDGetTokenT<HGCRecHitCollection> _recHitsEE;
  edm::EDGetTokenT<HGCRecHitCollection> _recHitsFH;
  edm::EDGetTokenT<HGCRecHitCollection> _recHitsBH;
  edm::EDGetTokenT<std::vector<TrackingVertex> > _vtx;
  edm::EDGetTokenT<std::vector<TrackingParticle> > _part;
  edm::EDGetTokenT<std::vector<CaloParticle> > _caloParticles;
  //  edm::EDGetTokenT<std::vector<reco::HGCalMultiCluster> > _multiClusters;

  std::string                detector;
  int                        algo;
  HGCalDepthPreClusterer     pre;
  bool                       rawRecHits;
  int                        CFDTimeCorrection;
  hgcal::RecHitTools         recHitTools;

  std::unique_ptr<hgcal::ClusterTools> clusterTools;

  RecHiTimeEstimator* myEstimator;


  std::vector<float> scaleCorrection;
  std::vector<float> weights;

  float keV2GeV;
  float keV2MeV;

  double keV2fC[2];
  double keV2MIP;

  double noisefC[2];
  double noiseMIP;
  //for cell type
  double fCPerMIP[3];

  TH1F* h_Vtx_x;
  TH1F* h_Vtx_y;
  TH1F* h_Vtx_z;

  TH1F* h_Vtx_dvx;
  TH1F* h_Vtx_dvy;
  TH1F* h_Vtx_dvz;

  TProfile2D* cellThick_Rvseta;
  TProfile2D* cellThick_RvsLayer;

  //cancella da qui
  TH2F* h2_YvsX_energy;
  TH2F* h_layer_energy;
  TH2F* h2_YvsX_charge;
  TH2F* h_layer_charge;
  TH2F* h2_YvsX_MIP;
  TH2F* h_layer_MIP;

  TH3F* h3_YvsXvsZ_mip;
  TH2F* h2_MIP_vsRadius;
  TH2F* h2_energy_vsRadius;


  TH1F* h_allRH_TimesOfRadius;

  TH2F* h2_dPhivsdEta_rhGen;
  TH2F* h2_dPhivsdEta_rhAxis;
  TH2F* h2_dPhivsdEta_GenAxis;

  //eta bins 1.65 - 1.85   1.85-2.05   2.05-2.25   2.25-2.45   2.45-2.65   2.65-2.85
  //eta bins 1.5-1.7  1.7-1.9   1.9-2.1 - 2.1-2.3   2.3-2.5 2.5-2.7   2.7-2.9   2.9-3.1
  //  eta bins 1.5-2.  2.-2.5  2.5-3.



  TH1F* hFractionEvents_PU_Eta_dRadius[6][4];
  TH1F* hFractionEvents_PUany_Eta_dRadius[6][4];
  TH1F* hFractionEvents_MainEvt_Eta_dRadius[6][4];
  TH1F* hFractionEvents_MainEvtany_Eta_dRadius[6][4];
  TH1F* hFractionEvents_Shared_Eta_dRadius[6][4];
  TH1F* hFractionEvents_SharedME_Eta_dRadius[6][4];
  TH1F* hFractionEvents_SharedPU_Eta_dRadius[6][4];

  TH1F* hFractionEvents_PU_Eta_dRadius_3MIP[6][4];
  TH1F* hFractionEvents_PUany_Eta_dRadius_3MIP[6][4];
  TH1F* hFractionEvents_MainEvt_Eta_dRadius_3MIP[6][4];
  TH1F* hFractionEvents_MainEvtany_Eta_dRadius_3MIP[6][4];
  TH1F* hFractionEvents_Shared_Eta_dRadius_3MIP[6][4];
  TH1F* hFractionEvents_SharedME_Eta_dRadius_3MIP[6][4];
  TH1F* hFractionEvents_SharedPU_Eta_dRadius_3MIP[6][4];

 
  TH1F* h_AllHits_noEcut;
  TH1F* h_AllHits_3MIPcut;

  TH1F* h_Energy;
  TH1F* h_Energy_PU;
  TH1F* h_Energy_MainEvt;
  TH1F* h_Energy_Shared;
  TH1F* h_Energy_SharedME;
  TH1F* h_Energy_SharedPU;

  TH1F* h_Fraction_Energy_PU;
  TH1F* h_Fraction_Energy_MainEvt;
  TH1F* h_Fraction_Energy_Shared;
  TH1F* h_Fraction_Energy_SharedME;
  TH1F* h_Fraction_Energy_SharedPU;

  TH1F* h_Energy_3MIP;
  TH1F* h_Energy_PU_3MIP;
  TH1F* h_Energy_MainEvt_3MIP;
  TH1F* h_Energy_Shared_3MIP;
  TH1F* h_Energy_SharedME_3MIP;
  TH1F* h_Energy_SharedPU_3MIP;

  TH1F* h_Fraction_Energy_PU_3MIP;
  TH1F* h_Fraction_Energy_MainEvt_3MIP;
  TH1F* h_Fraction_Energy_Shared_3MIP;
  TH1F* h_Fraction_Energy_SharedME_3MIP;
  TH1F* h_Fraction_Energy_SharedPU_3MIP;
 

  TH1F* h_NumberHits_Eta_dRadius[6][4];
  TH1F* h_NumberHits_PU_Eta_dRadius[6][4];
  TH1F* h_NumberHits_PUany_Eta_dRadius[6][4];
  TH1F* h_NumberHits_MainEvt_Eta_dRadius[6][4];
  TH1F* h_NumberHits_MainEvtany_Eta_dRadius[6][4];
  TH1F* h_NumberHits_Shared_Eta_dRadius[6][4];
  TH1F* h_NumberHits_SharedME_Eta_dRadius[6][4];
  TH1F* h_NumberHits_SharedPU_Eta_dRadius[6][4];
  TH1F* h_FractionHits_PU_Eta_dRadius[6][4];
  TH1F* h_FractionHits_PUany_Eta_dRadius[6][4];
  TH1F* h_FractionHits_MainEvt_Eta_dRadius[6][4];
  TH1F* h_FractionHits_MainEvtany_Eta_dRadius[6][4];
  TH1F* h_FractionHits_Shared_Eta_dRadius[6][4];
  TH1F* h_FractionHits_SharedME_Eta_dRadius[6][4];
  TH1F* h_FractionHits_SharedPU_Eta_dRadius[6][4];

  TH1F* h_Energy_Eta_dRadius[6][4];
  TH1F* h_Energy_PU_Eta_dRadius[6][4];
  TH1F* h_Energy_PUany_Eta_dRadius[6][4];
  TH1F* h_Energy_MainEvt_Eta_dRadius[6][4];
  TH1F* h_Energy_MainEvtany_Eta_dRadius[6][4];
  TH1F* h_Energy_Shared_Eta_dRadius[6][4];
  TH1F* h_Energy_SharedME_Eta_dRadius[6][4];
  TH1F* h_Energy_SharedPU_Eta_dRadius[6][4];
  TH1F* h_FractionEnergy_PU_Eta_dRadius[6][4];
  TH1F* h_FractionEnergy_PUany_Eta_dRadius[6][4];
  TH1F* h_FractionEnergy_MainEvt_Eta_dRadius[6][4];
  TH1F* h_FractionEnergy_MainEvtany_Eta_dRadius[6][4];
  TH1F* h_FractionEnergy_Shared_Eta_dRadius[6][4];
  TH1F* h_FractionEnergy_SharedME_Eta_dRadius[6][4];
  TH1F* h_FractionEnergy_SharedPU_Eta_dRadius[6][4];

  TH1F* h_NumberHits_Eta_dRadius_3MIP[6][4];
  TH1F* h_NumberHits_PU_Eta_dRadius_3MIP[6][4];
  TH1F* h_NumberHits_PUany_Eta_dRadius_3MIP[6][4];
  TH1F* h_NumberHits_MainEvt_Eta_dRadius_3MIP[6][4];
  TH1F* h_NumberHits_MainEvtany_Eta_dRadius_3MIP[6][4];
  TH1F* h_NumberHits_Shared_Eta_dRadius_3MIP[6][4];
  TH1F* h_NumberHits_SharedME_Eta_dRadius_3MIP[6][4];
  TH1F* h_NumberHits_SharedPU_Eta_dRadius_3MIP[6][4];
  TH1F* h_FractionHits_PU_Eta_dRadius_3MIP[6][4];
  TH1F* h_FractionHits_PUany_Eta_dRadius_3MIP[6][4];
  TH1F* h_FractionHits_MainEvt_Eta_dRadius_3MIP[6][4];
  TH1F* h_FractionHits_MainEvtany_Eta_dRadius_3MIP[6][4];
  TH1F* h_FractionHits_Shared_Eta_dRadius_3MIP[6][4];
  TH1F* h_FractionHits_SharedME_Eta_dRadius_3MIP[6][4];
  TH1F* h_FractionHits_SharedPU_Eta_dRadius_3MIP[6][4];

  TH1F* h_Energy_Eta_dRadius_3MIP[6][4];
  TH1F* h_Energy_PU_Eta_dRadius_3MIP[6][4];
  TH1F* h_Energy_PUany_Eta_dRadius_3MIP[6][4];
  TH1F* h_Energy_MainEvt_Eta_dRadius_3MIP[6][4]; 
  TH1F* h_Energy_MainEvtany_Eta_dRadius_3MIP[6][4]; 
  TH1F* h_Energy_Shared_Eta_dRadius_3MIP[6][4]; 
  TH1F* h_Energy_SharedME_Eta_dRadius_3MIP[6][4]; 
  TH1F* h_Energy_SharedPU_Eta_dRadius_3MIP[6][4]; 
  TH1F* h_FractionEnergy_PU_Eta_dRadius_3MIP[6][4];
  TH1F* h_FractionEnergy_PUany_Eta_dRadius_3MIP[6][4];
  TH1F* h_FractionEnergy_MainEvt_Eta_dRadius_3MIP[6][4];
  TH1F* h_FractionEnergy_MainEvtany_Eta_dRadius_3MIP[6][4];
  TH1F* h_FractionEnergy_Shared_Eta_dRadius_3MIP[6][4];
  TH1F* h_FractionEnergy_SharedME_Eta_dRadius_3MIP[6][4];
  TH1F* h_FractionEnergy_SharedPU_Eta_dRadius_3MIP[6][4];


  int totEvtsEtaRadius[6][4];
  int totEvtsPUEtaRadius[6][4];
  int totEvtsPUanyEtaRadius[6][4];
  int totEvtsMainEvtEtaRadius[6][4];
  int totEvtsMainEvtanyEtaRadius[6][4];
  int totEvtsSharedEtaRadius[6][4];
  int totEvtsSharedMEEtaRadius[6][4];
  int totEvtsSharedPUEtaRadius[6][4];

  int totEvtsEtaRadius_3MIP[6][4];
  int totEvtsPUEtaRadius_3MIP[6][4];
  int totEvtsPUanyEtaRadius_3MIP[6][4];
  int totEvtsMainEvtEtaRadius_3MIP[6][4];
  int totEvtsMainEvtanyEtaRadius_3MIP[6][4];
  int totEvtsSharedEtaRadius_3MIP[6][4];
  int totEvtsSharedMEEtaRadius_3MIP[6][4];
  int totEvtsSharedPUEtaRadius_3MIP[6][4];

  float radiusEtaRad[6][4];

  int nBinsEta;
  int nBinsRad;
  float binWidth;
  float binStart;
  float binEnd;

  bool debugCOUT;
  bool debugCOUT2;
  bool debugCOUT3;
  bool debugCOUT4;

  int cellType;
  float floorValue;
  int lifeAge;
  float absTrend;

};  

float HGCalTimingAnalyzerPU::getXmax(TH1F* histo, float& YMax){

  float yVal = 0.;
  int xBin = 2;
  for(int iB=2; iB<histo->GetNbinsX()-1; ++iB){
    float locValue = histo->GetBinContent(iB-1) + histo->GetBinContent(iB) + histo->GetBinContent(iB+1);
    if(locValue/3. > yVal){
      xBin = iB;
      yVal = locValue/3.;
      YMax = yVal;
    }
    if(yVal > 0 && locValue/3. < yVal/2.) break;
    
  }
  return histo->GetBinCenter(xBin);
}


bool HGCalTimingAnalyzerPU::comparePairs(const std::pair<float, float>& i, const std::pair<float, float>& j){
  return i.first < j.first;
}



HGCalTimingAnalyzerPU::HGCalTimingAnalyzerPU(const edm::ParameterSet& iConfig) :
  detector(iConfig.getParameter<std::string >("detector")),
  rawRecHits(iConfig.getParameter<bool>("rawRecHits")),
  CFDTimeCorrection(iConfig.getParameter<int>("CFDTimeCorrection")),
  cellType(iConfig.getParameter<int>("cellType")),
  floorValue(iConfig.getParameter<double>("floorValue")),
  lifeAge(iConfig.getParameter<int>("lifeAge")),
  absTrend(iConfig.getParameter<double>("absTrend"))
{

  debugCOUT = false;
  debugCOUT2 = false;
  debugCOUT3 = false;
  debugCOUT4 = true;
  //  if(debugCOUT) std::cout<< " >>> costruttore " << std::endl;

  /*
  nBinsEta = 8;
  binWidth = 0.2;
  //nBinsEta = 3;
  nBinsRad = 4;
  //  binWidth = 0.5;
  binStart = 1.4;
  */

  //with new indications
  nBinsEta = 6;
  binWidth = 0.2;
  binStart = 1.65;
  binEnd = 2.85;
  nBinsRad = 4;

  //now do what ever initialization is needed
  usesResource("TFileService");
  myEstimator = new RecHiTimeEstimator(iConfig);
  myEstimator->setOptions(cellType, floorValue, lifeAge, absTrend);

  if(detector=="all") {
    _recHitsEE = consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("HGCEEInput"));
    _recHitsFH = consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("HGCFHInput"));
    _recHitsBH = consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("HGCBHInput"));
    algo = 1;
  }else if(detector=="EM") {
    _recHitsEE = consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("HGCEEInput"));
    algo = 2;
  }else if(detector=="HAD") {
    _recHitsFH = consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("HGCFHInput"));
    _recHitsBH = consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("HGCBHInput"));
    algo = 3;
  }
  _vtx = consumes<std::vector<TrackingVertex> >(edm::InputTag("mix","MergedTrackTruth"));
  _part = consumes<std::vector<TrackingParticle> >(edm::InputTag("mix","MergedTrackTruth"));
  _caloParticles = consumes<std::vector<CaloParticle> >(edm::InputTag("mix","MergedCaloTruth"));
  //  _multiClusters = consumes<std::vector<reco::HGCalMultiCluster> >(edm::InputTag("hgcalLayerClusters"));


  keV2fC[0] =  iConfig.getParameter<double>("HGCEE_keV2fC");
  keV2fC[1] =  iConfig.getParameter<double>("HGCHEF_keV2fC");
  keV2MIP = iConfig.getParameter<double>("HGCHB_keV2MIP");
  fCPerMIP[0] =  (iConfig.getParameter<std::vector<double> >("HGCEE_fCPerMIP")).at(0);
  fCPerMIP[1] =  (iConfig.getParameter<std::vector<double> >("HGCEE_fCPerMIP")).at(1);
  fCPerMIP[2] =  (iConfig.getParameter<std::vector<double> >("HGCEE_fCPerMIP")).at(2);

  const auto& rcorr = iConfig.getParameter<std::vector<double> >("thicknessCorrection");
  scaleCorrection.clear();
  for( auto corr : rcorr ) {
    scaleCorrection.push_back(1.0/corr);
  }

  const auto& dweights = iConfig.getParameter<std::vector<double> >("dEdXweights");
  for( auto weight : dweights ) {
    weights.push_back(weight);
 }


  keV2GeV = 1e-6;
  keV2MeV = 1e-3;

  auto sumes = consumesCollector();
  clusterTools = std::make_unique<hgcal::ClusterTools>(iConfig,sumes);

  edm::Service<TFileService> fs;


  cellThick_Rvseta = fs->make<TProfile2D>("cellThick_Rvseta", "", 500, 1.5, 3., 200, 0., 200.);
  cellThick_RvsLayer = fs->make<TProfile2D>("cellThick_RvsLayer", "", 45, 0., 45., 200, 0., 200.);


  h_Vtx_x = fs->make<TH1F>("h_Vtx_x", "", 1000, -10., 10.);
  h_Vtx_y = fs->make<TH1F>("h_Vtx_y", "", 1000, -10., 10.);
  h_Vtx_z = fs->make<TH1F>("h_Vtx_z", "", 1000, -10., 10.);

  h_Vtx_dvx = fs->make<TH1F>("h_Vtx_dvx", "", 1000, -10., 10.);
  h_Vtx_dvy = fs->make<TH1F>("h_Vtx_dvy", "", 1000, -10., 10.);
  h_Vtx_dvz = fs->make<TH1F>("h_Vtx_dvz", "", 1000, -10., 10.);


  h_AllHits_noEcut = fs->make<TH1F>("h_AllHits_noEcut", "", 1000, 0, 3e6);
  h_AllHits_3MIPcut = fs->make<TH1F>("h_AllHits_3MIPcut", "", 1000, 0, 3e6);


  h_Energy = fs->make<TH1F>("h_Energy", "", 1000., 0., 200.);
  h_Energy_PU = fs->make<TH1F>("h_Energy_PU", "", 1000., 0., 200.);
  h_Energy_MainEvt = fs->make<TH1F>("h_Energy_MainEvt", "", 1000., 0., 200.);
  h_Energy_Shared = fs->make<TH1F>("h_Energy_Shared", "", 1000., 0., 200.);
  h_Energy_SharedME = fs->make<TH1F>("h_Energy_SharedME", "", 1000., 0., 200.);
  h_Energy_SharedPU = fs->make<TH1F>("h_Energy_SharedPU", "", 1000., 0., 200.);
  h_Fraction_Energy_PU = fs->make<TH1F>("h_Fraction_Energy_PU", "", 1000., 0., 1.);
  h_Fraction_Energy_MainEvt = fs->make<TH1F>("h_Fraction_Energy_MainEvt", "", 1000., 0., 1.);
  h_Fraction_Energy_Shared = fs->make<TH1F>("h_Fraction_Energy_Shared", "", 1000., 0., 1.);
  h_Fraction_Energy_SharedME = fs->make<TH1F>("h_Fraction_Energy_SharedME", "", 1000., 0., 1.);
  h_Fraction_Energy_SharedPU = fs->make<TH1F>("h_Fraction_Energy_SharedPU", "", 1000., 0., 1.);

  h_Energy_3MIP = fs->make<TH1F>("h_Energy_3MIP", "", 1000., 0., 200.);
  h_Energy_PU_3MIP = fs->make<TH1F>("h_Energy_PU_3MIP", "", 1000., 0., 200.);
  h_Energy_MainEvt_3MIP = fs->make<TH1F>("h_Energy_MainEvt_3MIP", "", 1000., 0., 200.);
  h_Energy_Shared_3MIP = fs->make<TH1F>("h_Energy_Shared_3MIP", "", 1000., 0., 200.);
  h_Energy_SharedME_3MIP = fs->make<TH1F>("h_Energy_SharedME_3MIP", "", 1000., 0., 200.);
  h_Energy_SharedPU_3MIP = fs->make<TH1F>("h_Energy_SharedPU_3MIP", "", 1000., 0., 200.);
  h_Fraction_Energy_PU_3MIP = fs->make<TH1F>("h_Fraction_Energy_PU_3MIP", "", 1000., 0., 1.);
  h_Fraction_Energy_MainEvt_3MIP = fs->make<TH1F>("h_Fraction_Energy_MainEvt_3MIP", "", 1000., 0., 1.);
  h_Fraction_Energy_Shared_3MIP = fs->make<TH1F>("h_Fraction_Energy_Shared_3MIP", "", 1000., 0., 1.);
  h_Fraction_Energy_SharedME_3MIP = fs->make<TH1F>("h_Fraction_Energy_SharedME_3MIP", "", 1000., 0., 1.);
  h_Fraction_Energy_SharedPU_3MIP = fs->make<TH1F>("h_Fraction_Energy_SharedPU_3MIP", "", 1000., 0., 1.);

  //finoQUIOK
  for(int ieta=0; ieta<nBinsEta; ++ieta){    
    if(debugCOUT) std::cout<< " ieta from = " << (binStart+ieta*binWidth) << " to " << binStart+binWidth+ieta*binWidth << std::endl;
    for(int iRad=0; iRad<nBinsRad; ++iRad){
      
      hFractionEvents_PU_Eta_dRadius[ieta][iRad] = fs->make<TH1F>(Form("hFractionEvents_PU_Eta_dRadius_Eta%.2f-%.2f_dRadius%d",
								       (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000, 0., 1.);
      hFractionEvents_PUany_Eta_dRadius[ieta][iRad] = fs->make<TH1F>(Form("hFractionEvents_PUany_Eta_dRadius_Eta%.2f-%.2f_dRadius%d",
									  (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000, 0., 1.);
      hFractionEvents_MainEvt_Eta_dRadius[ieta][iRad] = fs->make<TH1F>(Form("hFractionEvents_MainEvt_Eta_dRadius_Eta%.2f-%.2f_dRadius%d",
									    (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000, 0., 1.);
      hFractionEvents_MainEvtany_Eta_dRadius[ieta][iRad] = fs->make<TH1F>(Form("hFractionEvents_MainEvtany_Eta_dRadius_Eta%.2f-%.2f_dRadius%d",
									       (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000, 0., 1.);
      hFractionEvents_Shared_Eta_dRadius[ieta][iRad] = fs->make<TH1F>(Form("hFractionEvents_Shared_Eta_dRadius_Eta%.2f-%.2f_dRadius%d",
									   (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000., 0., 1.);
      hFractionEvents_SharedME_Eta_dRadius[ieta][iRad] = fs->make<TH1F>(Form("hFractionEvents_SharedME_Eta_dRadius_Eta%.2f-%.2f_dRadius%d",
									     (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000., 0., 1.);
      hFractionEvents_SharedPU_Eta_dRadius[ieta][iRad] = fs->make<TH1F>(Form("hFractionEvents_SharedPU_Eta_dRadius_Eta%.2f-%.2f_dRadius%d",
									     (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000., 0., 1.);

      hFractionEvents_PU_Eta_dRadius_3MIP[ieta][iRad] = fs->make<TH1F>(Form("hFractionEvents_PU_Eta_dRadius_Eta%.2f-%.2f_dRadius%d_3MIP",
									    (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000., 0., 1.);
      hFractionEvents_PUany_Eta_dRadius_3MIP[ieta][iRad] = fs->make<TH1F>(Form("hFractionEvents_PUany_Eta_dRadius_Eta%.2f-%.2f_dRadius%d_3MIP",
									       (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000., 0., 1.);
      hFractionEvents_MainEvt_Eta_dRadius_3MIP[ieta][iRad] = fs->make<TH1F>(Form("hFractionEvents_MainEvt_Eta_dRadius_Eta%.2f-%.2f_dRadius%d_3MIP",
										 (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000., 0., 1.);   
      hFractionEvents_MainEvtany_Eta_dRadius_3MIP[ieta][iRad] = fs->make<TH1F>(Form("hFractionEvents_MainEvtany_Eta_dRadius_Eta%.2f-%.2f_dRadius%d_3MIP",
										    (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000., 0., 1.);   
      hFractionEvents_Shared_Eta_dRadius_3MIP[ieta][iRad] = fs->make<TH1F>(Form("hFractionEvents_Shared_Eta_dRadius_Eta%.2f-%.2f_dRadius%d_3MIP",
										(binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000., 0., 1.);
      hFractionEvents_SharedME_Eta_dRadius_3MIP[ieta][iRad] = fs->make<TH1F>(Form("hFractionEvents_SharedME_Eta_dRadius_Eta%.2f-%.2f_dRadius%d_3MIP",
										  (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000., 0., 1.);
      hFractionEvents_SharedPU_Eta_dRadius_3MIP[ieta][iRad] = fs->make<TH1F>(Form("hFractionEvents_SharedPU_Eta_dRadius_Eta%.2f-%.2f_dRadius%d_3MIP",
										  (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000., 0., 1.);   
      
      h_NumberHits_Eta_dRadius[ieta][iRad] = fs->make<TH1F>(Form("h_NumberHits_Eta_dRadius_Eta%.2f-%.2f_dRadius%d", 
								 (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000, 0., 1000.);
      h_NumberHits_PU_Eta_dRadius[ieta][iRad] = fs->make<TH1F>(Form("h_NumberHits_PU_Eta_dRadius_Eta%.2f-%.2f_dRadius%d", 
								    (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000, 0., 1000.);
      h_NumberHits_PUany_Eta_dRadius[ieta][iRad] = fs->make<TH1F>(Form("h_NumberHits_PUany_Eta_dRadius_Eta%.2f-%.2f_dRadius%d", 
								       (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000, 0., 1000.);
      h_NumberHits_MainEvt_Eta_dRadius[ieta][iRad] = fs->make<TH1F>(Form("h_NumberHits_MainEvt_Eta_dRadius_Eta%.2f-%.2f_dRadius%d", 
									 (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000, 0., 1000.);
      h_NumberHits_MainEvtany_Eta_dRadius[ieta][iRad] = fs->make<TH1F>(Form("h_NumberHits_MainEvtany_Eta_dRadius_Eta%.2f-%.2f_dRadius%d", 
									    (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000, 0., 1000.);
      h_NumberHits_Shared_Eta_dRadius[ieta][iRad] = fs->make<TH1F>(Form("h_NumberHits_Shared_Eta_dRadius_Eta%.2f-%.2f_dRadius%d", 
									(binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000, 0., 1000.);
      h_NumberHits_SharedME_Eta_dRadius[ieta][iRad] = fs->make<TH1F>(Form("h_NumberHits_SharedME_Eta_dRadius_Eta%.2f-%.2f_dRadius%d", 
									  (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000, 0., 1000.);
      h_NumberHits_SharedPU_Eta_dRadius[ieta][iRad] = fs->make<TH1F>(Form("h_NumberHits_SharedPU_Eta_dRadius_Eta%.2f-%.2f_dRadius%d", 
									  (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000, 0., 1000.);
      h_FractionHits_PU_Eta_dRadius[ieta][iRad] = fs->make<TH1F>(Form("h_FractionHits_PU_Eta_dRadius_Eta%.2f-%.2f_dRadius%d", 
								      (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000., 0., 1.);
      h_FractionHits_PUany_Eta_dRadius[ieta][iRad] = fs->make<TH1F>(Form("h_FractionHits_PUany_Eta_dRadius_Eta%.2f-%.2f_dRadius%d", 
									 (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000., 0., 1.);
      h_FractionHits_MainEvt_Eta_dRadius[ieta][iRad] = fs->make<TH1F>(Form("h_FractionHits_MainEvt_Eta_dRadius_Eta%.2f-%.2f_dRadius%d", 
									   (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000., 0., 1.);
      h_FractionHits_MainEvtany_Eta_dRadius[ieta][iRad] = fs->make<TH1F>(Form("h_FractionHits_MainEvtany_Eta_dRadius_Eta%.2f-%.2f_dRadius%d", 
									      (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000., 0., 1.);
      h_FractionHits_Shared_Eta_dRadius[ieta][iRad] = fs->make<TH1F>(Form("h_FractionHits_Shared_Eta_dRadius_Eta%.2f-%.2f_dRadius%d", 
									  (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000., 0., 1.);
      h_FractionHits_SharedME_Eta_dRadius[ieta][iRad] = fs->make<TH1F>(Form("h_FractionHits_SharedME_Eta_dRadius_Eta%.2f-%.2f_dRadius%d", 
									  (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000., 0., 1.);
      h_FractionHits_SharedPU_Eta_dRadius[ieta][iRad] = fs->make<TH1F>(Form("h_FractionHits_SharedPU_Eta_dRadius_Eta%.2f-%.2f_dRadius%d", 
									  (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000., 0., 1.);

      h_Energy_Eta_dRadius[ieta][iRad] = fs->make<TH1F>(Form("h_Energy_Eta_dRadius_Eta%.2f-%.2f_dRadius%d", 
							     (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000, 0., 200.);
      h_Energy_PU_Eta_dRadius[ieta][iRad] = fs->make<TH1F>(Form("h_Energy_PU_Eta_dRadius_Eta%.2f-%.2f_dRadius%d", 
								(binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000, 0., 200.);
      h_Energy_PUany_Eta_dRadius[ieta][iRad] = fs->make<TH1F>(Form("h_Energy_PUany_Eta_dRadius_Eta%.2f-%.2f_dRadius%d", 
								   (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000, 0., 200.);
      h_Energy_MainEvt_Eta_dRadius[ieta][iRad] = fs->make<TH1F>(Form("h_Energy_MainEvt_Eta_dRadius_Eta%.2f-%.2f_dRadius%d", 
								     (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000, 0., 200.);
      h_Energy_MainEvtany_Eta_dRadius[ieta][iRad] = fs->make<TH1F>(Form("h_Energy_MainEvtany_Eta_dRadius_Eta%.2f-%.2f_dRadius%d", 
									(binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000, 0., 200.);
      h_Energy_Shared_Eta_dRadius[ieta][iRad] = fs->make<TH1F>(Form("h_Energy_Shared_Eta_dRadius_Eta%.2f-%.2f_dRadius%d", 
								    (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000, 0., 200.);
      h_Energy_SharedME_Eta_dRadius[ieta][iRad] = fs->make<TH1F>(Form("h_Energy_SharedME_Eta_dRadius_Eta%.2f-%.2f_dRadius%d",
								      (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000, 0., 200.);
      h_Energy_SharedPU_Eta_dRadius[ieta][iRad] = fs->make<TH1F>(Form("h_Energy_SharedPU_Eta_dRadius_Eta%.2f-%.2f_dRadius%d",
								      (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000, 0., 200.);
      h_FractionEnergy_PU_Eta_dRadius[ieta][iRad] = fs->make<TH1F>(Form("h_FractionEnergy_PU_Eta_dRadius_Eta%.2f-%.2f_dRadius%d", 
									(binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000., 0., 1.);
      h_FractionEnergy_PUany_Eta_dRadius[ieta][iRad] = fs->make<TH1F>(Form("h_FractionEnergy_PUany_Eta_dRadius_Eta%.2f-%.2f_dRadius%d", 
									   (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000., 0., 1.);
      h_FractionEnergy_MainEvt_Eta_dRadius[ieta][iRad] = fs->make<TH1F>(Form("h_FractionEnergy_MainEvt_Eta_dRadius_Eta%.2f-%.2f_dRadius%d", 
									     (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000., 0., 1.);
      h_FractionEnergy_MainEvtany_Eta_dRadius[ieta][iRad] = fs->make<TH1F>(Form("h_FractionEnergy_MainEvtany_Eta_dRadius_Eta%.2f-%.2f_dRadius%d", 
										(binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000., 0., 1.);
      h_FractionEnergy_Shared_Eta_dRadius[ieta][iRad] = fs->make<TH1F>(Form("h_FractionEnergy_Shared_Eta_dRadius_Eta%.2f-%.2f_dRadius%d", 
									    (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000., 0., 1.);
      h_FractionEnergy_SharedME_Eta_dRadius[ieta][iRad] = fs->make<TH1F>(Form("h_FractionEnergy_SharedME_Eta_dRadius_Eta%.2f-%.2f_dRadius%d", 
									    (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000., 0., 1.);
      h_FractionEnergy_SharedPU_Eta_dRadius[ieta][iRad] = fs->make<TH1F>(Form("h_FractionEnergy_SharedPU_Eta_dRadius_Eta%.2f-%.2f_dRadius%d", 
									    (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000., 0., 1.);
      
      
      h_NumberHits_Eta_dRadius_3MIP[ieta][iRad] = fs->make<TH1F>(Form("h_NumberHits_Eta_dRadius_Eta%.2f-%.2f_dRadius%d_3MIP",
								      (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000, 0., 1000.);
      h_NumberHits_PU_Eta_dRadius_3MIP[ieta][iRad] = fs->make<TH1F>(Form("h_NumberHits_PU_Eta_dRadius_Eta%.2f-%.2f_dRadius%d_3MIP",
									 (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000, 0., 1000.);
      h_NumberHits_PUany_Eta_dRadius_3MIP[ieta][iRad] = fs->make<TH1F>(Form("h_NumberHits_PUany_Eta_dRadius_Eta%.2f-%.2f_dRadius%d_3MIP",
									    (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000, 0., 1000.);
      h_NumberHits_MainEvt_Eta_dRadius_3MIP[ieta][iRad] = fs->make<TH1F>(Form("h_NumberHits_MainEvt_Eta_dRadius_Eta%.2f-%.2f_dRadius%d_3MIP",
									      (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000, 0., 1000.);
      h_NumberHits_MainEvtany_Eta_dRadius_3MIP[ieta][iRad] = fs->make<TH1F>(Form("h_NumberHits_MainEvtany_Eta_dRadius_Eta%.2f-%.2f_dRadius%d_3MIP",
										 (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000, 0., 1000.);
      h_NumberHits_Shared_Eta_dRadius_3MIP[ieta][iRad] = fs->make<TH1F>(Form("h_NumberHits_Shared_Eta_dRadius_Eta%.2f-%.2f_dRadius%d_3MIP", 
									     (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000, 0., 1000.);
      h_NumberHits_SharedME_Eta_dRadius_3MIP[ieta][iRad] = fs->make<TH1F>(Form("h_NumberHits_SharedME_Eta_dRadius_Eta%.2f-%.2f_dRadius%d_3MIP", 
									       (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000, 0., 1000.);
      h_NumberHits_SharedPU_Eta_dRadius_3MIP[ieta][iRad] = fs->make<TH1F>(Form("h_NumberHits_SharedPU_Eta_dRadius_Eta%.2f-%.2f_dRadius%d_3MIP", 
									       (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000, 0., 1000.);
      h_FractionHits_PU_Eta_dRadius_3MIP[ieta][iRad] = fs->make<TH1F>(Form("h_FractionHits_PU_Eta_dRadius_Eta%.2f-%.2f_dRadius%d_3MIP", 
									   (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000., 0., 1.);
      h_FractionHits_PUany_Eta_dRadius_3MIP[ieta][iRad] = fs->make<TH1F>(Form("h_FractionHits_PUany_Eta_dRadius_Eta%.2f-%.2f_dRadius%d_3MIP", 
									      (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000., 0., 1.);
      h_FractionHits_MainEvt_Eta_dRadius_3MIP[ieta][iRad] = fs->make<TH1F>(Form("h_FractionHits_MainEvt_Eta_dRadius_Eta%.2f-%.2f_dRadius%d_3MIP", 
										(binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000., 0., 1.);
      h_FractionHits_MainEvtany_Eta_dRadius_3MIP[ieta][iRad] = fs->make<TH1F>(Form("h_FractionHits_MainEvtany_Eta_dRadius_Eta%.2f-%.2f_dRadius%d_3MIP", 
										   (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000., 0., 1.);
      h_FractionHits_Shared_Eta_dRadius_3MIP[ieta][iRad] = fs->make<TH1F>(Form("h_FractionHits_Shared_Eta_dRadius_Eta%.2f-%.2f_dRadius%d_3MIP", 
									       (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000., 0., 1.);
      h_FractionHits_SharedME_Eta_dRadius_3MIP[ieta][iRad] = fs->make<TH1F>(Form("h_FractionHits_SharedME_Eta_dRadius_Eta%.2f-%.2f_dRadius%d_3MIP", 
									       (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000., 0., 1.);
      h_FractionHits_SharedPU_Eta_dRadius_3MIP[ieta][iRad] = fs->make<TH1F>(Form("h_FractionHits_SharedPU_Eta_dRadius_Eta%.2f-%.2f_dRadius%d_3MIP", 
									       (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000., 0., 1.);

      h_Energy_Eta_dRadius_3MIP[ieta][iRad] = fs->make<TH1F>(Form("h_Energy_Eta_dRadius_Eta%.2f-%.2f_dRadius%d_3MIP",
								  (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000, 0., 200.);
      h_Energy_PU_Eta_dRadius_3MIP[ieta][iRad] = fs->make<TH1F>(Form("h_Energy_PU_Eta_dRadius_Eta%.2f-%.2f_dRadius%d_3MIP",
								     (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000, 0., 200.);
      h_Energy_PUany_Eta_dRadius_3MIP[ieta][iRad] = fs->make<TH1F>(Form("h_Energy_PUany_Eta_dRadius_Eta%.2f-%.2f_dRadius%d_3MIP",
									(binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000, 0., 200.);
      h_Energy_MainEvt_Eta_dRadius_3MIP[ieta][iRad] = fs->make<TH1F>(Form("h_Energy_MainEvt_Eta_dRadius_Eta%.2f-%.2f_dRadius%d_3MIP",
									  (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000, 0., 200.);
      h_Energy_MainEvtany_Eta_dRadius_3MIP[ieta][iRad] = fs->make<TH1F>(Form("h_Energy_MainEvtany_Eta_dRadius_Eta%.2f-%.2f_dRadius%d_3MIP",
									     (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000, 0., 200.);
      h_Energy_Shared_Eta_dRadius_3MIP[ieta][iRad] = fs->make<TH1F>(Form("h_Energy_Shared_Eta_dRadius_Eta%.2f-%.2f_dRadius%d_3MIP",
									 (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000, 0., 200.);
      h_Energy_SharedME_Eta_dRadius_3MIP[ieta][iRad] = fs->make<TH1F>(Form("h_Energy_SharedME_Eta_dRadius_Eta%.2f-%.2f_dRadius%d_3MIP",
									   (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000, 0., 200.);
      h_Energy_SharedPU_Eta_dRadius_3MIP[ieta][iRad] = fs->make<TH1F>(Form("h_Energy_SharedPU_Eta_dRadius_Eta%.2f-%.2f_dRadius%d_3MIP",
									   (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000, 0., 200.);
      h_FractionEnergy_PU_Eta_dRadius_3MIP[ieta][iRad] = fs->make<TH1F>(Form("h_FractionEnergy_PU_Eta_dRadius_Eta%.2f-%.2f_dRadius%d_3MIP", 
									     (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000., 0., 1.);
      h_FractionEnergy_PUany_Eta_dRadius_3MIP[ieta][iRad] = fs->make<TH1F>(Form("h_FractionEnergy_PUany_Eta_dRadius_Eta%.2f-%.2f_dRadius%d_3MIP", 
										(binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000., 0., 1.);
      h_FractionEnergy_MainEvt_Eta_dRadius_3MIP[ieta][iRad] = fs->make<TH1F>(Form("h_FractionEnergy_MainEvt_Eta_dRadius_Eta%.2f-%.2f_dRadius%d_3MIP", 
										  (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000., 0., 1.);
      h_FractionEnergy_MainEvtany_Eta_dRadius_3MIP[ieta][iRad] = fs->make<TH1F>(Form("h_FractionEnergy_MainEvtany_Eta_dRadius_Eta%.2f-%.2f_dRadius%d_3MIP", 
										     (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000., 0., 1.);
      h_FractionEnergy_Shared_Eta_dRadius_3MIP[ieta][iRad] = fs->make<TH1F>(Form("h_FractionEnergy_Shared_Eta_dRadius_Eta%.2f-%.2f_dRadius%d_3MIP", 
										 (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000., 0., 1.);
      h_FractionEnergy_SharedME_Eta_dRadius_3MIP[ieta][iRad] = fs->make<TH1F>(Form("h_FractionEnergy_SharedME_Eta_dRadius_Eta%.2f-%.2f_dRadius%d_3MIP", 
										 (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000., 0., 1.);
      h_FractionEnergy_SharedPU_Eta_dRadius_3MIP[ieta][iRad] = fs->make<TH1F>(Form("h_FractionEnergy_SharedPU_Eta_dRadius_Eta%.2f-%.2f_dRadius%d_3MIP", 
										 (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000., 0., 1.);

    }
  }

  for(int iEta=0; iEta<nBinsEta; ++iEta){
    for(int iRad=0; iRad<nBinsRad; ++iRad){
	totEvtsEtaRadius[iEta][iRad] = 0;
	totEvtsPUEtaRadius[iEta][iRad] = 0;
	totEvtsPUanyEtaRadius[iEta][iRad] = 0;
	totEvtsMainEvtEtaRadius[iEta][iRad] = 0;
	totEvtsMainEvtanyEtaRadius[iEta][iRad] = 0;
	totEvtsSharedEtaRadius[iEta][iRad] = 0;
	totEvtsSharedMEEtaRadius[iEta][iRad] = 0;
	totEvtsSharedPUEtaRadius[iEta][iRad] = 0;

	totEvtsEtaRadius_3MIP[iEta][iRad] = 0;
	totEvtsPUEtaRadius_3MIP[iEta][iRad] = 0;
	totEvtsPUanyEtaRadius_3MIP[iEta][iRad] = 0;
	totEvtsMainEvtEtaRadius_3MIP[iEta][iRad] = 0;
	totEvtsMainEvtanyEtaRadius_3MIP[iEta][iRad] = 0;
	totEvtsSharedEtaRadius_3MIP[iEta][iRad] = 0;
	totEvtsSharedMEEtaRadius_3MIP[iEta][iRad] = 0;
	totEvtsSharedPUEtaRadius_3MIP[iEta][iRad] = 0;
    }
    radiusEtaRad[iEta][0] = 2.;
    radiusEtaRad[iEta][1] = 5.;
    radiusEtaRad[iEta][2] = 10.;
    radiusEtaRad[iEta][3] = 6000.;
  }
  
}

HGCalTimingAnalyzerPU::~HGCalTimingAnalyzerPU()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


void
HGCalTimingAnalyzerPU::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  if(debugCOUT) std::cout<< " >>> analyzer " << std::endl;
  using namespace edm;

  myEstimator->setEventSetup(iSetup);

  recHitTools.getEventSetup(iSetup);

  Handle<HGCRecHitCollection> recHitHandleEE;
  Handle<HGCRecHitCollection> recHitHandleFH;
  Handle<HGCRecHitCollection> recHitHandleBH;

  Handle<std::vector<TrackingVertex> > vtxHandle;
  Handle<std::vector<TrackingParticle> > partHandle;
  iEvent.getByToken(_vtx,vtxHandle);
  iEvent.getByToken(_part,partHandle);
  const std::vector<TrackingVertex>& vtxs = *vtxHandle;
  const std::vector<TrackingParticle>& part = *partHandle;

  Handle<std::vector<CaloParticle> > caloParticleHandle;
  iEvent.getByToken(_caloParticles, caloParticleHandle);
  const std::vector<CaloParticle>& caloParticles = *caloParticleHandle;

  // Handle<std::vector<reco::HGCalMultiCluster> > multiClusterHandle;
  // iEvent.getByToken(_multiClusters, multiClusterHandle);
  // const std::vector<reco::HGCalMultiCluster>& multiClusters = *multiClusterHandle;
  

  //  int reachedEE_posZ = 1;

  float vx = 0.;
  float vy = 0.;
  float vz = 0.;
  if(vtxs.size()!=0){
    vx = vtxs[0].position().x();
    vy = vtxs[0].position().y();
    vz = vtxs[0].position().z();
  }

  h_Vtx_x->Fill(vx);
  h_Vtx_y->Fill(vy);
  h_Vtx_z->Fill(vz);


  //  std::cout << " part.size() = " << part.size() << std::endl;
  for(unsigned int i=0;i<part.size();++i){
    //    std::cout << " part i " << i << std::endl;
    //if(part[i].parentVertex()->nGenVertices()>0){
    float dvx=0.;
    float dvy=0.;
    float dvz=0.;
    
    //    std::cout << " part i  decayVertices().size() = " << part[i].decayVertices().size() << " part[i].eta() = " << part[i].eta() << std::endl;
    if(part[i].decayVertices().size()>=1){   
      dvx=part[i].decayVertices()[0]->position().x();
      dvy=part[i].decayVertices()[0]->position().y();
      dvz=part[i].decayVertices()[0]->position().z();
      //      if( part[i].decayVertices()[0]->inVolume() && part[i].eta() >= 0.) reachedEE_posZ = 0;
      h_Vtx_dvx->Fill(dvx);
      h_Vtx_dvy->Fill(dvy);
      h_Vtx_dvz->Fill(dvz);
    }
  }
  

  HGCRecHitCollection NewrechitsEE;
  HGCRecHitCollection NewrechitsFH;

  //make a map detid-rechit
  std::map<DetId,const HGCRecHit*> hitmap;
  switch(algo){
  case 1:
    {
      iEvent.getByToken(_recHitsEE,recHitHandleEE);
      iEvent.getByToken(_recHitsFH,recHitHandleFH);
      iEvent.getByToken(_recHitsBH,recHitHandleBH);

      /*
      const HGCRecHitCollection& rechitsEEOld = *recHitHandleEE;
      const HGCRecHitCollection& rechitsFHOld = *recHitHandleFH;
      if(CFDTimeCorrection == 1){
	myEstimator->correctTime(rechitsEEOld, &NewrechitsEE);
	myEstimator->correctTime(rechitsFHOld, &NewrechitsFH);
      }
      else if(CFDTimeCorrection == 0){
	myEstimator->correctTimeFixThr(rechitsEEOld, &NewrechitsEE);
	myEstimator->correctTimeFixThr(rechitsFHOld, &NewrechitsFH);
      }
      else if(CFDTimeCorrection == -1){
	NewrechitsEE = *recHitHandleEE;
	NewrechitsFH = *recHitHandleFH;
      }
      */

      NewrechitsEE = *recHitHandleEE;
      NewrechitsFH = *recHitHandleFH;
      for(unsigned int i = 0; i < NewrechitsEE.size(); ++i){
	hitmap[NewrechitsEE[i].detid()] = &NewrechitsEE[i];
      }
      for(unsigned int i = 0; i < NewrechitsFH.size(); ++i){
	hitmap[NewrechitsFH[i].detid()] = &NewrechitsFH[i];
      }
      break;
    }
  case 2:
    {
      iEvent.getByToken(_recHitsEE,recHitHandleEE);
      /*
      const HGCRecHitCollection& rechitsEEOld = *recHitHandleEE;

      if(CFDTimeCorrection == 1) myEstimator->correctTime(rechitsEEOld, &NewrechitsEE);
      else if(CFDTimeCorrection == 0) myEstimator->correctTimeFixThr(rechitsEEOld, &NewrechitsEE);
      else if(CFDTimeCorrection == -1){
	NewrechitsEE = *recHitHandleEE;
      }
      */
      NewrechitsEE = *recHitHandleEE;
      for(unsigned int i = 0; i < NewrechitsEE.size(); i++){
	hitmap[NewrechitsEE[i].detid()] = &NewrechitsEE[i];
      }
      break;
    }
  case 3:
    {
      iEvent.getByToken(_recHitsFH,recHitHandleFH);
      iEvent.getByToken(_recHitsBH,recHitHandleBH);
      /*
      const HGCRecHitCollection& rechitsFHOld = *recHitHandleFH;
      const HGCRecHitCollection& rechitsBHOld = *recHitHandleBH;
      if(CFDTimeCorrection == 1){
	myEstimator->correctTime(rechitsFHOld, &NewrechitsFH);
      }
      else if(CFDTimeCorrection == 0){
	myEstimator->correctTimeFixThr(rechitsFHOld, &NewrechitsFH);
      }
      else if(CFDTimeCorrection == -1){
	NewrechitsFH = *recHitHandleFH;
      }
      */
      NewrechitsFH = *recHitHandleFH;
      for(unsigned int i = 0; i < NewrechitsFH.size(); i++){
	hitmap[NewrechitsFH[i].detid()] = &NewrechitsFH[i];
      }
      break;
    }
  default:
    break;
  }

  ////////////////////
  float etaGen = -1;
  float phiGen = -1;
  float xGen = -1;
  float yGen = -1;
  float zGen = -1;
  bool evtGood = false;

  SimClusterRefVector simClusterRefVectorChosen;
  std::map<DetId,float> hitsAndFractionChosen;

  if(debugCOUT) std::cout<< " >>> now caloparticles " << std::endl;

  //  std::array<double,3> vtx{ {vx, vy, vz}};
  std::array<double,3> vtx{ {0., 0., 0.}};
  int numbCaloPart = 0;
  if(debugCOUT4) std::cout << " >>> caloParticles.size() = " << caloParticles.size() << std::endl;
  // loop over caloParticles
  for (std::vector<CaloParticle>::const_iterator it_caloPart = caloParticles.begin(); it_caloPart != caloParticles.end(); ++it_caloPart){
    const SimClusterRefVector simClusterRefVector = it_caloPart->simClusters();

    if(debugCOUT4)    std::cout << " caloParticles loop 1 simClusterRefVector.size() = " << simClusterRefVector.size() << " eta = " << it_caloPart->eta() << " pdgID = " 
				<< it_caloPart->pdgId() << " energy = " << it_caloPart->pt() << " eventId().event() = " << it_caloPart->eventId().event() 
				<< " eventId().bunchCrossing() = " << it_caloPart->eventId().bunchCrossing() << std::endl; 
    if(CFDTimeCorrection == 22 && (simClusterRefVector.size() > 1 || it_caloPart->pdgId() != 22 || it_caloPart->pt() != 5.)) continue;
    //if(simClusterRefVector.size() > 1 || it_caloPart->pdgId() != 22 || it_caloPart->pt() != 60.) continue;
    if(CFDTimeCorrection == 130 && (simClusterRefVector.size() > 1 || it_caloPart->pdgId() != 130 || it_caloPart->pt() != 5.)) continue;
    if(CFDTimeCorrection == -22 && (simClusterRefVector.size() > 1 || it_caloPart->pdgId() != 22 || it_caloPart->pt() != 5.)) continue;
    if(CFDTimeCorrection == -130 && (simClusterRefVector.size() > 1 || it_caloPart->pdgId() != 130 || it_caloPart->pt() != 5.)) continue;
    if(CFDTimeCorrection == 14 && (simClusterRefVector.size() > 1 || 
				   (it_caloPart->pdgId() != 12 && it_caloPart->pdgId() != 14 && it_caloPart->pdgId() != 16) )) continue;
    if(CFDTimeCorrection == 12 && (simClusterRefVector.size() > 1 || it_caloPart->pdgId() != 12 || it_caloPart->pt() != 15. )) continue;

    etaGen = it_caloPart->eta();
    phiGen = it_caloPart->phi();
    xGen = it_caloPart->momentum().x();
    yGen = it_caloPart->momentum().y();
    zGen = it_caloPart->momentum().z();

    if(debugCOUT) std::cout<< " bau caloparticles inloop " << std::endl;

    if(debugCOUT2)    std::cout<< " >>> befor caloparticles eta = " << etaGen << std::endl;
    if(etaGen < 0) continue;
    if(etaGen < binStart || etaGen > binEnd) continue;
    ++numbCaloPart;
    if(numbCaloPart > 1) continue;

    evtGood = true;
    simClusterRefVectorChosen = simClusterRefVector;
    const SimCluster simClusterChosen = (*(*(simClusterRefVectorChosen.begin())));
    const std::vector<std::pair<uint32_t,float> > hits_and_fractions = simClusterChosen.hits_and_fractions();
    //loop over hits  
    //    for (std::vector<std::pair<uint32_t,float> >::const_iterator it_haf = hits_and_fractions.begin(); it_haf != hits_and_fractions.end(); ++it_haf) {
    //    }
    
    if(debugCOUT4) std::cout<< " bau caloparticles survived " << std::endl;

    
    if(debugCOUT2)    std::cout<< " >>> now caloparticles " << std::endl;


    //shower axis by recHits
    float axisX = 0;
    float axisY = 0;
    float axisZ = 0;
    float sumEnergyToNorm = 0;
    GlobalPoint showerAxis;
    
    if(debugCOUT) std::cout<< " bau before showerAxis " << std::endl;
  
    //loop on rechit - matched to gen => shower axis
    for (std::vector<std::pair<uint32_t,float> >::const_iterator it_haf = hits_and_fractions.begin(); it_haf != hits_and_fractions.end(); ++it_haf) {
      hitsAndFractionChosen[it_haf->first] = it_haf->second;
      DetId hitid = (it_haf->first);
      float rhEnergy = 0;
	
      std::map<DetId, const HGCRecHit*>::iterator trovatore = hitmap.find(hitid);
      if(trovatore == hitmap.end()){
	continue;
      }
      else if(recHitTools.getEta(hitid)*it_caloPart->eta() < 0){
	continue;
	//if(debugCOUT) std::cout<< " >>> MAJOR PROBLEM!!! " << std::endl;    
      }
      else{
	const HGCRecHit *hit = hitmap[hitid];

	rhEnergy = hit->energy();
	sumEnergyToNorm += rhEnergy*it_haf->second;
	axisX += rhEnergy*it_haf->second * recHitTools.getPosition(hitid).x();
	axisY += rhEnergy*it_haf->second * recHitTools.getPosition(hitid).y();
	axisZ += rhEnergy*it_haf->second * recHitTools.getPosition(hitid).z();
      }
    }
    showerAxis = GlobalPoint(axisX/sumEnergyToNorm, axisY/sumEnergyToNorm, (axisZ - vz)/sumEnergyToNorm);
    
    float axEta = showerAxis.eta();
    float axPhi = showerAxis.phi();
    float axX = showerAxis.x();
    float axY = showerAxis.y();
    float axZ = showerAxis.z();

    if(evtGood) break;
  }
  

  if(!evtGood) return;

  if(debugCOUT4){
    const SimCluster simClusterBeg = (*(*(simClusterRefVectorChosen.begin())));
    std::cout<< " evtGood = " << evtGood << " simClusterRefVectorChosen.size = " << simClusterRefVectorChosen.size() << " eta = " << simClusterBeg.eta() << std::endl;
  }


  
  if(debugCOUT) std::cout<< " bau after showerAxis " << std::endl;
  /////////////////////////////////////////////////////////////////
  UtilClasses utilsMet = UtilClasses(etaGen, phiGen);
  std::array<double,3> fromAxis;
  if(CFDTimeCorrection > 0) {
    fromAxis[0] = (xGen);
    fromAxis[1] = (yGen);
    fromAxis[2] = (zGen);
  }    
  if(CFDTimeCorrection < 0) {
    fromAxis[0] = (-1.*xGen);
    fromAxis[1] = (-1.*yGen);
    fromAxis[2] = (-1.*zGen);
  }

    //std::array<double,3> fromAxis{ {-1.*xGen, -1.*yGen, -1.*zGen} };
  //UtilClasses utilsMet = UtilClasses(axEta, axPhi);
  //std::array<double,3> fromAxis{ {axX, axY, axZ} };
  
  int totHitsAbove60 = 0;
  int totHitsAbove0 = 0;
  
  bool withTime_filled = false;
  bool allTime_filled = false;
  int totHitTime = 0.;
  
  if(debugCOUT) std::cout<< " now unmatched " << std::endl;
  
  if(debugCOUT) std::cout<< " bau 1 " << std::endl;
  
  float totEnergyPerEtaRadius[6][4];
  float totEnergyPUPerEtaRadius[6][4];
  float totEnergyPUanyPerEtaRadius[6][4];
  float totEnergyMainEvtPerEtaRadius[6][4];
  float totEnergyMainEvtanyPerEtaRadius[6][4];
  float totEnergySharedPerEtaRadius[6][4];
  float totEnergySharedMEPerEtaRadius[6][4];
  float totEnergySharedPUPerEtaRadius[6][4];

  int totRHPerEtaRadius[6][4];
  int totRHPUPerEtaRadius[6][4];
  int totRHPUanyPerEtaRadius[6][4];
  int totRHMainEvtPerEtaRadius[6][4];
  int totRHMainEvtanyPerEtaRadius[6][4];
  int totRHSharedPerEtaRadius[6][4];
  float totRHSharedMEPerEtaRadius[6][4];
  float totRHSharedPUPerEtaRadius[6][4];

  
  float totEnergyPerEtaRadius_3MIP[6][4];
  float totEnergyPUPerEtaRadius_3MIP[6][4];
  float totEnergyPUanyPerEtaRadius_3MIP[6][4];
  float totEnergyMainEvtPerEtaRadius_3MIP[6][4];
  float totEnergyMainEvtanyPerEtaRadius_3MIP[6][4];
  float totEnergySharedPerEtaRadius_3MIP[6][4];
  float totEnergySharedMEPerEtaRadius_3MIP[6][4];
  float totEnergySharedPUPerEtaRadius_3MIP[6][4];

  int totRHPerEtaRadius_3MIP[6][4];
  int totRHPUPerEtaRadius_3MIP[6][4];
  int totRHPUanyPerEtaRadius_3MIP[6][4];
  int totRHMainEvtPerEtaRadius_3MIP[6][4];
  int totRHMainEvtanyPerEtaRadius_3MIP[6][4];
  int totRHSharedPerEtaRadius_3MIP[6][4];
  float totRHSharedMEPerEtaRadius_3MIP[6][4];
  float totRHSharedPUPerEtaRadius_3MIP[6][4];


  float totEnergy = 0.;
  float totEnergy_PU = 0.;
  float totEnergy_MainEvt = 0.;
  float totEnergy_Shared = 0.;
  float totEnergy_SharedME = 0.;
  float totEnergy_SharedPU = 0.;

  float totEnergy_3MIP = 0.;
  float totEnergy_PU_3MIP = 0.;
  float totEnergy_MainEvt_3MIP = 0.;
  float totEnergy_Shared_3MIP = 0.;
  float totEnergy_SharedME_3MIP = 0.;
  float totEnergy_SharedPU_3MIP = 0.;


  for(int iet=0; iet<nBinsEta; ++iet){
    for(int irad=0; irad<nBinsRad; ++irad){
      totEnergyPerEtaRadius[iet][irad] = 0.;
      totEnergyPUPerEtaRadius[iet][irad] = 0.;
      totEnergyPUanyPerEtaRadius[iet][irad] = 0.;
      totEnergyMainEvtPerEtaRadius[iet][irad] = 0.;
      totEnergyMainEvtanyPerEtaRadius[iet][irad] = 0.;
      totEnergySharedPerEtaRadius[iet][irad] = 0.;
      totEnergySharedMEPerEtaRadius[iet][irad] = 0;
      totEnergySharedPUPerEtaRadius[iet][irad] = 0;
      
      totRHPerEtaRadius[iet][irad] = 0;
      totRHPUPerEtaRadius[iet][irad] = 0;
      totRHPUanyPerEtaRadius[iet][irad] = 0;
      totRHMainEvtPerEtaRadius[iet][irad] = 0;
      totRHMainEvtanyPerEtaRadius[iet][irad] = 0;
      totRHSharedPerEtaRadius[iet][irad] = 0;
      totRHSharedMEPerEtaRadius[iet][irad] = 0;
      totRHSharedPUPerEtaRadius[iet][irad] = 0;
      
      
      totEnergyPerEtaRadius_3MIP[iet][irad] = 0.;
      totEnergyPUPerEtaRadius_3MIP[iet][irad] = 0.;
      totEnergyPUanyPerEtaRadius_3MIP[iet][irad] = 0.;
      totEnergyMainEvtPerEtaRadius_3MIP[iet][irad] = 0.;
      totEnergyMainEvtanyPerEtaRadius_3MIP[iet][irad] = 0.;
      totEnergySharedPerEtaRadius_3MIP[iet][irad] = 0.;
      totEnergySharedMEPerEtaRadius_3MIP[iet][irad] = 0;
      totEnergySharedPUPerEtaRadius_3MIP[iet][irad] = 0;
      
      totRHPerEtaRadius_3MIP[iet][irad] = 0;
      totRHPUPerEtaRadius_3MIP[iet][irad] = 0;
      totRHPUanyPerEtaRadius_3MIP[iet][irad] = 0;
      totRHMainEvtPerEtaRadius_3MIP[iet][irad] = 0;
      totRHMainEvtanyPerEtaRadius_3MIP[iet][irad] = 0;
      totRHSharedPerEtaRadius_3MIP[iet][irad] = 0;
      totRHSharedMEPerEtaRadius_3MIP[iet][irad] = 0;
      totRHSharedPUPerEtaRadius_3MIP[iet][irad] = 0;
    }
  }
  
  int AllHitsNoEcut = 0.;
  int AllHits3MIPcut = 0.;
  
  /////////////////////////////
  for(std::map<DetId, const HGCRecHit*>::iterator iop=hitmap.begin(); iop != hitmap.end(); ++iop){
    const HGCalDetId hitid = iop->first;
    const HGCRecHit* hit = iop->second;
    
    bool found = false;
    float rhEnergy = hit->energy();
    float CPfraction = 0.;
    float rhX = recHitTools.getPosition(hitid).x();
    float rhY = recHitTools.getPosition(hitid).y();
    int rhL = recHitTools.getLayerWithOffset(hitid);
    float rhEta = recHitTools.getEta(recHitTools.getPosition(hitid));
    float rhZ = utilsMet.layerToZ(rhL, rhEta);
    float rhPt = rhEnergy/cosh(rhEta);
    //    bool skipEvt = false;

    if(debugCOUT3)    std::cout << " loop over hits eta = " << rhEta  << std::endl;
    //FIXME
    if(rhEta < 0 && CFDTimeCorrection > 0) continue;
    if(rhEta > 0 && CFDTimeCorrection < 0) continue;
    
    std::array<double,3> to{ {0., 0., rhZ} };
    utilsMet.layerIntersection(to, fromAxis, vtx);
    float deltaR = sqrt(pow(to[0] - rhX, 2) + pow(to[1] - rhY, 2));
    
    int etaBin = int((std::abs(etaGen) - binStart) / binWidth);
    int iRadBin = -1;
    //identify bin of radius around the gen axis
    for(int ir=0; ir<nBinsRad; ++ir){
      if(deltaR < radiusEtaRad[etaBin][ir]) {
	iRadBin = ir;
	break;
      }
    }
    // if(deltaR > 600) continue;
  
    // if(iRadBin != -1) {
    //   for(int ir=iRadBin; ir<nBinsRad; ++ir){
    // 	totRHPerEtaRadius[etaBin][ir] += 1;
    // 	totEnergyPerEtaRadius[etaBin][ir] += rhPt;
    //   }//loop over radii
    // }//radius ok


    if(debugCOUT4)    std::cout << " radius = " << deltaR  << std::endl;
    
    ///do matching with caloparticle
    // loop over caloParticles
    // for (std::vector<CaloParticle>::const_iterator it_caloPart = caloParticles.begin(); it_caloPart != caloParticles.end(); ++it_caloPart){
    //   if(it_caloPart->eta() < 0) continue;
    //   const SimClusterRefVector simClusterRefVector = it_caloPart->simClusters();
    //   if(simClusterRefVector.size() > 1) continue;
    // for (CaloParticle::sc_iterator it_sc = simClusterRefVectorChosen.begin(); it_sc != simClusterRefVectorChosen.end(); ++it_sc) {
    //   const SimCluster simCluster = (*(*it_sc));
    //   if(simCluster.eta() < 0.)continue;
    //   const std::vector<std::pair<uint32_t,float> > hits_and_fractions = simCluster.hits_and_fractions();     
    //   //loop over hits
    //   for (std::vector<std::pair<uint32_t,float> >::const_iterator it_haf = hits_and_fractions.begin(); it_haf != hits_and_fractions.end(); ++it_haf) {

    std::map<DetId, float>::iterator trovatoreCP = hitsAndFractionChosen.find(hitid);
    if(trovatoreCP != hitsAndFractionChosen.end()){
      found = true;
      CPfraction = float(hitsAndFractionChosen[hitid]);
    }

    //    if(debugCOUT3)    std::cout << " loop over caloparticle 2nd => skipEvt = " << skipEvt << " found = main event = " << found  << std::endl;
    // if(skipEvt == true) continue;
    

    if(debugCOUT4) std::cout << " >>> rhPt = " << rhPt << " rhEnergy = " << rhEnergy << " rhEta = " << rhEta << " found = " << found << std::endl;

    unsigned int layer = recHitTools.getLayerWithOffset(hitid);
    int thick = (hitid.det() != DetId::Forward) ? -1 : recHitTools.getSiThickness(hitid) / 100. - 1.;
    
    int sectionType = -1;
    if(hitid.subdetId() == HGCEE) sectionType = 0;
    else if(hitid.subdetId() == HGCHEF) sectionType = 1;
    else if(hitid.subdetId() == HGCHEB) sectionType = 2;
    
    int energyMIP = 0.;
    if(sectionType == 2) energyMIP = hit->energy()/keV2GeV * keV2MIP;
    else if(sectionType == 0 || sectionType == 1) energyMIP = hit->energy()/scaleCorrection.at(thick)/keV2GeV / (weights.at(layer)/keV2MeV);
    
    float energyCharge = 0.;
    if(sectionType == 2) energyCharge = energyMIP * 1.;
    else if(sectionType == 0 || sectionType == 1) energyCharge = energyMIP * fCPerMIP[thick];
    
    
    float charge = energyCharge;
    float MIP = energyMIP;
    

    cellThick_Rvseta->Fill(rhEta, sqrt(rhX*rhY+rhX*rhY), thick+1);
    cellThick_RvsLayer->Fill(rhL, sqrt(rhX*rhY+rhX*rhY), thick+1);

    if(debugCOUT4) std::cout << " >>> iRadBin = " << iRadBin << " MIP = " << MIP << " sectionType = " << sectionType << " thick = " << thick << " CPfraction = " << CPfraction << std::endl;

    //    std::cout << " CPfraction = " << CPfraction << " < 1.  = " << (CPfraction-1 < -0.001) << "  == 1. = " << (abs(CPfraction-1) < 0.001) << std::endl;

    if(iRadBin != -1) {
      for(int ir=iRadBin; ir<nBinsRad; ++ir){
	totRHPerEtaRadius[etaBin][ir] += 1;
	totEnergyPerEtaRadius[etaBin][ir] += rhPt;
	if(MIP > 3){
	  totRHPerEtaRadius_3MIP[etaBin][ir] += 1;
	  totEnergyPerEtaRadius_3MIP[etaBin][ir] += rhPt;
	}
	if(found){
	  if((CPfraction-1) < - 0.0001){
	    // std::cout << " first CPfraction = " << CPfraction << std::endl;
	    // std::cout << " >>> rhPt = " << rhPt << " rhEnergy = " << rhEnergy << " rhEta = " << rhEta << " found = " << found << std::endl;
	    totRHSharedPerEtaRadius[etaBin][ir] += 1;
	    totEnergySharedPerEtaRadius[etaBin][ir] += rhPt;

	    totRHSharedMEPerEtaRadius[etaBin][ir] += 1;
	    totEnergySharedMEPerEtaRadius[etaBin][ir] += rhPt*CPfraction;

	    totRHSharedPUPerEtaRadius[etaBin][ir] += 1;
	    totEnergySharedPUPerEtaRadius[etaBin][ir] += rhPt*(1. - CPfraction);

	    totRHMainEvtanyPerEtaRadius[etaBin][ir] += 1;
            totEnergyMainEvtanyPerEtaRadius[etaBin][ir] += rhPt*CPfraction;

	    totRHPUanyPerEtaRadius[etaBin][ir] += 1;
	    totEnergyPUanyPerEtaRadius[etaBin][ir] += rhPt*(1. - CPfraction);
	  }
	  else if(abs(CPfraction-1) < 0.0001){
	    // std::cout << " else if CPfraction = " << CPfraction << std::endl;
	    // std::cout << " >>> rhPt = " << rhPt << " rhEnergy = " << rhEnergy << " rhEta = " << rhEta << " found = " << found << std::endl;
	    totRHMainEvtPerEtaRadius[etaBin][ir] += 1;
	    totEnergyMainEvtPerEtaRadius[etaBin][ir] += rhPt;

	    totRHMainEvtanyPerEtaRadius[etaBin][ir] += 1;
            totEnergyMainEvtanyPerEtaRadius[etaBin][ir] += rhPt*CPfraction;
	  }
	  // else{
	  //   std::cout << " else CPfraction = " << CPfraction << std::endl;
	  //   std::cout << " >>> rhPt = " << rhPt << " rhEnergy = " << rhEnergy << " rhEta = " << rhEta << " found = " << found << std::endl;
	  //   // totRHPUPerEtaRadius_3MIP[etaBin][ir] += 1;
	  //   // totEnergyPUPerEtaRadius_3MIP[etaBin][ir] += rhPt;
	  // }
	  if(MIP > 3){
	    if( (CPfraction-1) < -0.00001){
	    // std::cout << "3MIP  first CPfraction = " << CPfraction << std::endl;
	    // std::cout << " >>> rhPt = " << rhPt << " rhEnergy = " << rhEnergy << " rhEta = " << rhEta << " found = " << found << std::endl;

	      totRHSharedPerEtaRadius_3MIP[etaBin][ir] += 1;
	      totEnergySharedPerEtaRadius_3MIP[etaBin][ir] += rhPt;	  

	      if(MIP*CPfraction > 3){
		totRHSharedMEPerEtaRadius_3MIP[etaBin][ir] += 1;
		totEnergySharedMEPerEtaRadius_3MIP[etaBin][ir] += rhPt*CPfraction;

		totRHMainEvtanyPerEtaRadius_3MIP[etaBin][ir] += 1;
		totEnergyMainEvtanyPerEtaRadius_3MIP[etaBin][ir] += rhPt*CPfraction;		
	      }
	      if(MIP * (1. - CPfraction) > 3){
		totRHSharedPUPerEtaRadius_3MIP[etaBin][ir] += 1;
		totEnergySharedPUPerEtaRadius_3MIP[etaBin][ir] += rhPt*(1. - CPfraction);

		totRHPUanyPerEtaRadius_3MIP[etaBin][ir] += 1;
		totEnergyPUanyPerEtaRadius_3MIP[etaBin][ir] += rhPt*(1. - CPfraction);
	      }	    
	    }
	    else if(abs(CPfraction-1) < 0.00001){
	      // std::cout << " elseif 3MIP CPfraction = " << CPfraction << std::endl;
	      // std::cout << " >>> rhPt = " << rhPt << " rhEnergy = " << rhEnergy << " rhEta = " << rhEta << " found = " << found << std::endl;
	      totRHMainEvtPerEtaRadius_3MIP[etaBin][ir] += 1;
	      totEnergyMainEvtPerEtaRadius_3MIP[etaBin][ir] += rhPt;

	      totRHMainEvtanyPerEtaRadius_3MIP[etaBin][ir] += 1;
	      totEnergyMainEvtanyPerEtaRadius_3MIP[etaBin][ir] += rhPt;		
	    }
	    // else{
	    //   std::cout << " else 3MIP CPfraction = " << CPfraction << std::endl;
	    //   std::cout << " >>> rhPt = " << rhPt << " rhEnergy = " << rhEnergy << " rhEta = " << rhEta << " found = " << found << std::endl;
	    //   // totRHPUPerEtaRadius_3MIP[etaBin][ir] += 1;
	    //   // totEnergyPUPerEtaRadius_3MIP[etaBin][ir] += rhPt;
	    // }
	  }
	}
	else{
	  totRHPUPerEtaRadius[etaBin][ir] += 1;
	  totEnergyPUPerEtaRadius[etaBin][ir] += rhPt;

	  totRHPUanyPerEtaRadius[etaBin][ir] += 1;
	  totEnergyPUanyPerEtaRadius[etaBin][ir] += rhPt;

	  if(MIP > 3){
	    totRHPUPerEtaRadius_3MIP[etaBin][ir] += 1;
	    totEnergyPUPerEtaRadius_3MIP[etaBin][ir] += rhPt;

	    totRHPUanyPerEtaRadius_3MIP[etaBin][ir] += 1;
	    totEnergyPUanyPerEtaRadius_3MIP[etaBin][ir] += rhPt;
	  }	    
	}
	
      }//loop over radii
    }//radius ok


    ++AllHitsNoEcut;
    if(MIP > 3) ++AllHits3MIPcut;

    //    if(iRadBin == 0)
    totEnergy += rhPt;
    if(MIP > 3) totEnergy_3MIP += rhPt;
    if(found){
      if(CPfraction == 1) totEnergy_MainEvt += rhPt;
      if(CPfraction < 1){
	totEnergy_Shared += rhPt;
	totEnergy_SharedME += rhPt*CPfraction;
	totEnergy_SharedPU += rhPt*(1. - CPfraction);
      }
      if(MIP > 3){
	if(CPfraction == 1) totEnergy_MainEvt_3MIP += rhPt;
	if(CPfraction < 1){
	  totEnergy_Shared_3MIP += rhPt;
	  if(MIP * CPfraction > 3) totEnergy_SharedME_3MIP += rhPt * CPfraction;
	  else if(MIP * (1. - CPfraction) > 3)totEnergy_SharedPU_3MIP += rhPt *(1. - CPfraction);
	}
      }
    }
    else{
      totEnergy_PU += rhPt;
      if(MIP > 3){
	totEnergy_PU_3MIP += rhPt;
      }
    } 
    
  }// loop over rechits

  
  if(debugCOUT)    std::cout << " bau 3 fine loop recHits " << std::endl;
  
  //fill average distribution per eta - rad bin


  if(debugCOUT4){
    std::cout << " >>> " << std::endl;

  }


  h_AllHits_noEcut->Fill(AllHitsNoEcut);
  h_AllHits_3MIPcut->Fill(AllHits3MIPcut);

  h_Energy->Fill(totEnergy);
  h_Energy_PU->Fill(totEnergy_PU);
  h_Energy_MainEvt->Fill(totEnergy_MainEvt);
  h_Energy_Shared->Fill(totEnergy_Shared);
  h_Energy_SharedME->Fill(totEnergy_SharedME);
  h_Energy_SharedPU->Fill(totEnergy_SharedPU);

  h_Energy_3MIP->Fill(totEnergy_3MIP);
  h_Energy_PU_3MIP->Fill(totEnergy_PU_3MIP);
  h_Energy_MainEvt_3MIP->Fill(totEnergy_MainEvt_3MIP);
  h_Energy_Shared_3MIP->Fill(totEnergy_Shared_3MIP);
  h_Energy_SharedME_3MIP->Fill(totEnergy_SharedME_3MIP);
  h_Energy_SharedPU_3MIP->Fill(totEnergy_SharedPU_3MIP);

  if(totEnergy != 0.){
    h_Fraction_Energy_PU->Fill(std::min(1.*totEnergy_PU/totEnergy, 0.999));
    h_Fraction_Energy_MainEvt->Fill(std::min(0.999, 1.*totEnergy_MainEvt/totEnergy));
    h_Fraction_Energy_Shared->Fill(std::min(0.999, 1.*totEnergy_Shared/totEnergy));
    if(totEnergy_Shared != 0.){
      h_Fraction_Energy_SharedME->Fill(std::min(0.999, 1.*totEnergy_SharedME/totEnergy_Shared));
      h_Fraction_Energy_SharedPU->Fill(std::min(0.999, 1.*totEnergy_SharedPU/totEnergy_Shared));
    }
    
    if(totEnergy_3MIP != 0.){
      h_Fraction_Energy_PU_3MIP->Fill(std::min(0.999, 1.*totEnergy_PU_3MIP/totEnergy_3MIP));
      h_Fraction_Energy_MainEvt_3MIP->Fill(std::min(0.999, 1.*totEnergy_MainEvt_3MIP/totEnergy_3MIP));
      h_Fraction_Energy_Shared_3MIP->Fill(std::min(0.999, 1.*totEnergy_Shared_3MIP/totEnergy_3MIP));
      if(totEnergy_Shared_3MIP != 0.){
	h_Fraction_Energy_SharedME_3MIP->Fill(std::min(0.999, 1.*totEnergy_SharedME_3MIP/totEnergy_Shared_3MIP));
	h_Fraction_Energy_SharedPU_3MIP->Fill(std::min(0.999, 1.*totEnergy_SharedPU_3MIP/totEnergy_Shared_3MIP));
      }
    }
  }
  

  if(debugCOUT3){
    std::cout << " totEnergy_PU = " << totEnergy_PU << " totEnergy " << totEnergy << std::endl;
    std::cout << " totEnergy_PU_3MIP = " << totEnergy_PU_3MIP << " totEnergy_3MIP " << totEnergy_3MIP << std::endl;
  }


  for(int iet=0; iet<nBinsEta; ++iet){
    for(int irad=0; irad<nBinsRad; ++irad){
      
      if(totRHPerEtaRadius[iet][irad] == 0.) continue;
      totEvtsEtaRadius[iet][irad] += 1;	

      if(totRHPUPerEtaRadius[iet][irad] > 2.)  totEvtsPUEtaRadius[iet][irad] += 1;	
      if(totRHPUanyPerEtaRadius[iet][irad] > 2.)  totEvtsPUanyEtaRadius[iet][irad] += 1;	
      if(totRHMainEvtPerEtaRadius[iet][irad] > 2.)  totEvtsMainEvtEtaRadius[iet][irad] += 1;	
      if(totRHMainEvtanyPerEtaRadius[iet][irad] > 2.)  totEvtsMainEvtanyEtaRadius[iet][irad] += 1;	
      if(totRHSharedPerEtaRadius[iet][irad] > 2.)  totEvtsSharedEtaRadius[iet][irad] += 1;	
      if(totRHSharedMEPerEtaRadius[iet][irad] > 2.)  totEvtsSharedMEEtaRadius[iet][irad] += 1;	
      if(totRHSharedPUPerEtaRadius[iet][irad] > 2.)  totEvtsSharedPUEtaRadius[iet][irad] += 1;	

      h_NumberHits_Eta_dRadius[iet][irad]->Fill(totRHPerEtaRadius[iet][irad]);
      h_NumberHits_PU_Eta_dRadius[iet][irad]->Fill(totRHPUPerEtaRadius[iet][irad]);
      h_NumberHits_PUany_Eta_dRadius[iet][irad]->Fill(totRHPUanyPerEtaRadius[iet][irad]);
      h_NumberHits_MainEvt_Eta_dRadius[iet][irad]->Fill(totRHMainEvtPerEtaRadius[iet][irad]);
      h_NumberHits_MainEvtany_Eta_dRadius[iet][irad]->Fill(totRHMainEvtanyPerEtaRadius[iet][irad]);
      h_NumberHits_Shared_Eta_dRadius[iet][irad]->Fill(totRHSharedPerEtaRadius[iet][irad]);
      h_NumberHits_SharedME_Eta_dRadius[iet][irad]->Fill(totRHSharedMEPerEtaRadius[iet][irad]);
      h_NumberHits_SharedPU_Eta_dRadius[iet][irad]->Fill(totRHSharedPUPerEtaRadius[iet][irad]);
      h_FractionHits_PU_Eta_dRadius[iet][irad]->Fill(std::min(0.999, 1.* totRHPUPerEtaRadius[iet][irad] / totRHPerEtaRadius[iet][irad]));
      h_FractionHits_PUany_Eta_dRadius[iet][irad]->Fill(std::min(0.999, 1.* totRHPUanyPerEtaRadius[iet][irad] / totRHPerEtaRadius[iet][irad]));
      h_FractionHits_MainEvt_Eta_dRadius[iet][irad]->Fill(std::min(0.999, 1.*totRHMainEvtPerEtaRadius[iet][irad] / totRHPerEtaRadius[iet][irad]));
      h_FractionHits_MainEvtany_Eta_dRadius[iet][irad]->Fill(std::min(0.999, 1.*totRHMainEvtanyPerEtaRadius[iet][irad] / totRHPerEtaRadius[iet][irad]));
      h_FractionHits_Shared_Eta_dRadius[iet][irad]->Fill(std::min(0.999, 1.*totRHSharedPerEtaRadius[iet][irad] / totRHPerEtaRadius[iet][irad]));
      if(totRHSharedPerEtaRadius[iet][irad] != 0.){
	h_FractionHits_SharedME_Eta_dRadius[iet][irad]->Fill(std::min(0.999, 1.*totRHSharedMEPerEtaRadius[iet][irad] / totRHSharedPerEtaRadius[iet][irad]));
	h_FractionHits_SharedPU_Eta_dRadius[iet][irad]->Fill(std::min(0.999, 1.*totRHSharedPUPerEtaRadius[iet][irad] / totRHSharedPerEtaRadius[iet][irad]));
      }

      h_Energy_Eta_dRadius[iet][irad]->Fill(1.*totEnergyPerEtaRadius[iet][irad]);
      h_Energy_PU_Eta_dRadius[iet][irad]->Fill(1.*totEnergyPUPerEtaRadius[iet][irad]);
      h_Energy_PUany_Eta_dRadius[iet][irad]->Fill(1.*totEnergyPUanyPerEtaRadius[iet][irad]);
      h_Energy_MainEvt_Eta_dRadius[iet][irad]->Fill(1.*totEnergyMainEvtPerEtaRadius[iet][irad]);
      h_Energy_Shared_Eta_dRadius[iet][irad]->Fill(1.*totEnergySharedPerEtaRadius[iet][irad]);
      h_Energy_SharedME_Eta_dRadius[iet][irad]->Fill(1.*totEnergySharedMEPerEtaRadius[iet][irad]);
      h_Energy_SharedPU_Eta_dRadius[iet][irad]->Fill(1.*totEnergySharedPUPerEtaRadius[iet][irad]);      
      h_FractionEnergy_PU_Eta_dRadius[iet][irad]->Fill(std::min(0.999, 1.*totEnergyPUPerEtaRadius[iet][irad] / totEnergyPerEtaRadius[iet][irad]));
      h_FractionEnergy_PUany_Eta_dRadius[iet][irad]->Fill(std::min(0.999, 1.*totEnergyPUanyPerEtaRadius[iet][irad] / totEnergyPerEtaRadius[iet][irad]));
      h_FractionEnergy_MainEvt_Eta_dRadius[iet][irad]->Fill(std::min(0.999, 1.*totEnergyMainEvtPerEtaRadius[iet][irad] / totEnergyPerEtaRadius[iet][irad]));
      h_FractionEnergy_MainEvtany_Eta_dRadius[iet][irad]->Fill(std::min(0.999, 1.*totEnergyMainEvtanyPerEtaRadius[iet][irad] / totEnergyPerEtaRadius[iet][irad]));
      h_FractionEnergy_Shared_Eta_dRadius[iet][irad]->Fill(std::min(0.999, 1.*totEnergySharedPerEtaRadius[iet][irad] / totEnergyPerEtaRadius[iet][irad]));
      if(totEnergySharedPerEtaRadius[iet][irad] != 0.){
	h_FractionEnergy_SharedME_Eta_dRadius[iet][irad]->Fill(std::min(0.999, 1.*totEnergySharedMEPerEtaRadius[iet][irad] / totEnergySharedPerEtaRadius[iet][irad]));
	h_FractionEnergy_SharedPU_Eta_dRadius[iet][irad]->Fill(std::min(0.999, 1.*totEnergySharedPUPerEtaRadius[iet][irad] / totEnergySharedPerEtaRadius[iet][irad]));
      }
      
      if(totRHPerEtaRadius_3MIP[iet][irad] == 0.) continue;
      totEvtsEtaRadius_3MIP[iet][irad] += 1;	
      
      if(totRHPUPerEtaRadius_3MIP[iet][irad] > 2.)  totEvtsPUEtaRadius_3MIP[iet][irad] += 1;	
      if(totRHPUanyPerEtaRadius_3MIP[iet][irad] > 2.)  totEvtsPUanyEtaRadius_3MIP[iet][irad] += 1;	
      if(totRHMainEvtPerEtaRadius_3MIP[iet][irad] > 2.)  totEvtsMainEvtEtaRadius_3MIP[iet][irad] += 1;	
      if(totRHMainEvtanyPerEtaRadius_3MIP[iet][irad] > 2.)  totEvtsMainEvtanyEtaRadius_3MIP[iet][irad] += 1;	
      if(totRHSharedPerEtaRadius_3MIP[iet][irad] > 2.)  totEvtsSharedEtaRadius_3MIP[iet][irad] += 1;	
      if(totRHSharedMEPerEtaRadius_3MIP[iet][irad] > 2.)  totEvtsSharedMEEtaRadius_3MIP[iet][irad] += 1;	
      if(totRHSharedPUPerEtaRadius_3MIP[iet][irad] > 2.)  totEvtsSharedPUEtaRadius_3MIP[iet][irad] += 1;	
      
      h_NumberHits_Eta_dRadius_3MIP[iet][irad]->Fill(totRHPerEtaRadius_3MIP[iet][irad]);
      h_NumberHits_PU_Eta_dRadius_3MIP[iet][irad]->Fill(totRHPUPerEtaRadius_3MIP[iet][irad]);
      h_NumberHits_PUany_Eta_dRadius_3MIP[iet][irad]->Fill(totRHPUanyPerEtaRadius_3MIP[iet][irad]);
      h_NumberHits_MainEvt_Eta_dRadius_3MIP[iet][irad]->Fill(totRHMainEvtPerEtaRadius_3MIP[iet][irad]);
      h_NumberHits_MainEvtany_Eta_dRadius_3MIP[iet][irad]->Fill(totRHMainEvtanyPerEtaRadius_3MIP[iet][irad]);
      h_NumberHits_Shared_Eta_dRadius_3MIP[iet][irad]->Fill(totRHSharedPerEtaRadius_3MIP[iet][irad]);
      h_NumberHits_SharedME_Eta_dRadius_3MIP[iet][irad]->Fill(totRHSharedMEPerEtaRadius_3MIP[iet][irad]);
      h_NumberHits_SharedPU_Eta_dRadius_3MIP[iet][irad]->Fill(totRHSharedPUPerEtaRadius_3MIP[iet][irad]);
      h_FractionHits_PU_Eta_dRadius_3MIP[iet][irad]->Fill(std::min(0.999, 1.* totRHPUPerEtaRadius_3MIP[iet][irad] / totRHPerEtaRadius_3MIP[iet][irad]));
      h_FractionHits_PUany_Eta_dRadius_3MIP[iet][irad]->Fill(std::min(0.999, 1.* totRHPUanyPerEtaRadius_3MIP[iet][irad] / totRHPerEtaRadius_3MIP[iet][irad]));
      h_FractionHits_MainEvt_Eta_dRadius_3MIP[iet][irad]->Fill(std::min(0.999, 1.*totRHMainEvtPerEtaRadius_3MIP[iet][irad] / totRHPerEtaRadius_3MIP[iet][irad]));
      h_FractionHits_MainEvtany_Eta_dRadius_3MIP[iet][irad]->Fill(std::min(0.999, 1.*totRHMainEvtanyPerEtaRadius_3MIP[iet][irad] / totRHPerEtaRadius_3MIP[iet][irad]));
      h_FractionHits_Shared_Eta_dRadius_3MIP[iet][irad]->Fill(std::min(0.999, 1.*totRHSharedPerEtaRadius_3MIP[iet][irad] / totRHPerEtaRadius_3MIP[iet][irad]));
      if(totRHSharedPerEtaRadius_3MIP[iet][irad] != 0.){
	h_FractionHits_SharedME_Eta_dRadius_3MIP[iet][irad]->Fill(std::min(0.999, 1.*totRHSharedMEPerEtaRadius_3MIP[iet][irad] / totRHSharedPerEtaRadius_3MIP[iet][irad]));
	h_FractionHits_SharedPU_Eta_dRadius_3MIP[iet][irad]->Fill(std::min(0.999, 1.*totRHSharedPUPerEtaRadius_3MIP[iet][irad] / totRHSharedPerEtaRadius_3MIP[iet][irad]));
      }
      
      h_Energy_Eta_dRadius_3MIP[iet][irad]->Fill(1.*totEnergyPerEtaRadius_3MIP[iet][irad]);
      h_Energy_PU_Eta_dRadius_3MIP[iet][irad]->Fill(1.*totEnergyPUPerEtaRadius_3MIP[iet][irad]);
      h_Energy_PUany_Eta_dRadius_3MIP[iet][irad]->Fill(1.*totEnergyPUanyPerEtaRadius_3MIP[iet][irad]);
      h_Energy_MainEvt_Eta_dRadius_3MIP[iet][irad]->Fill(1.*totEnergyMainEvtPerEtaRadius_3MIP[iet][irad]);
      h_Energy_MainEvtany_Eta_dRadius_3MIP[iet][irad]->Fill(1.*totEnergyMainEvtanyPerEtaRadius_3MIP[iet][irad]);
      h_Energy_Shared_Eta_dRadius_3MIP[iet][irad]->Fill(1.*totEnergySharedPerEtaRadius_3MIP[iet][irad]);
      h_Energy_SharedME_Eta_dRadius_3MIP[iet][irad]->Fill(1.*totEnergySharedMEPerEtaRadius_3MIP[iet][irad]);
      h_Energy_SharedPU_Eta_dRadius_3MIP[iet][irad]->Fill(1.*totEnergySharedPUPerEtaRadius_3MIP[iet][irad]);      
      h_FractionEnergy_PU_Eta_dRadius_3MIP[iet][irad]->Fill(std::min(0.999, 1.*totEnergyPUPerEtaRadius_3MIP[iet][irad] / totEnergyPerEtaRadius_3MIP[iet][irad]));
      h_FractionEnergy_PUany_Eta_dRadius_3MIP[iet][irad]->Fill(std::min(0.999, 1.*totEnergyPUanyPerEtaRadius_3MIP[iet][irad] / totEnergyPerEtaRadius_3MIP[iet][irad]));
      h_FractionEnergy_MainEvt_Eta_dRadius_3MIP[iet][irad]->Fill(std::min(0.999, 1.*totEnergyMainEvtPerEtaRadius_3MIP[iet][irad] / totEnergyPerEtaRadius_3MIP[iet][irad]));
      h_FractionEnergy_MainEvtany_Eta_dRadius_3MIP[iet][irad]->Fill(std::min(0.999, 1.*totEnergyMainEvtanyPerEtaRadius_3MIP[iet][irad] / totEnergyPerEtaRadius_3MIP[iet][irad]));
      h_FractionEnergy_Shared_Eta_dRadius_3MIP[iet][irad]->Fill(std::min(0.999, 1.*totEnergySharedPerEtaRadius_3MIP[iet][irad] / totEnergyPerEtaRadius_3MIP[iet][irad]));
   if(totEnergySharedPerEtaRadius_3MIP[iet][irad] != 0.){
   h_FractionEnergy_SharedME_Eta_dRadius_3MIP[iet][irad]->Fill(std::min(0.999, 1.*totEnergySharedMEPerEtaRadius_3MIP[iet][irad] / totEnergySharedPerEtaRadius_3MIP[iet][irad]));
   h_FractionEnergy_SharedPU_Eta_dRadius_3MIP[iet][irad]->Fill(std::min(0.999, 1.*totEnergySharedPUPerEtaRadius_3MIP[iet][irad] / totEnergySharedPerEtaRadius_3MIP[iet][irad]));
   }
    }
  }
  
  if(debugCOUT) std::cout<< " bau fine evento " << std::endl;
}

void
HGCalTimingAnalyzerPU::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
HGCalTimingAnalyzerPU::endJob()
{
  //  if(debugCOUT) std::cout<< " bau 7 " << std::endl;
  for(int iet=0; iet<nBinsEta; ++iet){
    for(int irad=0; irad<nBinsRad; ++irad){
      if(totEvtsEtaRadius[iet][irad] != 0){
	hFractionEvents_PU_Eta_dRadius[iet][irad]->Fill(std::min(0.999, 1.*totEvtsPUEtaRadius[iet][irad]/totEvtsEtaRadius[iet][irad]));
	hFractionEvents_PUany_Eta_dRadius[iet][irad]->Fill(std::min(0.999, 1.*totEvtsPUanyEtaRadius[iet][irad]/totEvtsEtaRadius[iet][irad]));
	hFractionEvents_MainEvt_Eta_dRadius[iet][irad]->Fill(std::min(0.999, 1.*totEvtsMainEvtEtaRadius[iet][irad]/totEvtsEtaRadius[iet][irad]));
	hFractionEvents_MainEvtany_Eta_dRadius[iet][irad]->Fill(std::min(0.999, 1.*totEvtsMainEvtanyEtaRadius[iet][irad]/totEvtsEtaRadius[iet][irad]));
	hFractionEvents_Shared_Eta_dRadius[iet][irad]->Fill(std::min(0.999, 1.*totEvtsSharedEtaRadius[iet][irad]/totEvtsEtaRadius[iet][irad]));
	if(totEvtsSharedEtaRadius[iet][irad] != 0){
	  hFractionEvents_SharedME_Eta_dRadius[iet][irad]->Fill(std::min(0.999, 1.*totEvtsSharedMEEtaRadius[iet][irad]/totEvtsSharedEtaRadius[iet][irad]));
	  hFractionEvents_SharedPU_Eta_dRadius[iet][irad]->Fill(std::min(0.999, 1.*totEvtsSharedPUEtaRadius[iet][irad]/totEvtsSharedEtaRadius[iet][irad]));
	}
      }
      if(totEvtsEtaRadius_3MIP[iet][irad] != 0){
	hFractionEvents_PU_Eta_dRadius_3MIP[iet][irad]->Fill(std::min(0.999, 1.*totEvtsPUEtaRadius_3MIP[iet][irad]/totEvtsEtaRadius_3MIP[iet][irad]));
	hFractionEvents_PUany_Eta_dRadius_3MIP[iet][irad]->Fill(std::min(0.999, 1.*totEvtsPUanyEtaRadius_3MIP[iet][irad]/totEvtsEtaRadius_3MIP[iet][irad]));
	hFractionEvents_MainEvt_Eta_dRadius_3MIP[iet][irad]->Fill(std::min(0.999, 1.*totEvtsMainEvtEtaRadius_3MIP[iet][irad]/totEvtsEtaRadius_3MIP[iet][irad]));
	hFractionEvents_MainEvtany_Eta_dRadius_3MIP[iet][irad]->Fill(std::min(0.999, 1.*totEvtsMainEvtanyEtaRadius_3MIP[iet][irad]/totEvtsEtaRadius_3MIP[iet][irad]));
	hFractionEvents_Shared_Eta_dRadius_3MIP[iet][irad]->Fill(std::min(0.999, 1.*totEvtsSharedEtaRadius_3MIP[iet][irad]/totEvtsEtaRadius_3MIP[iet][irad]));
	if(totEvtsSharedEtaRadius_3MIP[iet][irad] != 0){
	  hFractionEvents_SharedME_Eta_dRadius_3MIP[iet][irad]->Fill(std::min(0.999, 1.*totEvtsSharedMEEtaRadius_3MIP[iet][irad]/totEvtsSharedEtaRadius_3MIP[iet][irad]));
	  hFractionEvents_SharedPU_Eta_dRadius_3MIP[iet][irad]->Fill(std::min(0.999, 1.*totEvtsSharedPUEtaRadius_3MIP[iet][irad]/totEvtsSharedEtaRadius_3MIP[iet][irad]));
	}
      }
    }
  }
  
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HGCalTimingAnalyzerPU::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HGCalTimingAnalyzerPU);
