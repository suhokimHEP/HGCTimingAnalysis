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
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"


#include "RecoLocalCalo/HGCalRecAlgos/interface/ClusterTools.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/HGCalDepthPreClusterer.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include "HGCTimingAnalysis/HGCTiming/plugins/RecHiTimeEstimator.h"
#include "HGCTimingAnalysis/HGCTiming/interface/UtilClasses.h"



#include <string>
#include <map>

class HGCalTimingAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
//
// constructors and destructor
//
  explicit HGCalTimingAnalyzer(const edm::ParameterSet&);
  ~HGCalTimingAnalyzer();

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
  TH1F* h_allRH_TimesOfRadius;

  TH2F* h2_dPhivsdEta_rhGen;
  TH2F* h2_dPhivsdEta_rhAxis;
  TH2F* h2_dPhivsdEta_GenAxis;

  //eta bins 1.65 - 1.85   1.85-2.05   2.05-2.25   2.25-2.45   2.45-2.65   2.65-2.85
  //eta bins 1.5-1.7  1.7-1.9   1.9-2.1 - 2.1-2.3   2.3-2.5 2.5-2.7   2.7-2.9   2.9-3.1
  //  eta bins 1.5-2.  2.-2.5  2.5-3.

  TH1F* h_rhGen_Radius_Eta[6];
  TH1F* h_rhGen_RadiusDdeta_Eta[6];
  TH1F* h_rhGen_RadiusPdeta_Eta[6];
  TH1F* hAverageTime_Eta_dRadius[6][4];
  TH1F* hAverageTime_Eta_dRadius_AvgArm[6][4];
  TH1F* hAverageTime_Eta_dRadius_ResoWe[6][4];
  TH1F* hAverageTime_Eta_dRadius_AvgCutH[6][2];
  TH1F* hAverageTime_Eta_dRadius_Avg68[6][2];
  TH1F* hTime_Eta_dRadius[6][4];
  TH1F* hTimeCut_Eta_dRadius[6][4];
  TH1F* hEnergyWithTime_Eta_dRadius[6][4];
  TH1F* hNumberHitsWithTime_Eta_dRadius[6][4];

  TH1F* hEtaDistr_WithTime;
  TH1F* hEtaDistr_AllTime;
  TH1F* hEtaDistr_Fraction;

  ///////
  TH1F* hEtaDistrEvt3Hit;
  TH1F* hEtaDistrEvt;
  TH1F* hEtaDistrEvt3HitFraction;


  TH1F* energyChargeFraction;
  TH1F* energyCharge;
  TH1F* energyChargeAll;


  TH1F* hTotHits_Eta_dRadius[6][4];
  TH1F* hTotHitsWithTime_Eta_dRadius[6][4];
  TH1F* hFractionHitsWithTime_Eta_dRadius[6][4];

  TH1F* hFractionEvents_HitsWithTime_Eta_dRadius[6][4];

  TH1F* h_totEvtsEtaRadius_withTime_Eta_dRadius[6][4];
  TH1F* h_totEvtsEtaRadius_Eta_dRadius[6][4];


  int totEvtsEtaRadius[6][4];
  int totEvtsEtaRadius_withTime[6][4];


  float radiusEtaRad[6][4];

  int nBinsEta;
  int nBinsRad;
  float binWidth;
  float binStart;
  float binEnd;

  bool debugCOUT;
  bool debugCOUT2;


  int cellType;
  float floorValue;
  int lifeAge;
  float absTrend;

};  

float HGCalTimingAnalyzer::getXmax(TH1F* histo, float& YMax){

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

  //  std::cout << histo->GetName() << " " << histo->GetBinCenter(xBin) << std::endl; 
  return histo->GetBinCenter(xBin);
}


bool HGCalTimingAnalyzer::comparePairs(const std::pair<float, float>& i, const std::pair<float, float>& j){
  return i.first < j.first;
}



HGCalTimingAnalyzer::HGCalTimingAnalyzer(const edm::ParameterSet& iConfig) :
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

  h_Vtx_x = fs->make<TH1F>("h_Vtx_x", "", 1000, -10., 10.);
  h_Vtx_y = fs->make<TH1F>("h_Vtx_y", "", 1000, -10., 10.);
  h_Vtx_z = fs->make<TH1F>("h_Vtx_z", "", 1000, -10., 10.);

  h_Vtx_dvx = fs->make<TH1F>("h_Vtx_dvx", "", 1000, -10., 10.);
  h_Vtx_dvy = fs->make<TH1F>("h_Vtx_dvy", "", 1000, -10., 10.);
  h_Vtx_dvz = fs->make<TH1F>("h_Vtx_dvz", "", 1000, -10., 10.);

  h2_dPhivsdEta_rhGen = fs->make<TH2F>("h2_dPhivsdEta_rhGen", "", 500, -0.2, 0.2, 500, -0.2, 0.2);
  h2_dPhivsdEta_rhAxis = fs->make<TH2F>("h2_dPhivsdEta_rhAxis", "", 1000, -2., 2., 1000, -2., 2.);
  h2_dPhivsdEta_GenAxis = fs->make<TH2F>("h2_dPhivsdEta_GenAxis", "", 1000, -2., 2., 1000, -2., 2.);


  h_allRH_TimesOfRadius = fs->make<TH1F>("h_allRH_TimesOfRadius", "", 100., 0., 50.);


  hEtaDistr_WithTime = fs->make<TH1F>("hEtaDistr_WithTime", "", nBinsEta, binStart, binStart+nBinsEta*binWidth);
  hEtaDistr_AllTime = fs->make<TH1F>("hEtaDistr_AllTime", "", nBinsEta, binStart, binStart+nBinsEta*binWidth);
  hEtaDistr_Fraction = fs->make<TH1F>("hEtaDistr_Fraction", "", nBinsEta, binStart, binStart+nBinsEta*binWidth);


  hEtaDistrEvt3Hit = fs->make<TH1F>("hEtaDistrEvt3Hit", "", nBinsEta, binStart, binStart+nBinsEta*binWidth);
  hEtaDistrEvt = fs->make<TH1F>("hEtaDistrEvt", "", nBinsEta, binStart, binStart+nBinsEta*binWidth);
  hEtaDistrEvt3HitFraction = fs->make<TH1F>("hEtaDistrEvt3HitFraction", "", nBinsEta, binStart, binStart+nBinsEta*binWidth);
  energyChargeFraction = fs->make<TH1F>("energyChargeFraction", "", nBinsEta, binStart, binStart+nBinsEta*binWidth);
  energyCharge = fs->make<TH1F>("energyCharge", "", nBinsEta, binStart, binStart+nBinsEta*binWidth);
  energyChargeAll = fs->make<TH1F>("energyChargeAll", "", nBinsEta, binStart, binStart+nBinsEta*binWidth);

  for(int ieta=0; ieta<nBinsEta; ++ieta){    
    if(debugCOUT) std::cout<< " ieta from = " << (binStart+ieta*binWidth) << " to " << binStart+binWidth+ieta*binWidth << std::endl;
    h_rhGen_Radius_Eta[ieta] = fs->make<TH1F>(Form("h_rhGen_Radius_Eta_%.2f-%.2f", (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth), "", 1000, 0., 100.); 

    for(int iRad=0; iRad<nBinsRad; ++iRad){
      hTime_Eta_dRadius[ieta][iRad] = fs->make<TH1F>(Form("hTime_Eta_dRadius_Eta%.2f-%.2f_dRadius%d", (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000, -2., 48.);
      hTimeCut_Eta_dRadius[ieta][iRad] = fs->make<TH1F>(Form("hTimeCut_Eta_dRadius_Eta%.2f-%.2f_dRadius%d", (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000, -2., 48.);
      hEnergyWithTime_Eta_dRadius[ieta][iRad] = fs->make<TH1F>(Form("hEnergyWithTime_Eta_dRadius_Eta%.2f-%.2f_dRadius%d", (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000, 0., 1000.);
      hNumberHitsWithTime_Eta_dRadius[ieta][iRad] = fs->make<TH1F>(Form("hNumberHitsWithTime_Eta_dRadius_Eta%.2f-%.2f_dRadius%d", (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 2000, 0., 2000.);
      hAverageTime_Eta_dRadius[ieta][iRad] = fs->make<TH1F>(Form("hAverageTime_Eta%.2f-%.2f_dRadius%d", (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1500, -0.5, 1.);
      hAverageTime_Eta_dRadius_AvgArm[ieta][iRad] = fs->make<TH1F>(Form("hAverageTime_Eta%.2f-%.2f_dRadius%d_AvgArm", (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000, -0.2, 1.);
      hAverageTime_Eta_dRadius_ResoWe[ieta][iRad] = fs->make<TH1F>(Form("hAverageTime_Eta%.2f-%.2f_dRadius%d_ResoWe", (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1500, -0.5, 1.5);
      if(iRad < 2){
      hAverageTime_Eta_dRadius_AvgCutH[ieta][iRad] = fs->make<TH1F>(Form("hAverageTime_Eta%.2f-%.2f_dRadius%d_AvgCutH", (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000, -0.2, 1.);      
      hAverageTime_Eta_dRadius_Avg68[ieta][iRad] = fs->make<TH1F>(Form("hAverageTime_Eta%.2f-%.2f_dRadius%d_Avg68", (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000, -0.2, 1.);
      }

      hFractionHitsWithTime_Eta_dRadius[ieta][iRad] = fs->make<TH1F>(Form("hFractionHits_Eta%.2f-%.2f_dRadius%d", (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1100, 0, 1.1);

      hFractionEvents_HitsWithTime_Eta_dRadius[ieta][iRad] = fs->make<TH1F>(Form("hFractionEvents_Hits_Eta%.2f-%.2f_dRadius%d", (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1100, 0, 1.1);

      h_totEvtsEtaRadius_withTime_Eta_dRadius[ieta][iRad] =  fs->make<TH1F>(Form("h_totEvtsEtaRadius_withTime_Eta%.2f-%.2f_dRadius%d", (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1, 0., 1.);
      h_totEvtsEtaRadius_Eta_dRadius[ieta][iRad] =  fs->make<TH1F>(Form("h_totEvtsEtaRadius_Eta%.2f-%.2f_dRadius%d", (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1, 0., 1.);


      hTotHits_Eta_dRadius[ieta][iRad] = fs->make<TH1F>(Form("hTotHits_Eta_%.2f-%.2f_dRadius%d", (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000, 0, 100.);
      hTotHitsWithTime_Eta_dRadius[ieta][iRad] = fs->make<TH1F>(Form("hTotHitsWithTime_Eta_%.2f-%.2f_dRadius%d", (binStart+ieta*binWidth), binStart+binWidth+ieta*binWidth, iRad), "", 1000, 0, 1000.);
            
    }
  }


  for(int iEta=0; iEta<nBinsEta; ++iEta){
    for(int iRad=0; iRad<nBinsRad; ++iRad){
	totEvtsEtaRadius[iEta][iRad] = 0;
	totEvtsEtaRadius_withTime[iEta][iRad] = 0;
      }
    radiusEtaRad[iEta][0] = 2.;
    radiusEtaRad[iEta][1] = 5.;
    radiusEtaRad[iEta][2] = 10.;
    radiusEtaRad[iEta][3] = 600.;
  }

}

HGCalTimingAnalyzer::~HGCalTimingAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


void
HGCalTimingAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
  //  HGCRecHitCollection NewrechitsBH;

  //make a map detid-rechit
  std::map<DetId,const HGCRecHit*> hitmap;
  switch(algo){
  case 1:
    {
      iEvent.getByToken(_recHitsEE,recHitHandleEE);
      iEvent.getByToken(_recHitsFH,recHitHandleFH);
      iEvent.getByToken(_recHitsBH,recHitHandleBH);

      const HGCRecHitCollection& rechitsEEOld = *recHitHandleEE;
      const HGCRecHitCollection& rechitsFHOld = *recHitHandleFH;
      //      const HGCRecHitCollection& rechitsBHOld = *recHitHandleBH;
      if(CFDTimeCorrection == 1){
	myEstimator->correctTime(rechitsEEOld, &NewrechitsEE);
	myEstimator->correctTime(rechitsFHOld, &NewrechitsFH);
	//	myEstimator->correctTime(rechitsBHOld, &NewrechitsBH);
      }
      else if(CFDTimeCorrection == 0){
	myEstimator->correctTimeFixThr(rechitsEEOld, &NewrechitsEE);
	myEstimator->correctTimeFixThr(rechitsFHOld, &NewrechitsFH);
	//	myEstimator->correctTimeFixThr(rechitsBHOld, &NewrechitsBH);
      }
      else if(CFDTimeCorrection == -1){
	NewrechitsEE = *recHitHandleEE;
	NewrechitsFH = *recHitHandleFH;
	//	NewrechitsBH = *recHitHandleBH;
      }
      for(unsigned int i = 0; i < NewrechitsEE.size(); ++i){
	hitmap[NewrechitsEE[i].detid()] = &NewrechitsEE[i];
      }
      for(unsigned int i = 0; i < NewrechitsFH.size(); ++i){
	hitmap[NewrechitsFH[i].detid()] = &NewrechitsFH[i];
      }
      // for(unsigned int i = 0; i < NewrechitsBH.size(); ++i){
      // 	hitmap[NewrechitsBH[i].detid()] = &NewrechitsBH[i];
      // }
      break;
    }
  case 2:
    {
      iEvent.getByToken(_recHitsEE,recHitHandleEE);

      const HGCRecHitCollection& rechitsEEOld = *recHitHandleEE;

      if(CFDTimeCorrection == 1) myEstimator->correctTime(rechitsEEOld, &NewrechitsEE);
      else if(CFDTimeCorrection == 0) myEstimator->correctTimeFixThr(rechitsEEOld, &NewrechitsEE);
      else if(CFDTimeCorrection == -1){
	NewrechitsEE = *recHitHandleEE;
      }

      for(unsigned int i = 0; i < NewrechitsEE.size(); i++){
	hitmap[NewrechitsEE[i].detid()] = &NewrechitsEE[i];
      }
      break;
    }
  case 3:
    {
      iEvent.getByToken(_recHitsFH,recHitHandleFH);
      iEvent.getByToken(_recHitsBH,recHitHandleBH);

      const HGCRecHitCollection& rechitsFHOld = *recHitHandleFH;
      const HGCRecHitCollection& rechitsBHOld = *recHitHandleBH;
      if(CFDTimeCorrection == 1){
	myEstimator->correctTime(rechitsFHOld, &NewrechitsFH);
	//	myEstimator->correctTime(rechitsBHOld, &NewrechitsBH);
      }
      else if(CFDTimeCorrection == 0){
	myEstimator->correctTimeFixThr(rechitsFHOld, &NewrechitsFH);
	//	myEstimator->correctTimeFixThr(rechitsBHOld, &NewrechitsBH);
      }
      else if(CFDTimeCorrection == -1){
	NewrechitsFH = *recHitHandleFH;
	//	NewrechitsBH = *recHitHandleBH;
      }

      for(unsigned int i = 0; i < NewrechitsFH.size(); i++){
	hitmap[NewrechitsFH[i].detid()] = &NewrechitsFH[i];
      }
      // for(unsigned int i = 0; i < NewrechitsBH.size(); i++){
      // 	hitmap[NewrechitsBH[i].detid()] = &NewrechitsBH[i];
      // }
      break;
    }
  default:
    break;
  }

  ////////////////////

  if(debugCOUT) std::cout<< " >>> now caloparticles " << std::endl;

  std::array<double,3> vtx{ {vx, vy, vz}};
  int numbCaloPart = 0;
  // loop over caloParticles
  for (std::vector<CaloParticle>::const_iterator it_caloPart = caloParticles.begin(); it_caloPart != caloParticles.end(); ++it_caloPart){
    const SimClusterRefVector simClusterRefVector = it_caloPart->simClusters();

    //    std::cout << " simClusterRefVector.size() = " << simClusterRefVector.size() << std::endl; 
    if(simClusterRefVector.size() > 1) continue;

    float etaGen = it_caloPart->eta();
    float phiGen = it_caloPart->phi();
    float xGen = it_caloPart->momentum().x();
    float yGen = it_caloPart->momentum().y();
    float zGen = it_caloPart->momentum().z();

    if(debugCOUT) std::cout<< " bau caloparticles inloop " << std::endl;

    if(debugCOUT2)    std::cout<< " >>> befor caloparticles eta = " << etaGen << std::endl;
    if(etaGen < 0) continue;
    if(etaGen < binStart || etaGen > binEnd) continue;
    ++numbCaloPart;
    if(numbCaloPart > 1) continue;

    if(debugCOUT) std::cout<< " bau caloparticles survived " << std::endl;

    
    if(debugCOUT2)    std::cout<< " >>> now caloparticles " << std::endl;


    //shower axis by recHits
    float axisX = 0;
    float axisY = 0;
    float axisZ = 0;
    float sumEnergyToNorm = 0;
    GlobalPoint showerAxis;
    
    if(debugCOUT) std::cout<< " bau before showerAxis " << std::endl;
  
    //loop on rechit - matched to gen => shower axis
    for (CaloParticle::sc_iterator it_sc = simClusterRefVector.begin(); it_sc != simClusterRefVector.end(); ++it_sc) {
      const SimCluster simCluster = (*(*it_sc));
      const std::vector<std::pair<uint32_t,float> > hits_and_fractions = simCluster.hits_and_fractions();      
      for (std::vector<std::pair<uint32_t,float> >::const_iterator it_haf = hits_and_fractions.begin(); it_haf != hits_and_fractions.end(); ++it_haf) {
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
    }
    showerAxis = GlobalPoint(axisX/sumEnergyToNorm, axisY/sumEnergyToNorm, (axisZ - vz)/sumEnergyToNorm);
  
    float axEta = showerAxis.eta();
    float axPhi = showerAxis.phi();
    float axX = showerAxis.x();
    float axY = showerAxis.y();
    float axZ = showerAxis.z();
 

    std::vector<std::pair<float, float> > timePerEtaRadiusDistr[6][2];
    TH1F* timePerEtaRadiusHisto[6][2];
    for(int iet=0; iet<nBinsEta; ++iet){
      for(int irad=0; irad<2; ++irad){
	timePerEtaRadiusDistr[iet][irad].clear();
	timePerEtaRadiusHisto[iet][irad] = new TH1F(Form("timePerEtaRadiusHisto_eta%d_iR%d", iet, irad), "", 500, -0.2, 1.);
      }
    }
  
    if(debugCOUT) std::cout<< " bau after showerAxis " << std::endl;
    /////////////////////////////////////////////////////////////////
    UtilClasses utilsMet = UtilClasses(etaGen, phiGen);
    std::array<double,3> fromAxis{ {xGen, yGen, zGen} };
    //UtilClasses utilsMet = UtilClasses(axEta, axPhi);
    //std::array<double,3> fromAxis{ {axX, axY, axZ} };

    int totHitsAbove60 = 0;
    int totHitsAbove0 = 0;

    bool withTime_filled = false;
    bool allTime_filled = false;
    int totHitTime = 0.;

    if(debugCOUT) std::cout<< " now unmatched " << std::endl;

    if(debugCOUT) std::cout<< " bau 1 " << std::endl;

    float allRHSeedEnergy = 0;   
    float totEnergyPerEtaRadius[6][4];
    float totEnergyWithTimePerEtaRadius[6][4];
    float timePerEtaRadius[6][4];
    int totRHPerEtaRadius[6][4];
    float timePerEtaRadius_AvgArm[6][4];
    int totRHPerEtaRadius_AvgArm[6][4];
    float timePerEtaRadius_ResoWe[6][4];
    float totRHPerEtaRadius_ResoWe[6][4];
    float timePerEtaRadiusAvgCutH[6][2];
    int totRHPerEtaRadiusAvgCutH[6][2];
    float timePerEtaRadiusAvg68[6][2];
    int totRHPerEtaRadiusAvg68[6][2];

    int totRHPerEtaRadius_allTime[6][4];

    for(int iet=0; iet<nBinsEta; ++iet){
      for(int irad=0; irad<nBinsRad; ++irad){
	totEnergyPerEtaRadius[iet][irad] = 0.;
	totEnergyWithTimePerEtaRadius[iet][irad] = 0.;

	timePerEtaRadius[iet][irad] = 0.;
	totRHPerEtaRadius[iet][irad] = 0;
	timePerEtaRadius_AvgArm[iet][irad] = 0.;
	totRHPerEtaRadius_AvgArm[iet][irad] = 0;
	timePerEtaRadius_ResoWe[iet][irad] = 0.;
	totRHPerEtaRadius_ResoWe[iet][irad] = 0.;

	totRHPerEtaRadius_allTime[iet][irad] = 0;
	if(irad > 1) continue;
	timePerEtaRadiusAvgCutH[iet][irad] = 0.;
	timePerEtaRadiusAvg68[iet][irad] = 0.;
	totRHPerEtaRadiusAvgCutH[iet][irad] = 0;
	totRHPerEtaRadiusAvg68[iet][irad] = 0;
      }
    }
  
    if(debugCOUT) std::cout<< " bau 2 " << std::endl;   
    for (CaloParticle::sc_iterator it_sc = simClusterRefVector.begin(); it_sc != simClusterRefVector.end(); ++it_sc) {
      const SimCluster simCluster = (*(*it_sc));
      const std::vector<std::pair<uint32_t,float> > hits_and_fractions = simCluster.hits_and_fractions();     

      //loop over hits
      for (std::vector<std::pair<uint32_t,float> >::const_iterator it_haf = hits_and_fractions.begin(); it_haf != hits_and_fractions.end(); ++it_haf) {
	DetId hitid = (it_haf->first);
	//	if(recHitTools.getEta(hitid)*it_caloPart->eta() < 0) continue;

	bool found = false;
	float rhTime = -1;
	float rhEnergy = 0.;
	float rhPt = 0.;

	std::map<DetId, const HGCRecHit*>::iterator trovatore = hitmap.find(hitid);
	if(trovatore == hitmap.end()){
	  continue;
	}
	else if(recHitTools.getEta(hitid)*it_caloPart->eta() < 0) continue;
	else{
	  const HGCRecHit *hit = hitmap[hitid];
	  rhTime = hit->time() - 1.;
	  rhEnergy = hit->energy();
	  found = true;
	  // if(recHitTools.getEta(hitid)*it_caloPart->eta() < 0){
	  //   if(debugCOUT) std::cout<< " >>> MAJOR PROBLEM!!! " << std::endl;

	  // if(it_caloPart->eta() > 0){
	  //   h_dEta_pos->Fill(recHitTools.getEta(hitid) - it_caloPart->eta());
	  //   h_dPhi_pos->Fill(reco::deltaPhi(recHitTools.getPhi(hitid), it_caloPart->phi()) );
	  // }
	  // if(it_caloPart->eta() < 0){
	  //   h_dEta_neg->Fill(recHitTools.getEta(hitid) - it_caloPart->eta());
	  //   h_dPhi_neg->Fill(reco::deltaPhi(recHitTools.getPhi(hitid), it_caloPart->phi()) );
	  // }
	  if(hit->energy()*it_haf->second > allRHSeedEnergy){
	    allRHSeedEnergy = hit->energy()*it_haf->second;
	  }

	  if(allTime_filled == false){
	    hEtaDistr_AllTime->Fill(std::abs(etaGen));
	    allTime_filled = true;
	  }
	  if(rhTime > -2.) ++totHitTime;
	  if(totHitTime > 2 && withTime_filled == false){
	    hEtaDistr_WithTime->Fill(std::abs(etaGen));
	    withTime_filled = true;
	  }


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
	  

	  ++totHitsAbove0;
	  if(energyCharge > 60. && rhTime > -2.) ++totHitsAbove60; 
	}
	if(found){
	  if(debugCOUT) std::cout<< " bau 3 " << std::endl;
	  float rhEta = recHitTools.getEta(hitid);
	  float rhPhi = recHitTools.getPhi(hitid);

	  rhPt = rhEnergy/cosh(rhEta);
	  
	  h2_dPhivsdEta_rhGen->Fill(reco::deltaPhi(rhPhi, it_caloPart->phi()), rhEta - it_caloPart->eta());
	  h2_dPhivsdEta_rhAxis->Fill(reco::deltaPhi(rhPhi, axPhi), rhEta - axEta);
	  h2_dPhivsdEta_GenAxis->Fill(reco::deltaPhi(it_caloPart->phi(), axPhi), it_caloPart->eta() - axEta);

	  float rhX = recHitTools.getPosition(hitid).x();
	  float rhY = recHitTools.getPosition(hitid).y();
	  float rhL = recHitTools.getLayerWithOffset(hitid);
	  float rhZ = utilsMet.layerToZ(rhL, rhEta);

	  std::array<double,3> toRH{ {0., 0., rhZ} };
	  utilsMet.layerIntersection(toRH, fromAxis, vtx);
	  float Radius_RhGen = sqrt(pow(toRH[0]-rhX, 2) + pow(toRH[1]-rhY, 2));

	  //check which axis definition
	  //float Radius_RhGen = utilsMet.dsGenRecHit(axEta, axPhi, rhL, rhX, rhY);
	  //float Radius_RhGen = utilsMet.dsGenRecHit(etaGen, phiGen, rhL, rhX, rhY);

	  int etaBin = int((std::abs(etaGen) - binStart) / binWidth);
	  int iRadBin = -1;
	  //identify bin of radius around the gen axis
	  for(int ir=0; ir<nBinsRad; ++ir){
	    if(Radius_RhGen < radiusEtaRad[etaBin][ir]) {
	      iRadBin = ir;
	      break;
	    }
	  }
	  if(debugCOUT)	  std::cout << " bau 3 before call " << std::endl;

	  if(iRadBin != -1) {
	    for(int ir=iRadBin; ir<nBinsRad; ++ir){
	      totRHPerEtaRadius_allTime[etaBin][ir] += 1;
	      totEnergyPerEtaRadius[etaBin][ir] += rhPt;
	      if(rhTime > -2.){
		totEnergyWithTimePerEtaRadius[etaBin][ir] += rhPt;
		timePerEtaRadius[etaBin][ir] += rhTime; 
		totRHPerEtaRadius[etaBin][ir] += 1;

		timePerEtaRadius_AvgArm[etaBin][ir] += 1./(rhTime+1.); 
		totRHPerEtaRadius_AvgArm[etaBin][ir] += 1;
	
		hTime_Eta_dRadius[etaBin][ir]->Fill(rhTime);
	
		if(debugCOUT)   std::cout << " bau 3 before call " << std::endl;
		float resoWeigh = (float) myEstimator->getExpectedReso(hitid, rhEnergy);
		if(debugCOUT)   std::cout << " bau 3 before call resoWeigh = " << resoWeigh << std::endl;
		if(resoWeigh != 0.){
		  timePerEtaRadius_ResoWe[etaBin][ir] += rhTime * 1./(resoWeigh*resoWeigh); 
		  totRHPerEtaRadius_ResoWe[etaBin][ir] += 1./(resoWeigh*resoWeigh);
		}
		else{
		  timePerEtaRadius_ResoWe[etaBin][ir] += rhTime;
		  totRHPerEtaRadius_ResoWe[etaBin][ir] += 1;
		}
		if(ir < 2){
		  const std::pair<float, float> myPair(rhTime, resoWeigh);
		  timePerEtaRadiusDistr[etaBin][ir].push_back(myPair); 
		  timePerEtaRadiusHisto[etaBin][ir]->Fill(rhTime); 
		}
		if(debugCOUT)    std::cout << " bau 3 riempio raggio etaBin = " << etaBin << std::endl;
		h_rhGen_Radius_Eta[etaBin]->Fill(Radius_RhGen);
	      }
	      if(debugCOUT)    std::cout << " bau 3 fine hit with time " << std::endl;
	    }
	    if(debugCOUT)    std::cout << " bau 3 loop over rad " << std::endl;
	  }
	  if(debugCOUT)    std::cout << " bau 3 irad is reasonable " << std::endl;
	  h_allRH_TimesOfRadius->Fill(Radius_RhGen);
	}
	if(debugCOUT)    std::cout << " bau 3 fine found " << std::endl;
      }
      if(debugCOUT)    std::cout << " bau 3 hitsAndFractions " << std::endl;
    }// first loop over rechits

    if(debugCOUT)    std::cout << " bau 3 fine loop recHits " << std::endl;
  

    //compute Avg cutting biggest 32% and 68% around most probable
    for(int iet=0; iet<nBinsEta; ++iet){
      for(int irad=0; irad<2; ++irad){

	std::sort(timePerEtaRadiusDistr[iet][irad].begin(), timePerEtaRadiusDistr[iet][irad].end(), comparePairs);
	int totSize = timePerEtaRadiusDistr[iet][irad].size();
	float Ymax = 0.;
	float mpv = getXmax(timePerEtaRadiusHisto[iet][irad], Ymax);
	float sigma = timePerEtaRadiusHisto[iet][irad]->GetRMS();
	float sumCut = 0.;
	float sum68 = 0.;
	int numCut = 0;
	int num68 = 0;

	// if(debugCOUT)	std::cout << "ora problems size = " << totSize << std::endl;
	// for(int ij=0; ij<totSize; ++ij){
	//   if(float(ij+1) < totSize*0.3 || ij < 2){
	//     sumCut += timePerEtaRadiusDistr[iet][irad].at(ij);
	//     numCut += 1;
	//   }
	//   if(timePerEtaRadiusDistr[iet][irad].at(ij) < mpv + 2.*sigma || ij < 2){
	//     sum68 += timePerEtaRadiusDistr[iet][irad].at(ij);
	//     hTimeCut_Eta_dRadius[iet][irad]->Fill(timePerEtaRadiusDistr[iet][irad].at(ij));
        //     num68 += 1;
	//   }
	// }
	if(totSize > 0){
	  //	  std::cout << " >>> totSize = " << totSize << " max bin = " << int(totSize*0.32)  << " diff 0.68 = " << int(totSize*0.68) <<std::endl;
	float diffTimeValues = 999;
	float startTBin = 0;
	for(int ij=0; ij<int(totSize*0.32); ++ij){
	  float localDiff = fabs(timePerEtaRadiusDistr[iet][irad].at(ij).first - fabs(timePerEtaRadiusDistr[iet][irad].at(int(ij+totSize*0.68)).first));
	  //	  std::cout << " localDiff = " << localDiff << " jj = " << ij << " diffTimeValues = " << diffTimeValues << std::endl;
	  if(localDiff < diffTimeValues){
	    diffTimeValues = localDiff;
	    startTBin = ij;
	  }
	}



	int startBin = 0;
	int endBin = totSize;
	if(fabs(timePerEtaRadiusDistr[iet][irad].at(0).first - timePerEtaRadiusDistr[iet][irad].at(totSize-1).first) > 0.08){
	  startBin = startTBin;
	  endBin = int(startBin+totSize*0.68);
	  float HalfTimeDiff = std::abs(timePerEtaRadiusDistr[iet][irad].at(startBin).first - timePerEtaRadiusDistr[iet][irad].at(endBin).first) / 2.;
	  for(int ij=0; ij<startBin; ++ij){
	    if((timePerEtaRadiusDistr[iet][irad].at(ij).first) > (timePerEtaRadiusDistr[iet][irad].at(startBin).first - HalfTimeDiff) ){
	      startBin = ij;
	      break;
	    }
	  }
	  for(int ij=endBin; ij<totSize; ++ij){
	    if( (timePerEtaRadiusDistr[iet][irad].at(ij).first) > (timePerEtaRadiusDistr[iet][irad].at(endBin).first + HalfTimeDiff) ){
	      endBin = ij-1;
	      break;
	    }
	  }
	  
	}

	//	std::cout << " final bin start = " << startBin << " end bin = " << endBin << std::endl;

	for(int ij=startBin; ij<endBin; ++ij){
	  float resW = timePerEtaRadiusDistr[iet][irad].at(ij).second;
	  if(resW == 0) resW = 1;
	  sum68 += timePerEtaRadiusDistr[iet][irad].at(ij).first / (resW*resW);
	  hTimeCut_Eta_dRadius[iet][irad]->Fill(timePerEtaRadiusDistr[iet][irad].at(ij).first);
	  num68 += 1/(resW*resW);

	  sumCut += timePerEtaRadiusDistr[iet][irad].at(ij).first;
	  numCut += 1;
	}
	}
	timePerEtaRadiusAvgCutH[iet][irad] = sumCut;
	totRHPerEtaRadiusAvgCutH[iet][irad] = numCut;
	
	timePerEtaRadiusAvg68[iet][irad] = sum68;
	totRHPerEtaRadiusAvg68[iet][irad] = num68;


	timePerEtaRadiusHisto[iet][irad]->Delete();
	timePerEtaRadiusDistr[iet][irad].clear();
      }
    }
  
    if(debugCOUT) std::cout<< " bau 4 " << std::endl;


    //    if(debugCOUT) std::cout<< " bau 6 " << std::endl;
    if(debugCOUT) std::cout<< " >>>>>>>>> caloparticle "<< std::endl;
     
    //fill average distribution per eta - rad bin
    for(int iet=0; iet<nBinsEta; ++iet){
      for(int irad=0; irad<nBinsRad; ++irad){

	if(totRHPerEtaRadius_allTime[iet][irad] == 0.) continue;
	totEvtsEtaRadius[iet][irad] += 1;	

	hNumberHitsWithTime_Eta_dRadius[iet][irad]->Fill(totRHPerEtaRadius[iet][irad]);

	h_totEvtsEtaRadius_Eta_dRadius[iet][irad]->SetBinContent(1, h_totEvtsEtaRadius_Eta_dRadius[iet][irad]->GetBinContent(1)+1);

	//	hTime_Eta_dRadius[iet][irad]->Fill(1.*totEnergyPerEtaRadius[iet][irad]);

	if(totRHPerEtaRadius[iet][irad] > 2){
	  totEvtsEtaRadius_withTime[iet][irad] += 1;
	  h_totEvtsEtaRadius_withTime_Eta_dRadius[iet][irad]->SetBinContent(1, h_totEvtsEtaRadius_withTime_Eta_dRadius[iet][irad]->GetBinContent(1) + 1);

	  if(debugCOUT) std::cout<< " bau 6 v1 " << std::endl;

	  hEnergyWithTime_Eta_dRadius[iet][irad]->Fill(1.*totEnergyWithTimePerEtaRadius[iet][irad]);
	  hAverageTime_Eta_dRadius[iet][irad]->Fill(1.*timePerEtaRadius[iet][irad]/totRHPerEtaRadius[iet][irad]);
	  hAverageTime_Eta_dRadius_AvgArm[iet][irad]->Fill((1./(1.*timePerEtaRadius_AvgArm[iet][irad]/totRHPerEtaRadius_AvgArm[iet][irad])) - 1.);
	  hAverageTime_Eta_dRadius_ResoWe[iet][irad]->Fill(1.*timePerEtaRadius_ResoWe[iet][irad]/totRHPerEtaRadius_ResoWe[iet][irad]);
	  if(irad < 2){
	    hAverageTime_Eta_dRadius_AvgCutH[iet][irad]->Fill(1.*timePerEtaRadiusAvgCutH[iet][irad]/totRHPerEtaRadiusAvgCutH[iet][irad]);
	    hAverageTime_Eta_dRadius_Avg68[iet][irad]->Fill(1.*timePerEtaRadiusAvg68[iet][irad]/totRHPerEtaRadiusAvg68[iet][irad]);
	  }
	//	if(debugCOUT) std::cout<< " bau 6 v2 " << std::endl;
	  hFractionHitsWithTime_Eta_dRadius[iet][irad]->Fill(1.*totRHPerEtaRadius[iet][irad]/totRHPerEtaRadius_allTime[iet][irad]);	  
	}

	hTotHits_Eta_dRadius[iet][irad]->Fill(totRHPerEtaRadius_allTime[iet][irad]);
	hTotHitsWithTime_Eta_dRadius[iet][irad]->Fill(totRHPerEtaRadius[iet][irad]);
      }
    }

  
    if(debugCOUT2)    std::cout<< " >>> caloparticles fatta " << std::endl;
    hEtaDistrEvt->Fill(std::abs(etaGen));

    if(totHitsAbove60 > 2) hEtaDistrEvt3Hit->Fill(std::abs(etaGen));
    if(totHitsAbove0 > 0){
      energyChargeAll->Fill(std::abs(etaGen));
      if(totHitsAbove60 > 2) energyCharge->Fill(std::abs(etaGen));
    }

    if(debugCOUT) std::cout<< " bau fine caloparticle " << std::endl;

  }//caloparticle

  if(debugCOUT) std::cout<< " bau fine evento " << std::endl;
}

void
HGCalTimingAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
HGCalTimingAnalyzer::endJob()
{
  //  if(debugCOUT) std::cout<< " bau 7 " << std::endl;
  for(int iet=0; iet<nBinsEta; ++iet){
    for(int irad=0; irad<nBinsRad; ++irad){
      if(totEvtsEtaRadius[iet][irad] != 0) 
	hFractionEvents_HitsWithTime_Eta_dRadius[iet][irad]->Fill(1.*totEvtsEtaRadius_withTime[iet][irad]/totEvtsEtaRadius[iet][irad]);
    }
  }

  hEtaDistr_Fraction->Reset();
  hEtaDistr_Fraction->Divide(hEtaDistr_WithTime, hEtaDistr_AllTime, 1, 1);

  energyChargeFraction->Reset();
  energyChargeFraction->Divide(energyCharge, energyChargeAll);


  hEtaDistrEvt3HitFraction->Reset();
  hEtaDistrEvt3HitFraction->Divide(hEtaDistrEvt3Hit, hEtaDistrEvt);

  //  if(debugCOUT) std::cout<< " bau FINE " << std::endl;

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HGCalTimingAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HGCalTimingAnalyzer);
