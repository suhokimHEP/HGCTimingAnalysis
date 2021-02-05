// user include files

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/TruncatedPyramid.h"
#include "Geometry/HGCalCommonData/interface/HGCalDDDConstants.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCHEDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCScintillatorDetId.h"
#include "DataFormats/GeometrySurface/interface/BoundDisk.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/TrajectorySeed/interface/PropagationDirection.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/ClusterTools.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/HGCalDepthPreClusterer.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include "RecoLocalCalo/HGCalRecProducers/interface/ComputeClusterTime.h"

#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"

#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"

#include "HGCTimingAnalysis/HGCTiming/interface/UtilClasses.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"

#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TH2F.h"
#include "TH3F.h"

#include <vector>
#include <string>
#include <map>


class HGCalTimingAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources, edm::one::WatchRuns> {
  //class HGCalTimingAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
//
// constructors and destructor
//
  explicit HGCalTimingAnalyzer(const edm::ParameterSet&);
  ~HGCalTimingAnalyzer();
  

  void initialize(const edm::EventSetup &es);
  void buildLayers();
  void propagateTrack(const edm::Event &ev, const edm::EventSetup &es, const reco::Track &tk, int idx);

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  void beginRun(edm::Run const& iEvent, edm::EventSetup const& es) override;
  //void beginJob() override;
  void analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) override;
  void endJob() override;
  void endRun(edm::Run const& iEvent, edm::EventSetup const&) override {};

private:
  int  nEvents;
  int  nEventsGood;

 // ----------member data ---------------------------

  edm::EDGetTokenT<HGCRecHitCollection> _recHitsEE;
  edm::EDGetTokenT<HGCRecHitCollection> _recHitsFH;
  edm::EDGetTokenT<HGCRecHitCollection> _recHitsBH;
  edm::EDGetTokenT<std::vector<TrackingVertex> > _vtx;
  edm::EDGetTokenT<std::vector<TrackingParticle> > _part;
  edm::EDGetTokenT<std::vector<CaloParticle> > _caloParticles;
  edm::EDGetTokenT<std::vector<pat::Muon>> _muonSrcM;
  edm::EDGetTokenT<std::vector<reco::Muon>> _muonSrc;
  const HGCalDDDConstants* hgcons_;
  inline static const std::string detectorName_ = "HGCalEESensitive";
  inline static const std::string propName_ = "PropagatorWithMaterial";
  edm::ESHandle<MagneticField> bfield_;
  edm::ESGetToken<HGCalDDDConstants, IdealGeometryRecord> hdc_token_;
  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> bfield_token_;
  edm::ESGetToken<Propagator, TrackingComponentsRecord> propagator_token_;
  edm::ESHandle<Propagator> propagator_;
  //  const Propagator& prop

  std::string                detector;
  int                        algo;
  HGCalDepthPreClusterer     pre;
  bool                       rawRecHits;
  float                      particleGenPt;
  int                        CaloPartPDGID;
  float                      timeOffset;
  hgcal::RecHitTools         recHitTools;



  std::vector<int> layersToSkip;

  std::unique_ptr<hgcal::ClusterTools> clusterTools;

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
  TH1F* h_VtxSurvived_z;

  TH1F* h_Vtx_dvx;
  TH1F* h_Vtx_dvy;
  TH1F* h_Vtx_dvz;


  TH2F* h_minR_pos;
  TH2F* h_minR_neg;
  TH2F* h_minR;
  TH2F* h_maxR_pos;
  TH2F* h_maxR_neg;
  TH2F* h_maxR;

  TH1F* h_minL;

  TProfile* h_layersZ_pos;
  TProfile* h_layersZ_neg;
  TProfile* h_layersZ;

  inline static const int nSciLayers = 14;
  inline static const int firstSciLayer = 37;

  std::unique_ptr<GeomDet> firstDisk_[2][nSciLayers];

  TH2F* h2_occupancy[nSciLayers];
  TH2F* energyDistr[nSciLayers];
  TH2F* MIPDistr[nSciLayers];
  TH2F* SoNDistr[nSciLayers];
  TProfile2D* h2_RvsL;

  std::map<int, std::vector<float>> xposOnlayer[2];
  std::map<int, std::vector<float>> yposOnlayer[2];
  std::map<int, std::vector<float>> momOnlayer[2];
  std::map<int, std::vector<int>> muIdxOnlayer[2];


  TTree* newT;
  std::vector<float> muonP;
  std::vector<float> muonEta;
  std::vector<float> muonPhi;
  std::vector<float> crossX;
  std::vector<float> crossY;
  std::vector<float> crossZ;
  std::vector<float> crossL;
  std::vector<float> crossM;
  std::vector<float> recHitX;
  std::vector<float> recHitY;
  std::vector<float> recHitZ;
  std::vector<int> recHitiR;
  std::vector<int> recHitiPhi;
  std::vector<int> recHitiEta;
  std::vector<int> recHitL;
  std::vector<float> recHitEne;
  std::vector<int> recHitMip;
  std::vector<float> recHitNoise;

  bool debugCOUT;
};  



HGCalTimingAnalyzer::HGCalTimingAnalyzer(const edm::ParameterSet& iConfig) :
  detector(iConfig.getParameter<std::string >("detector")),
  rawRecHits(iConfig.getParameter<bool>("rawRecHits")),
  particleGenPt(iConfig.getParameter<double>("particleGENPT")),
  CaloPartPDGID(iConfig.getParameter<int>("CaloPartPDGID")),
  timeOffset(iConfig.getParameter<double>("timeOffset"))
{
  nEvents = 0;
  nEventsGood = 0;

  debugCOUT = false;

  //now do what ever initialization is needed
  usesResource("TFileService");
  auto sumes = consumesCollector();
  clusterTools = std::make_unique<hgcal::ClusterTools>(iConfig,sumes);


  if(detector=="all") {
    _recHitsEE = consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("HGCEEInput"));
    _recHitsFH = consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("HGCFHInput"));
    _recHitsBH = consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("HGCBHInput"));
    algo = 1;
  }
  _vtx = consumes<std::vector<TrackingVertex> >(edm::InputTag("mix","MergedTrackTruth"));
  _part = consumes<std::vector<TrackingParticle> >(edm::InputTag("mix","MergedTrackTruth"));
  _caloParticles = consumes<std::vector<CaloParticle> >(edm::InputTag("mix","MergedCaloTruth"));
  _muonSrcM = consumes<std::vector<pat::Muon> >(edm::InputTag("slimmedMuons"));
  _muonSrc = consumes<std::vector<reco::Muon> >(edm::InputTag("muons"));
  //propagation
  hdc_token_ = sumes.esConsumes<HGCalDDDConstants, IdealGeometryRecord, edm::Transition::BeginRun>(edm::ESInputTag("", detectorName_));
  bfield_token_ = sumes.esConsumes<MagneticField, IdealMagneticFieldRecord, edm::Transition::BeginRun>();
  propagator_token_ = sumes.esConsumes<Propagator, TrackingComponentsRecord, edm::Transition::BeginRun>(edm::ESInputTag("", propName_));
  //propagator_token_ = sumes.esConsumes<Propagator, TrackingComponentsRecord>(edm::ESInputTag("", propName_));

  //parameters to provide conversion GeV - MIP
  keV2fC[0] =  iConfig.getParameter<double>("HGCEE_keV2fC");
  keV2fC[1] =  iConfig.getParameter<double>("HGCHEF_keV2fC");
  keV2MIP = iConfig.getParameter<double>("HGCHB_keV2MIP");
  noisefC[0] = (iConfig.getParameter<std::vector<double> >("HGCEE_noisefC")).at(0);
  noisefC[1] = (iConfig.getParameter<std::vector<double> >("HGCEE_noisefC")).at(1);
  noisefC[2] = (iConfig.getParameter<std::vector<double> >("HGCEE_noisefC")).at(2);
  noiseMIP = iConfig.getParameter<double>("HGCBH_noiseMIP");
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
  //end param conversion

  edm::Service<TFileService> fs;
  newT = fs->make<TTree>("newT", "");
  newT->Branch("muonP", &muonP);
  newT->Branch("muonEta", &muonEta);
  newT->Branch("muonPhi", &muonPhi);
  newT->Branch("crossX", &crossX);
  newT->Branch("crossY", &crossY);
  newT->Branch("crossZ", &crossZ);
  newT->Branch("crossL", &crossL);
  newT->Branch("crossM", &crossM);
  newT->Branch("recHitX", &recHitX);
  newT->Branch("recHitY", &recHitY);
  newT->Branch("recHitZ", &recHitZ);
  newT->Branch("recHitiR", &recHitiR);
  newT->Branch("recHitiPhi", &recHitiPhi);
  newT->Branch("recHitiEta", &recHitiEta);
  newT->Branch("recHitL", &recHitL);
  newT->Branch("recHitEne", &recHitEne);
  newT->Branch("recHitMip", &recHitMip);
  newT->Branch("recHitNoise", &recHitNoise);


  h_minR_pos = fs->make<TH2F>("h_minR_pos", "", nSciLayers, 0., nSciLayers, 1000, 100., 600.);
  h_minR_neg = fs->make<TH2F>("h_minR_neg", "", nSciLayers, 0., nSciLayers, 1000, 100., 600.);
  h_minR = fs->make<TH2F>("h_minR", "", nSciLayers, 0., nSciLayers, 1000, 100., 600.);
  h_maxR_pos = fs->make<TH2F>("h_maxR_pos", "", nSciLayers, 0., nSciLayers, 1000, 100., 600.);
  h_maxR_neg = fs->make<TH2F>("h_maxR_neg", "", nSciLayers, 0., nSciLayers, 1000, 100., 600.);
  h_maxR = fs->make<TH2F>("h_maxR", "", nSciLayers, 0., nSciLayers, 1000, 100., 600.);

  h_minL = fs->make<TH1F>("h_minL", "", 50., 0., 50.);

  h_layersZ_pos = fs->make<TProfile>("h_layersZ_pos", "", 50., 0., 50.);
  h_layersZ_neg = fs->make<TProfile>("h_layersZ_neg", "", 50., 0., 50.);
  h_layersZ = fs->make<TProfile>("h_layersZ", "", 50., 0., 50.);


  for(int ij=0; ij<nSciLayers; ++ij){
    h2_occupancy[ij] = fs->make<TH2F>(Form("h2_occupancy_L%d", ij+firstSciLayer), "", 600, -300., 300, 600, -300., 300.);
    energyDistr[ij] = fs->make<TH2F>(Form("energyDistr_L%d", ij+firstSciLayer), "", 100, 0., 100., 100., 0., 100.);
    MIPDistr[ij] = fs->make<TH2F>(Form("MIPDistr_L%d", ij+firstSciLayer), "", 100, 0., 100., 100., 0., 100.);
    SoNDistr[ij] = fs->make<TH2F>(Form("SoNDistr_L%d", ij+firstSciLayer), "", 100, 0., 100., 100., 0., 100.);
  }
  h2_RvsL = fs->make<TProfile2D>("h2_RvsL", "", 50, 0., 50., 300, 0., 300);



  h_Vtx_x = fs->make<TH1F>("h_Vtx_x", "", 1000, -15., 15.);
  h_Vtx_y = fs->make<TH1F>("h_Vtx_y", "", 1000, -15., 15.);
  h_Vtx_z = fs->make<TH1F>("h_Vtx_z", "", 1000, -15., 15.);
  h_VtxSurvived_z = fs->make<TH1F>("h_VtxSurvived_z", "", 1000, -15., 15.);

  h_Vtx_dvx = fs->make<TH1F>("h_Vtx_dvx", "", 1000, -10., 10.);
  h_Vtx_dvy = fs->make<TH1F>("h_Vtx_dvy", "", 1000, -10., 10.);
  h_Vtx_dvz = fs->make<TH1F>("h_Vtx_dvz", "", 1000, -10., 10.);



  
}

HGCalTimingAnalyzer::~HGCalTimingAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

void HGCalTimingAnalyzer::beginRun(edm::Run const& iEvent, edm::EventSetup const& es) { 
  initialize(es); 
}


void HGCalTimingAnalyzer::initialize(const edm::EventSetup &es) {
  edm::ESHandle<HGCalDDDConstants> hdc = es.getHandle(hdc_token_);
  hgcons_ = hdc.product();
  
  buildLayers();

  bfield_ = es.getHandle(bfield_token_);
  propagator_ =  es.get<TrackingComponentsRecord>().getHandle(propagator_token_);
}

void HGCalTimingAnalyzer::buildLayers() {
  float minR[nSciLayers] = {153., 153., 153., 153., 137., 137., 119., 119., 119., 119., 104., 104., 104., 104};
  float maxR[nSciLayers] = {199., 203., 213., 218., 229., 240., 252., 252., 252., 252., 252., 252., 252., 252.};
  float zVal[nSciLayers] = {411.29, 416.739, 422.187, 427.636, 436.172, 444.722, 453.263, 461.817, 470.371, 478.925, 487.47, 496.024, 504.577, 513.127};


  for (int iSide = 0; iSide < 2; ++iSide) {
    for (int iL = 0; iL < nSciLayers; ++iL) {
      float zSide = (iSide == 0) ? (-1. * zVal[iL]) : zVal[iL];    
      
      firstDisk_[iSide][iL] = std::make_unique<GeomDet>(Disk::build(Disk::PositionType(0, 0, zSide),
								    Disk::RotationType(),
								    SimpleDiskBounds(minR[iL], maxR[iL], zSide - 0.5, zSide + 0.5)).get());
    }
  }
  
}


void HGCalTimingAnalyzer::propagateTrack(const edm::Event &ev,
					 const edm::EventSetup &es, const reco::Track& tk, int idx){
  auto bFieldProd = bfield_.product();
  const Propagator &prop = (*propagator_);

  if(debugCOUT)
  std::cout << " propagateTrack p = " << tk.p() << " pt = " << tk.pt() << " eta = " << tk.eta() << std::endl;

  FreeTrajectoryState fts = trajectoryStateTransform::outerFreeState((tk), bFieldProd);
  int iSide = int(tk.eta() > 0);
  for(int ij=0; ij<nSciLayers; ++ij){
    TrajectoryStateOnSurface tsos = prop.propagate(fts, firstDisk_[iSide][ij]->surface());

    if (tsos.isValid()) {
      auto position = tsos.globalPosition();
      math::XYZPoint pos(position.x(), position.y(), position.z()); 
      xposOnlayer[iSide][ij].push_back(pos.x());
      yposOnlayer[iSide][ij].push_back(pos.y());
      momOnlayer[iSide][ij].push_back(tk.p());
      muIdxOnlayer[iSide][ij].push_back(idx);
      if(debugCOUT)
	std::cout << " valid layer = " << ij << " trkX = " << pos.x() << " ttrkY = " << pos.y() << std::endl;

      muonP.push_back(tk.p());
      muonEta.push_back(tk.eta());
      muonPhi.push_back(tk.phi());
      crossX.push_back(pos.x());
      crossY.push_back(pos.y());
      crossZ.push_back(pos.z());
      crossL.push_back(ij+firstSciLayer);
      crossM.push_back(idx);
    }
  }
}


void
HGCalTimingAnalyzer::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup)
{

  ++nEvents;
  //if(nEvents != 14) return;
  // if(nEvents == 383) debugCOUT4 = true;
  if(debugCOUT) std::cout<< " >>> analyzer evt = " << nEvents << std::endl;
  using namespace edm;

  for(int ij=0; ij<2; ++ij){
    xposOnlayer[ij].clear();
    yposOnlayer[ij].clear();
    momOnlayer[ij].clear();
    muIdxOnlayer[ij].clear();
  }
  recHitTools.getEventSetup(iSetup);

  muonP.clear();
  muonEta.clear();
  muonPhi.clear();
  crossX.clear();
  crossY.clear();
  crossZ.clear();
  crossL.clear();
  crossM.clear();
  recHitX.clear();
  recHitY.clear();
  recHitZ.clear();
  recHitiR.clear();
  recHitiPhi.clear();
  recHitiEta.clear();
  recHitL.clear();
  recHitEne.clear();
  recHitMip.clear();
  recHitNoise.clear();
  

  Handle<HGCRecHitCollection> recHitHandleEE;
  Handle<HGCRecHitCollection> recHitHandleFH;
  Handle<HGCRecHitCollection> recHitHandleBH;

  /*
  Handle<std::vector<TrackingVertex> > vtxHandle;
  Handle<std::vector<TrackingParticle> > partHandle;
  iEvent.getByToken(_vtx,vtxHandle);
  iEvent.getByToken(_part,partHandle);
  const std::vector<TrackingVertex>& vtxs = *vtxHandle;
  const std::vector<TrackingParticle>& part = *partHandle;

  Handle<std::vector<CaloParticle> > caloParticleHandle;
  iEvent.getByToken(_caloParticles, caloParticleHandle);
  const std::vector<CaloParticle>& caloParticles = *caloParticleHandle;
  */


  /*
  edm::Handle<std::vector<pat::Muon>> muonsM;
  iEvent.getByToken(_muonSrcM, muonsM);
  */

  edm::Handle<std::vector<reco::Muon>> muons;
  iEvent.getByToken(_muonSrc, muons);

  if(debugCOUT)  std::cout << " nMuons = " << muons.product()->size() << std::endl;

  
  int mCount = -1;
  //  for(const pat::Muon & muon : *muons){
  for(const reco::Muon & muon : *muons){
    ++mCount;
    if (debugCOUT) std::cout << " recoMu = " << muon.p() << " " << muon.eta() << " " << muon.phi() << " isGlobal = " << muon.isGlobalMuon() << std::endl;
    if(!muon.isGlobalMuon() ) continue;
    //if(std::abs(muon.eta()) < 1.7 || std::abs(muon.eta()) > 2.8) continue;

    // if (!muon.combinedMuon().isNull()){
    //   reco::TrackRef muonTrk = muon.combinedMuon();
    //   if(debugCOUT)
    // 	std::cout << " combinedMuon innerPosition() = " << muonTrk->innerPosition() << " innerMomentum() = " << muonTrk->innerMomentum() 
    // 		  << " outerPosition() = " << muonTrk->outerPosition() << " outerMomentum() = " << muonTrk->outerMomentum() << std::endl;
    // }
    // else if (!muon.globalTrack().isNull()){
    //   reco::TrackRef muonTrk = muon.globalTrack();
    //   if(debugCOUT)
    // 	std::cout << "globalTrk  innerPosition() = " << muonTrk->innerPosition() << " innerMomentum() = " << muonTrk->innerMomentum() 
    // 		  << " outerPosition() = " << muonTrk->outerPosition() << " outerMomentum() = " << muonTrk->outerMomentum() << std::endl;
    // }
    // else if(!muon.standAloneMuon().isNull()){
    //   reco::TrackRef muonTrk = muon.standAloneMuon();
    //   if(debugCOUT)
    // 	std::cout << "standaloneMuon  innerPosition() = " << muonTrk->innerPosition() << " innerMomentum() = " << muonTrk->innerMomentum()
    // 		  << " outerPosition() = " << muonTrk->outerPosition() << " outerMomentum() = " << muonTrk->outerMomentum() << std::endl;
    // }
    // else 

    // if (!muon.globalTrack().isNull()){
    //   reco::Track muonTrk = *(muon.globalTrack());

    if(!muon.track().isNull()){
      reco::Track muonTrk = *(muon.track());
      if(debugCOUT)
	std::cout << "track  innerPosition() = " << muonTrk.innerPosition() << " innerMomentum() = " << muonTrk.innerMomentum()
		  << " outerPosition() = " << muonTrk.outerPosition() << " outerMomentum() = " << muonTrk.outerMomentum() << " mCount = " << mCount << std::endl;

      propagateTrack(iEvent, iSetup, muonTrk, mCount);
    }
  }
  
  if(debugCOUT){
  for(int iS=0; iS<2; ++iS){
    for(auto ipos : xposOnlayer[iS]){
      auto layer = ipos.first;   
      unsigned int iSize = ipos.second.size();
      for(unsigned int ij=0; ij<iSize; ++ij){
	std::cout << " iS = " << iS << " layer = " << layer 
		  << " x = " << xposOnlayer[iS][layer][ij] << " y = " << yposOnlayer[iS][layer][ij] 
		  << " momentum = " << momOnlayer[iS][layer][ij] 
		  << " muIdxOnlayer = " << muIdxOnlayer[iS][layer][ij] << " xCheck muon P = " << (*muons)[muIdxOnlayer[iS][layer][ij]].p() << std::endl;
      }
    }
  }
  }

  if(crossM.size() == 0) return;
  
  /*
  float vx = 0.;
  float vy = 0.;
  float vz = 0.;
  if(vtxs.size()!=0){
    vx = vtxs[0].position().x();
    vy = vtxs[0].position().y();
    vz = vtxs[0].position().z();
  }

  vz = 0.;

  h_Vtx_x->Fill(vx);
  h_Vtx_y->Fill(vy);
  h_Vtx_z->Fill(vz);


  for(unsigned int i=0;i<part.size();++i){
    float dvx=0.;
    float dvy=0.;
    float dvz=0.;
    
    if(part[i].decayVertices().size()>=1){   
      dvx=part[i].decayVertices()[0]->position().x();
      dvy=part[i].decayVertices()[0]->position().y();
      dvz=part[i].decayVertices()[0]->position().z();

      h_Vtx_dvx->Fill(dvx);
      h_Vtx_dvy->Fill(dvy);
      h_Vtx_dvz->Fill(dvz);
    }
  }
  */

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
      
      for (auto const& it : *recHitHandleEE)
        hitmap[it.detid().rawId()] = &(it);
      for (auto const& it : *recHitHandleFH)
        hitmap[it.detid().rawId()] = &(it);
      for (auto const& it : *recHitHandleBH)
        hitmap[it.detid().rawId()] = &(it);
      break;
    }
  default:
    break;
  }


  ///////////////////
  float minL = 99;
  float minRvalues[nSciLayers][2];
  float maxRvalues[nSciLayers][2];
  for(int ij=0; ij<nSciLayers; ++ij){
    minRvalues[ij][0] = 0.;
    minRvalues[ij][1] = 0.;
    maxRvalues[ij][0] = 0.;
    maxRvalues[ij][1] = 0.;
  }

  for(std::map<DetId, const HGCRecHit*>::iterator iop=hitmap.begin(); iop != hitmap.end(); ++iop){
    const HGCalDetId hitid = iop->first;
    const HGCRecHit* hit = iop->second;

    bool found = false;
    float rhEnergy = hit->energy();
    float rhTime = hit->time() - timeOffset;
    float CPfraction = 0.;
    float rhX = recHitTools.getPosition(hitid).x();
    float rhY = recHitTools.getPosition(hitid).y();
    float rhZ = recHitTools.getPosition(hitid).z();
    int rhL = recHitTools.getLayerWithOffset(hitid);
    float rhEta = recHitTools.getEta(recHitTools.getPosition(hitid));
    float rhPt = rhEnergy/cosh(rhEta);

    //to extract some conversion factors 
    unsigned int layer = recHitTools.getLayerWithOffset(hitid);
    int thick = int(recHitTools.getSiThickness(hitid)) / 100 - 1;

    int sectionType = -1;
    if (hitid.det() == DetId::HGCalEE) sectionType = 0; 
    else if (hitid.det() == DetId::HGCalHSi) sectionType = 1;
    else if (hitid.det() == DetId::HGCalHSc) sectionType = 2;

    if(rhZ > 0) {
      h_layersZ_pos->Fill(rhL, rhZ);
      h_layersZ->Fill(rhL, rhZ);
    }
    else{ 
      h_layersZ_neg->Fill(rhL, -1.*rhZ);
      h_layersZ->Fill(rhL, -1.*rhZ);
    }
    float rhR = sqrt(rhX * rhX + rhY * rhY);
    h2_RvsL->Fill(rhL, rhR, 1.*sectionType);


    if(sectionType != 2) continue;
    HGCScintillatorDetId sciId(hitid);
    int rhiEta = sciId.ietaAbs();
    int rhiPhi = sciId.iradiusphi().second;
    int rhiR = sciId.iradiusphi().first;


    if(rhL < minL){
      minL = rhL;
    }
 
    //    std::cout << " rhL = " << rhL << " firstSciLayer = " << firstSciLayer << " diff = " << rhL - firstSciLayer << std::endl;
    //    if((rhL - firstSciLayer) < 0 || (rhL - firstSciLayer) >= nSciLayers) continue;
    h2_occupancy[rhL - firstSciLayer]->Fill(rhX, rhY);

    int iSide = int(rhZ) > 0;

    minRvalues[rhL - firstSciLayer][iSide] = (minRvalues[rhL - firstSciLayer][iSide] == 0) ? rhR : std::min(rhR, minRvalues[rhL - firstSciLayer][iSide]);
    maxRvalues[rhL - firstSciLayer][iSide] = (maxRvalues[rhL - firstSciLayer][iSide] == 0) ? rhR : std::max(rhR, maxRvalues[rhL - firstSciLayer][iSide]);

    if(debugCOUT)
      std::cout << " thick = " << thick << " sectionType = " << sectionType 
		<< " rhL = " << rhL << " rhX = " << rhX << " rhY = " << rhY << " rhZ = " << rhZ << std::endl;

    int energyMIP = 0.;
    if(sectionType == 2) energyMIP = hit->energy()/keV2GeV * keV2MIP;
    else if(sectionType == 0 || sectionType == 1) energyMIP = hit->energy()/scaleCorrection.at(thick)/keV2GeV / (weights.at(layer)/keV2MeV);

    float energyCharge = 0.;
    if(sectionType == 2) energyCharge = energyMIP * 1.;
    else if(sectionType == 0 || sectionType == 1) energyCharge = energyMIP * fCPerMIP[thick];

    double sigmaNoiseMIP = 1.;
    if(sectionType == 2) sigmaNoiseMIP = noiseMIP;
    else if(sectionType == 0 || sectionType == 1) sigmaNoiseMIP = noisefC[thick]/fCPerMIP[thick];

    float charge = energyCharge;
    float MIP = energyMIP;
    float SoverN = energyMIP / sigmaNoiseMIP;

    recHitX.push_back(rhX);
    recHitY.push_back(rhY);
    recHitZ.push_back(rhZ);
    recHitiR.push_back(rhiR);
    recHitiPhi.push_back(rhiPhi);
    recHitiEta.push_back(rhiEta);
    recHitL.push_back(rhL);
    recHitEne.push_back(rhEnergy);
    recHitMip.push_back(MIP);
    recHitNoise.push_back(sigmaNoiseMIP);



    // std::cout << " rh iS = " << iSide << " rhL = " << rhL << " rhX = " << rhX << " rhY = " << rhY << std::endl; 

    //    std::cout << " xposOnlayer[rhL].size() = " << xposOnlayer[iSide][rhL- firstSciLayer].size() << " layer = " << rhL- firstSciLayer << std::endl;

    if(debugCOUT){
      int ic = -1;
      for(auto icount : xposOnlayer[iSide][rhL- firstSciLayer]){
	++ic;
	// std::cout << " ic = " << ic << " rhL = " << rhL- firstSciLayer << " rhX = " << rhX << " icount = " << icount << std::endl;
	// std::cout << " rhY = " << rhY << " icount = " << yposOnlayer[iSide][rhL- firstSciLayer][ic] << std::endl;
	if(std::abs(icount - rhX) < 3.){
	  if(std::abs(rhY - yposOnlayer[iSide][rhL- firstSciLayer][ic]) < 3.){
	    energyDistr[rhL - firstSciLayer]->Fill(momOnlayer[iSide][rhL- firstSciLayer][ic], rhEnergy);
	    MIPDistr[rhL - firstSciLayer]->Fill(momOnlayer[iSide][rhL- firstSciLayer][ic], MIP);
	    SoNDistr[rhL - firstSciLayer]->Fill(momOnlayer[iSide][rhL- firstSciLayer][ic], SoverN);
	    std::cout << " trovato layer = " << rhL- firstSciLayer << " energy = " << rhEnergy << " p = " << momOnlayer[iSide][rhL- firstSciLayer][ic]
		      << " rhX = " << rhX << " trkX = " << icount << " rhY = " << rhY << " trkY = " <<  yposOnlayer[iSide][rhL- firstSciLayer][ic]  << std::endl;
	  }
	}
      }
    }



    if(debugCOUT)
      std::cout << " MIP = " << MIP << " charge =  " << charge << " SoN = " << SoverN << std::endl;
  }

  h_minL->Fill(minL);

  for(int ij = 0; ij < nSciLayers; ++ ij){
    h_minR_pos->Fill(ij, minRvalues[ij][1]);
    h_minR_neg->Fill(ij, minRvalues[ij][0]);
    h_minR->Fill(ij, minRvalues[ij][0]);
    h_minR->Fill(ij, minRvalues[ij][1]);
    h_maxR_pos->Fill(ij, maxRvalues[ij][1]);
    h_maxR_neg->Fill(ij, maxRvalues[ij][0]);
    h_maxR->Fill(ij, maxRvalues[ij][0]);
    h_maxR->Fill(ij, maxRvalues[ij][1]);
  }

  newT->Fill();

  if(debugCOUT) std::cout<< " end of event process " << std::endl;
}

// void
// HGCalTimingAnalyzer::beginJob()
// {
// }

// ------------ method called once each job just after ending the event loop  ------------
void
HGCalTimingAnalyzer::endJob()
{

  // for(int ij=1; ij <= h_layersZ->GetNbinsX()+1; ++ij){
  //   std::cout << " layer = " << ij << " Z = " << h_layersZ->GetBinContent(ij) << std::endl;
  // }

  std::cout << " totEvents = " << nEvents << " events good = " << nEventsGood << " fraction non interacting = " << 1.*nEventsGood/nEvents << std::endl;

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
