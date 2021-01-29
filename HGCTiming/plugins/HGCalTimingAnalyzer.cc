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
#include "DataFormats/GeometrySurface/interface/BoundDisk.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

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
  void propagateTrack(const edm::Event &ev, const edm::EventSetup &es, const reco::Track &tk);

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  void beginRun(edm::Run const& iEvent, edm::EventSetup const& es) override;
  //  virtual void beginJob() override;
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
  edm::EDGetTokenT<std::vector<pat::Muon>> _muonSrc;
  const HGCalDDDConstants* hgcons_;
  inline static const std::string detectorName_ = "HGCalEESensitive";
  inline static const std::string propName_ = "PropagatorWithMaterial";
  edm::ESHandle<MagneticField> bfield_;
  std::unique_ptr<GeomDet> firstDisk_[2];
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
  _muonSrc = consumes<std::vector<pat::Muon> >(edm::InputTag("slimmedMuons"));
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

  h_minR_pos = fs->make<TH2F>("h_minR_pos", "", 11, 0., 11, 1000, 100., 600.);
  h_minR_neg = fs->make<TH2F>("h_minR_neg", "", 11, 0., 11, 1000, 100., 600.);
  h_minR = fs->make<TH2F>("h_minR", "", 11, 0., 11, 1000, 100., 600.);
  h_maxR_pos = fs->make<TH2F>("h_maxR_pos", "", 11, 0., 11, 1000, 100., 600.);
  h_maxR_neg = fs->make<TH2F>("h_maxR_neg", "", 11, 0., 11, 1000, 100., 600.);
  h_maxR = fs->make<TH2F>("h_maxR", "", 11, 0., 11, 1000, 100., 600.);



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
  //  initialize(es); 
}


void HGCalTimingAnalyzer::initialize(const edm::EventSetup &es) {
  edm::ESHandle<HGCalDDDConstants> hdc = es.getHandle(hdc_token_);
  hgcons_ = hdc.product();
  
  buildLayers();

  bfield_ = es.getHandle(bfield_token_);
  //  es.get<TrackingComponentsRecord>().get(propagator_token_, propagator_);
  // const auto& propagator = es.getData(propagator_token_);
  propagator_ =  es.get<TrackingComponentsRecord>().getHandle(propagator_token_);

  //  edm::ESHandle<Propagator> propH;
  //es.get<TrackingComponentsRecord>().get(propName_, propagator_);
  //const Propagator* prop = propH.product();

}

void HGCalTimingAnalyzer::buildLayers() {
  float zVal = hgcons_->waferZ(1, true);
  std::pair<double, double> rMinMax = hgcons_->rangeR(zVal, true);

  for (int iSide = 0; iSide < 2; ++iSide) {
    float zSide = (iSide == 0) ? (-1. * zVal) : zVal;
    firstDisk_[iSide] =
      std::make_unique<GeomDet>(Disk::build(Disk::PositionType(0, 0, zSide),
					    Disk::RotationType(),
					    SimpleDiskBounds(rMinMax.first, rMinMax.second, zSide - 0.5, zSide + 0.5))
				.get());
  }
}


void HGCalTimingAnalyzer::propagateTrack(const edm::Event &ev,
					 const edm::EventSetup &es, const reco::Track& tk){
  auto bFieldProd = bfield_.product();
  const Propagator &prop = (*propagator_);

  FreeTrajectoryState fts = trajectoryStateTransform::outerFreeState((tk), bFieldProd);
  int iSide = int(tk.eta() > 0);
  TrajectoryStateOnSurface tsos = prop.propagate(fts, firstDisk_[iSide]->surface());
  if (tsos.isValid()) {
    //      result.emplace_back(tsos.globalPosition(), tsos.globalMomentum(), iSide, i, trkId);
  }
  
  // sorting seeding region by descending momentum
  // std::sort(result.begin(), result.end(), [](const TICLSeedingRegion &a, const TICLSeedingRegion &b) {
  //     return a.directionAtOrigin.perp2() > b.directionAtOrigin.perp2();
  //   });

    /*
GlobalPoint v = tscbl.trackStateAtPCA().position();                                                                                                        
  math::XYZPoint pos(v.x(), v.y(), v.z());                                                                                                                   
  GlobalVector p = tscbl.trackStateAtPCA().momentum();                                                                                                       
  math::XYZVector mom(p.x(), p.y(), p.z());
     */
}

/*
bool HGCalTimingAnalyzer::getTrajectoryStateClosestToBeamLine(const Trajectory& traj,
					 const reco::BeamSpot& bs,
					 const Propagator* thePropagator,
					 TrajectoryStateClosestToBeamLine& tscbl) {
  // get the state closest to the beamline
    TrajectoryStateOnSurface stateForProjectionToBeamLineOnSurface =
      traj.closestMeasurement(GlobalPoint(bs.x0(), bs.y0(), bs.z0())).updatedState();

    if (!stateForProjectionToBeamLineOnSurface.isValid()) {
      edm::LogError("CannotPropagateToBeamLine") << "the state on the closest measurement isnot valid. skipping track.";
      return false;
    }

    const FreeTrajectoryState& stateForProjectionToBeamLine = *stateForProjectionToBeamLineOnSurface.freeState();

    TSCBLBuilderWithPropagator tscblBuilder(*thePropagator);
    tscbl = tscblBuilder(stateForProjectionToBeamLine, bs);

    return tscbl.isValid();
}

void HGCalTimingAnalyzer::getPositionOnLayer(){

  TrajectoryStateClosestToBeamLine tscbl;
  bool tsbcl_status = getTrajectoryStateClosestToBeamLine(traj, bs, thePropagator, tscbl);

  if (!tsbcl_status)
    return reco::Track();

  GlobalPoint v = tscbl.trackStateAtPCA().position();
  math::XYZPoint pos(v.x(), v.y(), v.z());
  GlobalVector p = tscbl.trackStateAtPCA().momentum();
  math::XYZVector mom(p.x(), p.y(), p.z());

}
*/

void
HGCalTimingAnalyzer::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup)
{

  ++nEvents;
  //if(nEvents != 14) return;
  // if(nEvents == 383) debugCOUT4 = true;
  if(debugCOUT) std::cout<< " >>> analyzer evt = " << nEvents << std::endl;
  using namespace edm;


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

  edm::Handle<std::vector<pat::Muon>> muons;
  iEvent.getByToken(_muonSrc, muons);

  if(debugCOUT)  std::cout << " nMuons = " << muons.product()->size() << std::endl;
  for(const pat::Muon & muon : *muons){
    if (muon.pt() < 10. && debugCOUT) std::cout << " recoMu = " << muon.pt() << " " << muon.eta() << " " << muon.phi() << std::endl;
    if (!muon.combinedMuon().isNull()){
      reco::TrackRef muonTrk = muon.combinedMuon();
      if(debugCOUT)
	std::cout << " combinedMuon innerPosition() = " << muonTrk->innerPosition() << " innerMomentum() = " << muonTrk->innerMomentum() 
		  << " outerPosition() = " << muonTrk->outerPosition() << " outerMomentum() = " << muonTrk->outerMomentum() << std::endl;
    }
    else if (!muon.globalTrack().isNull()){
      reco::TrackRef muonTrk = muon.globalTrack();
      if(debugCOUT)
	std::cout << "globalTrk  innerPosition() = " << muonTrk->innerPosition() << " innerMomentum() = " << muonTrk->innerMomentum() 
		  << " outerPosition() = " << muonTrk->outerPosition() << " outerMomentum() = " << muonTrk->outerMomentum() << std::endl;
    }
    else if(!muon.standAloneMuon().isNull()){
      reco::TrackRef muonTrk = muon.standAloneMuon();
      if(debugCOUT)
	std::cout << "standaloneMuon  innerPosition() = " << muonTrk->innerPosition() << " innerMomentum() = " << muonTrk->innerMomentum()
		  << " outerPosition() = " << muonTrk->outerPosition() << " outerMomentum() = " << muonTrk->outerMomentum() << std::endl;
    }
    else if(!muon.track().isNull()){
      reco::TrackRef muonTrk = muon.track();
      if(debugCOUT)
	std::cout << "track  innerPosition() = " << muonTrk->innerPosition() << " innerMomentum() = " << muonTrk->innerMomentum()
		  << " outerPosition() = " << muonTrk->outerPosition() << " outerMomentum() = " << muonTrk->outerMomentum() << std::endl;
    }
    
    //propagateTrack(iEvent, iSetup, muonTrk);

  }



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
  float minRvalues[11][2];
  float maxRvalues[11][2];
  for(int ij=0; ij<11; ++ij){
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

    if(sectionType != 2) continue;

    int iSide = rhZ > 0 ? 0 : 1;
    float rhR = sqrt(rhX * rhX + rhY * rhY);
    minRvalues[rhL - 39][iSide] = (minRvalues[rhL - 39][iSide] == 0) ? rhR : std::min(rhR, minRvalues[rhL - 39][iSide]);
    maxRvalues[rhL - 39][iSide] = (maxRvalues[rhL - 39][iSide] == 0) ? rhR : std::min(rhR, maxRvalues[rhL - 39][iSide]);

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

    if(debugCOUT)
      std::cout << " MIP = " << MIP << " charge =  " << charge << " SoN = " << SoverN << std::endl;
  }


  for(int ij = 0; ij < 11; ++ ij){
    h_minR_pos->Fill(ij, minRvalues[ij][1]);
    h_minR_neg->Fill(ij, minRvalues[ij][0]);
    h_minR->Fill(ij, minRvalues[ij][0]);
    h_minR->Fill(ij, minRvalues[ij][1]);
    h_maxR_pos->Fill(ij, maxRvalues[ij][1]);
    h_maxR_neg->Fill(ij, maxRvalues[ij][0]);
    h_maxR->Fill(ij, maxRvalues[ij][0]);
    h_maxR->Fill(ij, maxRvalues[ij][1]);
  }
  return;

  



  
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
