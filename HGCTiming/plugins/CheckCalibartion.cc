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
#include "SimDataFormats/EncodedEventId/interface/EncodedEventId.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/GeometryVector/interface/GlobalTag.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/TruncatedPyramid.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TH2F.h"
#include "TH3F.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Math/interface/Point3D.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/ClusterTools.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/HGCalDepthPreClusterer.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include "HGCTimingAnalysis/HGCTiming/interface/UtilClasses.h"


#include <vector>
#include <string>
#include <map>

typedef Point3DBase<float, GlobalTag> Global3DPoint;
typedef Global3DPoint GlobalPoint;
typedef ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float> > XYZPointF;

class CheckCalibartion : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
//
// constructors and destructor
//
  explicit CheckCalibartion(const edm::ParameterSet&);
  ~CheckCalibartion();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;


  int  nEvents;
  int  nEventsGood;

 // ----------member data ---------------------------

  edm::EDGetTokenT<HGCRecHitCollection> _recHitsEE;
  edm::EDGetTokenT<HGCRecHitCollection> _recHitsFH;
  edm::EDGetTokenT<HGCRecHitCollection> _recHitsBH;
  edm::EDGetTokenT<std::vector<TrackingVertex> > _vtx;
  edm::EDGetTokenT<std::vector<reco::Vertex> > _recoVtx;
  edm::EDGetTokenT<std::vector<TrackingParticle> > _part;
  edm::EDGetTokenT<std::vector<CaloParticle> > _caloParticles;

  //  edm::EDGetTokenT<GlobalPoint> genVtxPositionToken_;
  edm::EDGetTokenT<XYZPointF> genVtxPositionToken_;
  edm::EDGetTokenT<float> genVtxTimeToken_;

  std::string                detector;
  int                        algo;
  HGCalDepthPreClusterer     pre;
  bool                       rawRecHits;
  float                      particleGenPt;
  int                        CaloPartPDGID;
  float                      timeOffset;
  hgcal::RecHitTools         recHitTools;

  std::unique_ptr<hgcal::ClusterTools> clusterTools;


  TH1F* h_Vtx_x;
  TH1F* h_Vtx_y;
  TH1F* h_Vtx_z;
  TH1F* h_Vtx_t;
  TH1F* h_VtxSurvived_z;

  TH1F* h_Vtx_dvx;
  TH1F* h_Vtx_dvy;
  TH1F* h_Vtx_dvz;

  //as function of the layer
  //as function of the eta wafer
  //raising the threshold
  TH1F* h_timeHits_layer[50][2];
  TH1F* h_timeHits_layer_E1p2[50][2];
  TH2F* h2_timeHitsEnergy_layer[50][2];

  TH2F* h2_XvsY_Uwafer_layer16_eta2;
  TH2F* h2_XvsY_Vwafer_layer16_eta2;

  TH2F* h_timeHits_dRvsLayer_gamma;
  TH2F* h_timeHits_dRvsLayer_ele;
  TH2F* h_timeHits_dRvsLayer_pion;
  TH2F* h_timeHits_dRvsLayer_kaon;

  bool debugCOUT4;


  //  float cSpeed = (1.e-7 * CLHEP::c_light);
  float cSpeed = 29.9792458;   //ns/cm

};  



CheckCalibartion::CheckCalibartion(const edm::ParameterSet& iConfig) :
  detector(iConfig.getParameter<std::string >("detector")),
  rawRecHits(iConfig.getParameter<bool>("rawRecHits")),
  particleGenPt(iConfig.getParameter<double>("particleGENPT")),
  CaloPartPDGID(iConfig.getParameter<int>("CaloPartPDGID")),
  timeOffset(iConfig.getParameter<double>("timeOffset"))
{
  nEvents = 0;
  nEventsGood = 0;

  debugCOUT4 = true;



  //now do what ever initialization is needed
  usesResource("TFileService");

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
  _recoVtx = consumes<std::vector<reco::Vertex> >(edm::InputTag("offlinePrimaryVertices4D"));
  _part = consumes<std::vector<TrackingParticle> >(edm::InputTag("mix","MergedTrackTruth"));
  _caloParticles = consumes<std::vector<CaloParticle> >(edm::InputTag("mix","MergedCaloTruth"));

  //  genVtxPositionToken_ = consumes<GlobalPoint>(edm::InputTag("genParticles:xyz0"));
  genVtxPositionToken_ = consumes<XYZPointF>(edm::InputTag("genParticles:xyz0"));
  genVtxTimeToken_ = consumes<float>(edm::InputTag("genParticles:t0"));


  auto sumes = consumesCollector();
  clusterTools = std::make_unique<hgcal::ClusterTools>(iConfig,sumes);

  edm::Service<TFileService> fs;

  h_Vtx_x = fs->make<TH1F>("h_Vtx_x", "", 1000, -15., 15.);
  h_Vtx_y = fs->make<TH1F>("h_Vtx_y", "", 1000, -15., 15.);
  h_Vtx_z = fs->make<TH1F>("h_Vtx_z", "", 1000, -15., 15.);
  h_Vtx_t = fs->make<TH1F>("h_Vtx_t", "", 1000, -15., 15.);
  h_VtxSurvived_z = fs->make<TH1F>("h_VtxSurvived_z", "", 1000, -15., 15.);

  h_Vtx_dvx = fs->make<TH1F>("h_Vtx_dvx", "", 1000, -10., 10.);
  h_Vtx_dvy = fs->make<TH1F>("h_Vtx_dvy", "", 1000, -10., 10.);
  h_Vtx_dvz = fs->make<TH1F>("h_Vtx_dvz", "", 1000, -10., 10.);

  for(int ij=0; ij<50; ++ij){
    for(int kl=0; kl<2; ++kl){
      h_timeHits_layer[ij][kl] = fs->make<TH1F>(Form("h_timeHits_layer%d_iS%d", ij, kl), "", 1000, -1., 9.);
      h_timeHits_layer_E1p2[ij][kl] = fs->make<TH1F>(Form("h_timeHits_layer_E1p2%d_iS%d", ij, kl), "", 1000, -1., 9.);
      h2_timeHitsEnergy_layer[ij][kl] = fs->make<TH2F>(Form("h2_timeHitsEnergy_layer%d_iS%d", ij, kl), "", 2000, 0., 20., 1000, -1., 9.);
      //h2_timeHitsEnergy_layer[ij][kl] = fs->make<TH2F>(Form("h2_timeHitsEnergy_layer%d_iS%d", ij, kl), "", 2000, 0., 20., 5000, -25., 25.);
    }
  }
  h2_XvsY_Uwafer_layer16_eta2 = fs->make<TH2F>("h2_XvsY_Uwafer_layer16_eta2", "", 600, -300, 300, 600, -300, 300);
  h2_XvsY_Vwafer_layer16_eta2 = fs->make<TH2F>("h2_XvsY_Vwafer_layer16_eta2", "", 600, -300, 300, 600, -300, 300);

  h_timeHits_dRvsLayer_gamma = fs->make<TH2F>("h_timeHits_dRvsLayer_gamma", "", 50, 0, 50, 300, -10., 10.);
  h_timeHits_dRvsLayer_ele = fs->make<TH2F>("h_timeHits_dRvsLayer_ele", "", 50, 0, 50, 300, -10., 10.);
  h_timeHits_dRvsLayer_pion = fs->make<TH2F>("h_timeHits_dRvsLayer_pion", "", 50, 0, 50, 300, -10., 10.);
  h_timeHits_dRvsLayer_kaon = fs->make<TH2F>("h_timeHits_dRvsLayer_kaon", "", 50, 0, 50, 300, -10., 10.);

}

CheckCalibartion::~CheckCalibartion()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


void
CheckCalibartion::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  ++nEvents;
  // if(nEvents != 383) return;
  // if(nEvents == 383) debugCOUT4 = true;
  if(debugCOUT4) std::cout<< " >>> analyzer " << std::endl;
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

  Handle<std::vector<reco::Vertex> > recoVtxHandle;
  iEvent.getByToken(_recoVtx,recoVtxHandle);
  const std::vector<reco::Vertex>& recoVtxs = *recoVtxHandle;

  Handle<std::vector<CaloParticle> > caloParticleHandle;
  iEvent.getByToken(_caloParticles, caloParticleHandle);
  const std::vector<CaloParticle>& caloParticles = *caloParticleHandle;

  //RA
  const auto& genVtxPositionHandle = iEvent.getHandle(genVtxPositionToken_);
  const auto& genVtxTimeHandle = iEvent.getHandle(genVtxTimeToken_);
  // std::unique_ptr<math::XYZTLorentzVectorF> genPV(nullptr);
  // genPV = std::make_unique<math::XYZTLorentzVectorF>(genVtxPositionHandle->x(), 
  // 						     genVtxPositionHandle->y(), genVtxPositionHandle->z(), *(genVtxTimeHandle));


  float vx = 0.;
  float vy = 0.;
  float vz = 0.;
  float vt = 0.;

  //for TTBAR
  /*
  for(auto ij : vtxs){
    if(vtxs[0].eventId().bunchCrossing() == 0 && vtxs[0].eventId().event() == 0){
      //if(vtxs.size()!=0){
      vx = ij.position().x();
      vy = ij.position().y();
      vz = ij.position().z();
      vt = ij.position().t() * 1.e9;
      std::cout << " vt (ps) = " << vt << " in s = " << vtxs[0].position().t() << std::endl;
      //std::cout << " vtx bunchCrossing() = " << vtxs[0].eventId().bunchCrossing() << " event = " << vtxs[0].eventId().event() << std::endl;
    }
  }
  std::cout << " vtxs.size() = " << vtxs.size() << std::endl;
  */

  //using reco highest pT
  /*
  for(auto ij : recoVtxs){
    if(vtxs.size()!=0){
      vx = vtxs[0].position().x();
      vy = vtxs[0].position().y();
      vz = vtxs[0].position().z();
      vt = vtxs[0].position().t();
      std::cout << " vt (ns) = " << vt << std::endl;
      //std::cout << " vtx bunchCrossing() = " << vtxs[0].eventId().bunchCrossing() << " event = " << vtxs[0].eventId().event() << std::endl;
    }
  }
  */


  //RA

  vx = genVtxPositionHandle->x();
  vy = genVtxPositionHandle->y();
  vz = genVtxPositionHandle->z();
  vt = *(genVtxTimeHandle);
  //  std::cout << " vt (ps) = " << vt << " in s = " << vtxs[0].position().t() << std::endl;
  //  std::cout << " vtxs.size() = " << vtxs.size() << std::endl;


  std::array<double,3> vtx{ {vx, vy, vz}};


  h_Vtx_x->Fill(vx);
  h_Vtx_y->Fill(vy);
  h_Vtx_z->Fill(vz);
  h_Vtx_t->Fill(vt);


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

      NewrechitsEE = *recHitHandleEE;
      NewrechitsFH = *recHitHandleFH;
      for(unsigned int i = 0; i < NewrechitsEE.size(); ++i){
	
	const HGCalDetId hitid = NewrechitsEE[i].detid();
	int rhL = recHitTools.getLayerWithOffset(hitid);
	hitmap[NewrechitsEE[i].detid()] = &NewrechitsEE[i];
      }
      for(unsigned int i = 0; i < NewrechitsFH.size(); ++i){
	
	const HGCalDetId hitid = NewrechitsFH[i].detid();
	int rhL = recHitTools.getLayerWithOffset(hitid);
	hitmap[NewrechitsFH[i].detid()] = &NewrechitsFH[i];
      }
      break;
    }
  case 2:
    {
      iEvent.getByToken(_recHitsEE,recHitHandleEE);
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
      NewrechitsFH = *recHitHandleFH;
      for(unsigned int i = 0; i < NewrechitsFH.size(); i++){
	hitmap[NewrechitsFH[i].detid()] = &NewrechitsFH[i];
      }
      break;
    }
  default:
    break;
  }

  if(debugCOUT4)    std::cout << " now loop over hits " << std::endl;

  //check all the hits in the eta-layer region
  /////////////////////////////
  for(std::map<DetId, const HGCRecHit*>::iterator iop=hitmap.begin(); iop != hitmap.end(); ++iop){
    const HGCalDetId hitid = iop->first;
    const HGCRecHit* hit = iop->second;
    
    if(hit->time() == -1) continue;

    float rhEnergy = hit->energy();
    float rhTime = hit->time() - timeOffset;

    GlobalPoint rhPos = recHitTools.getPosition(hitid);

    float rhX = rhPos.x();
    float rhY = rhPos.y();
    float rhZ = rhPos.z();

    int rhL = recHitTools.getLayerWithOffset(hitid);

    float rhEta = recHitTools.getEta(rhPos);
    float rhPhi = recHitTools.getPhi(rhPos);
    float rhPt = rhEnergy/cosh(rhEta);

    int iSide = rhEta > 0 ? 1 : 0;

    //selection of layer
    //selection on energy
    
    if(debugCOUT4) std::cout << " >>> rhPt = " << rhPt << " rhEnergy = " << rhEnergy << " rhEta = " << rhEta << " time = " << rhTime << std::endl;


     float dist000 = sqrt(rhX*rhX + rhY*rhY + rhZ*rhZ);
     float distVtx = sqrt((rhX - vx)*(rhX - vx) + (rhY - vy)*(rhY - vy) + (rhZ - vz)*(rhZ - vz));
     //float deltaT = (dist000 - distVtx) / cSpeed - vt;
     float deltaT = (dist000 - distVtx) / cSpeed; // for nugun... problem with vtxT


     // float dist000 = rhZ;                                          
     // float distVtx = (rhZ - vz);
     //    float deltaT = vz / cSpeed;

    if(debugCOUT4) std::cout << " >>> deltaT = " << deltaT << std::endl;

    if(rhL < 49){
      h_timeHits_layer[rhL][iSide]->Fill(rhTime + deltaT); 
      if(rhEnergy > 1.2) h_timeHits_layer_E1p2[rhL][iSide]->Fill(rhTime + deltaT);
      h2_timeHitsEnergy_layer[rhL][iSide]->Fill(rhEnergy, (rhTime + deltaT)); 
    }

    int waferU = recHitTools.getWafer(hitid).first;
    int waferV = recHitTools.getWafer(hitid).second;

    if(rhL == 16 && std::abs(rhEta - 2.) < 0.5 && rhEta > 0 && waferV == 3){
      // std::cout << " getWaferU = " << recHitTools.getWafer(hitid).first 
      //  		<< " getWaferV = " << recHitTools.getWafer(hitid).second 
      // 	 	<< " X = " << rhX << " Y = " << rhY << " eta = " << rhEta << std::endl;
      h2_XvsY_Uwafer_layer16_eta2->Fill(rhY, rhX, waferU);
      h2_XvsY_Vwafer_layer16_eta2->Fill(rhY, rhX, waferV);
    }

  }// loop over rechits
  if(debugCOUT4) std::cout<< " end of calibration event process " << std::endl;
  float ptGen = -1;
  float eGen = -1;
  float etaGen = -1;
  float phiGen = -1;
  float xGen = -1;
  float yGen = -1;
  float zGen = -1;
  int numbCaloPart = 0;
  if(debugCOUT4) std::cout << " >>> caloParticles.size() = " << caloParticles.size() << std::endl;
  // loop over caloParticles
  for(std::vector<CaloParticle>::const_iterator it_caloPart = caloParticles.begin(); it_caloPart != caloParticles.end(); ++it_caloPart){
    ptGen = it_caloPart->pt();
    eGen = it_caloPart->energy();
    etaGen = it_caloPart->eta();
    phiGen = it_caloPart->phi();
    xGen = it_caloPart->momentum().x();
    yGen = it_caloPart->momentum().y();
    zGen = it_caloPart->momentum().z();
    if(debugCOUT4) std::cout << " ptGen = " << ptGen << " eGen = " << eGen  << " simCl size = " << it_caloPart->simClusters().size() 
			     << " pdgId = " << it_caloPart->pdgId() << std::endl;
    if(it_caloPart->simClusters().size() > 1) continue;

    std::vector<float> rHits_time;
    std::vector<float> rHits_x;
    std::vector<float> rHits_y;
    std::vector<float> rHits_z;
    std::vector<float> rHits_L;
    
    float axX = -1;
    float axY = -1;
    float axZ = -1;
    float sumEnergyToNorm = 0.;
    GlobalPoint showerAxis;

    for(CaloParticle::sc_iterator scIt=it_caloPart->simCluster_begin(); scIt!=it_caloPart->simCluster_end(); ++scIt){
      //all hits and energy fractions at sim level                                      
    std::vector<std::pair<uint32_t,float> > hits = (*scIt)->hits_and_fractions();
    if(debugCOUT4)  std::cout << " >>> hits.size() = " << hits.size() << std::endl;
    for(auto ij : hits) {
      DetId hitid = (ij.first); 
      auto finder = hitmap.find(hitid);
      if(finder == hitmap.end()) continue;
      else{
	if(debugCOUT4) std::cout << " trovata sim reco match hit " << std::endl;
	const HGCRecHit *hit = hitmap[hitid];
	float rhEnergy = hit->energy();
	sumEnergyToNorm += rhEnergy*ij.second;
	axX += rhEnergy*ij.second * recHitTools.getPosition(hitid).x();
	axY += rhEnergy*ij.second * recHitTools.getPosition(hitid).y(); 
	axZ += rhEnergy*ij.second * recHitTools.getPosition(hitid).z();

	if(debugCOUT4) std::cout << " time = " << hit->time() << " error = " << hit->timeError() << std::endl;
	if(hit->timeError() != -1){
	  float rhTime = hit->time() - timeOffset;
	  float rhX = recHitTools.getPosition(hitid).x();
	  float rhY = recHitTools.getPosition(hitid).y();
	  float rhZ = recHitTools.getPosition(hitid).z();
	  int rhL = recHitTools.getLayerWithOffset(hitid);
	  float rhEta = recHitTools.getEta(recHitTools.getPosition(hitid));
	  float rhPt = rhEnergy/cosh(rhEta);

	  rHits_time.push_back(rhTime);
	  rHits_x.push_back(rhX);
	  rHits_y.push_back(rhY);
	  rHits_z.push_back(rhZ);
	  rHits_L.push_back(rhL);

	}
      }
    }// loop over hits of simcluster
    }//loop over caloparticles simclusters

    showerAxis = GlobalPoint(axX/sumEnergyToNorm, axY/sumEnergyToNorm, (axZ - vz)/sumEnergyToNorm);  
    axX = showerAxis.x();                                                                                  
    axY = showerAxis.y();                                                                                  
    axZ = showerAxis.z();

    if(debugCOUT4) std::cout<< " after showerAxis " << std::endl;

    //    UtilClasses utilsMet = UtilClasses(scIt->eta(), scIt->phi());
    UtilClasses utilsMet = UtilClasses(etaGen, phiGen);
    std::array<double,3> fromAxis;
    /*
      fromAxis[0] = (axX);
      fromAxis[1] = (axY);
      fromAxis[2] = (axZ);
    */
      fromAxis[0] = (xGen);
      fromAxis[1] = (yGen);
      fromAxis[2] = (zGen);

      if(debugCOUT4) std::cout<< " rHits_time.size() = " << rHits_time.size() << std::endl;
     //now loop over recHits for the given simCluster
      for(unsigned int ij=0; ij<rHits_time.size(); ++ij){

      float rhX = rHits_x.at(ij);
      float rhY = rHits_y.at(ij);
      float rhZ = rHits_z.at(ij);
 
      std::array<double,3> to{ {0., 0., rhZ} };
      utilsMet.layerIntersection(to, fromAxis, vtx);
      float deltaR = sqrt(pow(to[0] - rhX, 2) + pow(to[1] - rhY, 2));

     
      float dist000 = sqrt(rhX*rhX + rhY*rhY + rhZ*rhZ);
      float distVtx = sqrt((rhX - vx)*(rhX - vx) + (rhY - vy)*(rhY - vy) + (rhZ - vz)*(rhZ - vz));
      float deltaT = (dist000 - distVtx) / cSpeed - vt;
      //float deltaT = (dist000 - distVtx) / cSpeed; // for nugun... problem with vtxT 

      if(it_caloPart->pdgId() == 22) h_timeHits_dRvsLayer_gamma->Fill(deltaR, rHits_L.at(ij), rHits_time.at(ij)+deltaT);
      if(std::abs(it_caloPart->pdgId()) == 11) h_timeHits_dRvsLayer_ele->Fill(deltaR, rHits_L.at(ij), rHits_time.at(ij)+deltaT);
      if(std::abs(it_caloPart->pdgId()) == 211) h_timeHits_dRvsLayer_pion->Fill(deltaR, rHits_L.at(ij), rHits_time.at(ij)+deltaT);
      if(std::abs(it_caloPart->pdgId()) == 130) h_timeHits_dRvsLayer_kaon->Fill(deltaR, rHits_L.at(ij), rHits_time.at(ij)+deltaT);
    }
  }//caloparticle
}

void
CheckCalibartion::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
CheckCalibartion::endJob()
{

}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
CheckCalibartion::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(CheckCalibartion);
