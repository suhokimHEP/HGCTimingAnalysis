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
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "Geometry/Records/interface/CaloGeometryRecord.h"
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


#include "RecoLocalCalo/HGCalRecAlgos/interface/ClusterTools.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/HGCalDepthPreClusterer.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

//#include "RecoParticleFlow/PFClusterProducer/plugins/SimMappers/ComputeClusterTime.h"
#include "HGCTimingAnalysis/HGCTiming/interface/UtilClasses.h"


#include "DataFormats/Common/interface/ValueMap.h"

#include <vector>
#include <string>
#include <map>

class HGCalTiming2DClustersAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
//
// constructors and destructor
//
  explicit HGCalTiming2DClustersAnalyzer(const edm::ParameterSet&);
  ~HGCalTiming2DClustersAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;



 // ----------member data ---------------------------

  edm::EDGetTokenT<HGCRecHitCollection> _recHitsEE;
  edm::EDGetTokenT<HGCRecHitCollection> _recHitsFH;
  edm::EDGetTokenT<HGCRecHitCollection> _recHitsBH;
  edm::EDGetTokenT<std::vector<TrackingVertex> > _vtxT;
  edm::EDGetTokenT<std::vector<SimVertex> > _vtx;
  edm::EDGetTokenT<std::vector<TrackingParticle> > _part;
  edm::EDGetTokenT<std::vector<CaloParticle> > _caloParticles;
  edm::EDGetTokenT<reco::CaloClusterCollection> _caloClusters;
  edm::EDGetTokenT<edm::ValueMap<float>> _timeMap;
  edm::EDGetTokenT<reco::GenParticleCollection> _genToken;




  std::string                detector;
  int                        algo;
  float                      particleGenPt;
  int                        CaloPartPDGID;
  int                        timeOffset;
  hgcal::RecHitTools         recHitTools;


  //  std::unique_ptr<hgcal::ClusterTools> clusterTools;

  //define the histograms/objects to fill and save
  TH1F* h_Vtx_x;
  TH1F* h_Vtx_y;
  TH1F* h_Vtx_z;
  TH1F* h_Vtx_t;

  TH1F* h_eta;
  TH1F* h_phi;
  TH1F* h_Recoeta;
  TH1F* h_Recophi;

  TH1F* nSimClusters;
  TH1F* zVtx_SimClustersTracks;
  TH1F* zVtx_SimClustersTracks_EE;
  TH1F* nCaloParticles;

  TH2F* h2_deltaR_vs_layer;
  TH2F* h2_Z_vs_layer;
  TH2F* h2_Z_vs_phi;
  TH2F* h2_Z_vs_eta;
  TH2F* h2_eta_vs_phi;
  TH2F* h2_Geta_vs_phi;

  TH2F* h2_dX_vs_layer;
  TH2F* h2_dY_vs_layer;

  TH1F* h_dX_Sci;
  TH1F* h_dY_Sci;
  TH2F* h2_dX_vs_dY_Sci;

  TH1F* h_nhits_energyWeighted[12];


  int nEvents;
  int nEventsGood;
  bool debugCOUT;

  int nParticles;
  int nParticlesGood;
  int noGOODParticle;

};  



HGCalTiming2DClustersAnalyzer::HGCalTiming2DClustersAnalyzer(const edm::ParameterSet& iConfig) :
  detector(iConfig.getParameter<std::string >("detector")),
  particleGenPt(iConfig.getParameter<double>("particleGENPT")),
  CaloPartPDGID(iConfig.getParameter<int>("CaloPartPDGID")),
  timeOffset(iConfig.getParameter<int>("timeOffset"))
{
  nEvents = 0;
  nEventsGood = 0;

  nParticles = 0;
  nParticlesGood = 0;
  noGOODParticle = 0;

  debugCOUT = true;
  

  //retrieve collections 
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
  _vtxT = consumes<std::vector<TrackingVertex> >(edm::InputTag("mix","MergedTrackTruth"));
  _vtx = consumes<std::vector<SimVertex> >(edm::InputTag("g4SimHits"));
  _part = consumes<std::vector<TrackingParticle> >(edm::InputTag("mix","MergedTrackTruth"));
  _caloParticles = consumes<std::vector<CaloParticle> >(edm::InputTag("mix","MergedCaloTruth"));
  _caloClusters = consumes<reco::CaloClusterCollection>(iConfig.getParameter<edm::InputTag>("caloClusterInput"));
  _timeMap = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("timeMapInput"));
  _genToken = consumes<reco::GenParticleCollection>( edm::InputTag("genParticles"));



  //  auto sumes = consumesCollector();
  //  clusterTools = std::make_unique<hgcal::ClusterTools>(iConfig,sumes);

  edm::Service<TFileService> fs;


  //initialize the histograms
  h_Vtx_x = fs->make<TH1F>("h_Vtx_x", "", 1000, -500., 500.);
  h_Vtx_y = fs->make<TH1F>("h_Vtx_y", "", 1000, -500., 500.);
  h_Vtx_z = fs->make<TH1F>("h_Vtx_z", "", 1000, -500., 500.);
  h_Vtx_t = fs->make<TH1F>("h_Vtx_t", "", 1000, -0.5, 0.5);

  h_eta = fs->make<TH1F>("h_eta", "", 600, -3., 3.);
  h_phi = fs->make<TH1F>("h_phi", "", 640, -3.2, 3.2);
  h_Recoeta = fs->make<TH1F>("h_Recoeta", "", 600, -3., 3.);
  h_Recophi = fs->make<TH1F>("h_Recophi", "", 640, -3.2, 3.2);

  nSimClusters = fs->make<TH1F>("nSimClusters", "", 1000, 0., 100.);
  zVtx_SimClustersTracks = fs->make<TH1F>("zVtx_SimClustersTracks", "", 1000, 0., 500.);
  zVtx_SimClustersTracks_EE = fs->make<TH1F>("zVtx_SimClustersTracks_EE", "", 1000, 0., 500.);
  nCaloParticles = fs->make<TH1F>("nCaloParticles", "", 100, 0., 100.);

  h2_deltaR_vs_layer = fs->make<TH2F>("h2_deltaR_vs_layer", "", 60, 0., 60., 300, 0., 300.);
  h2_Z_vs_layer = fs->make<TH2F>("h2_Z_vs_layer", "", 60, 0., 60., 1000, -500., 500.);

  h2_Z_vs_phi = fs->make<TH2F>("h2_Z_vs_phi", "", 640, -3.2, 3.2, 1000, -500., 500.);
  h2_Z_vs_eta = fs->make<TH2F>("h2_Z_vs_eta", "", 600, -3., 3., 300, -500., 500.);
  h2_eta_vs_phi = fs->make<TH2F>("h2_eta_vs_phi", "", 600, -3., 3., 640, -3.2, 3.2);
  h2_Geta_vs_phi = fs->make<TH2F>("h2_Geta_vs_phi", "", 600, -3., 3., 640, -3.2, 3.2);

  h2_Geta_vs_phi->Sumw2();

  h2_dX_vs_layer = fs->make<TH2F>("h2_dX_vs_layer", "", 60, 0., 60., 600, -300., 300.);
  h2_dY_vs_layer = fs->make<TH2F>("h2_dY_vs_layer", "", 60, 0., 60., 600, -300., 300.);

  h_dX_Sci = fs->make<TH1F>("h_dX_Sci", "", 6000, -300., 300.);
  h_dY_Sci = fs->make<TH1F>("h_dY_Sci", "", 6000, -300., 300.);
  h2_dX_vs_dY_Sci = fs->make<TH2F>("h2_dX_vs_dY_Sci", "", 6000, -300., 300., 6000, -300., 300.);

  for(int ij=0; ij<12; ++ij) h_nhits_energyWeighted[ij] = fs->make<TH1F>(Form("h_nhits_energyWeighted_layer%d", ij), "", 10, 0., 10.);

}

HGCalTiming2DClustersAnalyzer::~HGCalTiming2DClustersAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


void
HGCalTiming2DClustersAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  ++nEvents;
  if(debugCOUT) std::cout<< " >>> analyzer " << std::endl;
  using namespace edm;

  //  if(nEvents != 73) return;
  recHitTools.getEventSetup(iSetup);

  Handle<HGCRecHitCollection> recHitHandleEE;
  Handle<HGCRecHitCollection> recHitHandleFH;
  Handle<HGCRecHitCollection> recHitHandleBH;

  Handle<std::vector<TrackingVertex> > vtxHandleT;
  Handle<std::vector<SimVertex> > vtxHandle;
  Handle<std::vector<TrackingParticle> > partHandle;
  iEvent.getByToken(_vtxT,vtxHandleT);
  iEvent.getByToken(_vtx,vtxHandle);
  iEvent.getByToken(_part,partHandle);
  const std::vector<TrackingVertex>& vtxsT = *vtxHandleT;
  const std::vector<SimVertex>& vtxs = *vtxHandle;
  const std::vector<TrackingParticle>& part = *partHandle;

  Handle<std::vector<CaloParticle> > caloParticleHandle;
  iEvent.getByToken(_caloParticles, caloParticleHandle);
  const std::vector<CaloParticle>& caloParticles = *caloParticleHandle;

  Handle<reco::CaloClusterCollection> caloClusters;
  iEvent.getByToken(_caloClusters, caloClusters);

  edm::Handle<edm::ValueMap<float> > time2DMap;
  iEvent.getByToken(_timeMap, time2DMap);
  
  edm::Handle<reco::GenParticleCollection> genPart;
  iEvent.getByToken(_genToken, genPart);
  //  unsigned int gensize = genParticles->size();

  
  float vx = 0.;
  float vy = 0.;
  float vz = 0.;
  float vt = 0.;
  if(vtxs.size()!=0){
    vx = vtxs[0].position().x();
    vy = vtxs[0].position().y();
    vz = vtxs[0].position().z();
    vt = vtxs[0].position().t();
  }

  h_Vtx_x->Fill(vx);
  h_Vtx_y->Fill(vy);
  h_Vtx_z->Fill(vz);
  h_Vtx_t->Fill(vt);

  if(debugCOUT)
    std::cout << " simVtx = " << vx << " " << vy << " " << vz << " " << vt 
	      << " trckingVtx = " << vtxsT[0].position().x() << " " << vtxsT[0].position().y() << " " << vtxsT[0].position().z() << " " << vtxsT[0].position().t() << std::endl;

  /*
  const  CaloGeometry *geom_;
  edm::ESHandle<CaloGeometry> geom;
  iSetup.get<CaloGeometryRecord>().get(geom);
  geom_ = geom.product();
  */

  HGCRecHitCollection NewrechitsEE;
  HGCRecHitCollection NewrechitsFH;
  HGCRecHitCollection NewrechitsBH;

  //make a map detid-rechit
  std::map<DetId,const HGCRecHit*> hitmap;
  std::map<DetId,const HGCRecHit*> hitmapSci;
  switch(algo){
  case 1:
    {
      iEvent.getByToken(_recHitsEE,recHitHandleEE);
      iEvent.getByToken(_recHitsFH,recHitHandleFH);
      iEvent.getByToken(_recHitsBH,recHitHandleBH);

      NewrechitsEE = *recHitHandleEE;
      NewrechitsFH = *recHitHandleFH;
      NewrechitsBH = *recHitHandleBH;
      for(unsigned int i = 0; i < NewrechitsEE.size(); ++i){

        const HGCalDetId hitid = NewrechitsEE[i].detid();
        int rhL = recHitTools.getLayerWithOffset(hitid);
        hitmap[NewrechitsEE[i].detid()] = &NewrechitsEE[i];

	/*
	int subdet = hitid.det();
	//auto geomEE = static_cast<const HGCalGeometry*>(geom_->getSubdetectorGeometry(DetId::HGCalEE,ForwardSubdetector::ForwardEmpty));
	const CaloSubdetectorGeometry *subGeom = static_cast<const HGCalGeometry*>(geom_->getSubdetectorGeometry(DetId::Detector(subdet), ForwardSubdetector::ForwardEmpty));
	const std::vector<DetId> & ids = subGeom->getValidDetIds();
	std::cout << "on subdet " << subdet << " I got a total of " << ids.size() << " det ids " << std::endl;
	int zside = recHitTools.zside(hitid);
	std::cout << " >> EE det hitid.detector = " << hitid.det()  << " subdet = " << hitid.subdetId() << " zside = " << zside << " rhL = " << rhL << std::endl;

	HGCalDetId pippo(HGCalDetId(ForwardSubdetector::ForwardEmpty, zside,rhL,0,0,0)); 
	HGCSiliconDetId pippo2(DetId::Detector(subdet), zside, 0, rhL, 0, 0, 0, 0);
	std::cout << " >> EE det pippo .detector = " << pippo.det()  << " subdet = " << pippo.subdetId() << std::endl;
	std::cout << " >> EE det pippo2 .detector = " << pippo2.det()  << " subdet = " << pippo2.subdetId() 
		  << pippo2.zside() << " lay = " << pippo2.layer() << std::endl;
	int newZside = recHitTools.zside(pippo2);
	int newL = recHitTools.getLayerWithOffset(pippo2);
	std::cout << " zside = " << newZside << " newrhL = " << newL << std::endl;
	*/
      }
      for(unsigned int i = 0; i < NewrechitsFH.size(); ++i){

        const HGCalDetId hitid = NewrechitsFH[i].detid();
        int rhL = recHitTools.getLayerWithOffset(hitid);
        hitmap[NewrechitsFH[i].detid()] = &NewrechitsFH[i];

	/*
	int zside = recHitTools.zside(hitid);

	std::cout << " >> FH det hitid.detector = " << hitid.det() << " subdet = " << hitid.subdetId() << " zside = " << zside << " rhL = " << rhL << std::endl;
	int subdet = hitid.det();
	const CaloSubdetectorGeometry *subGeom = geom_->getSubdetectorGeometry(DetId::Detector(subdet), ForwardSubdetector::ForwardEmpty);
	const std::vector<DetId> & ids = subGeom->getValidDetIds();
	std::cout << "on subdet " << subdet << " I got a total of " << ids.size() << " det ids " << std::endl;

	HGCalDetId pippo(HGCalDetId(ForwardSubdetector::ForwardEmpty, zside,rhL,0,0,0)); 
	std::cout << " >> FH det pippo .detector = " << pippo.det()  << " subdet = " << pippo.subdetId() << std::endl;

	HGCSiliconDetId pippo2(DetId::Detector(subdet), zside, 1, rhL, 0, 0, 0, 0);
	std::cout << " >> FH det pippo2 .detector = " << pippo2.det()  << " subdet = " << pippo2.subdetId() << " zside = " 
		  << pippo2.zside() << " lay = " << pippo2.layer() << std::endl;
	int newZside = recHitTools.zside(pippo2);
	int newL = recHitTools.getLayerWithOffset(pippo2);
	std::cout << " zside = " << newZside << " newrhL = " << newL << std::endl;
	*/
      }
      for(unsigned int i = 0; i < NewrechitsBH.size(); ++i){
        hitmapSci[NewrechitsBH[i].detid()] = &NewrechitsBH[i];

	/*
	const HGCalDetId hitid = NewrechitsBH[i].detid();
        int rhL = recHitTools.getLayerWithOffset(hitid);
	int zside = recHitTools.zside(hitid);
	std::cout << " >> BH det hitid.detector = " << hitid.det() << " subdet = " << hitid.subdetId() << " zside = " << zside << " rhL = " << rhL << std::endl;
	int subdet = hitid.det();
	const CaloSubdetectorGeometry *subGeom = geom_->getSubdetectorGeometry(DetId::Detector(subdet), ForwardSubdetector::ForwardEmpty);
	const std::vector<DetId> & ids = subGeom->getValidDetIds();
	std::cout << "on subdet " << subdet << " I got a total of " << ids.size() << " det ids " << std::endl;

	HGCalDetId pippo(HGCalDetId(ForwardSubdetector::ForwardEmpty, zside,rhL,0,0,0)); 
	std::cout << " >> BH det pippo .detector = " << pippo.det()  << " subdet = " << pippo.subdetId() << std::endl;
	HGCScintillatorDetId pippo2(DetId::Detector(subdet), rhL, 0, 0);
	std::cout << " >> BH det pippo2 .detector = " << pippo2.det()  << " subdet = " << pippo2.subdetId() << std::endl;
	int newZside = recHitTools.zside(pippo2);
	int newL = recHitTools.getLayerWithOffset(pippo2);
	std::cout << " zside = " << newZside << " newrhL = " << newL << std::endl;

	float rhX = recHitTools.getPosition(hitid).x();
	float rhY = recHitTools.getPosition(hitid).y();
	h2_dX_vs_layer->Fill(rhL, rhX);
	h2_dY_vs_layer->Fill(rhL, rhY);
	h_dX_Sci->Fill(rhX);
	h_dY_Sci->Fill(rhY);
	if( rhL == 43) h2_dX_vs_dY_Sci->Fill(rhX, rhY);
	*/
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
      for(unsigned int i = 0; i < NewrechitsBH.size(); ++i){
        hitmapSci[NewrechitsBH[i].detid()] = &NewrechitsBH[i];
      }

      break;
    }
  default:
    break;
  }

   std::cout << "hitmapSci.size = " << hitmapSci.size() << std::endl;
  // std::cout << "NewrechitsBH.size() = " << NewrechitsBH.size() << std::endl;


  ////////////////////
  float etaGen = -1;
  float phiGen = -1;
  float xGen = -1;
  float yGen = -1;
  float zGen = -1;
  bool evtGood = false;
  // float axX = -1;
  // float axY = -1;
  // float axZ = -1;


  SimClusterRefVector simClusterRefVectorChosen;
  std::map<DetId,float> hitsAndFractionChosen;

  if(debugCOUT) std::cout<< " >>> now caloparticles " << std::endl;

  std::array<double,3> vtx{ {vx, vy, vz}};

  float nHitsLayer[12] = {0.};

  int numbCaloPart = 0;
  if(debugCOUT) std::cout << " >>> caloParticles.size() = " << caloParticles.size() << std::endl;
  nCaloParticles->Fill(caloParticles.size());
  // loop over caloParticles
  for (std::vector<CaloParticle>::const_iterator it_caloPart = caloParticles.begin(); it_caloPart != caloParticles.end(); ++it_caloPart){
    ++nParticles;

    //for every caloparticle check if 
    //- non interacting or interacting within HGCAL
    //- if the particle is not from PU
    //reinitialize to false
    bool isParticlesGood = false;

    //vector of simClusters from the current caloparticle
    const SimClusterRefVector simClusterRefVector = it_caloPart->simClusters();

    nSimClusters->Fill(simClusterRefVector.size());
    //look at event with genParticle non interacting 
    if(debugCOUT)    std::cout << " caloParticles loop 1 simClusterRefVector.size() = " << simClusterRefVector.size() << " eta = " << it_caloPart->eta() << " pdgID = " 
				<< it_caloPart->pdgId() << " energy = " << it_caloPart->pt() << " eventId().event() = " << it_caloPart->eventId().event() 
				<< " eventId().bunchCrossing() = " << it_caloPart->eventId().bunchCrossing() << " Z vtx = " << vz << std::endl; 


    // in presence of pileup need selections to identify hard scattering particles
    //if(CaloPartPDGID == 22 && (simClusterRefVector.size() > 2 || it_caloPart->pdgId() != 22 || it_caloPart->pt() != particleGenPt)) continue;
      
    etaGen = it_caloPart->eta();
    phiGen = it_caloPart->phi();
    xGen = it_caloPart->momentum().x();
    yGen = it_caloPart->momentum().y();
    zGen = it_caloPart->momentum().z();


    h_eta->Fill(etaGen);
    h_phi->Fill(phiGen);

    // in presence of pileup need selections to identify hard scattering particles
    //events produced at eta > 0, reject if in other direction 
    //if(etaGen < 0) continue;

    ++numbCaloPart;
    //expect 1! caloparticle per event in single particle gun
    //    if(numbCaloPart > 1) continue;

    bool notPassedParticle = true;
    bool checkVtx = false;
    //loop over simclusters of the caloparticle
    for(const auto& scl : simClusterRefVector){
      if(debugCOUT)
	std::cout << " simcl energy = " << (*(scl)).energy() << " position " << (*(scl)).p4() << " pt = " << (*(scl)).pt() << " eta = " << (*(scl)).eta()<< std::endl;

      checkVtx = false;
      //check vtxZ of the tracks associated to the simclusters
      if(simClusterRefVector.size() == 1) {isParticlesGood = true; evtGood = true;}
      else{
      const std::vector<SimTrack>& myTRacks = (*(scl)).g4Tracks();
      if(debugCOUT) std::cout << " myTRacks.size() = " << myTRacks.size() << std::endl;
      for(const auto& itT : myTRacks){
	int vertexIndex = itT.vertIndex();
	int genPartIndex = itT.genpartIndex();
	unsigned int genPartIndexNosign = itT.genpartIndex();
	if(vertexIndex != -1 && vertexIndex < int(vtxs.size())){
	  notPassedParticle = false;
	  if(debugCOUT) std::cout << " ZsimTrack= " << vtxs[itT.vertIndex()].position().z() << std::endl;
	  zVtx_SimClustersTracks->Fill(std::abs(vtxs[itT.vertIndex()].position().z()));
	  checkVtx = true;
	  if(simClusterRefVector.size() == 2) zVtx_SimClustersTracks_EE->Fill(std::abs(vtxs[itT.vertIndex()].position().z()));
	  if(std::abs(vtxs[itT.vertIndex()].position().z()) >= 319.805) {
	    evtGood = true;	    
	    isParticlesGood = true;
	  }
	}
	else{std::cout << " >>> vertexIndex == -1 itT.pdgId() = " << itT.type()
		       << " eta = " << itT.momentum().eta() << " n simclusters = " <<  simClusterRefVector.size()
		       << std::endl;

	}
	if(checkVtx) break;
      }
      //if(checkVtx) break;
      }
    
      //get the collection of hits from every simcluster
      const SimCluster simClusterChosen = (*(scl));
      const std::vector<std::pair<uint32_t,float> > hits_and_fractions = simClusterChosen.hits_and_fractions();

      std::cout << " >> nSimHits = " << hits_and_fractions.size() << std::endl;
      for(std::vector<std::pair<uint32_t,float> >::const_iterator it_haf = hits_and_fractions.begin(); it_haf != hits_and_fractions.end(); ++it_haf) {
	DetId hitid = (it_haf->first);
	float rhEnergy = 0;
	std::map<DetId, const HGCRecHit*>::iterator trovatore = hitmapSci.find(hitid);
	if(trovatore == hitmapSci.end()){
	  continue;
	}
	else{
	  const HGCRecHit *hit = hitmapSci[hitid];

	  float rhEnergy = hit->energy();
	  float rhX = recHitTools.getPosition(hitid).x();
	  float rhY = recHitTools.getPosition(hitid).y();
	  int rhL = recHitTools.getLayerWithOffset(hitid);

	  if(rhL > 40) ++nHitsLayer[rhL - 40];
	}
      }   

      //loop over simHits associated to every simCluster
      //     {
      //...
      // }
    }//loop over simclusters


    if(isParticlesGood) {++nParticlesGood;
     std::cout << " particle good " << std::endl;
    }
    if(notPassedParticle) ++ noGOODParticle;
  }////loop over caloparticles

  for(int ij=0; ij<12; ++ij){
    if(nHitsLayer[ij] > 0)
      h_nhits_energyWeighted[ij]->Fill(nHitsLayer[ij]);
  }

  if(!evtGood) return;
  ++nEventsGood;

  return;

  if(debugCOUT) std::cout<< " caloparticles survived zGen = " << zGen << " xGen = " << xGen << " yGen = " << yGen << " etaGen = " << etaGen << std::endl;

  UtilClasses utilsMet = UtilClasses(etaGen, phiGen);
  std::array<double,3> fromAxis;
  if(CaloPartPDGID > 0) {
    fromAxis[0] = (xGen);
    fromAxis[1] = (yGen);
    fromAxis[2] = (zGen);
  }

  for(const reco::CaloCluster &sCl : *caloClusters){
    float clX = sCl.x();
    float clY = sCl.y();
    float clZ = sCl.z();
    float clEta = sCl.eta();
    float clPhi = sCl.phi();

    if(sCl.energy() == 0.) continue;
   
    if(clEta * etaGen < 0){
      if(debugCOUT) std::cout << " >>> opposite to genEta " << std::endl;
      continue;
    }
    

    const HGCalDetId clSeedId = sCl.hitsAndFractions()[0].first;
    int rhSeedL = recHitTools.getLayerWithOffset(clSeedId);

    std::array<double,3> to{ {0., 0., clZ} };
    utilsMet.layerIntersection(to, fromAxis, vtx);

    float deltaR = sqrt(pow(to[0] - clX, 2) + pow(to[1] - clY, 2));
    float clR = sqrt(clX*clX + clY*clY);

    h2_Z_vs_layer->Fill(rhSeedL, clZ);
    h2_deltaR_vs_layer->Fill(rhSeedL, clR, sCl.energy());

    h_Recoeta->Fill(clEta);
    h_Recophi->Fill(clPhi);


    h2_Z_vs_phi->Fill(clPhi, clZ);
    h2_Z_vs_eta->Fill(clEta, clZ);
    h2_eta_vs_phi->Fill(clPhi, clEta);
    //if(rhSeedL > 40 && clZ > 420.)
    h2_Geta_vs_phi->Fill(phiGen, etaGen);
  }








  std::cout << " nEvents = " << nEvents << std::endl;
  
  if(debugCOUT) std::cout<< " end of event process " << std::endl;
}

void
HGCalTiming2DClustersAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
HGCalTiming2DClustersAnalyzer::endJob()
{
  std::cout << " total events = " << nEvents << " goodEvents = " << nEventsGood << std::endl;
  std::cout << " total particles = " << nParticles << " goodParticles = " << nParticlesGood << "  noGOODParticle = " << noGOODParticle << std::endl;

}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HGCalTiming2DClustersAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HGCalTiming2DClustersAnalyzer);
