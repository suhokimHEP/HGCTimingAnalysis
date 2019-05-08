import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

process = cms.Process('DEMO',eras.Phase2C4)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryExtended2023D28Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')

process.source = cms.Source("PoolSource",
                            #fileNames = readFiles,
                            fileNames = cms.untracked.vstring(
'file:/eos/cms/store/group/dpg_hgcal/comm_hgcal/amartell/HGCAL_Timing_vtxStudies/RECO_22Pt10_VtxZ320_LowEta_noTSmear_106/RECO_22Pt10_VtxZ320_LowEta_noTSmear_106_106.root'
        ), 
                            secondaryFileNames = cms.untracked.vstring(),
                            duplicateCheckMode = cms.untracked.string("noDuplicateCheck")
                            )



process.maxEvents = cms.untracked.PSet( 
    input = cms.untracked.int32(10) 
)


process.ana = cms.EDAnalyzer('HGCalTiming2DClustersAnalyzer',
                             detector = cms.string("all"),
                             particleGENPT = cms.double(50),
                             sigmaTIME = cms.double(0.03),
                             CaloPartPDGID = cms.int32(22),
                             timeOffset = cms.int32(5),
                             HGCEEInput = cms.InputTag('HGCalRecHit:HGCEERecHits'),
                             HGCFHInput = cms.InputTag('HGCalRecHit:HGCHEFRecHits'),
                             HGCBHInput = cms.InputTag('HGCalRecHit:HGCHEBRecHits'),
                             caloClusterInput = cms.InputTag('hgcalLayerClusters'),
                             timeMapInput = cms.InputTag('hgcalLayerClusters:timeLayerCluster'),
                             )



process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("file:pippo.root")
                                   )
process.p = cms.Path(process.ana)
