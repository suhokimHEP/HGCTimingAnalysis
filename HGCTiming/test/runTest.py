import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Phase2C9_cff import Phase2C9
process = cms.Process('Demo',Phase2C9)

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')

process.maxEvents = cms.untracked.PSet( 
    input = cms.untracked.int32(15) 
)


#import FWCore.Utilities.FileUtils as FileUtils
#readFiles = cms.untracked.vstring()
#readFiles.extend(FileUtils.loadListFromFile ('INPUTFILELIST') )


process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
                            #fileNames = cms.untracked.vstring('file:/eos/cms/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/NoPU_111X_mcRun4_realistic_T15_v1-v1/100000/66CDE37C-EE86-0A48-BC19-063DDD5BC953.root'),
                            fileNames = cms.untracked.vstring('file:/eos/cms/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/QCD_Pt-15to3000_TuneCP5_Flat_14TeV-pythia8/FEVT/NoPU_castor_all_pt_tracks_111X_mcRun4_realistic_T15_v1-v1/100000/01307FF4-AB5E-FE49-BA98-5D1D4A60A59F.root'),
                            #fileNames = readFiles,
                            duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
                            inputCommands=cms.untracked.vstring('keep *',
                                                  #'drop EcalEBTriggerPrimitiveDigisSorted_simEcalEBTriggerPrimitiveDigis_*_HLT',
                                                  'drop l1tEMTFHit2016Extras_simEmtfDigis_CSC_HLT',
                                                  'drop l1tEMTFHit2016Extras_simEmtfDigis_RPC_HLT',
                                                  'drop l1tEMTFHit2016s_simEmtfDigis__HLT',
                                                  'drop l1tEMTFTrack2016Extras_simEmtfDigis__HLT',
                                                  'drop l1tEMTFTrack2016s_simEmtfDigis__HLT'
                                                  )
                            )
        
# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '110X_mcRun4_realistic_v3', '')

from HGCTimingAnalysis.HGCTiming.timeRecHitEstimator_cfi import HGCalTimeEstimator

process.ana = cms.EDAnalyzer('HGCalTimingAnalyzer',
                             detector = cms.string("all"),
                             rawRecHits = cms.bool(True),                              
                             particleGENPT = cms.double(5),
                             CaloPartPDGID = cms.int32(22),
                             timeOffset = HGCalTimeEstimator.timeOFFSET,
                             HGCEEInput = cms.InputTag('HGCalRecHit:HGCEERecHits'),
                             HGCFHInput = cms.InputTag('HGCalRecHit:HGCHEFRecHits'),
                             HGCBHInput = cms.InputTag('HGCalRecHit:HGCHEBRecHits'),
                             caloClusterInput = cms.InputTag('hgcalLayerClusters'),
                             dEdXweights = HGCalTimeEstimator.dEdXweights,
                             thicknessCorrection = HGCalTimeEstimator.thicknessCorrection,
                             HGCEE_fCPerMIP = HGCalTimeEstimator.HGCEE_fCPerMIP,
                             HGCEE_noisefC = HGCalTimeEstimator.HGCEE_noisefC,
                             HGCEF_noisefC = HGCalTimeEstimator.HGCEF_noisefC,
                             HGCBH_noiseMIP = HGCalTimeEstimator.HGCBH_noiseMIP,
                             nonAgedNoises = HGCalTimeEstimator.nonAgedNoises,
                             endOfLifeNoises = HGCalTimeEstimator.endOfLifeNoises,
                             HGCEE_keV2fC  = HGCalTimeEstimator.HGCEE_keV2fC,
                             HGCHEF_keV2fC = HGCalTimeEstimator.HGCHEF_keV2fC,
                             HGCHB_keV2MIP = HGCalTimeEstimator.HGCHB_keV2MIP
                             )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("file:testMinBias.root")
                                   )
process.p = cms.Path(process.ana)
