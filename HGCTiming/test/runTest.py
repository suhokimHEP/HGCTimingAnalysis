import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Phase2C9_cff import Phase2C9
process = cms.Process('Demo',Phase2C9)

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.register('inputFile', None, VarParsing.multiplicity.singleton, VarParsing.varType.string, "input list to ntuplize")
options.parseArguments()

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')

process.maxEvents = cms.untracked.PSet( 
    input = cms.untracked.int32(-1) 
)
import FWCore.Utilities.FileUtils as FileUtils
readFiles = cms.untracked.vstring()
readFiles.extend(FileUtils.loadListFromFile ('lists/%s.list'%options.inputFile) )
#readFiles.extend(FileUtils.loadListFromFile ('EOL.txt') )

process.MessageLogger.cerr.FwkReport.reportEvery = 20
process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
#fileNames = cms.untracked.vstring(
#'root://cms-xrd-global.cern.ch//store/group/dpg_hgcal/comm_hgcal/suhokim/old/RECO/startup_sn2.0/GEN_13Pt10_Vtx0_flatEta_1p5_1p8_26D49_499_startup_sn2.0_step3.root',
#'file:GEN_13Pt10_Vtx0_flatEta_1p5_1p8_26D49_50mu_startup_sn2.0_step3.root',
#'file:GEN_13Pt10_Vtx0_flatEta_1p5_1p8_26D49_None_startup_step3.root',
#                            ),
                            fileNames = readFiles,
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
#process.GlobalTag = GlobalTag(process.GlobalTag, '110X_mcRun4_realistic_v3', '')
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T15', '')


from HGCTimingAnalysis.HGCTiming.timeRecHitEstimator_cfi import HGCalTimeEstimator

process.ana = cms.EDAnalyzer('HGCalTimingAnalyzer',
                             detectorSciName = cms.string("HGCalHEScintillatorSensitive"), 
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
                             sciThicknessCorrection = HGCalTimeEstimator.sciThicknessCorrection,
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
                                   #fileName = cms.string("file:MinBias140PU_0.root")
                                   fileName = cms.string("file:%s.root"%options.inputFile)
                                   )
process.p = cms.Path(process.ana)
