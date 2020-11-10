import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

process = cms.Process("Demo")
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )


import FWCore.Utilities.FileUtils as FileUtils
readFiles = cms.untracked.vstring()
readFiles.extend(FileUtils.loadListFromFile ('./listTTBar_TTbar_0PU_11_2_0_pre6.txt') )


process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
                            fileNames = readFiles,
                            #fileNames = cms.untracked.vstring(
                            #    '/store/relval/CMSSW_11_0_0/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_110X_mcRun4_realistic_v3_2026D49PU200-v1/20000/EF473A21-7CDA-054D-B119-47A1F92A3EC4.root'
                            #),
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
        

#from HGCTimingAnalysis.HGCTiming.timeRecHitEstimator_cfi import HGCalTimeEstimator

process.ana = cms.EDAnalyzer('CheckCalibartion',
                             detector = cms.string("all"),
                             rawRecHits = cms.bool(True),                              
                             particleGENPT = cms.double(5.0),
                             CaloPartPDGID = cms.int32(22),
                             timeOffset = cms.double(5.0),
                             HGCEEInput = cms.InputTag('HGCalRecHit:HGCEERecHits'),
                             HGCFHInput = cms.InputTag('HGCalRecHit:HGCHEFRecHits'),
                             HGCBHInput = cms.InputTag('HGCalRecHit:HGCHEBRecHits'),
                             #caloClusterInput = cms.InputTag('hgcalLayerClusters'),
                             #dEdXweights = HGCalTimeEstimator.dEdXweights,
                             #thicknessCorrection = HGCalTimeEstimator.thicknessCorrection,
                             #HGCEE_fCPerMIP = HGCalTimeEstimator.HGCEE_fCPerMIP,
                             #HGCEE_noisefC = HGCalTimeEstimator.HGCEE_noisefC,
                             #HGCEF_noisefC = HGCalTimeEstimator.HGCEF_noisefC,
                             #HGCBH_noiseMIP = HGCalTimeEstimator.HGCBH_noiseMIP,
                             #nonAgedNoises = HGCalTimeEstimator.nonAgedNoises,
                             #endOfLifeNoises = HGCalTimeEstimator.endOfLifeNoises,
                             #HGCEE_keV2fC  = HGCalTimeEstimator.HGCEE_keV2fC,
                             #HGCHEF_keV2fC = HGCalTimeEstimator.HGCHEF_keV2fC,
                             #HGCHB_keV2MIP = HGCalTimeEstimator.HGCHB_keV2MIP
                             )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("file:calibrationTime_gen_nov2020.root")
                                   )
process.p = cms.Path(process.ana)
