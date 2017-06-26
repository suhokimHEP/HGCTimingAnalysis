import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

process = cms.Process("Demo")
process.load('Configuration.Geometry.GeometryExtended2023D13Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")
#process.load("RecoLocalCalo.HGCalRecProducers.hgcalLayerClusters_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


import FWCore.Utilities.FileUtils as FileUtils
readFiles = cms.untracked.vstring()
readFiles.extend(FileUtils.loadListFromFile ('INPUTFILELIST') )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
                            fileNames = readFiles,
                            duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
                            inputCommands=cms.untracked.vstring('keep *',
                                                                'drop EcalEBTriggerPrimitiveDigisSorted_simEcalEBTriggerPrimitiveDigis_*_HLT'
                                                                )
                            )
        

from HGCalCalibration.HitValidation.timeRecHitEstimator_cfi import HGCalTimeEstimator

process.ana = cms.EDAnalyzer('HGCalTimingAnalyzer',
                             detector = cms.string("all"),
                             rawRecHits = cms.bool(True),                              
                             CFDTimeCorrection = cms.int32(CFDVal),
                             cellType = cms.int32(CELLT),
                             floorValue = cms.double(FLOORV),
                             lifeAge = cms.int32(LIFEA),
                             absTrend = cms.double(ABSTREND),
                             HGCEEInput = cms.InputTag('HGCalRecHit:HGCEERecHits'),
                             HGCFHInput = cms.InputTag('HGCalRecHit:HGCHEFRecHits'),
                             HGCBHInput = cms.InputTag('HGCalRecHit:HGCHEBRecHits'),
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
                                   fileName = cms.string("file:OUTFILE.root")
                                   )
process.p = cms.Path(process.ana)
