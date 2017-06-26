import FWCore.ParameterSet.Config as cms


from RecoLocalCalo.HGCalRecProducers.HGCalRecHit_cfi import dEdX_weights, HGCalRecHit
from RecoLocalCalo.HGCalRecProducers.HGCalUncalibRecHit_cfi import HGCalUncalibRecHit
from SimCalorimetry.HGCalSimProducers.hgcalDigitizer_cfi import hgceeDigitizer, hgchefrontDigitizer, hgchebackDigitizer, nonAgedNoises, endOfLifeNoises

HGCalTimeEstimator = cms.PSet(
    dEdXweights = cms.vdouble(dEdX_weights),
    thicknessCorrection = cms.vdouble(HGCalRecHit.thicknessCorrection),
    HGCEE_fCPerMIP = cms.vdouble(HGCalUncalibRecHit.HGCEEConfig.fCPerMIP),
    
    HGCEE_noisefC = cms.vdouble(hgceeDigitizer.digiCfg.noise_fC),
    HGCEF_noisefC = cms.vdouble(hgchefrontDigitizer.digiCfg.noise_fC), ##same as EE
    HGCBH_noiseMIP = hgchebackDigitizer.digiCfg.noise_MIP, ##BH
    
    nonAgedNoises = cms.vdouble(nonAgedNoises),
    endOfLifeNoises = cms.vdouble(endOfLifeNoises),

    HGCEE_keV2fC  = hgceeDigitizer.digiCfg.keV2fC,
    HGCHEF_keV2fC = hgchefrontDigitizer.digiCfg.keV2fC,
    HGCHB_keV2MIP = hgchebackDigitizer.digiCfg.keV2MIP ##BH
    )



