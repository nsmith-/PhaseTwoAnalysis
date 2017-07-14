import FWCore.ParameterSet.Config as cms

photonfilter = cms.EDProducer('PatPhotonFilter',
        photons     = cms.InputTag("slimmedPhotons"),
        beamspot      = cms.InputTag("offlineBeamSpot"),
        conversions   = cms.InputTag("reducedEgamma", "reducedConversions", "PAT"),
)
