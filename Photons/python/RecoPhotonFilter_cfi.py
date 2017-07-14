import FWCore.ParameterSet.Config as cms

photonfilter = cms.EDProducer('RecoPhotonFilter',
        photons    = cms.InputTag("gedPhotons"),
        beamspot     = cms.InputTag("offlineBeamSpot"),
        conversions  = cms.InputTag("particleFlowEGamma"),
        genParts     = cms.InputTag("genParticles"),
        vertices     = cms.InputTag("offlinePrimaryVertices"),
        puppiIsolationChargedHadrons = cms.InputTag("egmPhotonIsolationAODPUPPI","h+-DR030-BarVeto000-EndVeto001"),
        puppiIsolationNeutralHadrons = cms.InputTag("egmPhotonIsolationAODPUPPI","h0-DR030-BarVeto000-EndVeto000"),
        puppiIsolationPhotons        = cms.InputTag("egmPhotonIsolationAODPUPPI","gamma-DR030-BarVeto000-EndVeto008"),
)


from RecoEgamma.EgammaIsolationAlgos.egmPhotonIsolationPUPPI_cff import egmPhotonIsolationAODPUPPI

