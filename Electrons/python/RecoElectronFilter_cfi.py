import FWCore.ParameterSet.Config as cms

electronfilter = cms.EDProducer('RecoElectronFilter',
        electrons    = cms.InputTag("gedGsfElectrons"),
        beamspot     = cms.InputTag("offlineBeamSpot"),
        conversions  = cms.InputTag("particleFlowEGamma"),
        genParts     = cms.InputTag("genParticles"),
        vertices     = cms.InputTag("offlinePrimaryVertices"),
        puppiIsolationChargedHadrons = cms.InputTag("egmElectronIsolationAODPUPPI","h+-DR030-BarVeto000-EndVeto001"),
        puppiIsolationNeutralHadrons = cms.InputTag("egmElectronIsolationAODPUPPI","h0-DR030-BarVeto000-EndVeto000"),
        puppiIsolationPhotons        = cms.InputTag("egmElectronIsolationAODPUPPI","gamma-DR030-BarVeto000-EndVeto008"),
        puppiNoLeptonsIsolationChargedHadrons = cms.InputTag("egmElectronIsolationAODPUPPINoLep","h+-DR030-BarVeto000-EndVeto001"),
        puppiNoLeptonsIsolationNeutralHadrons = cms.InputTag("egmElectronIsolationAODPUPPINoLep","h0-DR030-BarVeto000-EndVeto000"),
        puppiNoLeptonsIsolationPhotons = cms.InputTag("egmElectronIsolationAODPUPPINoLep","gamma-DR030-BarVeto000-EndVeto008"),
)


from RecoEgamma.EgammaIsolationAlgos.egmElectronIsolationPUPPI_cff import egmElectronIsolationAODPUPPI
egmElectronIsolationAODPUPPINoLep = egmElectronIsolationAODPUPPI.clone(usePUPPINoLepton = cms.bool(True))

