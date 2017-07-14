import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.register('outFilename', 'FilteredEvents.root',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "Output file name"
                 )
options.register('inputFormat', 'PAT',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "format of the input files (PAT or RECO)"
                 )
options.register('skim', True,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "skim events with one lepton and 2 jets"
                 )
options.parseArguments()

process = cms.Process("EDMFilter")

# Geometry, GT, and other standard sequences
process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '91X_upgrade2023_realistic_v3', '')

# Log settings
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('MyAna')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
        limit = cms.untracked.int32(-1)
)

# Input
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) ) 

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(*(
        '/store/mc/PhaseIITDRSpring17MiniAOD/TTToSemiLepton_TuneCUETP8M1_14TeV-powheg-pythia8/MINIAODSIM/PU200_91X_upgrade2023_realistic_v3-v1/120000/008BFDF2-285E-E711-8055-001E674FC887.root',
    ))
)
if (options.inputFormat.lower() == "reco"):
    process.source.fileNames = cms.untracked.vstring(*(
        '/store/mc/PhaseIITDRSpring17DR/TTToSemiLepton_TuneCUETP8M1_14TeV-powheg-pythia8/AODSIM/PU200_91X_upgrade2023_realistic_v3-v1/120000/000CD008-7A58-E711-82DB-1CB72C0A3A61.root',
    ))
process.source.inputCommands = cms.untracked.vstring("keep *")

# run Puppi 
process.load('CommonTools/PileupAlgos/Puppi_cff')
process.load('CommonTools/PileupAlgos/PhotonPuppi_cff')
process.load('CommonTools/PileupAlgos/softKiller_cfi')
from CommonTools.PileupAlgos.PhotonPuppi_cff        import setupPuppiPhoton
from PhysicsTools.PatAlgos.slimming.puppiForMET_cff import makePuppies
makePuppies(process)

# recluster jets
process.load('RecoJets/Configuration/RecoPFJets_cff')
process.ak4PUPPIJets  = process.ak4PFJets.clone(rParam=0.4, src = cms.InputTag('puppi'))

# recompute MET
process.load('RecoMET.METProducers.PFMET_cfi')
process.puppiMet = process.pfMet.clone()
process.puppiMet.src = cms.InputTag('puppi')

process.puSequence = cms.Sequence(process.pfNoLepPUPPI * process.puppi * process.puppiNoLep * process.ak4PUPPIJets * process.puppiMet)


# electron producer
moduleElecName = "PatElectronFilter"    
if (options.inputFormat.lower() == "reco"):
    moduleElecName = "RecoElectronFilter"
process.electronfilter = cms.EDProducer(moduleElecName)
process.load("PhaseTwoAnalysis.Electrons."+moduleElecName+"_cfi")

# muon producer
moduleMuonName = "PatMuonFilter"    
if (options.inputFormat.lower() == "reco"):
    moduleMuonName = "RecoMuonFilter"
process.muonfilter = cms.EDProducer(moduleMuonName)
process.load("PhaseTwoAnalysis.Muons."+moduleMuonName+"_cfi")
if (options.inputFormat.lower() == "reco"):
    process.muonfilter.pfCandsNoLep = cms.InputTag("puppiNoLep")

# producer
moduleJetName = "PatJetFilter"    
if (options.inputFormat.lower() == "reco"):
    moduleJetName = "RecoJetFilter"
process.jetfilter = cms.EDProducer(moduleJetName)
process.load("PhaseTwoAnalysis.Jets."+moduleJetName+"_cfi")
if (options.inputFormat.lower() == "reco"):
    process.jetfilter.jets = "ak4PUPPIJets"    
        
# output
process.out = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('keep *_*_*_*',
                                           'drop patElectrons_slimmedElectrons_*_*',
                                           'drop recoGsfElectrons_gedGsfElectrons_*_*',
                                           'drop patMuons_slimmedMuons_*_*',
                                           'drop recoMuons_muons_*_*',
                                           'drop patJets_slimmedJetsPuppi_*_*',
                                           'drop reco*_ak4*Jets*_*_*',
                                           'drop recoPFMETs_pfMet_*_*',
                                           ),    
    fileName = cms.untracked.string(options.outFilename)
)

# run
if (options.inputFormat.lower() == "reco"):
    process.p = cms.Path(process.puSequence * process.electronfilter * process.muonfilter * process.jetfilter)
else:
    process.p = cms.Path(process.electronfilter * process.muonfilter * process.jetfilter)

process.e = cms.EndPath(process.out)
    
