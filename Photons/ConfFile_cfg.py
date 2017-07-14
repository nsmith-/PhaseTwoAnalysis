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
options.parseArguments()


process = cms.Process("PhotonFilter")

# Log settings
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1

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

# run Puppi 
process.load('CommonTools/PileupAlgos/Puppi_cff')
process.load('CommonTools/PileupAlgos/PhotonPuppi_cff')
process.load('CommonTools/PileupAlgos/softKiller_cfi')
from CommonTools.PileupAlgos.PhotonPuppi_cff        import setupPuppiPhoton
from PhysicsTools.PatAlgos.slimming.puppiForMET_cff import makePuppies
makePuppies(process)
process.puSequence = cms.Sequence(process.pfNoLepPUPPI * process.puppiNoLep)

# producer
moduleName = "PatPhotonFilter"    
if (options.inputFormat.lower() == "reco"):
    moduleName = "RecoPhotonFilter"
process.photonfilter = cms.EDProducer(moduleName)
process.load("PhaseTwoAnalysis.Photons."+moduleName+"_cfi")

process.out = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('keep *_*_*_*',
                                           'drop patPhotons_slimmedPhotons_*_*',
                                           'drop recoPhotons_gedPhotons_*_*'),
    fileName = cms.untracked.string(options.outFilename)
)
  
if (options.inputFormat.lower() == "reco"):
    process.p = cms.Path(process.puSequence * process.photonfilter)
else:
    process.p = cms.Path(process.photonfilter)

process.e = cms.EndPath(process.out)
