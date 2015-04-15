
# Auto generated configuration file
# using: 
# Revision: 1.20 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step3_diPhotonP10 --conditions auto:upgradePLS3 -n 100 --eventcontent FEVTDEBUGHLT -s RAW2DIGI,L1Reco,RECO --datatier GEN-SIM-RECO --customise SLHCUpgradeSimulations/Configuration/combinedCustoms.cust_2023HGCalMuon --geometry Extended2023HGCalMuon,Extended2023HGCalMuonReco --magField 38T_PostLS1 --filein file:step2_diPhotonPt10_100evts_1.root --fileout file:step3_diPhotonPt10_100evts_1.root --no_exec
import FWCore.ParameterSet.Config as cms

process = cms.Process('arborRECO')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2023HGCalMuonReco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring(
    'file:/afs/cern.ch/user/r/ruan/workspace/CMSSW_6_2_0_SLHC23/sample/ReReco/300GeV/step3_300GeV.root',
)
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.20 $'),
    annotation = cms.untracked.string('step3_diPhotonP10 nevts:100'),
    name = cms.untracked.string('Applications')
)

# Output definition

process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
    fileName = cms.untracked.string('file:step3_Arbor_CHS.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-RECO')
    )
)


# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)

#Arbor PFlow Candidate
process.particleFlowArbor = cms.EDProducer('HGCArborProducer')

process.FEVTDEBUGHLTEventContent.outputCommands.extend(['keep *_particleFlowArbor*_*_*'])

#ak4 ArborPFJet

process.ak4ArborPFJets=process.ak4PFJets.clone()
process.ak4ArborPFJets.src=cms.InputTag("particleFlowArbor")


#Arbor PFlowPtr for CHS

#Ptr
process.particleFlowArborPtrs=cms.EDProducer("PFCandidateFwdPtrProducer",
    src = cms.InputTag("particleFlowArbor")
)

#Identify PU particles
process.pfPileUpJMEArbor=process.pfPileUpJME.clone()
process.pfPileUpJMEArbor.PFCandidates= cms.InputTag("particleFlowArborPtrs")

#Idenfity noPU particles
process.pfNoPileUpJMEArbor=cms.EDProducer("TPPFCandidatesOnPFCandidates",
    bottomCollection = cms.InputTag("particleFlowArborPtrs"),
    enable = cms.bool(True),
    topCollection = cms.InputTag("pfPileUpJMEArbor"),
    name = cms.untracked.string('pileUpOnPFCandidates'),
    verbose = cms.untracked.bool(False)
)

process.ak4ArborPFJetsCHS=process.ak4PFJetsCHS.clone()
process.ak4ArborPFJetsCHS.src=cms.InputTag("pfNoPileUpJMEArbor")

process.FEVTDEBUGHLTEventContent.outputCommands.extend(['keep *_ak4ArborPFJets*_*_*'])

process.reconstruction_step = cms.Path(process.particleFlowArbor*process.ak4ArborPFJets*process.particleFlowArborPtrs*process.pfPileUpJMEArbor*process.pfNoPileUpJMEArbor*process.ak4ArborPFJetsCHS)

process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)

# Schedule definition
process.schedule = cms.Schedule(process.reconstruction_step,process.endjob_step,process.FEVTDEBUGHLToutput_step)

# customisation of the process.

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.combinedCustoms
from SLHCUpgradeSimulations.Configuration.combinedCustoms import cust_2023HGCalMuon 

#call to customisation function cust_2023HGCalMuon imported from SLHCUpgradeSimulations.Configuration.combinedCustoms
process = cust_2023HGCalMuon(process)

# End of customisation functions
