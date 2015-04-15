import FWCore.ParameterSet.Config as cms

#from RecoParticleFlow.PFClusterProducer.particleFlowRecHitHGCEE_cfi import *
#from RecoParticleFlow.PFClusterProducer.particleFlowRecHitHGCHEF_cfi import *
#from RecoParticleFlow.PFClusterProducer.particleFlowRecHitHGCHEB_cfi import *

process = cms.Process("HGCPFClustersAnalysis")

process.load('RecoParticleFlow.PFClusterProducer.particleFlowRecHitHGC_cff')

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')    
process.load('FWCore.MessageService.MessageLogger_cfi') 
process.load('Configuration.Geometry.GeometryExtended2023HGCalMuonReco_cff')
process.load('Configuration.Geometry.GeometryExtended2023HGCalMuon_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')

ffile=0
step=-1
preFix='SingleMuon_v15'

#configure from command line
#import sys
#if(len(sys.argv)>2):
#    preFix=sys.argv[2]
#    if(len(sys.argv)>3):
#        if(sys.argv[3].isdigit()) : ffile=int(sys.argv[3])
#    if(len(sys.argv)>4):
#        if(sys.argv[4].isdigit()) : step=int(sys.argv[4])


#process.particleFlowRecHitHGC = cms.Sequence( process.particleFlowRecHitHGCEE  +
#                                      process.particleFlowRecHitHGCHEF +
#                                      process.particleFlowRecHitHGCHEB   )

process.TFileService = cms.Service("TFileService", fileName = cms.string('Tau_HTT.root'))#/tmp/psilva/%s_SimHits_%d.root'%(preFix,ffile)))

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) ) 

#from UserCode.HGCanalysis.storeTools_cff import fillFromStore

process.source = cms.Source("PoolSource",                            
                            fileNames=cms.untracked.vstring()
                            )

#if preFix.find('RelVal')>=0 :
#    process.source.fileNames=fillFromStore('/store/relval/CMSSW_6_2_0_SLHC13/%s/GEN-SIM/DES23_62_V1_UPGHGCalMuon-v1/00000/'%preFix,ffile,step)
#else :
#    process.source.fileNames=fillFromStore('/store/cmst3/group/hgcal/CMSSW/%s'%preFix,ffile,step)

process.source.fileNames=cms.untracked.vstring(
# 'file:/afs/cern.ch/user/r/ruan/workspace/CMSSW_6_2_0_SLHC23/sample/step3.root '
# ' /store/cmst3/group/hgcal/CMSSW/Single211_CMSSW_6_2_0_SLHC23_patch1/mqRECO/Events_1.root '
# 'file:/afs/cern.ch/user/r/ruan/workspace/CMSSW_6_2_0_SLHC23/sample/PionScan/step3_drop.root',
# 'file:/afs/cern.ch/user/r/ruan/workspace/CMSSW_6_2_0_SLHC23/sample/ReReco/50GeV/step3_50GeV.root',
  '/store/cmst3/group/hgcal/CMSSW/RelValQQH1352T_Tauola_14TeV/PU0/Events_1.root',
  '/store/cmst3/group/hgcal/CMSSW/RelValQQH1352T_Tauola_14TeV/PU0/Events_2.root',
  '/store/cmst3/group/hgcal/CMSSW/RelValQQH1352T_Tauola_14TeV/PU0/Events_3.root',
  '/store/cmst3/group/hgcal/CMSSW/RelValQQH1352T_Tauola_14TeV/PU0/Events_4.root',
#  '/store/relval/CMSSW_6_2_0_SLHC20/RelValSingleGammaPt35Extended/GEN-SIM-RECO/PU_DES23_62_V1_HGCalV5PU140split-v1/00000/004AC59D-D96A-E411-AE50-02163E00FEA9.root',
#  '/store/relval/CMSSW_6_2_0_SLHC20/RelValSinglePiE50HCAL/GEN-SIM-RECO/DES23_62_V1_UPGHGCalV5-v1/00000/0A797444-AB5F-E411-98ED-002618943885.root',
 #'root://cms-xrd-global.cern.ch//store/path/to/file.root',
 #  'root://cms-xrd-global.cern.ch//store/user/mhaytmyr/HGCMC_Samples/step3_SingleGamma_E100p0_eta1p75_RAW2DIGI_L1Reco_RECO.root',
 # 'root://cms-xrd-global.cern.ch//store/relval/CMSSW_6_2_0_SLHC20/RelValSingleGammaPt35Extended/GEN-SIM-RECO/PU_DES23_62_V1_HGCalV5PU140split-v1/00000/004AC59D-D96A-E411-AE50-02163E00FEA9.root',
  #'root://cms-xrd-global.cern.ch//store/user/skalafut/HGCal/chgdPionGunFixedEAndEta_SLHC16/step3_chgdPionGun_E100_Eta1.75_2_200evts.root',
  #'root://cms-xrd-global.cern.ch//store/user/skalafut/HGCal/chgdPionGunFixedEAndEta_SLHC16/step3_chgdPionGun_E100_Eta1.75_3_200evts.root',
  #'root://cms-xrd-global.cern.ch//store/user/skalafut/HGCal/chgdPionGunFixedEAndEta_SLHC16/step3_chgdPionGun_E100_Eta1.75_4_200evts.root',
  #'root://cms-xrd-global.cern.ch//store/user/skalafut/HGCal/chgdPionGunFixedEAndEta_SLHC16/step3_chgdPionGun_E100_Eta1.75_5_200evts.root',

 )

process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

print process.source.fileNames

process.Taus = cms.EDAnalyzer("TauAnalysis")

process.p = cms.Path(process.Taus)

