#ifndef _HGCPFJetFDAnalyzer_h_
#define _HGCPFJetFDAnalyzer_h_

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Candidate/interface/Candidate.h"

//PFBlock and PFBlockElement
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElement.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockFwd.h"

#include "DataFormats/ParticleFlowReco/interface/PFBlockElementCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementTrack.h"
//Track
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
// vertex
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
//PFCluster
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
//PFRecHit
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
//PFCanidate
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "DataFormats/JetReco/interface/GenJetCollection.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"

#include <string>

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//#include "RecoParticleFlow/PFClusterProducer/src/SimpleArborClusterizer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
//#include "RecoParticleFlow/PFClusterProducer/interface/Arbor.hh"
#include "TMath.h"

/**
   @class PFAnalyzer
   @author T. Cheng (CERN : IHEP, CAS)
*/

class HGCPFJetFDAnalyzer : public edm::EDAnalyzer 
{
  
 public:
  
  explicit HGCPFJetFDAnalyzer( const edm::ParameterSet& );
  ~HGCPFJetFDAnalyzer();
  virtual void analyze( const edm::Event&, const edm::EventSetup& );

  typedef math::XYZPoint Point;

 private:
    
  long int Index_Convertor( float PosX, float PosY, float PosZ, int CellSize );

    // cluster

     // cluster sorting
     
    double FD_PFJet(std::vector<reco::PFCandidate>);

    // geometry

    double _cellSizeEE, _layerThicknessEE;
    double _cellSizeHF, _layerThicknessHF;
    double _cellSizeHB, _layerThicknessHB;

    // debug level

    bool debug_;  
    void setDebug(bool isDebug){ debug_= isDebug; };
   
    // final outputs

    TTree *t_;

    // 


    int run, lumi, event;
 
    int ngen;

    float RecoJetEnergy, RecoJetPhi, RecoJetEta, RecoJetPt; 
    float GenJetEnergy, GenJetPhi, GenJetEta, GenJetPt; 
    int NRecoJet, NGenJet; 
    float GenEmE, GenHadE, RecoChE, RecoNeE, DeltaR;
    float SumGenJetEn, SumRecoJetEn;
    float EFrac[5]; 

    edm::Service<TFileService> *fs_;
    std::map<std::string,TH2F*> histContainer2D_;

};
 

#endif


