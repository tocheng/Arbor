#ifndef _TauAnalysis_h_
#define _TauAnalysis_h_

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Candidate/interface/Candidate.h"

/*
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
*/

#include "TH1F.h"
#include "TTree.h"

#include <string>

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// Arbor

//#include "FWCore/MessageLogger/interface/MessageLogger.h"
//#include "TMath.h"

class TauAnalysis : public edm::EDAnalyzer {
	public:
		explicit TauAnalysis(const edm::ParameterSet&);
		~TauAnalysis();

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


	private:
		virtual void beginJob() override;
		virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
		virtual void endJob() override;

		//virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
		//virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
		//virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
		//virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

		// ----------member data ---------------------------

		TTree *t_, *t_1;

		int run, lumi, event, nClu, nHEClu, nTau, TauDecayType, OneDType, NDaughter, NNeutrino;

		double TauEn, TauVisEn, CluEn, DisCluTau, TauEta, TauPhi; 
		double TauP[3], TauVisP[3], SumCluEn[7], DisClu[7];
		double CluP[3];

//		double PFOEn[4], PCluEn[4];
//		double PEta, PPhi;
//		int Charge;

		edm::Service<TFileService> *fs_;

};



#endif


