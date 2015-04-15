#include "Arbor/TauAnalysis/plugins/HGCPFJetFDAnalyzer.h"

#include "DataFormats/ParticleFlowReco/interface/PFBlockFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/Common/interface/RefToPtr.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/JetReco/interface/PFJet.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DetectorDescription/OfflineDBLoader/interface/GeometryInfoDump.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "SimG4CMS/Calo/interface/CaloHitID.h"

#include "DetectorDescription/Core/interface/DDFilter.h"
#include "DetectorDescription/Core/interface/DDFilteredView.h"
#include "DetectorDescription/Core/interface/DDSolid.h"

#include "DataFormats/GeometrySurface/interface/TrapezoidalPlaneBounds.h"

#include "DataFormats/GeometryVector/interface/Basic3DVector.h"

#include "CLHEP/Units/GlobalSystemOfUnits.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Geometry/Plane3D.h"
#include "CLHEP/Geometry/Vector3D.h"

#include "TVector2.h"

#include <iostream>

#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "Fireworks/Core/interface/fwLog.h"

#include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCHEDetId.h"

// vertex
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/ForwardDetId/interface/ForwardSubdetector.h"

#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

using namespace std;

//
HGCPFJetFDAnalyzer::HGCPFJetFDAnalyzer( const edm::ParameterSet &iConfig )
{
/*
  _cellSizeEE = iConfig.getParameter<double>("cellSizeEE");
  _layerThicknessEE = iConfig.getParameter<double>("layerThicknessEE");

  //////

  _cellSizeHF = iConfig.getParameter<double>("cellSizeHF");
  _layerThicknessHF = iConfig.getParameter<double>("layerThicknessHF");

  //////

  _cellSizeHB = iConfig.getParameter<double>("cellSizeHB");
  _layerThicknessHB = iConfig.getParameter<double>("layerThicknessHB");
*/

  edm::Service<TFileService> fs;
  fs_=&fs;

  t_=fs->make<TTree>("HGCJetFDAnalyzer","Event Summary of Jet FD info");

  // event info
    
  t_->Branch("run",       &run,        "run/I");
  t_->Branch("lumi",      &lumi,       "lumi/I");
  t_->Branch("event",     &event,      "event/I");

  // gen P info

  t_->Branch("GenJetEta", &GenJetEta, "GenJetEta/F");
  t_->Branch("GenJetPhi", &GenJetPhi, "GenJetPhi/F");
  t_->Branch("GenJetEnergy", &GenJetEnergy, "GenJetEnergy/F");
  t_->Branch("GenJetPt", &GenJetPt, "GenJetPt/F");  
  t_->Branch("GenEmE", &GenEmE, "GenEmE/F");
  t_->Branch("GenHadE", &GenHadE, "GenHadE/F");
  t_->Branch("RecoJetEta", &RecoJetEta, "RecoJetEta/F");
  t_->Branch("RecoJetPhi", &RecoJetPhi, "RecoJetPhi/F");
  t_->Branch("RecoJetEnergy", &RecoJetEnergy, "RecoJetEnergy/F");
	
  t_->Branch("SumGenJetEn", &SumGenJetEn, "SumGenJetEn/F");
  t_->Branch("SumRecoJetEn", &SumRecoJetEn, "SumRecoJetEn/F");
  t_->Branch("RecoJetPt", &RecoJetPt, "RecoJetPt/F");  
  t_->Branch("RecoChE", &RecoChE, "RecoChE/F");
  t_->Branch("EFrac", EFrac, "EFrac[5]/F");
  t_->Branch("RecoNeE", &RecoNeE, "RecoNeE/F");

  t_->Branch("NRecoJet", &NRecoJet, "NRecoJet/I");
  t_->Branch("NGenJet", &NGenJet, "NGenJet/I");

  t_->Branch("DeltaR", &DeltaR, "DeltaR/F");
}

//
HGCPFJetFDAnalyzer::~HGCPFJetFDAnalyzer()
{

}


//
void HGCPFJetFDAnalyzer::analyze( const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
    bool isDebug = false;

    setDebug(isDebug);

    if(debug_) cout << "[HGCPFAnalyzer::analyze] start" << endl;
  
    //event header
    run    = iEvent.id().run();
    lumi   = iEvent.luminosityBlock();
    event  = iEvent.id().event();

    //get gen level information

    // PFJets
    edm::Handle<reco::PFJetCollection> AK4PFJetsCHS;
    iEvent.getByLabel("ak4ArborPFJets", "", AK4PFJetsCHS);    
    
    edm::Handle<reco::GenJetCollection> AK4GenJets;
    iEvent.getByLabel("ak4GenJets", "", AK4GenJets); 

    //Leading Gen Jet; Matched PF Jet

    GenJetEnergy = 0;
    GenJetEta = -1000;
    GenJetPhi = -1000;
    GenJetPt = -1000; 
    GenEmE = 0; 
    GenHadE = 0;  
    NRecoJet = 0;
    NGenJet = 0;
    RecoJetEta = -1000;
    RecoJetPhi = -1000;
    RecoJetEnergy = -1000;
    RecoJetPt = -1000; 
    RecoChE = 0; 
    RecoNeE = 0; 
    DeltaR = -1000; 
    SumGenJetEn = 0; 
    SumRecoJetEn = 0; 
    for(int p = 0; p < 5; p++)
	{
		EFrac[p] = 0;
	}

    for(reco::GenJetCollection::const_iterator jet_it = AK4GenJets->begin(); jet_it != AK4GenJets->end(); ++ jet_it){

	if( abs(jet_it->eta()) < 2.8 && abs(jet_it->eta()) > 1.6 )
	{
		NGenJet ++;
		if(jet_it->energy() > GenJetEnergy)
		{
			GenJetEta = jet_it->eta();
			GenJetPhi = jet_it->phi();
			GenJetEnergy = jet_it->energy();
			GenJetPt = jet_it->pt();
			GenEmE = jet_it->emEnergy();
			GenHadE = jet_it->hadEnergy(); 
		}
		SumGenJetEn += jet_it->energy();
	}
    }
	
    for(reco::PFJetCollection::const_iterator jet_it = AK4PFJetsCHS->begin(); jet_it != AK4PFJetsCHS->end(); ++ jet_it){

        if( abs(jet_it->eta()) < 2.8 && abs(jet_it->eta()) > 1.6 )
        {
                NRecoJet ++;

		if(jet_it->energy() > RecoJetEnergy)
		{
			RecoJetEta = jet_it->eta();
			RecoJetPhi = jet_it->phi();
			RecoJetPt = jet_it->pt();
			RecoJetEnergy = jet_it->energy();
			EFrac[0] = jet_it->chargedHadronEnergy();
			EFrac[1] = jet_it->neutralHadronEnergy();
			EFrac[2] = jet_it->photonEnergy();
			EFrac[3] = jet_it->electronEnergy();
			EFrac[4] = jet_it->muonEnergy();
		}

		SumRecoJetEn += jet_it->energy();

	}
    }

    DeltaR = sqrt( (RecoJetEta - GenJetEta)*(RecoJetEta - GenJetEta) + (RecoJetPhi - GenJetPhi)*(RecoJetPhi - GenJetPhi) );

    if(debug_) std::cout<<"start filling tree"<<endl;

    t_->Fill();

    if(debug_) std::cout<<"end filling tree"<<endl;

}

double HGCPFJetFDAnalyzer::FD_PFJet(std::vector<reco::PFCandidate> hits)
{

	std::map <long int, float> idmap;   //id map to Energy      

	int NHitsOri = hits.size();

	float FD = 0;

	int Scales[10] = {2, 3, 4, 5, 6, 7, 8, 9, 10, 15};

	for(int i = 0; i < 10; i++)
	{
		idmap.clear();
		for( int i_hit = 0; i_hit < NHitsOri; i_hit++  )
		{

			int currhitIndex = HGCPFJetFDAnalyzer::Index_Convertor( hits[i_hit].eta(), hits[i_hit].phi(), 0.0, Scales[i] ); // map to eta-phi plane

			if(idmap.find(currhitIndex) == idmap.end())
			{
				idmap[currhitIndex] = hits[i_hit].energy();
			}
			else
			{
				idmap[currhitIndex] += hits[i_hit].energy();
			}

		}

		if(debug_)  std::cout<<NHitsOri<<" NHits at Scale "<<Scales[i]<<" is "<<idmap.size()<<std::endl;
		FD += 0.1 * TMath::Log(float(NHitsOri)/idmap.size())/TMath::Log(Scales[i]);
	}

	//if(debug_) std::cout<<"FD: "<<FD<<std::endl<<std::endl;

	return FD;

}

long int HGCPFJetFDAnalyzer::Index_Convertor( float PosX, float PosY, float PosZ, int CellSize )
{
	int Index_X = int( (PosX + 1000) ) / CellSize;
	int Index_Y = int( (PosY + 1000) ) / CellSize;
	int Index_Z = int( (PosZ + 1000) );

	long int TotalIndex = Index_Z*100000000 + Index_X*10000 + Index_Y;

	return TotalIndex;
}

DEFINE_FWK_MODULE(HGCPFJetFDAnalyzer);
