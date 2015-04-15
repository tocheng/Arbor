// -*- C++ -*-
//
// Package:    HGCArborProducer
// Class:      HGCArborProducer
// 
/**\class HGCArborProducer HGCArborProducer.cc Arbor/HGCArborProducer/plugins/HGCArborProducer.cc
  
 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Manqi Ruan
//         Created:  Mon, 23 Feb 2015 18:05:43 GMT
// $Id$
//
//


// system include files
#include <memory>
#include <unordered_map>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementTrack.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementCluster.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/TruncatedPyramid.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "Geometry/EcalAlgo/interface/EcalEndcapGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalBarrelGeometry.h"
#include "Geometry/CaloTopology/interface/EcalEndcapTopology.h"
#include "Geometry/CaloTopology/interface/EcalBarrelTopology.h"

#include "Geometry/HcalTowerAlgo/interface/HcalGeometry.h"

// Addition for HGC geometry
 #include "Geometry/Records/interface/IdealGeometryRecord.h"
 #include "Geometry/CaloGeometry/interface/FlatTrd.h"
 #include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"
 #include "DataFormats/ForwardDetId/interface/HGCHEDetId.h"
 #include "Geometry/FCalGeometry/interface/HGCalGeometry.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

//#include "FastSimulation/BaseParticlePropagator/interface/BaseParticlePropagator.h"
//#include "MagneticField/Engine/interface/MagneticField.h"
//#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"


#include "RecoParticleFlow/PFClusterProducer/interface/InitialClusteringStepBase.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFraction.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "RecoParticleFlow/PFClusterProducer/interface/PFCPositionCalculatorBase.h"

#include "DataFormats/ForwardDetId/interface/ForwardSubdetector.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/TrackReco/interface/Track.h"

#include "DataFormats/TrajectorySeed/interface/PropagationDirection.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/MaterialEffects/interface/PropagatorWithMaterial.h"
#include "DataFormats/GeometrySurface/interface/BoundDisk.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "Geometry/FCalGeometry/interface/HGCalGeometry.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/ForwardDetId/interface/ForwardSubdetector.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/RootAutoLibraryLoader/interface/RootAutoLibraryLoader.h"

#include "TGClient.h"
#include "TVirtualX.h"
#include "TROOT.h"
#include "TVector3.h"
#include "TRint.h"
#include "TTree.h"

#include <fstream>
#include <iostream>
#include <cstdlib>

//
// class declaration
//


class HGCArborProducer : public edm::EDProducer {
   public:
      explicit HGCArborProducer(const edm::ParameterSet&);
      ~HGCArborProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      
      double CalibHit(double Z);

      int run, lumi, event;
      edm::Service<TFileService> *fs_;
      int NChPFO[3], NNeutralHad, NFrag, NPhoton, NPFNeutral, NLeakClu; 
      float ChEn[3], PhEn, FragEn, NeutralHadEn, PFNeutralEn, LeakCluEn; 
      TTree *t_;

      // debug
      float TrackEn[100]; float ClusterEn[100];

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
HGCArborProducer::HGCArborProducer(const edm::ParameterSet& iConfig)
{
   //register your products
/* Examples
   produces<ExampleData2>();

   //if do put with a label
   produces<ExampleData2>("label");
 
   //if you want to put into the Run
   produces<ExampleData2,InRun>();
*/
   //now do what ever other initialization is needed

   produces<reco::PFCandidateCollection>();  
//   produces<reco::PFBlockCollection>();
    
   edm::Service<TFileService> fs;
   fs_=&fs;

   run = 0; 
   lumi  = 0;
   event = 0;
   for(int s = 0; s < 3; s++)
   {
	   NChPFO[s] = 0;
	   ChEn[s] = 0;
   }
   NNeutralHad = 0;  NeutralHadEn = 0;
   NPhoton = 0; PhEn = 0; 
   NFrag = 0; FragEn = 0;
   NPFNeutral = 0; PFNeutralEn = 0;
   NLeakClu = 0; LeakCluEn = 0; 

   t_=fs->make<TTree>("HGCJet","Event Summary of Jet FD info");

   // event info

   t_->Branch("run",       &run,        "run/I");
   t_->Branch("lumi",      &lumi,       "lumi/I");
   t_->Branch("event",     &event,      "event/I");
   t_->Branch("NChPFO",    NChPFO,      "NChPFO[3]/I") ;
   t_->Branch("ChEn",      ChEn,        "ChEn[3]/F") ;

   t_->Branch("TrackEn",      TrackEn,        "TrackEn[100]/F") ;
   t_->Branch("ClusterEn",      ClusterEn,        "ClusterEn[100]/F") ;

   t_->Branch("NPhoton",   &NPhoton,    "NPhoton/I") ;
   t_->Branch("PhEn",      &PhEn,       "PhEn/F") ;
   t_->Branch("NFrag",     &NFrag,      "NFrag/I") ;
   t_->Branch("FragEn",    &FragEn,     "FragEn/F") ;
   t_->Branch("NPFNeutral",&NPFNeutral, "NPFNeutral/I");
   t_->Branch("PFNeutralEn", &PFNeutralEn, "PFNeutralEn/F");
   t_->Branch("NNeutralHad", &NNeutralHad, "NNeutralHad/I");
   t_->Branch("NeutralHadEn", &NeutralHadEn, "NeutralHadEn/F");
   t_->Branch("NLeakClu",  &NLeakClu, "NLeakClu/I");
   t_->Branch("LeakCluEn", &LeakCluEn,"LeakCluEn/F");
}




HGCArborProducer::~HGCArborProducer()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

double HGCArborProducer::CalibHit(double Z)
{
	double CalibConst = 0;

	if(Z < 320)
	{
		std::cout<<"Barrel Hit?"<<std::endl;
		CalibConst = 0;
	}
	else if(Z < 310)
	{
		CalibConst = 1.2*43.1;//1.0/55e-6 * 0.01 * 0.237
	}
	else if(Z < 329)
	{
		CalibConst = 1.2*155.1;
	}
	else if(Z < 339)
	{
		CalibConst = 1.2*185.3;
	}
	else if(Z < 350)
	{
		CalibConst = 1.2*241.3;
	}
	else if(Z < 356)
	{
		CalibConst = 1.2*765.4;
	}
	else if(Z < 410)
	{
		CalibConst = 1.2*618.2;
	}
	else
	{
		CalibConst = 1.2*103.6;
	}

	//return 1.2*CalibConst;

        return CalibConst;
}


// ------------ method called to produce the data  ------------
void HGCArborProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;
	using namespace std;

	run    = iEvent.id().run();
	lumi   = iEvent.luminosityBlock();
	event  = iEvent.id().event();

	for(int s = 0; s < 3; s++)
	{
		NChPFO[s] = 0;
		ChEn[s] = 0;
	}
	NNeutralHad = 0;  NeutralHadEn = 0; 
	NPhoton = 0; PhEn = 0; 
	NFrag = 0; FragEn = 0;
        NPFNeutral = 0; PFNeutralEn = 0;
	NLeakClu = 0; LeakCluEn = 0;

	std::auto_ptr<reco::PFCandidateCollection> pfCandidates_(new reco::PFCandidateCollection);
	//	std::auto_ptr<reco::PFBlockCollection> pfblocks(new reco::PFBlockCollection);

        // input tracks
	edm::Handle<reco::TrackCollection> tracks;
	iEvent.getByLabel("generalTracks", tracks);
	const reco::TrackCollection ptracks = *tracks;
	int NTrk = tracks->size();

        // input clusters
	edm::Handle<reco::PFClusterCollection> pfclustersHGCEE;
	iEvent.getByLabel("particleFlowClusterHGCEE", pfclustersHGCEE);
	const reco::PFClusterCollection pclusters = *pfclustersHGCEE;
	int NClu = pfclustersHGCEE->size();

        // rechits
	edm::Handle<reco::PFRecHitCollection> pfRecHitsHGC;
	iEvent.getByLabel("particleFlowRecHitHGCEE", pfRecHitsHGC);
	const reco::PFRecHitCollection phits = *pfRecHitsHGC;

	std::map<DetId, TVector3> DetID_2_Position;
	std::map<DetId, double> DetID_2_Energy;
	DetID_2_Position.clear();
	DetID_2_Energy.clear();

	std::map<int, double> TrkEn;
	std::map<int, double> CluEn;
	std::map<int, int> CluType; 
	TrkEn.clear();
	CluEn.clear();
	CluType.clear();

	for(unsigned i = 0; i<pfRecHitsHGC->size(); i++)
	{
		const reco::PFRecHit& a_hit = phits[i];
		DetId a_Id = a_hit.detId();
		const math::XYZPoint pos = a_hit.position();
		TVector3 a_pos;
		a_pos.SetXYZ(pos.x(), pos.y(), pos.z());
		DetID_2_Position[a_Id] = a_pos;

		double CalibConstSSS = HGCArborProducer::CalibHit(double(fabs(pos.z())));

		// DetID_2_Energy[a_Id] = a_hit.energy()*100;

		if( fabs(pos.z()) < 350 || (fabs(pos.z()) < 410 && a_hit.energy()*CalibConstSSS < 2) || a_hit.energy()*CalibConstSSS < 3 )
		{
			DetID_2_Energy[a_Id] = a_hit.energy()*CalibConstSSS;
		}
		else
		{
			DetID_2_Energy[a_Id] = 0.0000111;
		}
	}

	//std::cout<<NTrk<<std::endl; 
	/*
	   int pid = 22;
	   double charge = 0;
	   */

	std::vector<int> _usable_tracks;
	_usable_tracks.clear();
	_usable_tracks.reserve(NTrk);

	for(int i = 0; i < NTrk; i++)
	{
		_usable_tracks.push_back(i);
                if(i<100) TrackEn[i] = ptracks[i].pt();
	}

	std::sort(_usable_tracks.begin(),_usable_tracks.end(), [&](const unsigned i, const unsigned j)
			{
			return ptracks[i].pt() > ptracks[j].pt();
			});

	/*
	   for(const unsigned i : _usable_tracks)
	   {
	   const reco::Track& tk = ptracks[i];
	   cout<<tk.pt()<<endl; 
	   }
	   */

	cout<<"Track Prepared"<<endl;

	std::vector<int> _usable_clusters;
	_usable_clusters.clear();
	_usable_clusters.reserve(NClu);
	//double SumCluHitEn = 0;

	double EcalE = 0;
	double HFE = 0; 
	double HBE = 0;
	double HEndE = 0; 

	for(int i = 0; i < NClu; i++)
	{
		_usable_clusters.push_back(i);
		const reco::PFCluster& b_clu = pclusters[i];
		const std::vector<std::pair<DetId, float> >& clusterDetIds_b = b_clu.hitsAndFractions();
		double SumCluHitEn = 0;
		for(int t1 = 0; t1 < int( clusterDetIds_b.size() ); t1 ++)
		{
			TVector3 Hit_b = DetID_2_Position[clusterDetIds_b[t1].first];
			double currHitEn = DetID_2_Energy[clusterDetIds_b[t1].first];
			
			EcalE = 0; HFE = 0; HBE = 0; HEndE = 0; 

			// SumCluHitEn += currHitEn;
		
			if(fabs(Hit_b.Z()) < 350)
			{
				SumCluHitEn += currHitEn*0.8;
				EcalE += currHitEn*0.8;
			}
			else if(fabs(Hit_b.Z()) < 410 )
			{
				SumCluHitEn += currHitEn*1.2;	
				HFE += currHitEn*1.2;
			}
			else
			{
				SumCluHitEn += currHitEn*1.0;
				HBE += currHitEn*1.2;
				if( fabs(Hit_b.Z()) > 450 )
				{
					HEndE += currHitEn*1.2;
				}
			}
		}
		CluEn[i] = SumCluHitEn;

		if(EcalE/CluEn[i] > 0.9)	//Photon or Frag
		{
			CluType[i] = 22; 
		}
		else if(HEndE/CluEn[i] > 0.2)
		{
			CluType[i] = 101;	//Leakage
		}
		else 
		{
			CluType[i] = 211; 
		}

                if(i<100) ClusterEn[i] = SumCluHitEn;

	}

	std::sort(_usable_clusters.begin(),_usable_clusters.end(), [&](const unsigned i, const unsigned j)
			{
			return pclusters[i].energy() > pclusters[j].energy();
			});


	cout<<"Cluster Prepared"<<endl; 


	edm::ESHandle<MagneticField> _bField;
	edm::ESHandle<TrackerGeometry> _tkGeom;
	std::unique_ptr<PropagatorWithMaterial> _mat_prop;

	iSetup.get<IdealMagneticFieldRecord>().get(_bField);
	iSetup.get<TrackerDigiGeometryRecord>().get(_tkGeom);
	_mat_prop.reset( new PropagatorWithMaterial(alongMomentum, 0.1396, _bField.product()) );

	ReferenceCountingPointer<BoundDisk> boundECAL, boundHCALF, boundHCALB, _plusSurface_ECAL, _minusSurface_ECAL, _plusSurface_FHCAL, _minusSurface_FHCAL, _plusSurface_BHCAL, _minusSurface_BHCAL;

	Surface::RotationType rot;      // Hmmmm

	_minusSurface_ECAL = ReferenceCountingPointer<BoundDisk> ( new BoundDisk( Surface::PositionType(0,0,-320.4), rot, new SimpleDiskBounds( 0, 200, -0.001, 0.001)));
	_plusSurface_ECAL = ReferenceCountingPointer<BoundDisk> ( new BoundDisk( Surface::PositionType(0,0,+320.4), rot, new SimpleDiskBounds( 0, 200, -0.001, 0.001)));
	_minusSurface_FHCAL = ReferenceCountingPointer<BoundDisk> ( new BoundDisk( Surface::PositionType(0,0,-354.8), rot, new SimpleDiskBounds( 0, 250, -0.001, 0.001)));
	_plusSurface_FHCAL = ReferenceCountingPointer<BoundDisk> ( new BoundDisk( Surface::PositionType(0,0,+354.8), rot, new SimpleDiskBounds( 0, 250, -0.001, 0.001)));
	_minusSurface_BHCAL = ReferenceCountingPointer<BoundDisk> ( new BoundDisk( Surface::PositionType(0,0,-425.7), rot, new SimpleDiskBounds( 0, 300, -0.001, 0.001)));
	_plusSurface_BHCAL = ReferenceCountingPointer<BoundDisk> ( new BoundDisk( Surface::PositionType(0,0,+425.7), rot, new SimpleDiskBounds( 0, 300, -0.001, 0.001)));

	cout<<"Trajectory Prepared"<<endl; 

	TVector3 CluPos, TrackSecECAL, TrackSecHF, TrackSecHB, DiffCluSec;
	float DistanceTrajectory = 0;

	std::map<std::pair<unsigned, unsigned>, float> DisCluClu; 
	std::map<std::pair<unsigned, unsigned>, float> DisTrkClu;	
	DisCluClu.clear();
	DisTrkClu.clear();

	//int NHit_clu_a = 0;
	//int NHit_clu_b = 0;
	double MinDisBushes = 1.0E10;
	// double tmpDisBushes = 1.0E10;
	TVector3 Hit_a, Hit_b, diffAB;
	int count = 0; 

	for(const unsigned s : _usable_clusters)
	{
		const reco::PFCluster& a_clu = pclusters[s];
		//const std::vector<std::pair<DetId, float> >& clusterDetIds_a = a_clu.hitsAndFractions();
		//NHit_clu_a = clusterDetIds_a.size();
		count ++; 
		if(count%10 == 0)
			cout<<"Counts: "<<count<<endl; 

		for(const unsigned t : _usable_clusters)
		{
			// DisCluClu[s][t] = -1;
			if(s > t )
			{
				const reco::PFCluster& b_clu = pclusters[t];
				//const std::vector<std::pair<DetId, float> >& clusterDetIds_b = b_clu.hitsAndFractions();
				//NHit_clu_b = clusterDetIds_b.size();

				if( (a_clu.position().z() * b_clu.position().z()) > 0)
				{
					// MinDisBushes = 1.0E10;

					/*
					   for(int s1 = 0; s1 < NHit_clu_a; s1 ++)
					   {
					   Hit_a = DetID_2_Position[clusterDetIds_a[s1].first];           

					   for(int t1 = 0; t1 < NHit_clu_b; t1 ++)
					   {
					   Hit_b = DetID_2_Position[clusterDetIds_b[t1].first];
					   tmpDisBushes = (Hit_a - Hit_b).Mag();

					   if(tmpDisBushes < MinDisBushes)
					   {
					   MinDisBushes = tmpDisBushes;
					   }
					   }
					   }
					   */

					MinDisBushes = (a_clu.position().x() - b_clu.position().x())*(a_clu.position().x() - b_clu.position().x()) + (a_clu.position().y() - b_clu.position().y())*(a_clu.position().y() - b_clu.position().y());	

					if(MinDisBushes < 225)	//cm
					{
						std::pair<unsigned, unsigned> clupair, cluinvpair; 
						clupair.first = s; 
						clupair.second = t; 
						cluinvpair.first = t; 
						cluinvpair.second = s; 
						DisCluClu[clupair] = MinDisBushes;
						DisCluClu[cluinvpair] = MinDisBushes;
					}	
				}
			}
		}
	}

	for(const unsigned i : _usable_tracks)
	{
		const reco::Track& tk = ptracks[i];
		// /MinDis = 1.0E9;

		const TrajectoryStateOnSurface myTSOS = trajectoryStateTransform::outerStateOnSurface(tk, *(_tkGeom.product()),_bField.product());

		double currTrkEn = 0.5*(exp(tk.outerEta()) + exp(-1 * tk.outerEta()))*tk.pt();
		TrkEn[int(i)] = currTrkEn;

		if(myTSOS.globalPosition().z() > 0)
		{
			boundECAL = _plusSurface_ECAL;
			boundHCALF = _plusSurface_FHCAL;
			boundHCALB = _plusSurface_BHCAL;
		}
		else
		{
			boundECAL = _minusSurface_ECAL; //Only look at the entrance 
			boundHCALF = _minusSurface_FHCAL;
			boundHCALB = _minusSurface_BHCAL;
		}

		TrajectoryStateOnSurface piStateAtSurface_E = _mat_prop->propagate(myTSOS, *boundECAL);
		TrajectoryStateOnSurface piStateAtSurface_HF = _mat_prop->propagate(myTSOS, *boundHCALF);
		TrajectoryStateOnSurface piStateAtSurface_HB = _mat_prop->propagate(myTSOS, *boundHCALB);

		for(const unsigned j : _usable_clusters)
		{
			const reco::PFCluster& clu = pclusters[j];
			CluPos.SetXYZ(clu.position().x(), clu.position().y(), clu.position().z());

			if(piStateAtSurface_E.isValid())        //other 3 should be similar??
			{
				GlobalPoint pt = piStateAtSurface_E.globalPosition();
				TrackSecECAL.SetXYZ(pt.x(), pt.y(), pt.z());
			}
			if(piStateAtSurface_HF.isValid())
			{
				GlobalPoint ptHF = piStateAtSurface_HF.globalPosition();
				TrackSecHF.SetXYZ(ptHF.x(), ptHF.y(), ptHF.z());
			}
			if(piStateAtSurface_HB.isValid())
			{
				GlobalPoint ptHB = piStateAtSurface_HB.globalPosition();
				TrackSecHB.SetXYZ(ptHB.x(), ptHB.y(), ptHB.z());
			}

			if( abs(CluPos.Z()) > 425 )     //HBHits
			{
				DiffCluSec = CluPos - TrackSecHB;
			}
			else if(abs(CluPos.Z()) > 350 )
			{
				DiffCluSec = CluPos - TrackSecHF;
			}
			else
			{
				DiffCluSec = CluPos - TrackSecECAL;
			}

			if( piStateAtSurface_E.isValid() || piStateAtSurface_HF.isValid() || piStateAtSurface_HB.isValid() )
			{
				if(myTSOS.globalPosition().z() * clu.position().z() > 0)
				{
					DistanceTrajectory = DiffCluSec.Perp();
					if(DistanceTrajectory < 25)
					{
						std::pair<unsigned, unsigned> trkclupair; 
						trkclupair.first = i;
						trkclupair.second = j;
						DisTrkClu[trkclupair] = DistanceTrajectory;
					}
				}
			}
		}
	}

	cout<<"DisCluClu Matrix "<<DisCluClu.size()<<" : "<<DisTrkClu.size()<<endl;

	std::map<unsigned int, int> TypeFlag;           
	std::map<unsigned int, double> CluDepth;
	CluDepth.clear();

	int NType[3] = {0, 0, 0};

	std::map<int, int> TouchFlag;
	TouchFlag.clear();

	for(const unsigned i1 : _usable_clusters)
	{
		double MinDisToTrk = 1.0E9;
		double Depth = -1;
		double energy = CluEn[i1];

		if(pclusters[i1].position().z() > 0)
			Depth = pclusters[i1].position().z() - 320;
		else
			Depth = -320 - pclusters[i1].position().z();

		CluDepth[i1] = Depth;

		for(const unsigned p : _usable_tracks)
		{
			std::pair<unsigned, unsigned> qtrkclu; 
			qtrkclu.first = p; 
			qtrkclu.second = i1; 
			if(DisTrkClu.find(qtrkclu) != DisTrkClu.end())
			{
				if(DisTrkClu[qtrkclu] < MinDisToTrk )
				{
					MinDisToTrk = DisTrkClu[qtrkclu];
				}
			}
		}

		if(Depth < 2 && MinDisToTrk > 4 && energy > 2.0)        // To be protected
		{
			cout<<"A potential EM Cluster Candidate "<<endl;
			TypeFlag[i1] = 1;
			NType[0] ++;
		}
		else if(energy < 0.2 && Depth > 5.0)    //&& HitSize , Depth, FD
		{
			TypeFlag[i1] = 10;      //Fragments
			NType[1] ++;
		}
		else
		{
			TypeFlag[i1] = 12;      //Likely to be linked...
			NType[2] ++;
		}

		// cout<<"CluEn Com "<<CluEn[i1]<<" =?= "<<pclusters[i1].energy()<<endl; 
	}

	// cout<<endl<<"NTypes: "<<NType[0]<<" : "<<NType[1]<<" : "<<NType[2]<<" : "<<NClu<<endl<<endl;

	std::vector<std::pair<double,int>> ClusortTodistance;   //Find leading 1 - 2 clusters...
	std::vector<int> Index;
	std::vector<int> LinkedCluIndex;
	double currCluEn = 0;   //all clusterzed energy for given trk
	double a_cluEn = 0;

	/*
	   double ChEn = 0; 
	   double VetoCluEn = 0; 
	   double NeEn = 0; 
	   */

	for(const unsigned i : _usable_tracks)          //      Tag Core. 
	{
		ClusortTodistance.clear();
		Index.clear();
		LinkedCluIndex.clear();
		currCluEn = 0;

		const reco::Track& tk = ptracks[i];

		reco::TrackRef currTrkRef(tracks, i);

		// cout<<"Com Trk Ref "<<tk.pt()<<" =? "<<currTrkRef->pt()<<endl; 

		for(const unsigned j : _usable_clusters)
		{
			std::pair<unsigned, unsigned> ptrkclu; 
			ptrkclu.first = i;
			ptrkclu.second = j; 

			if( DisTrkClu.find(ptrkclu) != DisTrkClu.end() && TouchFlag.find(int(j)) == TouchFlag.end() && TypeFlag[j] != 1 ) // && CluDepth[j] < 15 ) // && TouchFlag.find(int(j)) == TouchFlag.end() )
			{
				std::pair<double, int> currPair;
				currPair.first = DisTrkClu[ptrkclu];
				currPair.second = j;
				ClusortTodistance.push_back(currPair); //*DisEn);
				Index.push_back(j);
			}
		}	

		std::sort(ClusortTodistance.begin(), ClusortTodistance.end(), [&](const std::pair<double,int> i1, const std::pair<double,int> i2)
				{
				return i1.first < i2.first;
				});

		for(unsigned int k = 0; k <  ClusortTodistance.size(); k++)     //Core Tagging
		{
			int index_tmp = ClusortTodistance[k].second;
			a_cluEn = CluEn[index_tmp];

			if( k == 0 && CluDepth[index_tmp] < 3 ) // && a_cluEn < TrkEn[i] + 3*sqrt(TrkEn[i]) + 2 ) //&& ClusortTodistance[k].first < 3 ) // || ClusortTodistance[k].first < 4 + 0.2*CluDepth[index_tmp]    )     
			{
				currCluEn += a_cluEn;
				LinkedCluIndex.push_back(index_tmp);
			}
			else if( k == 1 && currCluEn < 0.5*TrkEn[i] && CluDepth[index_tmp] < 4 )
			{
				currCluEn += a_cluEn;
				LinkedCluIndex.push_back(index_tmp);
			}
			else if( ClusortTodistance[k].first < 2 + 0.2*CluDepth[index_tmp] && currCluEn < TrkEn[i] + 2*sqrt(TrkEn[i]) )  //Energy Protection...
			{
				currCluEn += a_cluEn;
				LinkedCluIndex.push_back(index_tmp);
			}
		}	

		if( currCluEn < TrkEn[i] + 1.*sqrt(TrkEn[i]) + 2 )
		{
			for(unsigned int t = 0; t < LinkedCluIndex.size(); t++)
			{
				int CoreCluIndex = LinkedCluIndex[t];

				for(const unsigned u : _usable_clusters)        //should sort according to distance as well.
				{

					std::pair<unsigned, unsigned> currpair; 
					currpair.first = CoreCluIndex; 
					currpair.second = u; 

					if( DisCluClu.find(currpair) != DisCluClu.end())
					{
						if( DisCluClu[currpair] < 2 + 0.2*CluDepth[u] && TouchFlag.find(int(u)) == TouchFlag.end() && TypeFlag[u] != 1 && CluDepth[u] > 2  )
						{

							if( currCluEn + CluEn[u] < TrkEn[i] + 2*sqrt(TrkEn[i]) + 2 && find(LinkedCluIndex.begin(), LinkedCluIndex.end(), u) == LinkedCluIndex.end())
							{
								LinkedCluIndex.push_back(u);
								currCluEn += CluEn[u];
							}
						}
					}
				}
			}
		}

		//reco::PFCandidate ChargedCandidate; 
                
		if( fabs(currCluEn - TrkEn[i]) < 1.*sqrt(TrkEn[i]) + 2 ) //&& LinkedCluIndex.size() > 0 )   //&& 1 == 0 )
		{
			// charge = 1;	// Depend on Omega...

			NChPFO[0] ++;
			math::XYZTLorentzVector pf4Momentum( tk.px(), tk.py(), tk.pz(), sqrt( tk.px()*tk.px() + tk.py()*tk.py() + tk.pz()*tk.pz() ) );
			reco::PFCandidate thisPFCandidate(tk.charge(), pf4Momentum, reco::PFCandidate::h);
			//ChargedCandidate = thisPFCandidate;
			thisPFCandidate.setTrackRef(currTrkRef);
			pfCandidates_->push_back(thisPFCandidate);
			ChEn[0] += sqrt( tk.px()*tk.px() + tk.py()*tk.py() + tk.pz()*tk.pz() );

                        //if(NChPFO[0]-1<100) { TrackEn[NChPFO[0]-1] = tk.pt(); ClusterEn[NChPFO[0]-1] = currCluEn; }

		}
		else if ( currCluEn <= TrkEn[i] - sqrt(TrkEn[i]) - 2 && currCluEn > 0 )		// MIP, or crazy track, MIP need to be added.
		{
			NChPFO[1] ++;
			double wEn = (currCluEn + 0.5*LinkedCluIndex.size())/sqrt( tk.px()*tk.px() + tk.py()*tk.py() + tk.pz()*tk.pz() );
			math::XYZTLorentzVector pf4Momentum( wEn*tk.px(), wEn*tk.py(), wEn*tk.pz(), currCluEn + 0.5*LinkedCluIndex.size() );
			reco::PFCandidate thisPFCandidate(tk.charge(), pf4Momentum, reco::PFCandidate::h);
			//ChargedCandidate = thisPFCandidate;
			thisPFCandidate.setTrackRef(currTrkRef);
			pfCandidates_->push_back(thisPFCandidate);			
			ChEn[1] += currCluEn + 0.5*LinkedCluIndex.size();
		}
		else if ( currCluEn >= TrkEn[i] + sqrt(TrkEn[i]) + 2 )				//Enable PF
		{
			NChPFO[2] ++;
			math::XYZTLorentzVector pf4Momentum( tk.px(), tk.py(), tk.pz(), sqrt( tk.px()*tk.px() + tk.py()*tk.py() + tk.pz()*tk.pz() ) );
			reco::PFCandidate thisPFCandidate(tk.charge(), pf4Momentum, reco::PFCandidate::h);
			//ChargedCandidate = thisPFCandidate;
			thisPFCandidate.setTrackRef(currTrkRef);
			pfCandidates_->push_back(thisPFCandidate);
			ChEn[2] += sqrt( tk.px()*tk.px() + tk.py()*tk.py() + tk.pz()*tk.pz() );	// == TrkEn?

			double wEn = (currCluEn - TrkEn[i] + 0.5*LinkedCluIndex.size() )/sqrt( tk.px()*tk.px() + tk.py()*tk.py() + tk.pz()*tk.pz() );
			math::XYZTLorentzVector pf4MomentumNeutral( wEn*tk.px(), wEn*tk.py(), wEn*tk.pz(), currCluEn - TrkEn[i]+ 0.5*LinkedCluIndex.size() );
			reco::PFCandidate AOverlapNeutralCandidate(0, pf4MomentumNeutral, reco::PFCandidate::h0);
			pfCandidates_->push_back(AOverlapNeutralCandidate);
			NPFNeutral ++;
			PFNeutralEn += currCluEn - TrkEn[i] + 0.5*LinkedCluIndex.size();
		}

		//ChargedCandidate.setTrackRef(currTrkRef);
		for(unsigned s = 0; s < LinkedCluIndex.size(); s++)
		{
			TouchFlag[int(LinkedCluIndex[s])] = 3;
		}
		//pfCandidates_->push_back(ChargedCandidate);

	}

	int NLowEPhoton = 0; 
	int NHighEPhoton = 0;
	double SumLowEPhotonEn = 0; 
	double SumHighEPhotonEn = 0; 

	for(const unsigned i : _usable_clusters)
	{

		if( TouchFlag.find(i) == TouchFlag.end() ) 
		{
			// double wEn = (CluEn[i] + 0.5)
			// cout<<"Terrible strange "<<CluEn[i]<<" vs "<<pclusters[i].energy()<<endl; 

			double wEn = (CluEn[i] + 0.5)/sqrt(pclusters[i].position().x() * pclusters[i].position().x() + pclusters[i].position().y() * pclusters[i].position().y() + pclusters[i].position().z() * pclusters[i].position().z() );

			math::XYZTLorentzVector pf4Momentum( wEn*pclusters[i].position().x(), wEn*pclusters[i].position().y(), wEn*pclusters[i].position().z(), CluEn[i] + 0.5 );
			reco::PFCandidate thisPFCandidate(0, pf4Momentum, reco::PFCandidate::gamma); // reco::PFCandidate::X);
			pfCandidates_->push_back(thisPFCandidate);

			if(CluType[i] == 101)
			{
				NLeakClu++;
				LeakCluEn += CluEn[i] + 0.5; 
			}
			else if(CluEn[i] > 1)
			{
				NHighEPhoton++;
				SumHighEPhotonEn += CluEn[i] + 0.5;

				if(CluType[i] == 22)
				{
					NPhoton ++;
					PhEn += CluEn[i] + 0.5;
				}
				else
				{
					NNeutralHad++;
					NeutralHadEn += CluEn[i] + 0.5;
				}
			}
			else
			{
				NLowEPhoton++;
				SumLowEPhotonEn += CluEn[i] + 0.5;
				NFrag ++;
				FragEn += CluEn[i] + 0.5;
			}
			// NeEn += (CluEn[i] + 0.5);
		}
	}

	cout<<"Splitttttttttttttt???????????? "<<NLowEPhoton<<" vs "<<NHighEPhoton<<"  and En: L/H "<<SumLowEPhotonEn<<" vs "<<SumHighEPhotonEn<<endl; 

	// cout<<"Energy Partition, Veto vs Ch vs Ne "<<VetoCluEn<<" : "<<ChEn<<" : "<<NeEn<<endl; 

	iEvent.put(pfCandidates_);

	t_->Fill();

	/* This is an event example
	//Read 'ExampleData' from the Event
	Handle<ExampleData> pIn;
	iEvent.getByLabel("example",pIn);

	//Use the ExampleData to create an ExampleData2 which 
	// is put into the Event
	std::auto_ptr<ExampleData2> pOut(new ExampleData2(*pIn));
	iEvent.put(pOut);
	*/

	/* this is an EventSetup example
	//Read SetupData from the SetupRecord in the EventSetup
	ESHandle<SetupData> pSetup;
	iSetup.get<SetupRecord>().get(pSetup);
	*/

}

// ------------ method called once each job just before starting event loop  ------------
	void 
HGCArborProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HGCArborProducer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
/*
   void
   HGCArborProducer::beginRun(edm::Run const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when ending the processing of a run  ------------
/*
   void
   HGCArborProducer::endRun(edm::Run const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when starting to processes a luminosity block  ------------
/*
   void
   HGCArborProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when ending the processing of a luminosity block  ------------
/*
   void
   HGCArborProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HGCArborProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HGCArborProducer);
