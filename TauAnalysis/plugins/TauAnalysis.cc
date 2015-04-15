// -*- C++ -*-
//
// Package:    TauAnalysis
// Class:      TauAnalysis
// 
/**\class TauAnalysis TauAnalysis.cc Arbor/TauAnalysis/plugins/TauAnalysis.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Manqi Ruan
//         Created:  Wed, 25 Feb 2015 11:02:08 GMT
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "Arbor/TauAnalysis/plugins/TauAnalysis.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"



//
// class declaration
//

/*
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
};
*/

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
TauAnalysis::TauAnalysis(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

	edm::Service<TFileService> fs;
	fs_=&fs;

	// int run, lumi, event, nCluEE[2], nTau; 

	t_ = fs->make<TTree>("TauAnalyzer","Event Summary of PF info");

	t_->Branch("run",       &run,        "run/I");
        t_->Branch("lumi",      &lumi,       "lumi/I");
        t_->Branch("event",     &event,      "event/I");
        t_->Branch("nClu",    	&nClu,	     "nClu/I");	
	t_->Branch("nHEClu",    &nHEClu,     "nHEClu/I");
	t_->Branch("nTau",	&nTau, 	"nTau/I");	// At most 2 taus... always 2 taus. 
	t_->Branch("TauEn",     &TauEn ,  "TauEn/D");
	t_->Branch("TauDecayType", &TauDecayType, "TauDecayType/I");	//1, evv, 2, muvv, 3, 3prong, 4,h,  5, h+pi0, 5 h+...
	t_->Branch("NDaughter",     &NDaughter ,  "NDaughter/I");
	t_->Branch("NNeutrino",     &NNeutrino ,  "NNeutrino/I");
	t_->Branch("OneDType",     &OneDType ,  "OneDType/I");
	t_->Branch("TauVisEn",     &TauVisEn ,  "TauVisEn/D");
	t_->Branch("TauP",     TauP ,  "TauP[3]/D");
	t_->Branch("TauEta",     &TauEta ,  "TauEta/D");
	t_->Branch("TauPhi",     &TauPhi ,  "TauPhi/D");
	t_->Branch("TauVisP",   TauVisP ,  "TauVisP[3]/D");
	t_->Branch("CluEn",     &CluEn ,  "CluEn/D");
	t_->Branch("CluP",     CluP ,  "CluP[3]/D");
	t_->Branch("DisCluTau",     &DisCluTau ,  "DisCluTau/D");
	t_->Branch("SumCluEn",	SumCluEn, "SumCluEn[7]/D");
	t_->Branch("DisClu",  DisClu, "DisClu[7]/D");
}


TauAnalysis::~TauAnalysis()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
	void
TauAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;
	using namespace std;

	run    = iEvent.id().run();
        lumi   = iEvent.luminosityBlock();
        event  = iEvent.id().event();
	nHEClu = 0;
	nClu = 0; 
	nTau = 0; 
	TauEn = 0; 
	TauDecayType = -1; 
	TauVisEn = 0; 
	CluEn = 0;
	OneDType = 0;
	DisCluTau = 1.0E10;		//
	NNeutrino = 0;
	for(int i = 0; i < 3; i++)
	{
		CluP[i] = 0;
		TauP[i] = 0;
		TauVisP[i] = 0;
	}

/*
#ifdef THIS_IS_AN_EVENT_EXAMPLE
	Handle<ExampleData> pIn;
	iEvent.getByLabel("example",pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
	ESHandle<SetupData> pSetup;
	iSetup.get<SetupRecord>().get(pSetup);
#endif
*/
	edm::Handle<reco::PFClusterCollection> pfclustersHGCEE;
	iEvent.getByLabel("particleFlowClusterHGCEE", pfclustersHGCEE);
	const reco::PFClusterCollection pclusters = *pfclustersHGCEE;
	int nClu = pclusters.size();
	std::vector<int> HighEnergyClusters;
	HighEnergyClusters.clear();

	for(int i = 0; i < nClu; i++)
        {
		if( pclusters[i].energy() > 20 + 15*abs(pclusters[i].eta()) )
		{
			nHEClu++;
			HighEnergyClusters.push_back(i);
		}
	}

	int EventFlag = 0; 

	edm::Handle<edm::View<reco::Candidate> > genParticles;
	iEvent.getByLabel("genParticles", genParticles);

	int ChNeutrinoFl = 0; 
	int NCharge = 0; 
	double neutrinoEnergy = 0; 
	double neutrinoP[3] = {0, 0, 0}; 

	for(size_t i = 0; i < genParticles->size(); ++ i) 
	{
		const reco::GenParticle & p = dynamic_cast<const reco::GenParticle &>( (*genParticles)[i] );
		neutrinoEnergy = 0;

		if(abs(p.pdgId()) == 15 && p.status() == 2 && abs(p.eta()) < 3 && abs(p.eta()) > 1.5) // && p.status() == 2 )
		{
			nTau++;
			TauEn = p.energy();
			TauP[0]	= p.px();
			TauP[1] = p.py();
			TauP[2] = p.pz();
			TauEta = p.eta();
			TauPhi = p.phi();

			NNeutrino = 0;
			NDaughter = p.numberOfDaughters();

			for(int s = 0; s < NDaughter; s++)
			{
				const reco::Candidate* daug = p.daughter(s);
				cout<<s<<"th Daughter type: "<<daug->pdgId()<<" : "<<daug->energy()<<endl; 
				if(abs(daug->pdgId()) == 12 || abs(daug->pdgId()) == 14 || abs(daug->pdgId()) == 16  )
				{
					neutrinoEnergy += daug->energy();
					neutrinoP[0] += daug->px();
					neutrinoP[1] += daug->py();
					neutrinoP[2] += daug->pz();
					NNeutrino ++;
					if(  abs(daug->pdgId()) != 16 )
					{
						ChNeutrinoFl = abs(daug->pdgId());
					}
				}
				else
				{
					OneDType = daug->pdgId();
				}
				if(p.charge() != 0)
				{
					NCharge++;
				}
			}
			EventFlag++;
			TauVisEn = p.energy() - neutrinoEnergy;
			TauVisP[0] = TauP[0] - neutrinoP[0];			
			TauVisP[1] = TauP[1] - neutrinoP[1];
			TauVisP[2] = TauP[2] - neutrinoP[2];
	
			if(NCharge == 3)	//3 prong
			{
				TauDecayType = 3;
			}
			else if(NCharge > 3)	//Really??
			{
				TauDecayType = 10 * NCharge; 
			}
			else if(NCharge == 1)
			{
				if(NNeutrino == 2)
				{
					if(ChNeutrinoFl == 12)
					{
						TauDecayType = 1;	//eVV
					}
					else if(ChNeutrinoFl == 14)
					{
						TauDecayType = 2;	//mVV
					}
				}
				else if(NNeutrino == 1)	//
				{
					TauDecayType = 1 + NDaughter*10;	//Number of Final Sta...	
				}
			}

			DisCluTau = 1.0E10; 
			CluEn = -1; 
			double CluPositionMod = 0; 

			for(int j = 0; j < nHEClu; j++)		
			{
				const reco::PFCluster& clu = pclusters[  HighEnergyClusters[j] ];
				double dis = (clu.eta() - p.eta())*(clu.eta() - p.eta()) + (clu.phi() - p.phi())*(clu.phi() - p.phi()); 
				if(dis < DisCluTau)
				{
					DisCluTau = dis; 
					CluEn = clu.energy();
					CluPositionMod = sqrt(clu.position().x()*clu.position().x() + clu.position().y()*clu.position().y() + clu.position().z()*clu.position().z());
					CluP[0] = CluEn * clu.position().x()/CluPositionMod;
					CluP[1] = CluEn * clu.position().y()/CluPositionMod;
					CluP[2] = CluEn * clu.position().z()/CluPositionMod;
				}	
			}

			std::vector<std::pair<double,int>> ClusortTodistance;
			ClusortTodistance.clear();	

			for(int k = 0; k < nClu; k++)
			{
				const reco::PFCluster& clu = pclusters[k];
				if(clu.eta() * p.eta() > 0)
				{
					double dis = (clu.eta() - p.eta())*(clu.eta() - p.eta()) + (clu.phi() - p.phi())*(clu.phi() - p.phi());
					std::pair<double, int> currPair;
					currPair.first = dis;
					currPair.second = k;
					ClusortTodistance.push_back(currPair); //*DisEn);
				}
			}

			std::sort(ClusortTodistance.begin(), ClusortTodistance.end(), [&](const std::pair<double,int> i1, const std::pair<double,int> i2)
                                {
                                return i1.first < i2.first;
                                });

			double TmpSumCluEn = 0; 

			for(int s = 0; s < 7; s++)
			{
				SumCluEn[s] = 0;
				DisClu[s] = 0;
			}

			for(unsigned int t = 0; t < ClusortTodistance.size(); t++)
			{
				if(t > 6) break; 
				TmpSumCluEn += pclusters[ ClusortTodistance[t].second ].energy();
				SumCluEn[t] = TmpSumCluEn;
				DisClu[t] = ClusortTodistance[t].first; 
			}

			//Should tag the object within a small cone of the track... and calculate the mass. 
			//From matching to Cluster
			//Finding the closet Cluster to the Taudirection/TauVisDirection

			cout<<"Tau Tagged, Energy "<<p.energy()<<" visible Energy "<<TauVisEn <<endl;

			t_->Fill();

		}
	}

	cout<<iEvent.id().event()<<" has "<<nClu<<" clusters found, with Flag"<<EventFlag<<" visible Tau energy "<<TauVisEn<<endl;


}


// ------------ method called once each job just before starting event loop  ------------
	void 
TauAnalysis::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
	void 
TauAnalysis::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
   void 
   TauAnalysis::beginRun(edm::Run const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when ending the processing of a run  ------------
/*
   void 
   TauAnalysis::endRun(edm::Run const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when starting to processes a luminosity block  ------------
/*
   void 
   TauAnalysis::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when ending the processing of a luminosity block  ------------
/*
   void 
   TauAnalysis::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TauAnalysis::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TauAnalysis);
