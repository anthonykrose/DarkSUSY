// -*- C++ -*-
//
// Package:    FakeMuons
// Class:      FakeMuons
// 
/**\class FakeMuons FakeMuons.cc AnalysisAlgos/FakeMuons/src/FakeMuons.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  James Pivarski
//         Created:  Mon Dec 20 15:33:51 CST 2010
// $Id: FakeMuons.cc,v 1.1 2010/12/21 01:16:30 pivarski Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

//
// class declaration
//

class FakeMuons : public edm::EDProducer {
   public:
      explicit FakeMuons(const edm::ParameterSet&);
      ~FakeMuons();

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      bool sameTrack(const reco::Track *one, const reco::Track *two) const;

      // ----------member data ---------------------------
      edm::InputTag m_tracks;
      edm::InputTag m_muons;
      double m_minPt;
      double m_maxAbsEta;
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
FakeMuons::FakeMuons(const edm::ParameterSet& iConfig)
   : m_tracks(iConfig.getParameter<edm::InputTag>("tracks"))
   , m_muons(iConfig.getParameter<edm::InputTag>("muons"))
   , m_minPt(iConfig.getParameter<double>("minPt"))
   , m_maxAbsEta(iConfig.getParameter<double>("maxAbsEta"))
{
   //register your products
   produces<pat::MuonCollection>();

   //now do what ever other initialization is needed
}


FakeMuons::~FakeMuons()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

bool FakeMuons::sameTrack(const reco::Track *one, const reco::Track *two) const {
   return (fabs(one->px() - two->px()) < 1e-10  &&
	   fabs(one->py() - two->py()) < 1e-10  &&
	   fabs(one->pz() - two->pz()) < 1e-10  &&
	   fabs(one->vx() - two->vx()) < 1e-10  &&
	   fabs(one->vy() - two->vy()) < 1e-10  &&
	   fabs(one->vz() - two->vz()) < 1e-10);
}

// ------------ method called to produce the data  ------------
void
FakeMuons::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   edm::Handle<reco::TrackCollection> tracks;
   iEvent.getByLabel(m_tracks, tracks);

   edm::Handle<pat::MuonCollection> muons;
   iEvent.getByLabel(m_muons, muons);

   std::auto_ptr<pat::MuonCollection> fakeMuons(new pat::MuonCollection);

   for (unsigned int index = 0;  index < tracks->size();  index++) {
      reco::TrackRef trackRef(tracks, index);
      if (trackRef->pt() > m_minPt  &&  abs(trackRef->eta()) < m_maxAbsEta) {
	 const pat::Muon *theMuon = NULL;
	 for (pat::MuonCollection::const_iterator muon = muons->begin();  muon != muons->end();  ++muon) {
	    if (muon->innerTrack().isAvailable()  &&  sameTrack(&*(muon->innerTrack()), &*trackRef)) {
	       theMuon = &*muon;
	       break;
	    }
	 }

	 if (theMuon == NULL) {
	    reco::Muon newRecoMuon(trackRef->charge(), reco::Candidate::LorentzVector(trackRef->px(), trackRef->py(), trackRef->pz(), sqrt(pow(trackRef->p(), 2) + pow(0.105658367, 2))));
	    newRecoMuon.setInnerTrack(trackRef);
	    std::vector<reco::MuonChamberMatch> emptyMatches;
	    newRecoMuon.setMatches(emptyMatches);

	    fakeMuons->push_back(pat::Muon(newRecoMuon));
	    fakeMuons->back().embedTrack();
	 }
	 else {
	    fakeMuons->push_back(pat::Muon(*theMuon));
	    fakeMuons->back().embedTrack();
	 }
      }
   }

   iEvent.put(fakeMuons);
}

// ------------ method called once each job just before starting event loop  ------------
void 
FakeMuons::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
FakeMuons::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(FakeMuons);
