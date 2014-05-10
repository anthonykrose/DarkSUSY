// -*- C++ -*-
//
// Package:    MuJetPlusFakeMuon
// Class:      MuJetPlusFakeMuon
// 
/**\class MuJetPlusFakeMuon MuJetPlusFakeMuon.cc AnalysisAlgos/MuJetPlusFakeMuon/src/MuJetPlusFakeMuon.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Jim Pivarski,,,
//         Created:  Mon Dec  6 12:13:45 CST 2010
// $Id: MuJetPlusFakeMuon.cc,v 1.2 2010/12/08 14:50:12 pivarski Exp $
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

#include "AnalysisDataFormats/MuJetAnalysis/interface/MultiMuon.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "TRandom3.h"
#include "TMath.h"

//
// class declaration
//

class MuJetPlusFakeMuon : public edm::EDProducer {
   public:
      explicit MuJetPlusFakeMuon(const edm::ParameterSet&);
      ~MuJetPlusFakeMuon();

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      enum {
	 kGroupByDeltaR,
	 kGroupByMass,
	 kGroupByVertexProb,
	 kGroupByDeltaRAndVertexProb,
	 kGroupByMassAndVertexProb,
	 kGroupByDeltaROrMass,
	 kGroupByDeltaROrMassAndVertexProb,
	 kGroupByMassAndVertexProbOrDeltaR,

	 kGroupByAnyCharge,
	 kGroupByOppositeCharge,
	 kGroupBySameCharge
      };

      bool trackOkay(const reco::Track &track);

      // ----------member data ---------------------------

      unsigned int m_numberOfFakes;
      edm::InputTag m_muJets;
      edm::InputTag m_muons;
      edm::InputTag m_tracks;
      edm::InputTag m_caloTowers;
      double m_minPt;
      double m_minPmag;
      double m_maxAbsEta;
      int m_minTrackerHits;
      double m_maxTrackerNormChi2;
      double m_maxTrackerDxy;
      double m_maxTrackerDz;
      double m_maxQoverpError;
      double m_maxPhiError;
      double m_maxEtaError;
      double m_maxDxyError;
      double m_maxDzError;
      bool m_calculateVertex;
      bool m_calculateIsolation;
      std::string m_groupingMode_string;
      int m_groupingMode;
      double m_maxDeltaR;
      double m_maxMass;
      double m_minVertexProb;
      std::string m_groupByCharge_string;
      int m_groupByCharge;

      TRandom3 m_trandom3;
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
MuJetPlusFakeMuon::MuJetPlusFakeMuon(const edm::ParameterSet& iConfig)
   : m_numberOfFakes(iConfig.getParameter<unsigned int>("numberOfFakes"))
   , m_muJets(iConfig.getParameter<edm::InputTag>("muJets"))
   , m_muons(iConfig.getParameter<edm::InputTag>("muons"))
   , m_tracks(iConfig.getParameter<edm::InputTag>("tracks"))
   , m_caloTowers(iConfig.getParameter<edm::InputTag>("caloTowers"))
   , m_minPt(iConfig.getParameter<double>("minPt"))
   , m_minPmag(iConfig.getParameter<double>("minPmag"))
   , m_maxAbsEta(iConfig.getParameter<double>("maxAbsEta"))
   , m_minTrackerHits(iConfig.getParameter<int>("minTrackerHits"))
   , m_maxTrackerNormChi2(iConfig.getParameter<double>("maxTrackerNormChi2"))
   , m_maxTrackerDxy(iConfig.getParameter<double>("maxTrackerDxy"))
   , m_maxTrackerDz(iConfig.getParameter<double>("maxTrackerDz"))
   , m_maxQoverpError(iConfig.getParameter<double>("maxQoverpError"))
   , m_maxPhiError(iConfig.getParameter<double>("maxPhiError"))
   , m_maxEtaError(iConfig.getParameter<double>("maxEtaError"))
   , m_maxDxyError(iConfig.getParameter<double>("maxDxyError"))
   , m_maxDzError(iConfig.getParameter<double>("maxDzError"))
   , m_calculateVertex(iConfig.getParameter<bool>("calculateVertex"))
   , m_calculateIsolation(iConfig.getParameter<bool>("calculateIsolation"))
   , m_groupingMode_string(iConfig.getParameter<std::string>("groupingMode"))
   , m_maxDeltaR(iConfig.getParameter<double>("maxDeltaR"))
   , m_maxMass(iConfig.getParameter<double>("maxMass"))
   , m_minVertexProb(iConfig.getParameter<double>("minVertexProb"))
   , m_groupByCharge_string(iConfig.getParameter<std::string>("groupByCharge"))
{
   //register your products
   produces<pat::MuonCollection>("fakeMuons");
   produces<pat::MultiMuonCollection>();

   //now do what ever other initialization is needed

   if (m_groupingMode_string == "GroupByDeltaR") m_groupingMode = kGroupByDeltaR;
   else if (m_groupingMode_string == "GroupByMass") m_groupingMode = kGroupByMass;
   else if (m_groupingMode_string == "GroupByVertexProb") m_groupingMode = kGroupByVertexProb;
   else if (m_groupingMode_string == "GroupByDeltaRAndVertexProb") m_groupingMode = kGroupByDeltaRAndVertexProb;
   else if (m_groupingMode_string == "GroupByMassAndVertexProb") m_groupingMode = kGroupByMassAndVertexProb;
   else if (m_groupingMode_string == "GroupByDeltaROrMass") m_groupingMode = kGroupByDeltaROrMass;
   else if (m_groupingMode_string == "GroupByDeltaROrMassAndVertexProb") m_groupingMode = kGroupByDeltaROrMassAndVertexProb;
   else if (m_groupingMode_string == "GroupByMassAndVertexProbOrDeltaR") m_groupingMode = kGroupByMassAndVertexProbOrDeltaR;
   else {
      throw cms::Exception("BadConfig") << "groupingMode string not recognized (see MuJetProducer.cc)" << std::endl;
   }

   if (m_groupByCharge_string == "any") m_groupByCharge = kGroupByAnyCharge;
   else if (m_groupByCharge_string == "opposite") m_groupByCharge = kGroupByOppositeCharge;
   else if (m_groupByCharge_string == "same") m_groupByCharge = kGroupBySameCharge;
   else {
      throw cms::Exception("BadConfig") << "groupByCharge string not recognized (see MuJetProducer.cc)" << std::endl;
   }
}


MuJetPlusFakeMuon::~MuJetPlusFakeMuon()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

bool MuJetPlusFakeMuon::trackOkay(const reco::Track &track) {
   if (track.pt() < m_minPt  ||  track.p() < m_minPmag  ||  fabs(track.eta()) > m_maxAbsEta) return false;

   if (m_minTrackerHits > 0) {
      if (track.numberOfValidHits() < m_minTrackerHits) return false;
   }

   if (m_maxTrackerNormChi2 > 0.) {
      if (track.normalizedChi2() > m_maxTrackerNormChi2) return false;
   }

   if (m_maxTrackerDxy > 0.) {
      if (fabs(track.dxy()) > m_maxTrackerDxy) return false;
   }

   if (m_maxTrackerDz > 0.) {
      if (fabs(track.dsz()) > m_maxTrackerDz) return false;
   }

   if (m_maxQoverpError > 0.) {
      if (track.qoverpError() > m_maxQoverpError) return false;
   }

   if (m_maxPhiError > 0.) {
      if (track.phiError() > m_maxPhiError) return false;
   }

   if (m_maxEtaError > 0.) {
      if (track.etaError() > m_maxEtaError) return false;
   }

   if (m_maxDxyError > 0.) {
      if (track.dxyError() > m_maxDxyError) return false;
   }

   if (m_maxDzError > 0.) {
      if (track.dzError() > m_maxDzError) return false;
   }

   return true;
}

//
// member functions
//

// ------------ method called to produce the data  ------------
void
MuJetPlusFakeMuon::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   edm::Handle<pat::MultiMuonCollection> muJets;
   iEvent.getByLabel(m_muJets, muJets);

   edm::Handle<reco::TrackCollection> tracks;
   iEvent.getByLabel(m_tracks, tracks);

   edm::Handle<pat::MuonCollection> muons;
   edm::Handle<CaloTowerCollection> caloTowers;
   const pat::MuonCollection *muons_ptr = NULL;
   const CaloTowerCollection *caloTowers_ptr = NULL;
   if (m_calculateIsolation) {
      iEvent.getByLabel(m_muons, muons);
      iEvent.getByLabel(m_caloTowers, caloTowers);
      muons_ptr = &*muons;
      caloTowers_ptr = &*caloTowers;
   }

   edm::ESHandle<TransientTrackBuilder> transientTrackBuilder;
   const TransientTrackBuilder *transientTrackBuilder_ptr = NULL;
   if (m_calculateVertex) {
      iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", transientTrackBuilder);
      transientTrackBuilder_ptr = &*transientTrackBuilder;
   }

   KalmanVertexFitter vertexFitter;

   std::auto_ptr<pat::MuonCollection> fakeMuons(new pat::MuonCollection);
   std::auto_ptr<pat::MultiMuonCollection> newMuJets(new pat::MultiMuonCollection);

   for (pat::MultiMuonCollection::const_iterator muJet = muJets->begin();  muJet != muJets->end();  ++muJet) {
      double muonMass = -1.;
      std::vector<unsigned int> compatible;
      for (unsigned int muonIndex = 0;  muonIndex < muJet->numberOfDaughters();  muonIndex++) {
	 const pat::Muon *chosenMuon = muJet->muon(muonIndex);
	 muonMass = chosenMuon->mass();

	 unsigned int trackIndex = 0;
	 for (reco::TrackCollection::const_iterator track = tracks->begin();  track != tracks->end();  ++track, ++trackIndex) {
	    if (trackOkay(*track)) {
	       // 1. make sure it satisfies the charge conditions (get rid of half the candidates right away)
	       if (m_groupByCharge == kGroupByOppositeCharge  &&  track->charge() == chosenMuon->charge()) continue;
	       if (m_groupByCharge == kGroupBySameCharge  &&  track->charge() != chosenMuon->charge()) continue;

	       // 2. make sure this isn't one of the muons or one of the already-selected tracks
	       bool match = false;
	       for (unsigned int muonIndex2 = 0;  muonIndex2 < muJet->numberOfDaughters();  muonIndex2++) {
		  const pat::Muon *muon = muJet->muon(muonIndex2);
		  if (muon->innerTrack().isAvailable()  &&  muJet->sameTrack(&*(muon->innerTrack()), &*track)) match = true;
	       }
	       for (unsigned int trackIndex2 = 0;  trackIndex2 < compatible.size();  trackIndex2++) {
		  if (muJet->sameTrack(&*(reco::TrackRef(tracks, compatible[trackIndex2])), &*track)) match = true;
	       }
	       if (match) continue;

	       // 3. calculate vertex compatibility
	       double vertexProb = -1.;
	       GlobalPoint vertexPoint(0., 0., 0.);
	       GlobalVector trackMomentum(track->momentum().x(), track->momentum().y(), track->momentum().z());
	       GlobalVector muonMomentum(chosenMuon->momentum().x(), chosenMuon->momentum().y(), chosenMuon->momentum().z());

	       if (m_calculateVertex) {
		  std::vector<reco::TransientTrack> tracksToVertex;
		  tracksToVertex.push_back(transientTrackBuilder_ptr->build(*track));

		  if (chosenMuon->innerTrack().isAvailable()) {
		     tracksToVertex.push_back(transientTrackBuilder_ptr->build(chosenMuon->innerTrack()));
		  }
		  else if (chosenMuon->outerTrack().isAvailable()) {
		     tracksToVertex.push_back(transientTrackBuilder_ptr->build(chosenMuon->outerTrack()));
		  }
	  
		  CachingVertex<5> fittedVertex = vertexFitter.vertex(tracksToVertex);
	  
		  if (fittedVertex.isValid()  &&  fittedVertex.totalChiSquared() >= 0.) {
		     vertexProb = TMath::Prob(fittedVertex.totalChiSquared(), fittedVertex.degreesOfFreedom());
		     vertexPoint = GlobalPoint(fittedVertex.position().x(), fittedVertex.position().y(), fittedVertex.position().z());
		     TrajectoryStateClosestToPoint trackTSCTP = tracksToVertex[0].trajectoryStateClosestToPoint(vertexPoint);
		     trackMomentum = trackTSCTP.momentum();
		     TrajectoryStateClosestToPoint muonTSCTP = tracksToVertex[1].trajectoryStateClosestToPoint(vertexPoint);
		     muonMomentum = muonTSCTP.momentum();
		  }
	       }

	       // 4. calculate deltaR
	       double dphi = trackMomentum.phi() - muonMomentum.phi();
	       while (dphi > M_PI) dphi -= 2.*M_PI;
	       while (dphi < -M_PI) dphi += 2.*M_PI;
	       double deta = trackMomentum.eta() - muonMomentum.eta();
	       double dR = sqrt(dphi*dphi + deta*deta);

	       // 5. calculate invariant mass
	       double total_energy = sqrt(pow(chosenMuon->mass(), 2) + pow(trackMomentum.mag(), 2)) + sqrt(pow(chosenMuon->mass(), 2) + pow(muonMomentum.mag(), 2));
	       double total_px = trackMomentum.x() + muonMomentum.x();
	       double total_py = trackMomentum.y() + muonMomentum.y();
	       double total_pz = trackMomentum.z() + muonMomentum.z();
	       double mass = sqrt(total_energy*total_energy - total_px*total_px - total_py*total_py - total_pz*total_pz);

	       // 6. now see if it matches the mujet conditions
	       bool satisfied_deltaR = (dR < m_maxDeltaR);
	       bool satisfied_mass = (mass < m_maxMass);
	       bool satisfied_vertexProb = (vertexProb > m_minVertexProb);

	       bool satisfied = false;
	       if (m_groupingMode == kGroupByDeltaR) satisfied = satisfied_deltaR;
	       else if (m_groupingMode == kGroupByMass) satisfied = satisfied_mass;
	       else if (m_groupingMode == kGroupByVertexProb) satisfied = satisfied_vertexProb;
	       else if (m_groupingMode == kGroupByDeltaRAndVertexProb) satisfied = satisfied_deltaR  &&  satisfied_vertexProb;
	       else if (m_groupingMode == kGroupByMassAndVertexProb) satisfied = satisfied_mass  &&  satisfied_vertexProb;
	       else if (m_groupingMode == kGroupByDeltaROrMass) satisfied = satisfied_deltaR  ||  satisfied_mass;
	       else if (m_groupingMode == kGroupByDeltaROrMassAndVertexProb) satisfied = (satisfied_deltaR  ||  satisfied_mass)  &&  satisfied_vertexProb;
	       else if (m_groupingMode == kGroupByMassAndVertexProbOrDeltaR) satisfied = (satisfied_mass  &&  satisfied_vertexProb)  ||  satisfied_deltaR;
	       else assert(false);

	       if (satisfied) {
		  compatible.push_back(trackIndex);
	       }
	    }
	 }
      }
      assert(muonMass > 0.);

      if (compatible.size() > 0) {
	 unsigned int numberOfFakes = m_numberOfFakes;
	 if (compatible.size() < m_numberOfFakes) numberOfFakes = compatible.size();

	 std::map<unsigned int,bool> used;
	 for (unsigned int fake = 0;  fake < numberOfFakes; fake++) {
	    unsigned int chosen = compatible[m_trandom3.Integer(compatible.size())];
	    while (used.find(chosen) != used.end()) {
	       chosen = compatible[m_trandom3.Integer(compatible.size())];
	    }
	    used[chosen] = true;
	
	    reco::TrackRef trackRef(tracks, chosen);
	    reco::Muon newRecoMuon(trackRef->charge(), reco::Candidate::LorentzVector(trackRef->px(), trackRef->py(), trackRef->pz(), sqrt(pow(trackRef->p(), 2) + pow(muonMass, 2))));
	    newRecoMuon.setInnerTrack(trackRef);
	    std::vector<reco::MuonChamberMatch> emptyMatches;
	    newRecoMuon.setMatches(emptyMatches);

	    fakeMuons->push_back(pat::Muon(newRecoMuon));
	    fakeMuons->back().embedTrack();
	 }

	 std::vector<const pat::Muon*> muonsPlusTracks;
	 for (unsigned int i = 0;  i < muJet->numberOfDaughters();  i++) {
	    muonsPlusTracks.push_back(muJet->muon(i));
	 }
	 for (pat::MuonCollection::const_iterator fakeMuon = fakeMuons->begin();  fakeMuon != fakeMuons->end();  ++fakeMuon) {
	    muonsPlusTracks.push_back(&*fakeMuon);
	 }

	 pat::MultiMuon muJetWithTracks(muonsPlusTracks, transientTrackBuilder_ptr, &*tracks, muons_ptr, caloTowers_ptr, muJet->centralTrackIsolationCone(), muJet->unionTrackIsolationCone(), muJet->centralTrackThresholdPt(), muJet->unionTrackThresholdPt(), muJet->centralCaloIsolationCone(), muJet->unionCaloIsolationCone(), muJet->centralNumberAboveThresholdCone(), muJet->unionNumberAboveThresholdCone(), muJet->centralNumberAboveThresholdPt(), muJet->unionNumberAboveThresholdPt());
	 muJetWithTracks.addUserInt("numReal", muJet->numberOfDaughters());
	 muJetWithTracks.addUserInt("numFakes", fakeMuons->size());
	 newMuJets->push_back(muJetWithTracks);
      }
      else {
	 pat::MultiMuon muJetCopy(*muJet);
	 muJetCopy.addUserInt("numReal", muJet->numberOfDaughters());
	 muJetCopy.addUserInt("numFakes", 0);
	 newMuJets->push_back(muJetCopy);
      }
   } // end loop over muJets

   iEvent.put(fakeMuons, "fakeMuons");
   iEvent.put(newMuJets);
}

// ------------ method called once each job just before starting event loop  ------------
void 
MuJetPlusFakeMuon::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuJetPlusFakeMuon::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuJetPlusFakeMuon);
