// System include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// user include files
#include "TTree.h"
#include "TRandom3.h"
#include "TMath.h"


#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/TriggerReport.h"
#include "FWCore/Framework/interface/TriggerNamesService.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/GeometrySurface/interface/Cylinder.h"
#include "DataFormats/GeometrySurface/interface/Plane.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/MuonReco/interface/MuonChamberMatch.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSegmentMatch.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "AnalysisDataFormats/MuJetAnalysis/interface/MultiMuon.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "AnalysisDataFormats/MuJetAnalysis/src/eig3.cpp"
#include "AnalysisDataFormats/MuJetAnalysis/interface/eig3.h"

#include "FWCore/Framework/interface/EventSetup.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"



using namespace math;


//******************************************************************************
//              Auxiliary function: Order objects by pT                         
//******************************************************************************
bool PtOrder (const reco::GenParticle* p1, const reco::GenParticle* p2) { return (p1->pt() > p2->pt() ); }
bool PtOrderrec (const reco::Muon* p1, const reco::Muon* p2) { return (p1->pt() > p2->pt() ); }

//******************************************************************************
// Auxiliary function: Calculate difference between two angles: -PI < phi < PI  
//******************************************************************************
double My_dPhi (double phi1, double phi2) {
  double dPhi = phi1 - phi2;
  if (dPhi >  M_PI) dPhi -= 2.*M_PI;
  if (dPhi < -M_PI) dPhi += 2.*M_PI;
  return dPhi;
}


//******************************************************************************
// Auxiliary function: check if track is associated to dimuon
//******************************************************************************

bool sameTrack(const reco::Track *one, const reco::Track *two) {
   return (fabs(one->px() - two->px()) < 1e-10  &&
	   fabs(one->py() - two->py()) < 1e-10  &&
	   fabs(one->pz() - two->pz()) < 1e-10  &&
	   fabs(one->vx() - two->vx()) < 1e-10  &&
	   fabs(one->vy() - two->vy()) < 1e-10  &&
	   fabs(one->vz() - two->vz()) < 1e-10);
}


//******************************************************************************
// Cholesky decomposition of matrix A
//******************************************************************************

//Cholesky decomposition of matrix A
std::vector<std::vector<double> > cholesky(std::vector<std::vector<double> > A)
{
  int an = A.size();
  double sum1 = 0.0;
  double sum2 = 0.0;
  double sum3 = 0.0;

  std::vector<std::vector<double> > l(an, std::vector<double> (an));

  l[0][0] = sqrt(A[0][0]);
  for (int j = 1; j <= an-1; j++)
    l[j][0] = A[j][0]/l[0][0];
  for (int i = 1; i <= (an-2); i++)
    {
      for (int k = 0; k <= (i-1); k++)
	sum1 += pow(l[i][k], 2);
      l[i][i]= sqrt(A[i][i]-sum1);
      for (int j = (i+1); j <= (an-1); j++)
	{
	  for (int k = 0; k <= (i-1); k++)
	    sum2 += l[j][k]*l[i][k];
	  l[j][i]= (A[j][i]-sum2)/l[i][i];
	}
    }
  for (int k = 0; k <= (an-2); k++)
    sum3 += pow(l[an-1][k], 2);
  l[an-1][an-1] = sqrt(A[an-1][an-1]-sum3);
  return l;
}

struct EventContainer 
{

  double vtx_dzmj1;
  double vtx_dzmj2;
  double vtx_dxymj1;
  double vtx_dxymj2; 
  double vtx_massmj1;
  double vtx_massmj2;
  XYZTLorentzVector vtx_mj1;
  XYZTLorentzVector vtx_mj2;
  Global3DPoint vtx_1;
  Global3DPoint vtx_2;

  EventContainer(){
    vtx_dzmj1 = 0.;
    vtx_dzmj2 = 0.;
    vtx_dxymj1 = 0.;
    vtx_dxymj2 = 0.; 
    vtx_massmj1 = 0.;
    vtx_massmj2 = 0.;
    vtx_mj1 = XYZTLorentzVector(0.,0.,0.,0.);
    vtx_mj2 = XYZTLorentzVector(0.,0.,0.,0.);
    vtx_1 = Global3DPoint(0., 0., 0.);
    vtx_2 = Global3DPoint(0., 0., 0.);
  }

  void set_vtx_1(const Global3DPoint x){
    vtx_1 = x;
  }
  Global3DPoint get_vtx_1(){
    return vtx_1;
  }
  void set_vtx_2(const Global3DPoint x){
    vtx_2 = x;
  }
  Global3DPoint get_vtx_2(){
    return vtx_2;
  }
  void set_vtx_dzmj1(const double x){
    vtx_dzmj1 = x;
  }
  double get_vtx_dzmj1(){
    return vtx_dzmj1;
  }
  void set_vtx_dzmj2(const double x){
    vtx_dzmj2 = x;
  }
  double get_vtx_dzmj2(){
    return vtx_dzmj2;
  }
  void set_vtx_dxymj1(const double x){
    vtx_dxymj1 = x;
  }
  double get_vtx_dxymj1(){
    return vtx_dxymj1;
  }
  void set_vtx_dxymj2(const double x){
    vtx_dxymj2 = x;
  }
  double get_vtx_dxymj2(){
    return vtx_dxymj2;
  }
  void set_vtx_massmj1(const double x){
    vtx_massmj1 = x;
  }
  double get_vtx_massmj1(){
    return vtx_massmj1;
  }
  void set_vtx_massmj2(const double x){
    vtx_massmj2 = x;
  }
  double get_vtx_massmj2(){
    return vtx_massmj2;
  }
  void set_vtx_mj1(XYZTLorentzVector x){
    vtx_mj1 = x;
  }
  XYZTLorentzVector get_vtx_mj1(){
    return vtx_mj1;
  }
  void set_vtx_mj2(XYZTLorentzVector x){
    vtx_mj2 = x;
  }
  XYZTLorentzVector get_vtx_mj2(){
    return vtx_mj2;
  }
};



//******************************************************************************
//                           Class declaration                                  
//******************************************************************************


class Analysis : public edm::EDAnalyzer {
public:
  explicit Analysis(const edm::ParameterSet&);
  ~Analysis();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:

  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      
  // ---------- Member data ----------
  
  TTree * AnaTree;  // Pointer to Tree
  
  Int_t run;    // run number    
  Int_t lumi;   // lumi number   
  Int_t event;  // event number  
  Int_t trigger;
  Int_t isVtx;

  // ---------- Generator Level ----------
  
  edm::InputTag m_genParticles;   // Label to access generator particles
  

  // Gen branches in ROOT tree

  // ---------- GEN level ---------------

  Int_t genmuons;
  Int_t gengamD;
  Int_t genHiggs;

  Float_t ptgenMuons[10];
  Float_t etagenMuons[10];

  Float_t ptgenMuonsA0[5];
  Float_t ptgenMuonsA1[5];
  Float_t etagenMuonsA0[5];
  Float_t etagenMuonsA1[5];
  Float_t ptgengamD[5];
  Float_t ptgenHiggs[5];

  Float_t gen_match_dzmj1;
  Float_t gen_match_dzmj2;

  Float_t gen_match_x1diff;
  Float_t gen_match_y1diff;
  Float_t gen_match_z1diff;
  Float_t gen_match_x2diff;
  Float_t gen_match_y2diff;
  Float_t gen_match_z2diff;


  Float_t gen_match_x1diff_old;
  Float_t gen_match_y1diff_old;
  Float_t gen_match_z1diff_old;
  Float_t gen_match_x2diff_old;
  Float_t gen_match_y2diff_old;
  Float_t gen_match_z2diff_old;


  Float_t gen_match_massmj1;
  Float_t gen_match_massmj2;

  Float_t gen_mismatch_dzmj1;
  Float_t gen_mismatch_dzmj2;
  Float_t gen_mismatch_massmj1;
  Float_t gen_mismatch_massmj2;

  Float_t gen_match_chi2mj1;
  Float_t gen_match_chi2mj2;
  Float_t gen_mismatch_chi2mj1;
  Float_t gen_mismatch_chi2mj2;


  // ---------- RECO Level ----------
  
  Int_t recmuons;
  Float_t ptrecMuons[30];
  Float_t etarecMuons[30];

  Int_t recmujets;
  Int_t mj1muons;
  Int_t mj2muons;
  Float_t ptmj1muons[4];
  Float_t ptmj2muons[4];
  Float_t etamj1muons[4];
  Float_t etamj2muons[4];

  Float_t dzmj1;
  Float_t dzmj2;

  Int_t three_sig_x1;
  Int_t three_sig_y1;
  Int_t three_sig_z1;
  Int_t three_sig_1;
  Int_t three_sig_x2;
  Int_t three_sig_y2;
  Int_t three_sig_z2;
  Int_t three_sig_2;


  Float_t dxymj1;
  Float_t dxymj2;

  Float_t dzmj1_old;
  Float_t dzmj2_old;

  Float_t dphi1;
  Float_t dphi2;

  Float_t massmj1;
  Float_t massmj2;

  Float_t isomj1;
  Float_t isomj2;

  Float_t reco_match_dzmj1;
  Float_t reco_match_dzmj2;


  // Labels to access
  edm::InputTag m_muons;  // reconstructed muons
  edm::InputTag m_muJets; // muon jets built from reconstructed muons

  // Auxiliary variables
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
Analysis::Analysis(const edm::ParameterSet& iConfig)

{

  
  genmuons=0;
  gengamD=0;
  genHiggs=0;
  recmujets=0;
  mj1muons=0;
  mj2muons=0;
  genmuons=0;
  recmuons=0;
  dzmj1=0.;
  dzmj2=0.;

  dxymj1=0.;
  dxymj2=0.;


  three_sig_x1 = 0;
  three_sig_y1 = 0;
  three_sig_z1 = 0;
  three_sig_1 = 0;
  three_sig_x2 = 0;
  three_sig_y2 = 0;
  three_sig_z2 = 0;
  three_sig_2 = 0;


  massmj1=0.0;
  massmj2=0.0;
  isomj1=0.0;
  isomj2=0.0;
  isVtx=0;


  dzmj1_old=0.;
  dzmj2_old=0.;

  dphi1=0.;
  dphi2=0.;

  gen_match_dzmj1=0.;
  gen_match_dzmj2=0.;
  gen_match_massmj1=0.0;
  gen_match_massmj2=0.0;
  gen_match_chi2mj1=0.0;
  gen_match_chi2mj2=0.0;
  gen_mismatch_dzmj1=0.;
  gen_mismatch_dzmj2=0.;
  gen_mismatch_massmj1=0.0;
  gen_mismatch_massmj2=0.0;
  gen_mismatch_chi2mj1=0.0;
  gen_mismatch_chi2mj2=0.0;


  gen_match_x1diff = 0.0;
  gen_match_y1diff = 0.0;
  gen_match_z1diff = 0.0;
  gen_match_x2diff = 0.0;
  gen_match_y2diff = 0.0;
  gen_match_z2diff = 0.0;

  gen_match_x1diff_old = 0.0;
  gen_match_y1diff_old = 0.0;
  gen_match_z1diff_old = 0.0;
  gen_match_x2diff_old = 0.0;
  gen_match_y2diff_old = 0.0;
  gen_match_z2diff_old = 0.0;

  reco_match_dzmj1=0.;
  reco_match_dzmj2=0.;



  AnaTree  = NULL;    
  
  // ---------- Generator Level ----------
  m_genParticles = iConfig.getParameter<edm::InputTag>("genParticles");
  
  
  // ---------- RECO Level ----------
  m_muons = iConfig.getParameter<edm::InputTag>("muons");
  m_muJets = iConfig.getParameter<edm::InputTag>("muJets");


  //  m_maxIsoDiMuons  = iConfig.getParameter<double>("maxIsoDiMuons");

}


Analysis::~Analysis()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}

//
// member functions
//

// ------------ method called for each event  ------------
void
Analysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  // get the run, lumi and event number
  run   = iEvent.id().run();
  lumi  = iEvent.id().luminosityBlock();
  event = iEvent.id().event();


  // if (event != 24465) return;

  // if (event != 20484) return;

  // if (event != 20246)  return;

  // if ((event != 20246) && (event != 20444) && (event != 20802) && (event != 20849) && (event != 20918) && (event != 21269) && (event != 21385) && (event != 75056) && (event != 21646)) return;

  // if (event!=20768) return;

  // if (event!=39954) return;

  // if (event != 29143) return;

  // if (event != 29143 && event != 20349 && event != 20849 && event != 22124 && event != 22912 && event != 22924 && event != 23541 && event != 25289 && event != 26405 && event != 75522 && event != 26585 && event != 27817) return;

  // std::cout << "Event: " << event << std::endl;
  // if (event != 20349) return;


  // if (event != 20246 && event != 20392 && event != 20802 && event != 20918 && event != 20956 && event != 21151 && event != 21759 && event != 21998 && event != 78889 && event != 25002) return;


  // if (event != 39954 && event != 20802 && event != 22014 && event != 75151 && event != 25020 && event != 25578 && event != 28503 && event != 29802 && event != 12989 && event != 14753 && event != 18956 ) return;



  bool verbose = false;
  bool e39995 = false;
  bool e20775 = false;

  // if (event == 39995) e39995 = true;
  // if (event == 20775) e20775 = true;

  // std::cout << std::endl << "run: " << run << ", lumi: " << lumi << ", event: " << event << std::endl;


  //****************************************************************************
  //                              GEN LEVEL                                     
  //****************************************************************************
  
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByLabel(m_genParticles, genParticles);
  
  // Loop over all genParticles and save prompt muons from particles with codes 36 (a1) or 3000022 (gammaD) in vector genMuons
  std::vector<const reco::GenParticle*> genH;
  std::vector<const reco::GenParticle*> genA;
  std::vector<const reco::GenParticle*> genMuons;
  std::vector<const reco::Candidate*>   genMuonMothers;

  // Loop over all gen particles
  int counterGenParticle = 0;
  for(reco::GenParticleCollection::const_iterator iGenParticle = genParticles->begin();  iGenParticle != genParticles->end();  ++iGenParticle) {
    counterGenParticle++;
    
    if ( fabs( iGenParticle->pdgId() ) == 13 && iGenParticle->status() == 1 ) {
      // Mother of the muon can be muon. Find the last muon in this chain: genMuonCand
      // Example: a1 -> mu+ (status = 3) mu- (status = 3)
      //          mu- (status = 3) -> mu- (status = 2) -> mu- (status = 1)
      const reco::Candidate *genMuonCand = &(*iGenParticle);
      bool isMuonMother = true;
      while(isMuonMother) {
        isMuonMother = false;
        for ( size_t iMother = 0; iMother < genMuonCand->numberOfMothers(); iMother++ ) {
          if ( fabs( genMuonCand->mother(iMother)->pdgId() ) == 13 ) {
            isMuonMother = true;
            genMuonCand = genMuonCand->mother(iMother);
          }
        }
      }
      // Loop over all real (non-muon) mothers of the muon (here we use genMuonCand)
      for ( size_t iMother = 0; iMother < genMuonCand->numberOfMothers(); iMother++ ) {
        // Check if mother is CP-odd Higgs (PdgId = 36) or gamma_Dark (PdgId = 3000022)
        if ( genMuonCand->mother(iMother)->pdgId() == 36 || genMuonCand->mother(iMother)->pdgId() == 3000022 || genMuonCand->mother(iMother)->pdgId() == 443 ) {
          // Store the muon (stable, first in chain) into vector
          genMuons.push_back(&(*iGenParticle));
	  // 	  std::cout << "genMuon pt: " << genMuonCand->pt() << std::endl;
	  // 	  std::cout << "genMuon mother pt: " << genMuonCand->mother(iMother)->pt() << std::endl; 
          // Store mother of the muon into vector. We need this to group muons into dimuons later
          genMuonMothers.push_back(genMuonCand->mother(iMother));
        }
      }
    }

    // Check if gen particle is decaying (status = 3) CP-even Higgs (pdgId = +/-35)
    if ( iGenParticle->status() == 3 && iGenParticle->pdgId() == 35 ) {
      genH.push_back(&(*iGenParticle)); // Store the Higgs into vector
    }
    // Check if gen particle is
    if (    ( iGenParticle->status() == 3 && iGenParticle->pdgId() == 36      )     // decaying (status = 3) CP-odd Higgs (pdgId = 36)
	    || ( iGenParticle->status() == 3 && iGenParticle->pdgId() == 3000022 )     // decaying (status = 3) gamma_Dark (pdgId = 3000022)
	    || ( iGenParticle->status() == 2 && iGenParticle->pdgId() == 443     ) ) { // decaying (status = 2) J/psi (pdgId = 443)
      genA.push_back(&(*iGenParticle));
    }
  }
  

  genmuons = genMuons.size();
  gengamD  = genA.size();
  genHiggs = genH.size();


  
  // Group muons from the same mother into dimuons
  std::vector< std::vector<const reco::GenParticle*> > genMuonGroups;
  std::vector<const reco::GenParticle*> genMuonsTMP1       = genMuons;
  std::vector<const reco::Candidate*>   genMuonMothersTMP1 = genMuonMothers;
  unsigned int nMuonGroup = 0;
  while ( genMuonsTMP1.size() > 0 ) {
    std::vector<const reco::GenParticle*> genMuonsTMP2;
    std::vector<const reco::Candidate*>   genMuonMothersTMP2;
    std::vector<const reco::GenParticle*> genMuonsSameMother;
    for ( unsigned int j = 0; j < genMuonsTMP1.size(); j++ ) {
      // Check if mothers are the same particle
      if ( fabs( genMuonMothersTMP1[0]->pt() - genMuonMothersTMP1[j]->pt() ) < 0.00001 ) {
	genMuonsSameMother.push_back( genMuonsTMP1[j] );
      } else {
	genMuonsTMP2.push_back( genMuonsTMP1[j] );
	genMuonMothersTMP2.push_back( genMuonMothersTMP1[j] );
      }
    }
    genMuonGroups.push_back(genMuonsSameMother);
    genMuonsTMP1       = genMuonsTMP2;
    genMuonMothersTMP1 = genMuonMothersTMP2;
    nMuonGroup++;
  }

  std::sort( genMuons.begin(), genMuons.end(), PtOrder );

  for(int k=0;k<genmuons;k++){
    ptgenMuons[k] = genMuons[k]->pt();
    etagenMuons[k] = genMuons[k]->eta();
  }

  if ( genH.size() == 1 ) {
    ptgenHiggs[0] = genH[0]->pt();
  }
  
  if ( genA.size() >= 2 ) {
    // Sort genA by pT (leading pT first)
    std::sort (genA.begin(), genA.end(), PtOrder);
    ptgengamD[0] = genA[0]->pt();
    ptgengamD[1] = genA[1]->pt();
  } else {
    std::cout << "WARNING! genA.size() < 2" << std::endl;
  }
  
  if ( genMuonGroups.size() == 2 && genMuonGroups[0].size() == 2 && genMuonGroups[1].size() == 2 ) {
    std::sort( genMuonGroups[0].begin(), genMuonGroups[0].end(), PtOrder );
    std::sort( genMuonGroups[1].begin(), genMuonGroups[1].end(), PtOrder );

    ptgenMuonsA0[0] = genMuonGroups[0][0]->pt();
    ptgenMuonsA0[1] = genMuonGroups[0][1]->pt();

    ptgenMuonsA1[0] = genMuonGroups[1][0]->pt();
    ptgenMuonsA1[1] = genMuonGroups[1][1]->pt();

    etagenMuonsA0[0] = genMuonGroups[0][0]->eta();
    etagenMuonsA0[1] = genMuonGroups[0][1]->eta();

    etagenMuonsA1[0] = genMuonGroups[1][0]->eta();
    etagenMuonsA1[1] = genMuonGroups[1][1]->eta();
  }


  
  //================= Reco =======================//
  
  edm::Handle<pat::MuonCollection> muons;
  iEvent.getByLabel(m_muons, muons);

  std::vector<const reco::Muon*> recMuons;

  for (pat::MuonCollection::const_iterator iMuon = muons->begin();  iMuon != muons->end();  ++iMuon) {
    if ( fabs(iMuon->eta()) < 2.4 && ( iMuon->isTrackerMuon() || iMuon->isGlobalMuon() ) ) {
      recMuons.push_back( &(*iMuon));
    }
  }

  recmuons = recMuons.size();

  // Sort recMuons by pT (leading pT first)
  if ( recMuons.size() > 1 ) std::sort( recMuons.begin(), recMuons.end(), PtOrderrec );

  for(int k=0;k<recmuons;k++){
    ptrecMuons[k] = recMuons[k]->pt();
    etarecMuons[k] = recMuons[k]->eta();
  }
  

  //================ MuJets  ======================//

  edm::Handle<pat::MultiMuonCollection> muJets;
  iEvent.getByLabel(m_muJets, muJets);

  const pat::MultiMuon *muJet1 = NULL;
  const pat::MultiMuon *muJet2 = NULL;

  unsigned int nMuJets = muJets->size();
  recmujets = muJets->size();

  //do delta R matching here

  edm::Handle<reco::BeamSpot> beamSpot;
  iEvent.getByLabel("offlineBeamSpot",beamSpot);

  const pat::MultiMuon *muJet1_tmp = NULL;
  const pat::MultiMuon *muJet2_tmp = NULL;

  int rec11;
  int rec12;
  int rec21;
  int rec22;


  int matched = 0;
  if (nMuJets == 2 && genMuonGroups.size() == 2 && genMuonGroups[0].size() == 2 && genMuonGroups[1].size() == 2) {
    int mj1muons_tmp = (*muJets)[0].numberOfDaughters();
    int mj2muons_tmp = (*muJets)[1].numberOfDaughters();
    muJet1_tmp = &((*muJets)[0]);
    muJet2_tmp = &((*muJets)[1]);     

    if ( muJet1_tmp != NULL && muJet2_tmp != NULL ) {
	
      for (unsigned int q = 0; q < genMuonGroups.size(); ++q) { // loop over gen muon jets

	double min_dR1_1 = 1000;
	double min_dR1_2 = 1000;
	double min_dR2_1 = 1000;
	double min_dR2_2 = 1000;

	for(int m=0;m<mj1muons_tmp;m++){ 
	  double dPhi1_1 = My_dPhi((*muJets)[0].muon(m)->phi(), genMuonGroups[q][0]->phi());
	  double dEta1_1 = (*muJets)[0].muon(m)->eta() - genMuonGroups[q][0]->eta();
	  double dPhi1_2 = My_dPhi((*muJets)[0].muon(m)->phi(), genMuonGroups[q][1]->phi());
	  double dEta1_2 = (*muJets)[0].muon(m)->eta() - genMuonGroups[q][1]->eta();
	  double dR1_1 = sqrt( dPhi1_1*dPhi1_1 + dEta1_1*dEta1_1 ); 

	  if (dR1_1< min_dR1_1) {
	    min_dR1_1 = dR1_1;
	  }
	  double dR1_2 = sqrt( dPhi1_2*dPhi1_2 + dEta1_2*dEta1_2 ); 

	  if (dR1_2< min_dR1_2){
	    min_dR1_2 = dR1_2;
	  }

	}
	for(int m=0;m<mj2muons_tmp;m++){ 

	  double dPhi2_1 = My_dPhi((*muJets)[1].muon(m)->phi(), genMuonGroups[q][0]->phi());
	  double dEta2_1 = (*muJets)[1].muon(m)->eta() - genMuonGroups[q][0]->eta();
	  double dPhi2_2 = My_dPhi((*muJets)[1].muon(m)->phi(), genMuonGroups[q][1]->phi());
	  double dEta2_2 = (*muJets)[1].muon(m)->eta() - genMuonGroups[q][1]->eta();
	  double dR2_1 = sqrt( dPhi2_1*dPhi2_1 + dEta2_1*dEta2_1 ); 

	  if (dR2_1< min_dR2_1){
	    min_dR2_1 = dR2_1;
	  }
	  double dR2_2 = sqrt( dPhi2_2*dPhi2_2 + dEta2_2*dEta2_2 ); 

	  if (dR2_2< min_dR2_2){
	    min_dR2_2 = dR2_2;
	  }
	}
	if (min_dR1_1 < 0.1 && min_dR1_2 < 0.1 ){
	  matched++;

	  if (q==0){
	    rec11 = 11;
	    rec12 = 12;
	  }
	  if (q==1){
	    rec11 = 21;
	    rec12 = 22;
	  }


	  if (muJet1_tmp->vertexValid())  gen_match_chi2mj1= muJet1_tmp->vertexNormalizedChi2();

	  double gez = genMuonGroups[q][0]->vertex().z()-beamSpot->position().z();
	  double gex = genMuonGroups[q][0]->vertex().x()-beamSpot->position().x();
	  double gey = genMuonGroups[q][0]->vertex().y()-beamSpot->position().y();
	  double gepx = genMuonGroups[q][0]->p4().x() + genMuonGroups[q][1]->p4().x();
	  double gepy = genMuonGroups[q][0]->p4().y() + genMuonGroups[q][1]->p4().y();
	  double gepz = genMuonGroups[q][0]->p4().z() + genMuonGroups[q][1]->p4().z();
	  double gept =  pow(gepx*gepx + gepy*gepy,0.5);

	  gen_match_dzmj1 = (gez) - ((gex)*gepx+(gey)*gepy)/gept * gepz/gept;

	  XYZTLorentzVector p4Sum;
	  p4Sum += XYZTLorentzVector(gepx, gepy, gepz, genMuonGroups[q][0]->p4().e() + genMuonGroups[q][1]->p4().e());

	  gen_match_massmj1= p4Sum.mass();

	  gen_mismatch_massmj1= 0.;
	  gen_mismatch_chi2mj1= 0.;
	  gen_mismatch_dzmj1 = 0.;
	}
	if (min_dR2_1 < 0.1 && min_dR2_2 < 0.1 ){
	  matched++;

	  if (q==0){
	    rec21 = 11;
	    rec22 = 12;
	  }
	  if (q==1){
	    rec21 = 21;
	    rec22 = 22;
	  }


	  if (muJet2_tmp->vertexValid()) gen_match_chi2mj2= muJet2_tmp->vertexNormalizedChi2();
	  gen_mismatch_massmj2= 0.;
	  gen_mismatch_dzmj2 = 0.;
	  gen_mismatch_chi2mj2 = 0.;

	  double gez = genMuonGroups[q][0]->vertex().z()-beamSpot->position().z();
	  double gex = genMuonGroups[q][0]->vertex().x()-beamSpot->position().x();
	  double gey = genMuonGroups[q][0]->vertex().y()-beamSpot->position().y();

	  double gepx = genMuonGroups[q][0]->p4().x() + genMuonGroups[q][1]->p4().x();
	  double gepy = genMuonGroups[q][0]->p4().y() + genMuonGroups[q][1]->p4().y();
	  double gepz = genMuonGroups[q][0]->p4().z() + genMuonGroups[q][1]->p4().z();
	  double gept =  pow(gepx*gepx + gepy*gepy,0.5);
	  gen_match_dzmj2 = (gez) - ((gex)*gepx+(gey)*gepy)/gept * gepz/gept;

	  XYZTLorentzVector p4Sum;
	  p4Sum += XYZTLorentzVector(gepx, gepy, gepz, genMuonGroups[q][0]->p4().e() + genMuonGroups[q][1]->p4().e());
	  gen_match_massmj2= p4Sum.mass();

	}
	if (q == 1 && matched == 0){
	  gen_mismatch_massmj1= muJet1_tmp->mass();
	  gen_mismatch_massmj2= muJet2_tmp->mass();
	  gen_mismatch_dzmj1 = muJet1_tmp->dz(beamSpot->position());
	  gen_mismatch_dzmj2 = muJet2_tmp->dz(beamSpot->position());
	  if (muJet1_tmp->vertexValid()) gen_mismatch_chi2mj1= muJet1_tmp->vertexNormalizedChi2();
	  if (muJet2_tmp->vertexValid()) gen_mismatch_chi2mj2= muJet2_tmp->vertexNormalizedChi2();

	  gen_match_massmj1= 0;
	  gen_match_massmj2= 0;
	  gen_match_dzmj1 = 0;
	  gen_match_dzmj2 = 0;
	  gen_match_chi2mj1= 0;
	  gen_match_chi2mj2= 0;
	}
      }
    }
  }  


  if (nMuJets == 2) {

    mj1muons = (*muJets)[0].numberOfDaughters();
    mj2muons = (*muJets)[1].numberOfDaughters();
      
    muJet1 = &((*muJets)[0]);
    muJet2 = &((*muJets)[1]);

    for(int m=0;m<mj1muons;m++){
      ptmj1muons[m]  = (*muJets)[0].muon(m)->pt();
      etamj1muons[m]  = (*muJets)[0].muon(m)->eta();
    }

    for(int m=0;m<mj2muons;m++){
      ptmj2muons[m]  = (*muJets)[1].muon(m)->pt();
      etamj2muons[m]  = (*muJets)[1].muon(m)->eta();
    }


  }

  //standard:

  double temp_dz_diff;
  double temp_dzmj1;
  double temp_dzmj2;

  // double dphi1 = 0.;
  // double dphi2 = 0.;

  if ( muJet1 != NULL && muJet2 != NULL ) {
    dzmj1 = muJet1->dz(beamSpot->position());
    dzmj2 = muJet2->dz(beamSpot->position());
    temp_dz_diff = fabs(dzmj1-dzmj2);
    temp_dzmj1 = dzmj1;
    temp_dzmj2 = dzmj2;

    dzmj1_old = dzmj1;
    dzmj2_old = dzmj2;

    dphi1 = My_dPhi((*muJets)[0].muon(0)->phi(), (*muJets)[0].muon(1)->phi() );
    dphi2 = My_dPhi((*muJets)[1].muon(0)->phi(), (*muJets)[1].muon(1)->phi() );

    // std::cout <<"loop dphi1: " << dphi1 << std::endl;
    // std::cout << "loop dphi2: " << dphi2 << std::endl;


  }

  if ( muJet1!= NULL && muJet2!= NULL ) {
    massmj1= muJet1->mass();
    massmj2= muJet2->mass();
  }

  if (verbose){
    if ( muJet1!= NULL && muJet2!= NULL ) {
      if (nMuJets == 2) {
	if (((*muJets)[0].vertexValid()) && ((*muJets)[1].vertexValid())){
      
	  std::cout << "*************" << std::endl;
	  std::cout << "Event start: " << event << std::endl;
	  std::cout << "dzmj1: " << dzmj1 << ", dzmj2: " << dzmj2 << std::endl;
	  std::cout << "fabs(dzmj1-dzmj2): " << fabs(dzmj1-dzmj2) << std::endl;
	  std::cout << "massmj1: " << massmj1 << ", massmj2: " << massmj2 << std::endl;
  
	}
      }
    }
  }

  const unsigned int seed = 0;

  TRandom3 *rnd = new TRandom3(seed);

  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByLabel("generalTracks", tracks);

   std::vector<reco::TransientTrack> tracksToVertex1;
   std::vector<reco::TransientTrack> tracksToVertex2;
  
  if ( muJet1 != NULL && muJet2 != NULL ) {     

   edm::ESHandle<TransientTrackBuilder> transientTrackBuilder;
   iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", transientTrackBuilder);

   for (unsigned int i = 0;  i < (*muJets)[0].numberOfDaughters();  i++) {
     if ((*muJets)[0].muon(i)->innerTrack().isAvailable()) {
       tracksToVertex1.push_back(transientTrackBuilder->build((*muJets)[0].muon(i)->innerTrack()));
     }
     else if ((*muJets)[0].muon(i)->outerTrack().isAvailable()) {
       tracksToVertex1.push_back(transientTrackBuilder->build((*muJets)[0].muon(i)->outerTrack()));
     }
   }

   for (unsigned int i = 0;  i < (*muJets)[1].numberOfDaughters();  i++) {
     if ((*muJets)[1].muon(i)->innerTrack().isAvailable()) {
       tracksToVertex2.push_back(transientTrackBuilder->build((*muJets)[1].muon(i)->innerTrack()));
     }
     else if ((*muJets)[1].muon(i)->outerTrack().isAvailable()) {
       tracksToVertex2.push_back(transientTrackBuilder->build((*muJets)[1].muon(i)->outerTrack()));
     }
   }
  }

  double min_dz_n = 10000.;
  double min_dxy_n = 10000.;
  double min_mass_n = 10000.;
  double newx1 = 0;
  double newy1 = 0;
  double newz1 = 0;
  double newx2 = 0;
  double newy2 = 0;
  double newz2 = 0;
  std::vector<EventContainer> events_arr;



  XYZTLorentzVector diMuonTmp;


  if ( muJet1 != NULL && muJet2 != NULL ) {     

    if (nMuJets == 2) {
      if (((*muJets)[0].vertexValid()) && ((*muJets)[1].vertexValid())){

	double A1[3][3];
	double V1[3][3];
	double d1[3];
	d1[0]=0.;
	d1[1]=0.;
	d1[2]=0.;

	std::vector<std::vector<double> > cA1(3,std::vector<double>(3));
	std::vector<std::vector<double> > cA2(3,std::vector<double>(3));

	double A2[3][3];
	double V2[3][3];
	double d2[3];
	d2[0]=0.;
	d2[1]=0.;
	d2[2]=0.;

	verbose = false;

	std::cout << "Event: " << event << std::endl;
	for (unsigned int a = 0; a < 3; ++a){
	  for (unsigned int b = 0; b < 3; ++b){
	    A1[a][b] = (*muJets)[0].vertexCovariance(a,b);
	    A2[a][b] = (*muJets)[1].vertexCovariance(a,b);
	    cA1[a][b] = (*muJets)[0].vertexCovariance(a,b);
	    cA2[a][b] = (*muJets)[1].vertexCovariance(a,b);
	    // if (a==b){

	    if (verbose){
	      std::cout << "(*muJets)[0].vertexCovariance(" << a << "," << b << "): " << (*muJets)[0].vertexCovariance(a,b) << std::endl;
	      std::cout << "(*muJets)[1].vertexCovariance(" << a << "," << b << "): " << (*muJets)[1].vertexCovariance(a,b) << std::endl;

	      std::cout << "cA1[" << a << "][" << b << "]: " << cA1[a][b] << std::endl;
	      std::cout << "cA2[" << a << "][" << b << "]: " << cA2[a][b] << std::endl;
	    }
	    // }
	  }
	}


	std::vector<std::vector<double> > cfA1 = cholesky(cA1);
	std::vector<std::vector<double> > cfA2 = cholesky(cA2);

	verbose = false;
	if (verbose){
	  for (unsigned int a = 0; a < 3; ++a){
	    for (unsigned int b = 0; b < 3; ++b){
	      std::cout << "cfA1[" << a << "][" << b << "]: " << cfA1[a][b] << std::endl;
	    }
	  }
	  for (unsigned int a = 0; a < 3; ++a){
	    for (unsigned int b = 0; b < 3; ++b){
	      std::cout << "cfA2[" << a << "][" << b << "]: " << cfA2[a][b] << std::endl;
	    }
	  }
	}
	verbose = false;


	// Perform the decomposition       
	eigen_decomposition(A1, V1, d1);
	eigen_decomposition(A2, V2, d2);

	verbose = false;
	if (verbose)
	  {
	    for (unsigned int a = 0; a < 3; ++a){
	      for (unsigned int b = 0; b < 3; ++b){
		std::cout << "eigenvector V1[" << a << "][" << b <<"]: " << V1[a][b] << std::endl;   
		std::cout << "eigenvector V2[" << a << "][" << b <<"]: " << V2[a][b] << std::endl;   
	      }
	    }

	    for (unsigned int a = 0; a < 3; ++a){
	      std::cout << "eigenvalue d1[" << a << "]: " << d1[a] << std::endl;   
	      std::cout << "eigenvalue d2[" << a << "]: " << d2[a] << std::endl;       
	    }

	    std::cout << "muonjet1: " << std::endl;
	    std::cout << "lambda1*nu1 = (" << d1[2]*V1[0][2] << ", " << d1[2]*V1[1][2] << ", " << d1[2]*V1[2][2] << ")" << std::endl;
	    std::cout << "lambda2*nu2 = (" << d1[1]*V1[0][1] << ", " << d1[1]*V1[1][1] << ", " << d1[1]*V1[2][1] << ")" << std::endl;
	    std::cout << "lambda3*nu3 = (" << d1[0]*V1[0][0] << ", " << d1[0]*V1[1][0] << ", " << d1[0]*V1[2][0] << ")" << std::endl;

	    std::cout << "muonjet2: " << std::endl;
	    std::cout << "lambda1*nu1 = (" << d2[2]*V2[0][2] << ", " << d2[2]*V2[1][2] << ", " << d2[2]*V2[2][2] << ")" << std::endl;
	    std::cout << "lambda2*nu2 = (" << d2[1]*V2[0][1] << ", " << d2[1]*V2[1][1] << ", " << d2[1]*V2[2][1] << ")" << std::endl;
	    std::cout << "lambda3*nu3 = (" << d2[0]*V2[0][0] << ", " << d2[0]*V2[1][0] << ", " << d2[0]*V2[2][0] << ")" << std::endl;
	  }
	verbose = false;
	double elipiseX1Size = 0;
	double elipiseX2Size = 0;
	double elipiseY1Size = 0;
	double elipiseY2Size = 0;
	double elipiseZ1Size = 0;
	double elipiseZ2Size = 0;

	if (((*muJets)[0].vertexCovariance(0,0)) > ((*muJets)[0].vertexCovariance(1,1))){

	  if (((*muJets)[0].vertexCovariance(0,0)) > ((*muJets)[0].vertexCovariance(2,2))){

	    elipiseX1Size = sqrt(d1[2]);

	      if (((*muJets)[0].vertexCovariance(1,1)) > ((*muJets)[0].vertexCovariance(2,2))){
		elipiseY1Size = sqrt(d1[1]);
		elipiseZ1Size = sqrt(d1[0]);		
	      }
	      else {
		elipiseY1Size = sqrt(d1[0]);		
		elipiseZ1Size = sqrt(d1[1]);		
	      }
	  }
	  else {
	    elipiseX1Size = sqrt(d1[1]);
	    elipiseY1Size = sqrt(d1[0]);
	    elipiseZ1Size = sqrt(d1[2]);		
	  }
	}


	if (((*muJets)[0].vertexCovariance(1,1)) > ((*muJets)[0].vertexCovariance(0,0))){

	  if (((*muJets)[0].vertexCovariance(1,1)) > ((*muJets)[0].vertexCovariance(2,2))){

	    elipiseY1Size = sqrt(d1[2]);

	      if (((*muJets)[0].vertexCovariance(0,0)) > ((*muJets)[0].vertexCovariance(2,2))){
		elipiseX1Size = sqrt(d1[1]);
		elipiseZ1Size = sqrt(d1[0]);		
	      }
	      else {
		elipiseX1Size = sqrt(d1[0]);		
		elipiseZ1Size = sqrt(d1[1]);		
	      }
	  }
	  else {
	    elipiseX1Size = sqrt(d1[0]);
	    elipiseY1Size = sqrt(d1[1]);
	    elipiseZ1Size = sqrt(d1[2]);		
	  }
	}



	/////////////////////////

	if (((*muJets)[1].vertexCovariance(0,0)) > ((*muJets)[1].vertexCovariance(1,1))){

	  if (((*muJets)[1].vertexCovariance(0,0)) > ((*muJets)[1].vertexCovariance(2,2))){

	    elipiseX2Size = sqrt(d2[2]);

	      if (((*muJets)[1].vertexCovariance(1,1)) > ((*muJets)[1].vertexCovariance(2,2))){
		elipiseY2Size = sqrt(d2[1]);
		elipiseZ2Size = sqrt(d2[0]);		
	      }
	      else {
		elipiseY2Size = sqrt(d2[0]);		
		elipiseZ2Size = sqrt(d2[1]);		
	      }
	  }
	  else {
	    elipiseX2Size = sqrt(d2[1]);
	    elipiseY2Size = sqrt(d2[0]);
	    elipiseZ2Size = sqrt(d2[2]);		
	  }
	}


	if (((*muJets)[1].vertexCovariance(1,1)) > ((*muJets)[1].vertexCovariance(0,0))){

	  if (((*muJets)[1].vertexCovariance(1,1)) > ((*muJets)[1].vertexCovariance(2,2))){

	    elipiseY2Size = sqrt(d2[2]);

	      if (((*muJets)[1].vertexCovariance(0,0)) > ((*muJets)[1].vertexCovariance(2,2))){
		elipiseX2Size = sqrt(d2[1]);
		elipiseZ2Size = sqrt(d2[0]);		
	      }
	      else {
		elipiseX2Size = sqrt(d2[0]);		
		elipiseZ2Size = sqrt(d2[1]);		
	      }
	  }
	  else {
	    elipiseX2Size = sqrt(d2[0]);
	    elipiseY2Size = sqrt(d2[1]);
	    elipiseZ2Size = sqrt(d2[2]);		
	  }
	}


	// if (elipiseX1Size > elipiseY1Size) std::cout << "(elipiseX1Size > elipiseY1Size)" << std::endl;
	// if (elipiseX1Size > elipiseZ1Size) std::cout << "(elipiseX1Size > elipiseZ1Size)" << std::endl;
	// if (elipiseY1Size > elipiseZ1Size) std::cout << "(elipiseY1Size > elipiseZ1Size)" << std::endl;
	// if (elipiseX2Size > elipiseY2Size) std::cout << "(elipiseX2Size > elipiseY2Size)" << std::endl;
	// if (elipiseX2Size > elipiseZ2Size) std::cout << "(elipiseX2Size > elipiseZ2Size)" << std::endl;
	// if (elipiseY2Size > elipiseZ2Size) std::cout << "(elipiseY2Size > elipiseZ2Size)" << std::endl;

	if (elipiseX1Size == 0) std::cout << "elipiseX1Size == 0 (grep)" << std::endl;
	if (elipiseY1Size == 0) std::cout << "elipiseY1Size == 0 (grep)" << std::endl;
	if (elipiseZ1Size == 0) std::cout << "elipiseZ1Size == 0 (grep)" << std::endl;
	if (elipiseX2Size == 0) std::cout << "elipiseX2Size == 0 (grep)" << std::endl;
	if (elipiseY2Size == 0) std::cout << "elipiseY2Size == 0 (grep)" << std::endl;
	if (elipiseZ2Size == 0) std::cout << "elipiseZ2Size == 0 (grep)" << std::endl;





	// double elipiseX1Size = 3*sqrt(d1[1]);
	// double elipiseX2Size = 3*sqrt(d2[1]);
	// double elipiseY1Size = 3*sqrt(d1[0]);
	// double elipiseY2Size = 3*sqrt(d2[0]);

	// double elipiseX1Size = 0;
	// double elipiseX2Size = 0;
	// double elipiseY1Size = 0;
	// double elipiseY2Size = 0;




	/////////////////////////////////////////

	//     CODE TO USE MAX ELLIPSE SIZE

	/////////////////////////////////////////


	double max1 = 0;
	double max2 = 0;

	if (sqrt(d1[1]) > sqrt(d1[0])){
	  if (sqrt(d1[2]) > sqrt(d1[1])){
	    max1 = sqrt(d1[2]);
	  }
	  else{
	    max1 = sqrt(d1[1]);
	  }
	}
	else {
	  if (sqrt(d1[2]) > sqrt(d1[0])){
	    max1 = sqrt(d1[2]);
	  }
	  else{
	    max1 = sqrt(d1[0]);
	  }
	}

	if (sqrt(d2[1]) > sqrt(d2[0])){
	  if (sqrt(d2[2]) > sqrt(d2[1])){
	    max2 = sqrt(d2[2]);
	  }
	  else{
	    max2 = sqrt(d2[1]);
	  }
	}
	else {
	  if (sqrt(d2[2]) > sqrt(d2[0])){
	    max2 = sqrt(d2[2]);
	  }
	  else{
	    max2 = sqrt(d2[0]);
	  }
	}


	// elipiseX1Size = max1;
	// elipiseY1Size = max1;
	// elipiseZ1Size = max1;

	// elipiseX2Size = max2;
	// elipiseY2Size = max2;
	// elipiseZ2Size = max2;



	// elipiseX1Size = sqrt(d1[0]);
	// elipiseX2Size = sqrt(d2[0]);
	// elipiseY1Size = sqrt(d1[1]);
	// elipiseY2Size = sqrt(d2[1]);
	// elipiseZ1Size = sqrt(d1[2]);
	// elipiseZ2Size = sqrt(d2[2]);



	// elipiseX1Size = cfA1[0][0];
	// elipiseY1Size = cfA1[1][1];
	// elipiseZ1Size = cfA1[2][2];

	// elipiseX2Size = cfA2[0][0];
	// elipiseY2Size = cfA2[1][1];
	// elipiseZ2Size = cfA2[2][2];

	elipiseX1Size = cfA1[0][0];
	elipiseY1Size = cfA1[1][0];
	elipiseZ1Size = cfA1[2][0];

	elipiseX2Size = cfA2[0][0];
	elipiseY2Size = cfA2[1][0];
	elipiseZ2Size = cfA2[2][0];




	verbose = false;
	if (verbose){
	  std::cout << "elipiseX1Size: " << elipiseX1Size << std::endl;
	  std::cout << "elipiseY1Size: " << elipiseY1Size << std::endl;
	  std::cout << "elipiseZ1Size: " << elipiseZ1Size << std::endl;
	  std::cout << "elipiseX2Size: " << elipiseX2Size << std::endl;
	  std::cout << "elipiseY2Size: " << elipiseY2Size << std::endl;
	  std::cout << "elipiseZ2Size: " << elipiseZ2Size << std::endl;
	}
	verbose = false;
	Global3DPoint muonjet1v = muJet1->vertexPoint();
	Global3DPoint muonjet2v = muJet2->vertexPoint();
	Global3DVector muonjet1 = muJet1->vertexMomentum();
	Global3DVector muonjet2 = muJet2->vertexMomentum();

	int throws = 100000;
	int min_dice = throws;

	XYZTLorentzVector final_mj1_vtx;
	XYZTLorentzVector final_mj2_vtx;

	verbose = false;
	if (verbose){
	  std::cout <<"Event: " << event << std::endl;
	  std::cout << "gen1 xyz: " << genMuonGroups[0][0]->vertex().x() << ", " << genMuonGroups[0][0]->vertex().y() << ", " << genMuonGroups[0][0]->vertex().z() << std::endl;
	  std::cout << "gen2 xyz: " << genMuonGroups[1][0]->vertex().x() << ", " << genMuonGroups[1][0]->vertex().y() << ", " << genMuonGroups[1][0]->vertex().z() << std::endl;
	  std::cout << "reco vertex 1 xyz: " << muonjet1v.x() << ", " << muonjet1v.y() << ", " << muonjet1v.z() << std::endl;
	  std::cout << "reco vertex 2 xyz: " << muonjet2v.x() << ", " << muonjet2v.y() << ", " << muonjet2v.z() << std::endl;
	}
	verbose = false;

	// double dz_arr_mj1[throws];

	std::vector<double> dz_arr_mj1;
	std::vector<double> dz_arr_mj2;

	// double dz_arr_mj2[throws];
	double dz_sum_mj1 = 0;
	double dz_sum_mj2 = 0;


	std::vector<double> dzdiff_arr;
	std::vector<double> dxydiff_arr;
	std::vector<double> dzdiff_abs_arr;
	std::vector<double> dxydiff_abs_arr;


	// double dzdiff_arr[throws];
	double dzdiff_sum = 0;
	// double dxydiff_arr[throws];
	double dxydiff_sum = 0;

	// double dzdiff_abs_arr[throws];
	double dzdiff_abs_sum = 0;
	// double dxydiff_abs_arr[throws];
	double dxydiff_abs_sum = 0;


	std::vector<double> dxy_arr_mj1;
	std::vector<double> dxy_arr_mj2;


	// double dxy_arr_mj1[throws];
	// double dxy_arr_mj2[throws];
	double dxy_sum_mj1 = 0;
	double dxy_sum_mj2 = 0;


	double x1rnd = 0;
	double x2rnd = 0;
	double y1rnd = 0;
	double y2rnd = 0;
	double z1rnd = 0;
	double z2rnd = 0;


	bool goodGenDz = false;
	bool foundEvent = false;

	  if (rec11<20 && rec22 > 20){

	    x1rnd = fabs(muonjet1v.x() - genMuonGroups[0][0]->vertex().x());
	    x2rnd = fabs(muonjet2v.x() - genMuonGroups[1][0]->vertex().x());
	    y1rnd = fabs(muonjet1v.y() - genMuonGroups[0][0]->vertex().y());
	    y2rnd = fabs(muonjet2v.y() - genMuonGroups[1][0]->vertex().y());
	    z1rnd = fabs(muonjet1v.z() - genMuonGroups[0][0]->vertex().z());
	    z2rnd = fabs(muonjet2v.z() - genMuonGroups[1][0]->vertex().z());

	  }
	  else if (rec11 > 20 && rec22 < 20){

	    x1rnd = fabs(muonjet1v.x() - genMuonGroups[1][0]->vertex().x());
	    x2rnd = fabs(muonjet2v.x() - genMuonGroups[0][0]->vertex().x());
	    y1rnd = fabs(muonjet1v.y() - genMuonGroups[1][0]->vertex().y());
	    y2rnd = fabs(muonjet2v.y() - genMuonGroups[0][0]->vertex().y());
	    z1rnd = fabs(muonjet1v.z() - genMuonGroups[1][0]->vertex().z());
	    z2rnd = fabs(muonjet2v.z() - genMuonGroups[0][0]->vertex().z());
	     
	  }


	for (int dice = 0; dice < throws; dice++) {


	  ///////////////////////////////////////////////////////////

	  //            Standard ellipse setup

	  ///////////////////////////////////////////////////////////


	  double rnd_gZ1 = rnd->Gaus(0., fabs(elipiseZ1Size));
	  double rnd_gZ2 = rnd->Gaus(0., fabs(elipiseZ2Size));
	  double rnd_gX1 = rnd->Gaus(muonjet1v.x(), fabs(elipiseX1Size));
	  double rnd_gY1 = rnd->Gaus(muonjet1v.y(), fabs(elipiseY1Size));
	  double rnd_gX2 = rnd->Gaus(muonjet2v.x(), fabs(elipiseX2Size));
	  double rnd_gY2 = rnd->Gaus(muonjet2v.y(), fabs(elipiseY2Size));

	  if (verbose){
	    std::cout << "x1: " << rnd_gX1 << std::endl;
	    std::cout << "y1: " << rnd_gY1 << std::endl;
	    std::cout << "z1: " << muonjet1v.z()+rnd_gZ1 << std::endl;

	    std::cout << "x2: " << rnd_gX2 << std::endl;
	    std::cout << "y2: " << rnd_gY2 << std::endl;
	    std::cout << "z2: " << muonjet2v.z()+rnd_gZ2 << std::endl;
	  }
	  // double rnd_gZ1 = rnd->Gaus(0., z1rnd);
	  // double rnd_gZ2 = rnd->Gaus(0., z2rnd);
	  // double rnd_gX1 = rnd->Gaus(muonjet1v.x(), x1rnd);
	  // double rnd_gY1 = rnd->Gaus(muonjet1v.y(), y1rnd);
	  // double rnd_gX2 = rnd->Gaus(muonjet2v.x(), x2rnd);
	  // double rnd_gY2 = rnd->Gaus(muonjet2v.y(), y2rnd);

	  if (verbose){
	    std::cout << "gaus 1: " << rnd->Gaus(0., 1.) << std::endl;
	    std::cout << "gaus 2: " << rnd->Gaus(0., 1.) << std::endl;
	    std::cout << "gaus 3: " << rnd->Gaus(0., 1.) << std::endl;
	  }

	  // double rnd_gZ1 = 0.;
	  // double rnd_gZ2 = 0.;
	  // double rnd_gX1 = 0.;
	  // double rnd_gY1 = 0.;
	  // double rnd_gX2 = 0.;
	  // double rnd_gY2 = 0.;

	  XYZTLorentzVector gen_mj1;
	  XYZTLorentzVector gen_mj2;


	  /*
	  if (rec11<20 && rec22 > 20){
	    rnd_gX1 = rnd->Gaus(genMuonGroups[0][0]->vertex().x(), x1rnd);
	    rnd_gY1 = rnd->Gaus(genMuonGroups[0][0]->vertex().y(), y1rnd);
	    rnd_gZ1 = rnd->Gaus(0., z1rnd);
	    rnd_gX2 = rnd->Gaus(genMuonGroups[1][0]->vertex().x(), x2rnd);
	    rnd_gY2 = rnd->Gaus(genMuonGroups[1][0]->vertex().y(), y2rnd);
	    rnd_gZ2 = rnd->Gaus(0., z2rnd);


	    gen_mj1 += XYZTLorentzVector(genMuonGroups[0][0]->p4().x(), genMuonGroups[0][0]->p4().y(), genMuonGroups[0][0]->p4().z(), genMuonGroups[0][0]->p4().t());
	    gen_mj1 += XYZTLorentzVector(genMuonGroups[0][1]->p4().x(), genMuonGroups[0][1]->p4().y(), genMuonGroups[0][1]->p4().z(), genMuonGroups[0][1]->p4().t());

	    gen_mj2 += XYZTLorentzVector(genMuonGroups[1][0]->p4().x(), genMuonGroups[1][0]->p4().y(), genMuonGroups[1][0]->p4().z(), genMuonGroups[1][0]->p4().t());
	    gen_mj2 += XYZTLorentzVector(genMuonGroups[1][1]->p4().x(), genMuonGroups[1][1]->p4().y(), genMuonGroups[1][1]->p4().z(), genMuonGroups[1][1]->p4().t());



	  }
	  else if (rec11 > 20 && rec22 < 20){
	    rnd_gX1 = rnd->Gaus(genMuonGroups[1][0]->vertex().x(), x1rnd);
	    rnd_gY1 = rnd->Gaus(genMuonGroups[1][0]->vertex().y(), y1rnd);
	    rnd_gZ1 = rnd->Gaus(0., z1rnd);
	    rnd_gX2 = rnd->Gaus(genMuonGroups[0][0]->vertex().x(), x2rnd);
	    rnd_gY2 = rnd->Gaus(genMuonGroups[0][0]->vertex().y(), y2rnd);
	    rnd_gZ2 = rnd->Gaus(0., z2rnd);


	    gen_mj2 += XYZTLorentzVector(genMuonGroups[0][0]->p4().x(), genMuonGroups[0][0]->p4().y(), genMuonGroups[0][0]->p4().z(), genMuonGroups[0][0]->p4().t());
	    gen_mj2 += XYZTLorentzVector(genMuonGroups[0][1]->p4().x(), genMuonGroups[0][1]->p4().y(), genMuonGroups[0][1]->p4().z(), genMuonGroups[0][1]->p4().t());

	    gen_mj1 += XYZTLorentzVector(genMuonGroups[1][0]->p4().x(), genMuonGroups[1][0]->p4().y(), genMuonGroups[1][0]->p4().z(), genMuonGroups[1][0]->p4().t());
	    gen_mj1 += XYZTLorentzVector(genMuonGroups[1][1]->p4().x(), genMuonGroups[1][1]->p4().y(), genMuonGroups[1][1]->p4().z(), genMuonGroups[1][1]->p4().t());

	  }

	  */


	  //////////////////////////////////////////////////

	  //         Standard reco setup

	  //////////////////////////////////////////////////


	  double t1 = rnd_gZ1/muonjet1.z();
	  double t2 = rnd_gZ2/muonjet2.z();

	  double newX1 = rnd_gX1 + (t1*muonjet1.x());
	  double newY1 = rnd_gY1 + (t1*muonjet1.y());
	  double newX2 = rnd_gX2 + (t2*muonjet2.x());
	  double newY2 = rnd_gY2 + (t2*muonjet2.y());

	  double newZ1 = (muonjet1v.z()) + (t1*muonjet1.z());
	  double newZ2 = (muonjet2v.z()) + (t2*muonjet2.z());

	  verbose = false;	
	  if (verbose){
	    std::cout << "new1 xyz: " << newX1 << ", " << newY1 << ", " << newZ1 << std::endl;
	    std::cout << "reco vertex 1 xyz: " << muonjet1v.x() << ", " << muonjet1v.y() << ", " << muonjet1v.z() << std::endl;

	    std::cout << "new2 xyz: " << newX2 << ", " << newY2 << ", " << newZ2 << std::endl;
	    std::cout << "reco vertex 2 xyz: " << muonjet2v.x() << ", " << muonjet2v.y() << ", " << muonjet2v.z() << std::endl;

	    std::cout << "muonjet1.x(): " << muonjet1.x() << std::endl;
	    std::cout << "muonjet1.y(): " << muonjet1.y() << std::endl;
	    std::cout << "muonjet1.z(): " << muonjet1.z() << std::endl;

	    std::cout << "muonjet2.x(): " << muonjet2.x() << std::endl;
	    std::cout << "muonjet2.y(): " << muonjet2.y() << std::endl;
	    std::cout << "muonjet2.z(): " << muonjet2.z() << std::endl;
	  }
	  /*
	  double t1 = rnd_gZ1/gen_mj1.z();
	  double t2 = rnd_gZ2/gen_mj2.z();

	  double newX1 = 0;
	  double newY1 = 0;
	  double newX2 = 0;
	  double newY2 = 0;
	  double newZ1 = 0;
	  double newZ2 = 0;

	  if (rec11<20 && rec22 > 20){
	    newX1 = rnd_gX1 + (t1*gen_mj1.x());
	    newY1 = rnd_gY1 + (t1*gen_mj1.y());
	    newX2 = rnd_gX2 + (t2*gen_mj2.x());
	    newY2 = rnd_gY2 + (t2*gen_mj2.y());
	    newZ1 = (genMuonGroups[0][0]->vertex().z()) + (t1*gen_mj1.z());
	    newZ2 = (genMuonGroups[1][0]->vertex().z()) + (t2*gen_mj2.z());
	  }
	  else if (rec11 > 20 && rec22 < 20){
	    newX1 = rnd_gX1 + (t1*gen_mj1.x());
	    newY1 = rnd_gY1 + (t1*gen_mj1.y());
	    newX2 = rnd_gX2 + (t2*gen_mj2.x());
	    newY2 = rnd_gY2 + (t2*gen_mj2.y());
	    newZ1 = (genMuonGroups[1][0]->vertex().z()) + (t1*gen_mj1.z());
	    newZ2 = (genMuonGroups[0][0]->vertex().z()) + (t2*gen_mj2.z());
	  }

	  */


	  // double newX1 = muonjet1v.x() + (t1*muonjet1.x()) + rnd_gX1;
	  // double newY1 = muonjet1v.y() + (t1*muonjet1.y()) + rnd_gY1;
	  // double newZ1 = muonjet1v.z() + (t1*muonjet1.z());

	  // double newX2 = muonjet2v.x() + (t2*muonjet2.x()) + rnd_gX2;
	  // double newY2 = muonjet2v.y() + (t2*muonjet2.y()) + rnd_gY2;
	  // double newZ2 = muonjet2v.z() + (t2*muonjet2.z());



	  // double rnd_X1 = rnd->Gaus(muonjet1v.x(), elipiseX1Size);
	  // double rnd_Y1 = rnd->Gaus(muonjet1v.y(), elipiseY1Size);
	  // double rnd_Z1 = rnd->Gaus(muonjet1v.z(), elipiseZ1Size);
	  // double rnd_X2 = rnd->Gaus(muonjet2v.x(), elipiseX2Size);
	  // double rnd_Y2 = rnd->Gaus(muonjet2v.y(), elipiseY2Size);
	  // double rnd_Z2 = rnd->Gaus(muonjet2v.z(), elipiseZ2Size);


	  // std::cout << "new1 xyz: " << newX1 << ", " << newY1 << ", " << newZ1 << std::endl;
	  // std::cout << "reco vertex 1 xyz: " << muonjet1v.x() << ", " << muonjet1v.y() << ", " << muonjet1v.z() << std::endl;

	  // std::cout << "new2 xyz: " << newX2 << ", " << newY2 << ", " << newZ2 << std::endl;
	  // std::cout << "reco vertex 2 xyz: " << muonjet2v.x() << ", " << muonjet2v.y() << ", " << muonjet2v.z() << std::endl;


	  // double rnd_X1 = 0;
	  // double rnd_Y1 = 0;
	  // double rnd_Z1 = 0;
	  // double rnd_X2 = 0;
	  // double rnd_Y2 = 0;
	  // double rnd_Z2 = 0;



	  /*
	  if (rec11<20 && rec22 > 20){
	    newX1 = genMuonGroups[0][0]->vertex().x();
	    newX2 = genMuonGroups[1][0]->vertex().x();
	    newY1 = genMuonGroups[0][0]->vertex().y();
	    newY2 = genMuonGroups[1][0]->vertex().y();
	    newZ1 = genMuonGroups[0][0]->vertex().z();
	    newZ2 = genMuonGroups[1][0]->vertex().z();
	  }
	  else if (rec11 > 20 && rec22 < 20){
	    newX1 = genMuonGroups[1][0]->vertex().x();
	    newX2 = genMuonGroups[0][0]->vertex().x();
	    newY1 = genMuonGroups[1][0]->vertex().y();
	    newY2 = genMuonGroups[0][0]->vertex().y();
	    newZ1 = genMuonGroups[1][0]->vertex().z();
	    newZ2 = genMuonGroups[0][0]->vertex().z();
	  }
	  */



	  ///////////////////////////////////////////

	  //      STANDARD IMPLEMENTATION:

	  ///////////////////////////////////////////

	  // double rnd_X1 = newX1;
	  // double rnd_Y1 = newY1;
	  // double rnd_Z1 = newZ1;
	  // double rnd_X2 = newX2;
	  // double rnd_Y2 = newY2;
	  // double rnd_Z2 = newZ2;


	  // Sample Gaussian of mean 0 and width 1:

	  double random1x = rnd->Gaus(0., 1.);
	  double random1y = rnd->Gaus(0., 1.);
	  double random1z = rnd->Gaus(0., 1.);
	  double random2x = rnd->Gaus(0., 1.);
	  double random2y = rnd->Gaus(0., 1.);
	  double random2z = rnd->Gaus(0., 1.);


	  double choms_rnd1[3][3];
	  double choms_rnd2[3][3];

	  for (unsigned int a = 0; a < 3; ++a){
	    for (unsigned int b = 0; b < 3; ++b){
	      choms_rnd1[a][b]=0.;
	      choms_rnd2[a][b]=0.;
	    }
	  }

	  choms_rnd1[0][0] = random1x*cfA1[0][0];
	  choms_rnd1[0][1] = random1x*cfA1[0][1];
	  choms_rnd1[0][2] = random1x*cfA1[0][2];

	  choms_rnd1[1][0] = random1y*cfA1[1][0];
	  choms_rnd1[1][1] = random1y*cfA1[1][1];
	  choms_rnd1[1][2] = random1y*cfA1[1][2];

	  choms_rnd1[2][0] = random1z*cfA1[2][0];
	  choms_rnd1[2][1] = random1z*cfA1[2][1];
	  choms_rnd1[2][2] = random1z*cfA1[2][2];

	  choms_rnd2[0][0] = random2x*cfA2[0][0];
	  choms_rnd2[0][1] = random2x*cfA2[0][1];
	  choms_rnd2[0][2] = random2x*cfA2[0][2];

	  choms_rnd2[1][0] = random2y*cfA2[1][0];
	  choms_rnd2[1][1] = random2y*cfA2[1][1];
	  choms_rnd2[1][2] = random2y*cfA2[1][2];

	  choms_rnd2[2][0] = random2z*cfA2[2][0];
	  choms_rnd2[2][1] = random2z*cfA2[2][1];
	  choms_rnd2[2][2] = random2z*cfA2[2][2];

	  if (verbose){
	    for (unsigned int a = 0; a < 3; ++a){
	      for (unsigned int b = 0; b < 3; ++b){
		std::cout << "choms_rnd1[" << a << "][" << b << "]: " << choms_rnd1[a][b] << std::endl;
		std::cout << "choms_rnd2[" << a << "][" << b << "]: " << choms_rnd2[a][b] << std::endl;
	      }
	    }
	  }

	  double x1sum = choms_rnd1[0][0]+choms_rnd1[0][1]+choms_rnd1[0][2];
	  double y1sum = choms_rnd1[1][0]+choms_rnd1[1][1]+choms_rnd1[1][2];
	  double z1sum = choms_rnd1[2][0]+choms_rnd1[2][1]+choms_rnd1[2][2];

	  double x2sum = choms_rnd2[0][0]+choms_rnd2[0][1]+choms_rnd2[0][2];
	  double y2sum = choms_rnd2[1][0]+choms_rnd2[1][1]+choms_rnd2[1][2];
	  double z2sum = choms_rnd2[2][0]+choms_rnd2[2][1]+choms_rnd2[2][2];


	  x1sum+=muonjet1v.x();
	  y1sum+=muonjet1v.y();
	  z1sum+=muonjet1v.z();
	  x2sum+=muonjet2v.x();
	  y2sum+=muonjet2v.y();
	  z2sum+=muonjet2v.z();

	  double rnd_X1 = x1sum;
	  double rnd_Y1 = y1sum;
	  double rnd_Z1 = z1sum;
	  double rnd_X2 = x2sum;
	  double rnd_Y2 = y2sum;
	  double rnd_Z2 = z2sum;

	  // double rnd_X1 = muonjet1v.x()+(cfA1[0][0]*rnd->Gaus(0., 1.));
	  // double rnd_Y1 = muonjet1v.y()+(cfA1[1][1]*rnd->Gaus(0., 1.));
	  // double rnd_Z1 = muonjet1v.z()+(cfA1[2][2]*rnd->Gaus(0., 1.));

	  // double rnd_X2 = muonjet2v.x()+(cfA2[0][0]*rnd->Gaus(0., 1.));
	  // double rnd_Y2 = muonjet2v.y()+(cfA2[1][1]*rnd->Gaus(0., 1.));
	  // double rnd_Z2 = muonjet2v.z()+(cfA2[2][2]*rnd->Gaus(0., 1.));





	  verbose = false;
	  if (verbose){
	    std::cout << "rndx1->Fill(" << rnd_X1 << ");" << std::endl;
	    std::cout << "rndy1->Fill(" << rnd_Y1 << ");" << std::endl;
	    std::cout << "rndz1->Fill(" << rnd_Z1 << ");" << std::endl;
	    std::cout << "rndx2->Fill(" << rnd_X2 << ");" << std::endl;
	    std::cout << "rndy2->Fill(" << rnd_Y2 << ");" << std::endl;
	    std::cout << "rndz2->Fill(" << rnd_Z2 << ");" << std::endl;
	  }
	  verbose = false;

	  Global3DPoint newVtx1(rnd_X1, rnd_Y1, rnd_Z1);
	  Global3DPoint newVtx2(rnd_X2, rnd_Y2, rnd_Z2);


	  std::vector<XYZTLorentzVector> new_vertex1_P4;
	  std::vector<XYZTLorentzVector> new_vertex2_P4;

	  for (unsigned int i = 0;  i < tracksToVertex1.size();  i++) {
	    TrajectoryStateClosestToPoint TSCTP1 = tracksToVertex1[i].trajectoryStateClosestToPoint(newVtx1);
	    Global3DVector momentum1 = TSCTP1.momentum();      
	    new_vertex1_P4.push_back(XYZTLorentzVector(momentum1.x(), momentum1.y(), momentum1.z(), sqrt(momentum1.mag2() + pow((*muJets)[0].muon(i)->mass(), 2))));

	  }
	  for (unsigned int i = 0;  i < tracksToVertex2.size();  i++) {
	    TrajectoryStateClosestToPoint TSCTP2 = tracksToVertex2[i].trajectoryStateClosestToPoint(newVtx2);
	    Global3DVector momentum2 = TSCTP2.momentum();
	    new_vertex2_P4.push_back(XYZTLorentzVector(momentum2.x(), momentum2.y(), momentum2.z(), sqrt(momentum2.mag2() + pow((*muJets)[1].muon(i)->mass(), 2))));

	  }

	  XYZTLorentzVector new_mj1_vtx;
	  XYZTLorentzVector new_mj2_vtx;

	  for (unsigned int i = 0;  i < tracksToVertex1.size();  i++) {
	    new_mj1_vtx+= new_vertex1_P4[i];
	  }
	  for (unsigned int i = 0;  i < tracksToVertex2.size();  i++) {
	    new_mj2_vtx+= new_vertex2_P4[i];
	  }

	  double rez1 = rnd_Z1-beamSpot->position().z();
	  double rex1 = rnd_X1-beamSpot->position().x();
	  double rey1 = rnd_Y1-beamSpot->position().y();
	  double repx1 = new_mj1_vtx.x();
	  double repy1 = new_mj1_vtx.y();
	  double repz1 = new_mj1_vtx.z();
	  double rept1 = pow((new_mj1_vtx.x()*new_mj1_vtx.x())+(new_mj1_vtx.y()*new_mj1_vtx.y()),0.5);			

	  double re_dzmj1 = (rez1) - ((rex1)*repx1+(rey1)*repy1)/rept1 * repz1/rept1;

	  double rez2 = rnd_Z2-beamSpot->position().z();
	  double rex2 = rnd_X2-beamSpot->position().x();
	  double rey2 = rnd_Y2-beamSpot->position().y();
	  double repx2 = new_mj2_vtx.x();
	  double repy2 = new_mj2_vtx.y();
	  double repz2 = new_mj2_vtx.z();
	  double rept2 = pow((new_mj2_vtx.x()*new_mj2_vtx.x())+(new_mj2_vtx.y()*new_mj2_vtx.y()),0.5);
	    
	  double re_dzmj2 = (rez2) - ((rex2)*repx2+(rey2)*repy2)/rept2 * repz2/rept2;


 	  // double dxy1 =(  (rex1) * repy1 + (rey1) * repx1 ) / rept1; 
 	  // double dxy2 =(  (rex2) * repy2 + (rey2) * repx2 ) / rept2; 

 	  double dxy1 =( - (rex1) * repy1 + (rey1) * repx1 ) / rept1; 
 	  double dxy2 =( - (rex2) * repy2 + (rey2) * repx2 ) / rept2; 


	  if (rec11<20 && rec22 > 20){

	    if ((fabs(rnd_X1 - genMuonGroups[0][0]->vertex().x()) < 3.0) && (fabs(rnd_Y1 - genMuonGroups[0][0]->vertex().y()) < 3.0) && (fabs(rnd_Z1 - genMuonGroups[0][0]->vertex().z()) < 3.0)) three_sig_1 = 1;

	    if ((fabs(rnd_X2 - genMuonGroups[1][0]->vertex().x()) < 3.0) && (fabs(rnd_Y2 - genMuonGroups[1][0]->vertex().y()) < 3.0) && (fabs(rnd_Z2 - genMuonGroups[1][0]->vertex().z()) < 3.0)) three_sig_2 = 1;

	    if (fabs(rnd_X1 - genMuonGroups[0][0]->vertex().x()) < 3.0) three_sig_x1 = 1;
	    if (fabs(rnd_Y1 - genMuonGroups[0][0]->vertex().y()) < 3.0) three_sig_y1 = 1;
	    if (fabs(rnd_Z1 - genMuonGroups[0][0]->vertex().z()) < 3.0) three_sig_z1 = 1;
	    if (fabs(rnd_X2 - genMuonGroups[1][0]->vertex().x()) < 3.0) three_sig_x2 = 1;
	    if (fabs(rnd_Y2 - genMuonGroups[1][0]->vertex().y()) < 3.0) three_sig_y2 = 1;
	    if (fabs(rnd_Z2 - genMuonGroups[1][0]->vertex().z()) < 3.0) three_sig_z2 = 1;

	  }
	  else if (rec11 > 20 && rec22 < 20){


	    if ((fabs(rnd_X1 - genMuonGroups[1][0]->vertex().x()) < 3.0) && (fabs(rnd_Y1 - genMuonGroups[1][0]->vertex().y()) < 3.0) && (fabs(rnd_Z1 - genMuonGroups[1][0]->vertex().z()) < 3.0)) three_sig_1 = 1;

	    if ((fabs(rnd_X2 - genMuonGroups[0][0]->vertex().x()) < 3.0) && (fabs(rnd_Y2 - genMuonGroups[0][0]->vertex().y()) < 3.0) && (fabs(rnd_Z2 - genMuonGroups[0][0]->vertex().z()) < 3.0)) three_sig_2 = 1;

	    if (fabs(rnd_X1 - genMuonGroups[1][0]->vertex().x()) < 3.0) three_sig_x1 = 1;
	    if (fabs(rnd_Y1 - genMuonGroups[1][0]->vertex().y()) < 3.0) three_sig_y1 = 1;
	    if (fabs(rnd_Z1 - genMuonGroups[1][0]->vertex().z()) < 3.0) three_sig_z1 = 1;
	    if (fabs(rnd_X2 - genMuonGroups[0][0]->vertex().x()) < 3.0) three_sig_x2 = 1;
	    if (fabs(rnd_Y2 - genMuonGroups[0][0]->vertex().y()) < 3.0) three_sig_y2 = 1;
	    if (fabs(rnd_Z2 - genMuonGroups[0][0]->vertex().z()) < 3.0) three_sig_z2 = 1;
	     
	  }


	  if (dice == 0){
	    if (rec11<20 && rec22 > 20){

	      verbose = false;
	      if (verbose){
		std::cout << "genMuonGroups[0][0]->vertex().x()-rnd_X1: " << genMuonGroups[0][0]->vertex().x()-rnd_X1 << std::endl;
		std::cout << "genMuonGroups[0][0]->vertex().y()-rnd_Y1: " << genMuonGroups[0][0]->vertex().y()-rnd_Y1 << std::endl;
		std::cout << "genMuonGroups[0][0]->vertex().z()-rnd_Z1: " << genMuonGroups[0][0]->vertex().z()-rnd_Z1 << std::endl;
		std::cout << "genMuonGroups[1][0]->vertex().x()-rnd_X2: " << genMuonGroups[1][0]->vertex().x()-rnd_X2 << std::endl;
		std::cout << "genMuonGroups[1][0]->vertex().y()-rnd_Y2: " << genMuonGroups[1][0]->vertex().y()-rnd_Y2 << std::endl;
		std::cout << "genMuonGroups[1][0]->vertex().z()-rnd_Z2: " << genMuonGroups[1][0]->vertex().z()-rnd_Z2 << std::endl;
	      }
	      double gex1 = genMuonGroups[0][0]->vertex().x()-beamSpot->position().x();
	      double gey1 = genMuonGroups[0][0]->vertex().y()-beamSpot->position().y();
	      double gez1 = genMuonGroups[0][0]->vertex().z()-beamSpot->position().z();

	      double gex2 = genMuonGroups[1][0]->vertex().x()-beamSpot->position().x();
	      double gey2 = genMuonGroups[1][0]->vertex().y()-beamSpot->position().y();
	      double gez2 = genMuonGroups[1][0]->vertex().z()-beamSpot->position().z();
	      
	      double ge_dzmj1 = (gez1) - ((gex1)*repx1+(gey1)*repy1)/rept1 * repz1/rept1;
	      double ge_dzmj2 = (gez2) - ((gex2)*repx2+(gey2)*repy2)/rept2 * repz2/rept2;

	      double ge_dxymj1 =( - (gex1) * repy1 + (gey1) * repx1 ) / rept1; 
	      double ge_dxymj2 =( - (gex2) * repy2 + (gey2) * repx2 ) / rept2; 

	      if (fabs(ge_dzmj1-ge_dzmj2) < 0.3) goodGenDz = true;

	      verbose = false;
	      if (verbose){
		std::cout << "ge_dzmj1->Fill(" << ge_dzmj1 << ");" << std::endl;
		std::cout << "ge_dzmj2->Fill(" << ge_dzmj2 << ");" << std::endl;
		std::cout << "ge_dzdiff->Fill(" << ge_dzmj1-ge_dzmj2 << ");" << std::endl;
		std::cout << "ge_dzdiff_abs->Fill(" << fabs(ge_dzmj1-ge_dzmj2) << ");" << std::endl;
		std::cout << "ge_dxymj1->Fill(" << ge_dxymj1 << ");" << std::endl;
		std::cout << "ge_dxymj2->Fill(" << ge_dxymj2 << ");" << std::endl;
		std::cout << "ge_dxydiff->Fill(" << ge_dxymj1-ge_dxymj2 << ");" << std::endl;


		if (fabs(genMuonGroups[1][0]->vertex().x())>4 && fabs(genMuonGroups[0][0]->vertex().x())>4){
		  std::cout << "ge_dzdiff_bothGreater4_x->Fill(" << ge_dzmj1-ge_dzmj2 << ");" << std::endl;
		  std::cout << "ge_dxydiff_bothGreater4_x->Fill(" << ge_dxymj1-ge_dxymj2 << ");" << std::endl;
		}
		if (fabs(genMuonGroups[1][0]->vertex().x())>4 && fabs(genMuonGroups[0][0]->vertex().x())<4){
		  std::cout << "ge_dzdiff_oneGreater4_x->Fill(" << ge_dzmj1-ge_dzmj2 << ");" << std::endl;
		  std::cout << "ge_dxydiff_oneGreater4_x->Fill(" << ge_dxymj1-ge_dxymj2 << ");" << std::endl;
		}
		if (fabs(genMuonGroups[1][0]->vertex().x())<4 && fabs(genMuonGroups[0][0]->vertex().x())>4){
		  std::cout << "ge_dzdiff_oneGreater4_x->Fill(" << ge_dzmj1-ge_dzmj2 << ");" << std::endl;
		  std::cout << "ge_dxydiff_oneGreater4_x->Fill(" << ge_dxymj1-ge_dxymj2 << ");" << std::endl;
		}
		if (fabs(genMuonGroups[1][0]->vertex().x())<4 && fabs(genMuonGroups[0][0]->vertex().x())<4){
		  std::cout << "ge_dzdiff_bothLess4_x->Fill(" << ge_dzmj1-ge_dzmj2 << ");" << std::endl;
		  std::cout << "ge_dxydiff_bothLess4_x->Fill(" << ge_dxymj1-ge_dxymj2 << ");" << std::endl;
		}

	      }
	      verbose = false;
	    }
	    else if (rec11 > 20 && rec22 < 20){

	      if (verbose){
		std::cout << "genMuonGroups[1][0]->vertex().x()-rnd_X1: " << genMuonGroups[1][0]->vertex().x()-rnd_X1 << std::endl;
		std::cout << "genMuonGroups[1][0]->vertex().y()-rnd_Y1: " << genMuonGroups[1][0]->vertex().y()-rnd_Y1 << std::endl;
		std::cout << "genMuonGroups[1][0]->vertex().z()-rnd_Z1: " << genMuonGroups[1][0]->vertex().z()-rnd_Z1 << std::endl;
		std::cout << "genMuonGroups[0][0]->vertex().x()-rnd_X2: " << genMuonGroups[0][0]->vertex().x()-rnd_X2 << std::endl;
		std::cout << "genMuonGroups[0][0]->vertex().y()-rnd_Y2: " << genMuonGroups[0][0]->vertex().y()-rnd_Y2 << std::endl;
		std::cout << "genMuonGroups[0][0]->vertex().z()-rnd_Z2: " << genMuonGroups[0][0]->vertex().z()-rnd_Z2 << std::endl;
	      }
	      double gex1 = genMuonGroups[1][0]->vertex().x()-beamSpot->position().x();
	      double gey1 = genMuonGroups[1][0]->vertex().y()-beamSpot->position().y();
	      double gez1 = genMuonGroups[1][0]->vertex().z()-beamSpot->position().z();

	      double gex2 = genMuonGroups[0][0]->vertex().x()-beamSpot->position().x();
	      double gey2 = genMuonGroups[0][0]->vertex().y()-beamSpot->position().y();
	      double gez2 = genMuonGroups[0][0]->vertex().z()-beamSpot->position().z();
	      
	      double ge_dzmj1 = (gez1) - ((gex1)*repx1+(gey1)*repy1)/rept1 * repz1/rept1;
	      double ge_dzmj2 = (gez2) - ((gex2)*repx2+(gey2)*repy2)/rept2 * repz2/rept2;

	      double ge_dxymj1 =( - (gex1) * repy1 + (gey1) * repx1 ) / rept1; 
	      double ge_dxymj2 =( - (gex2) * repy2 + (gey2) * repx2 ) / rept2; 

	      verbose = false;
	      if (verbose){
		std::cout << "ge_dzmj1->Fill(" << ge_dzmj1 << ");" << std::endl;
		std::cout << "ge_dzmj2->Fill(" << ge_dzmj2 << ");" << std::endl;
		std::cout << "ge_dzdiff->Fill(" << ge_dzmj1-ge_dzmj2 << ");" << std::endl;
		std::cout << "ge_dzdiff_abs->Fill(" << fabs(ge_dzmj1-ge_dzmj2) << ");" << std::endl;
		std::cout << "ge_dxymj1->Fill(" << ge_dxymj1 << ");" << std::endl;
		std::cout << "ge_dxymj2->Fill(" << ge_dxymj2 << ");" << std::endl;
		std::cout << "ge_dxydiff->Fill(" << ge_dxymj1-ge_dxymj2 << ");" << std::endl;

		if (fabs(genMuonGroups[1][0]->vertex().x())>4 && fabs(genMuonGroups[0][0]->vertex().x())>4){
		  std::cout << "ge_dzdiff_bothGreater4_x->Fill(" << ge_dzmj1-ge_dzmj2 << ");" << std::endl;
		  std::cout << "ge_dxydiff_bothGreater4_x->Fill(" << ge_dxymj1-ge_dxymj2 << ");" << std::endl;
		}
		if (fabs(genMuonGroups[1][0]->vertex().x())>4 && fabs(genMuonGroups[0][0]->vertex().x())<4){
		  std::cout << "ge_dzdiff_oneGreater4_x->Fill(" << ge_dzmj1-ge_dzmj2 << ");" << std::endl;
		  std::cout << "ge_dxydiff_oneGreater4_x->Fill(" << ge_dxymj1-ge_dxymj2 << ");" << std::endl;
		}
		if (fabs(genMuonGroups[1][0]->vertex().x())<4 && fabs(genMuonGroups[0][0]->vertex().x())>4){
		  std::cout << "ge_dzdiff_oneGreater4_x->Fill(" << ge_dzmj1-ge_dzmj2 << ");" << std::endl;
		  std::cout << "ge_dxydiff_oneGreater4_x->Fill(" << ge_dxymj1-ge_dxymj2 << ");" << std::endl;
		}
		if (fabs(genMuonGroups[1][0]->vertex().x())<4 && fabs(genMuonGroups[0][0]->vertex().x())<4){
		  std::cout << "ge_dzdiff_bothLess4_x->Fill(" << ge_dzmj1-ge_dzmj2 << ");" << std::endl;
		  std::cout << "ge_dxydiff_bothLess4_x->Fill(" << ge_dxymj1-ge_dxymj2 << ");" << std::endl;
		}


	      }
	      verbose = false;
	    }
	  }



	  dz_arr_mj1.push_back(re_dzmj1);
	  dz_arr_mj2.push_back(re_dzmj2);

	  dz_sum_mj1 += re_dzmj1;
	  dz_sum_mj2 += re_dzmj2;

	  dzdiff_abs_arr.push_back(fabs(re_dzmj1-re_dzmj2));
	  dzdiff_abs_sum += fabs(re_dzmj1-re_dzmj2);
	  dxydiff_abs_arr.push_back(fabs(dxy1-dxy2));
	  dxydiff_abs_sum += fabs(dxy1-dxy2);

	  dzdiff_arr.push_back(re_dzmj1-re_dzmj2);
	  dzdiff_sum += re_dzmj1-re_dzmj2;
	  dxydiff_arr.push_back(dxy1-dxy2);
	  dxydiff_sum += dxy1-dxy2;

	  dxy_arr_mj1.push_back(dxy1);
	  dxy_arr_mj2.push_back(dxy2);
	  dxy_sum_mj1 += dxy1;
	  dxy_sum_mj2 += dxy2;


	  // std::cout << "dzmj1->Fill(" << re_dzmj1 << ");" << std::endl;
	  // std::cout << "dzmj2->Fill(" << re_dzmj2 << ");" << std::endl;

	  // std::cout << "dz_diff->Fill(" << re_dzmj1-re_dzmj2 << ");" << std::endl;
	  // std::cout << "dxy_diff->Fill(" << dxy1-dxy2 << ");" << std::endl;

	  // std::cout << "dz_diff_abs->Fill(" << fabs(re_dzmj1-re_dzmj2) << ");" << std::endl;
	  // std::cout << "dxy_diff_abs->Fill(" << fabs(dxy1-dxy2) << ");" << std::endl;
 

	  // std::cout << "dxymj1->Fill(" << dxy1 << ");" << std::endl;
	  // std::cout << "dxymj2->Fill(" << dxy2 << ");" << std::endl;


 	  // double dxy1 =( - (vx()-beamSpot->position().x()) * py() + (vy()-beamSpot->position().y()) * px() ) / pt(); 

	  // std::cout << "dxy1: " << dxy1 << std::endl;
	  // std::cout << "dxy2: " << dxy2 << std::endl;


	  // if (e39995) std::cout << "fabs(re_dzmj1-re_dzmj2): " << fabs(re_dzmj1-re_dzmj2) << std::endl;


	  EventContainer thisEvent;
	  thisEvent.set_vtx_dzmj1(re_dzmj1);
	  thisEvent.set_vtx_dzmj2(re_dzmj2);
	  thisEvent.set_vtx_dxymj1(dxy1);
	  thisEvent.set_vtx_dxymj2(dxy2);
	  thisEvent.set_vtx_massmj1(new_mj1_vtx.M());
	  thisEvent.set_vtx_massmj2(new_mj2_vtx.M());
	  thisEvent.set_vtx_mj1(new_mj1_vtx);
	  thisEvent.set_vtx_mj2(new_mj2_vtx);
	  thisEvent.set_vtx_1(newVtx1);
	  thisEvent.set_vtx_2(newVtx2);

	  events_arr.push_back(thisEvent);
	  

	  if (fabs(re_dzmj1-re_dzmj2) < min_dz_n){

	    if (verbose){
	    std::cout << "************************* " << std::endl;
	    std::cout << "dice: " << dice << std::endl;
	    std::cout << "new min_dz_n: " << min_dz_n << std::endl;
	    std::cout << "old xyz1: " << muJet1->vertexPoint().x() << ", " << muJet1->vertexPoint().y() << ", " << muJet1->vertexPoint().z() << std::endl;
	    std::cout << "new xyz1: " << newx1 << ", " << newy1 << ", " << newz1 << std::endl;
	    std::cout << "gen xyz1: " << genMuonGroups[0][0]->vertex().x() << ", " << genMuonGroups[0][0]->vertex().y() << ", " << genMuonGroups[0][0]->vertex().z() << std::endl ;
	    std::cout << "old xyz2: " << muJet2->vertexPoint().x() << ", " << muJet2->vertexPoint().y() << ", " << muJet2->vertexPoint().z() << std::endl;
	    std::cout << "new xyz2: " << newx2 << ", " << newy2 << ", " << newz2 << std::endl;
	    std::cout << "gen xyz2: " << genMuonGroups[1][0]->vertex().x() << ", " << genMuonGroups[1][0]->vertex().y() << ", " << genMuonGroups[1][0]->vertex().z() << std::endl;
	    }

	    min_dice = dice;
	    // newx1 = rnd_X1;
	    // newy1 = rnd_Y1;
	    // newz1 = rnd_Z1;
	    // newx2 = rnd_X2;
	    // newy2 = rnd_Y2;
	    // newz2 = rnd_Z2;

	    // massmj1 = new_mj1_vtx.M();
	    // massmj2 = new_mj2_vtx.M();
	    // dzmj1 =  re_dzmj1;
	    // dzmj2 =  re_dzmj2;

	    // dxymj1 = dxy1;
	    // dxymj2 = dxy2;

	    min_dz_n = fabs(re_dzmj1-re_dzmj2);	    
	    min_mass_n = fabs(massmj1-massmj2);	    

	    // final_mj1_vtx = new_mj1_vtx;
	    // final_mj2_vtx = new_mj2_vtx;

	    // if (rec11<20 && rec22 > 20){

	    //   gen_match_x1diff = fabs(rnd_X1 - genMuonGroups[0][0]->vertex().x()); 
	    //   if (verbose) std::cout << "gen_match_x1diff: " << gen_match_x1diff << ", rnd_X1: " << rnd_X1 << ", genMuonGroups[0][0]->vertex().x(): " << genMuonGroups[0][0]->vertex().x() << std::endl;
	    //   gen_match_y1diff = fabs(rnd_Y1 - genMuonGroups[0][0]->vertex().y()); 
	    //   if (verbose) std::cout << "gen_match_y1diff: " << gen_match_y1diff << ", rnd_Y1: " << rnd_Y1 << ", genMuonGroups[0][0]->vertex().y(): " << genMuonGroups[0][0]->vertex().y() << std::endl;
	    //   gen_match_z1diff = fabs(rnd_Z1 - genMuonGroups[0][0]->vertex().z()); 
	    //   if (verbose) std::cout << "gen_match_z1diff: " << gen_match_z1diff << ", rnd_Z1: " << rnd_Z1 << ", genMuonGroups[0][0]->vertex().z(): " << genMuonGroups[0][0]->vertex().z() << std::endl;
	    //   gen_match_x2diff = fabs(rnd_X2 - genMuonGroups[1][0]->vertex().x()); 
	    //   if (verbose) std::cout << "gen_match_x2diff: " << gen_match_x2diff << ", rnd_X2: " << rnd_X2 << ", genMuonGroups[1][0]->vertex().x(): " << genMuonGroups[1][0]->vertex().x() << std::endl;
	    //   gen_match_y2diff = fabs(rnd_Y2 - genMuonGroups[1][0]->vertex().y()); 
	    //   if (verbose) std::cout << "gen_match_y2diff: " << gen_match_y2diff << ", rnd_Y2: " << rnd_Y2 << ", genMuonGroups[1][0]->vertex().y(): " << genMuonGroups[1][0]->vertex().y() << std::endl;
	    //   gen_match_z2diff = fabs(rnd_Z2 - genMuonGroups[1][0]->vertex().z()); 
	    //   if (verbose) std::cout << "gen_match_z2diff: " << gen_match_z2diff << ", rnd_Z2: " << rnd_Z2 << ", genMuonGroups[1][0]->vertex().z(): " << genMuonGroups[1][0]->vertex().z() << std::endl;
	    // }
	    // else if (rec11 > 20 && rec22 < 20){
	    //   gen_match_x1diff = fabs(rnd_X1 - genMuonGroups[1][0]->vertex().x()); 
	    //   if (verbose) std::cout << "gen_match_x1diff: " << gen_match_x1diff << ", rnd_X1: " << rnd_X1 << ", genMuonGroups[1][0]->vertex().x(): " << genMuonGroups[1][0]->vertex().x() << std::endl;
	    //   gen_match_y1diff = fabs(rnd_Y1 - genMuonGroups[1][0]->vertex().y()); 
	    //   if (verbose) std::cout << "gen_match_y1diff: " << gen_match_y1diff << ", rnd_Y1: " << rnd_Y1 << ", genMuonGroups[1][0]->vertex().y(): " << genMuonGroups[1][0]->vertex().y() << std::endl;
	    //   gen_match_z1diff = fabs(rnd_Z1 - genMuonGroups[1][0]->vertex().z()); 
	    //   if (verbose) std::cout << "gen_match_z1diff: " << gen_match_z1diff << ", rnd_Z1: " << rnd_Z1 << ", genMuonGroups[1][0]->vertex().z(): " << genMuonGroups[1][0]->vertex().z() << std::endl;
	    //   gen_match_x2diff = fabs(rnd_X2 - genMuonGroups[0][0]->vertex().x()); 
	    //   if (verbose) std::cout << "gen_match_x2diff: " << gen_match_x2diff << ", rnd_X2: " << rnd_X2 << ", genMuonGroups[0][0]->vertex().x(): " << genMuonGroups[0][0]->vertex().x() << std::endl;
	    //   gen_match_y2diff = fabs(rnd_Y2 - genMuonGroups[0][0]->vertex().y()); 
	    //   if (verbose) std::cout << "gen_match_y2diff: " << gen_match_y2diff << ", rnd_Y2: " << rnd_Y2 << ", genMuonGroups[0][0]->vertex().y(): " << genMuonGroups[0][0]->vertex().y() << std::endl;
	    //   gen_match_z2diff = fabs(rnd_Z2 - genMuonGroups[0][0]->vertex().z()); 
	    //   if (verbose) std::cout << "gen_match_z2diff: " << gen_match_z2diff << ", rnd_Z2: " << rnd_Z2 << ", genMuonGroups[0][0]->vertex().z(): " << genMuonGroups[0][0]->vertex().z() << std::endl;
	    // }


	    // if (min_dz_n < 0.5){

	    //   if (verbose){
	    // 	std::cout << " **************************************** " << std::endl;
	    // 	std::cout << " **************************************** " << std::endl;
	    // 	std::cout << " **************************************** " << std::endl;
	    // 	std::cout << " **************************************** " << std::endl;
	    // 	std::cout << " **************************************** " << std::endl;
	    //   }
	    //   for ( unsigned int i = 0; i <= 1; i++ ) { 
	    // 	if (verbose) std::cout << " ---------------------------------------- " << std::endl;

	    // 	double isoTmp = 0.0;
	    // 	double n_dz;
	    // 	if ( i == 0 ){
	    // 	  diMuonTmp = final_mj1_vtx;
	    // 	  n_dz = dzmj1;
	    // 	if (verbose) std::cout << "i = 0, mj1" << std::endl;
	    // 	}
	    // 	if ( i == 1 ){
	    // 	  diMuonTmp = final_mj2_vtx;
	    // 	  n_dz = dzmj2;
	    // 	  if (verbose) std::cout << "i = 0, mj2" << std::endl;
	    // 	}
	      
	    // 	for (reco::TrackCollection::const_iterator track = tracks->begin(); track != tracks->end(); ++track) {
	    // 	  bool trackIsMuon = false;
	    // 	  if ( sameTrack( &*track, &*((*muJets)[i].muon(0)->innerTrack()) ) || sameTrack( &*track, &*((*muJets)[i].muon(1)->innerTrack()) ) ) trackIsMuon = true;
	    // 	  if (!trackIsMuon) {
	    // 	    double dPhi = My_dPhi( diMuonTmp.phi(), track->phi() );
	    // 	    double dEta = diMuonTmp.eta() - track->eta();
	    // 	    double dR = sqrt( dPhi*dPhi + dEta*dEta ); 

	    // 	    if ( dR < 0.4 && track->pt() > 0.5 ) {
	    // 	      double dz = fabs( track->dz(beamSpot->position()) - n_dz);
	    // 	      if (verbose) std::cout << "u i: " << i << ", track->pt(): " << track->pt() << ", track->eta(): " << track->eta() << ", track->phi(): " << track->phi() << ", dPhi: " << dPhi << ", dEta: " << dEta << ", dR: " << dR << ", dz " << dz << std::endl;

	    // 	      if (verbose) std::cout << "u i: " << i << ", diMuonTmp.pt(): " << diMuonTmp.pt() << ", diMuonTmp.eta(): " << diMuonTmp.eta() << ", diMuonTmp.phi(): " << diMuonTmp.phi() << ", dPhi: " << dPhi << ", dEta: " << dEta << ", dR: " << dR << ", dz " << dz << ", n_dz: " << n_dz << std::endl;
		      
	    // 	      if ( dz < 0.1 ){
	    // 		isoTmp += track->pt();
	    // 		if (verbose) std::cout << "u isoTmp: " << isoTmp << std::endl;
	    // 	      }
	    // 	    }    
	    // 	  }
	    // 	}
	    // 	if ( i == 0 ) isomj1 = isoTmp;
	    // 	if ( i == 1 ) isomj2 = isoTmp;
	    //   }
	    // }
	  

	  }

	}

	if (verbose){
	  std::cout << "dz_sum_mj1: " << dz_sum_mj1 << std::endl;
	  std::cout << "dz_sum_mj2: " << dz_sum_mj2 << std::endl;
	}
	double dz_ave_mj1 = dz_sum_mj1/static_cast<double>(throws);
	double dz_ave_mj2 = dz_sum_mj2/static_cast<double>(throws);

	double dzdiff_ave = dzdiff_sum/static_cast<double>(throws);
	double dxydiff_ave = dxydiff_sum/static_cast<double>(throws);

	double dzdiff_abs_ave = dzdiff_abs_sum/static_cast<double>(throws);
	double dxydiff_abs_ave = dxydiff_abs_sum/static_cast<double>(throws);

	double dz_var_mj1 = 0.;
	double dz_var_mj2 = 0.;

	double dzdiff_var = 0.;
	double dxydiff_var = 0.;

	double dzdiff_abs_var = 0.;
	double dxydiff_abs_var = 0.;

	double dxy_ave_mj1 = dxy_sum_mj1/static_cast<double>(throws);
	double dxy_ave_mj2 = dxy_sum_mj2/static_cast<double>(throws);

	double dxy_var_mj1 = 0.;
	double dxy_var_mj2 = 0.;


	double chi2_dz1 = 0.;
	double chi2_dz2 = 0.;
	double chi2_dxy1 = 0.;
	double chi2_dxy2 = 0.;

	double chi2_dz = 0.;
	double chi2_dxy = 0.;
	double final_chi2_dz = 0.;
	double final_chi2_dxy = 0.;

	double final_chi2_dz_val = 0.;
	double final_chi2_dxy_val = 0.;


	double final_chi2_dz1 = 0.;
	double final_chi2_dz2 = 0.;
	double final_chi2_dxy1 = 0.;
	double final_chi2_dxy2 = 0.;



	for(int i = 0; i<throws; ++i) 
	  {
	    dz_var_mj1 += ((dz_arr_mj1[i]-dz_ave_mj1)*(dz_arr_mj1[i]-dz_ave_mj1))/throws;
	    dz_var_mj2 += ((dz_arr_mj2[i]-dz_ave_mj2)*(dz_arr_mj2[i]-dz_ave_mj2))/throws;

	    dzdiff_var += ((dzdiff_arr[i]-dzdiff_ave)*(dzdiff_arr[i]-dzdiff_ave))/throws;
	    dxydiff_var += ((dxydiff_arr[i]-dxydiff_ave)*(dxydiff_arr[i]-dxydiff_ave))/throws;

	    dzdiff_abs_var += ((dzdiff_abs_arr[i]-dzdiff_abs_ave)*(dzdiff_abs_arr[i]-dzdiff_abs_ave))/throws;
	    dxydiff_abs_var += ((dxydiff_abs_arr[i]-dxydiff_abs_ave)*(dxydiff_abs_arr[i]-dxydiff_abs_ave))/throws;

	    dxy_var_mj1 += ((dxy_arr_mj1[i]-dxy_ave_mj1)*(dxy_arr_mj1[i]-dxy_ave_mj1))/throws;
	    dxy_var_mj2 += ((dxy_arr_mj2[i]-dxy_ave_mj2)*(dxy_arr_mj2[i]-dxy_ave_mj2))/throws;


	  }

	double threshold1 = 10000.0;
	double threshold2 = 10000.0;


	for(int i = 0; i<throws; ++i) 
	  {

	    chi2_dz1 = ((dz_arr_mj1[i]-dz_ave_mj1)/(sqrt(dz_var_mj1)))*((dz_arr_mj1[i]-dz_ave_mj1)/(sqrt(dz_var_mj1)));
	    chi2_dz2 = ((dz_arr_mj2[i]-dz_ave_mj2)/(sqrt(dz_var_mj2)))*((dz_arr_mj2[i]-dz_ave_mj2)/(sqrt(dz_var_mj2)));

	    //	    chi2_dz = ((dzdiff_arr[i])/(0.0257))*((dzdiff_arr[i])/(0.0257));
	    //	    chi2_dxy = ((dxydiff_arr[i])/(0.00819594))*((dxydiff_arr[i])/(0.00819594));

	    if (fabs(events_arr[i].get_vtx_1().x())>4 && fabs(events_arr[i].get_vtx_2().x())>4){
	      chi2_dz = ((dzdiff_arr[i])/(0.035))*((dzdiff_arr[i])/(0.035));
	      chi2_dxy = ((dxydiff_arr[i])/(0.016))*((dxydiff_arr[i])/(0.016));
	    }
	    if (fabs(events_arr[i].get_vtx_1().x())>4 && fabs(events_arr[i].get_vtx_2().x())<4){
	      chi2_dz = ((dzdiff_arr[i])/(0.02513))*((dzdiff_arr[i])/(0.02513));
	      chi2_dxy = ((dxydiff_arr[i])/(0.007702))*((dxydiff_arr[i])/(0.007702));
	    }
	    if (fabs(events_arr[i].get_vtx_1().x())<4 && fabs(events_arr[i].get_vtx_2().x())>4){
	      chi2_dz = ((dzdiff_arr[i])/(0.02513))*((dzdiff_arr[i])/(0.02513));
	      chi2_dxy = ((dxydiff_arr[i])/(0.007702))*((dxydiff_arr[i])/(0.007702));
	    }
	    if (fabs(events_arr[i].get_vtx_1().x())<4 && fabs(events_arr[i].get_vtx_2().x())<4){
	      chi2_dz = ((dzdiff_arr[i])/(0.008912))*((dzdiff_arr[i])/(0.008912));
	      chi2_dxy = ((dxydiff_arr[i])/(0.004607))*((dxydiff_arr[i])/(0.004607));
	    }
	    chi2_dxy1 = ((dxy_arr_mj1[i])/(0.00284177))*((dxy_arr_mj1[i])/(0.00284177));
	    chi2_dxy2 = ((dxy_arr_mj2[i])/(0.00377043))*((dxy_arr_mj2[i])/(0.00377043));



	    bool singleMatch = false;

	    if (verbose)
	      {
		
	    if (rec11<20 && rec22 > 20){

	      if (singleMatch){
		if ((fabs(events_arr[i].get_vtx_1().x() - genMuonGroups[0][0]->vertex().x())<0.1) && (fabs(events_arr[i].get_vtx_1().y() - genMuonGroups[0][0]->vertex().y())<0.1) && (fabs(events_arr[i].get_vtx_1().z() - genMuonGroups[0][0]->vertex().z())<0.1)){

		  std::cout << "************MUJET1 CLOSE TO GEN************" << std::endl;

		  std::cout << "events_arr[i].get_vtx_1().x(): " << events_arr[i].get_vtx_1().x() << std::endl;
		  std::cout << "events_arr[i].get_vtx_1().y(): " << events_arr[i].get_vtx_1().y() << std::endl;
		  std::cout << "events_arr[i].get_vtx_1().z(): " << events_arr[i].get_vtx_1().z() << std::endl;
		  std::cout << "events_arr[i].get_vtx_2().x(): " << events_arr[i].get_vtx_2().x() << std::endl;
		  std::cout << "events_arr[i].get_vtx_2().y(): " << events_arr[i].get_vtx_2().y() << std::endl;
		  std::cout << "events_arr[i].get_vtx_2().z(): " << events_arr[i].get_vtx_2().z() << std::endl;

		  std::cout << "dz_arr_mj1[i]: " << dz_arr_mj1[i] << std::endl;
		  std::cout << "dz_arr_mj2[i]: " << dz_arr_mj2[i] << std::endl;
		  std::cout << "dzdiff_arr[i]: " << dzdiff_arr[i] << std::endl;
		  std::cout << "chi2_dz: " << chi2_dz << std::endl <<std::endl;

		  std::cout << "dxydiff_arr[i]: " << dxydiff_arr[i] << std::endl;
		  std::cout << "chi2_dxy: " << chi2_dxy << std::endl;
		  std::cout << "dxy1: " << dxy_arr_mj1[i] << std::endl;
		  std::cout << "dxy2: " << dxy_arr_mj2[i] << std::endl;
		  std::cout << "chi2_dxy1: " << chi2_dxy1 << std::endl;
		  std::cout << "chi2_dxy2: " << chi2_dxy2 << std::endl;
		}

		if ((fabs(events_arr[i].get_vtx_2().x() - genMuonGroups[1][0]->vertex().x())<0.1) && (fabs(events_arr[i].get_vtx_2().y() - genMuonGroups[1][0]->vertex().y())<0.1) && (fabs(events_arr[i].get_vtx_2().z() - genMuonGroups[1][0]->vertex().z())<0.1)){

		  std::cout << "************MUJET2 CLOSE TO GEN************" << std::endl;

		  std::cout << "events_arr[i].get_vtx_1().x(): " << events_arr[i].get_vtx_1().x() << std::endl;
		  std::cout << "events_arr[i].get_vtx_1().y(): " << events_arr[i].get_vtx_1().y() << std::endl;
		  std::cout << "events_arr[i].get_vtx_1().z(): " << events_arr[i].get_vtx_1().z() << std::endl;
		  std::cout << "events_arr[i].get_vtx_2().x(): " << events_arr[i].get_vtx_2().x() << std::endl;
		  std::cout << "events_arr[i].get_vtx_2().y(): " << events_arr[i].get_vtx_2().y() << std::endl;
		  std::cout << "events_arr[i].get_vtx_2().z(): " << events_arr[i].get_vtx_2().z() << std::endl;
		  std::cout << "dz_arr_mj1[i]: " << dz_arr_mj1[i] << std::endl;
		  std::cout << "dz_arr_mj2[i]: " << dz_arr_mj2[i] << std::endl;
		  std::cout << "dzdiff_arr[i]: " << dzdiff_arr[i] << std::endl;
		  std::cout << "chi2_dz: " << chi2_dz << std::endl <<std::endl;

		  std::cout << "dxydiff_arr[i]: " << dxydiff_arr[i] << std::endl;
		  std::cout << "chi2_dxy: " << chi2_dxy << std::endl;
		  std::cout << "dxy1: " << dxy_arr_mj1[i] << std::endl;
		  std::cout << "dxy2: " << dxy_arr_mj2[i] << std::endl;
		  std::cout << "chi2_dxy1: " << chi2_dxy1 << std::endl;
		  std::cout << "chi2_dxy2: " << chi2_dxy2 << std::endl;

		}
	      }

	      // if ((fabs(events_arr[i].get_vtx_2().x() - genMuonGroups[1][0]->vertex().x())<0.3) && (fabs(events_arr[i].get_vtx_2().y() - genMuonGroups[1][0]->vertex().y())<0.3) && (fabs(events_arr[i].get_vtx_2().z() - genMuonGroups[1][0]->vertex().z())<0.3) && (fabs(events_arr[i].get_vtx_1().x() - genMuonGroups[0][0]->vertex().x())<0.3) && (fabs(events_arr[i].get_vtx_1().y() - genMuonGroups[0][0]->vertex().y())<0.3) && (fabs(events_arr[i].get_vtx_1().z() - genMuonGroups[0][0]->vertex().z())<0.3) ){
	      if ((fabs(events_arr[i].get_vtx_2().x() - genMuonGroups[1][0]->vertex().x())<0.3) && (fabs(events_arr[i].get_vtx_2().y() - genMuonGroups[1][0]->vertex().y())<0.3) && (fabs(events_arr[i].get_vtx_2().z() - genMuonGroups[1][0]->vertex().z())<0.3) && (fabs(events_arr[i].get_vtx_1().x() - genMuonGroups[0][0]->vertex().x())<0.3) && (fabs(events_arr[i].get_vtx_1().y() - genMuonGroups[0][0]->vertex().y())<0.3) && (fabs(events_arr[i].get_vtx_1().z() - genMuonGroups[0][0]->vertex().z())<0.3) && goodGenDz && !foundEvent ){

		std::cout << "************BOTH JETS CLOSE TO GEN************" << std::endl;

		if (goodGenDz) std::cout << "gen dz < 0.3" << std::endl;
		if (!goodGenDz) std::cout << "gen dz > 0.3" << std::endl;

		std::cout << "events_arr[i].get_vtx_1().x(): " << events_arr[i].get_vtx_1().x() << std::endl;
		std::cout << "events_arr[i].get_vtx_1().y(): " << events_arr[i].get_vtx_1().y() << std::endl;
		std::cout << "events_arr[i].get_vtx_1().z(): " << events_arr[i].get_vtx_1().z() << std::endl;
		std::cout << "events_arr[i].get_vtx_2().x(): " << events_arr[i].get_vtx_2().x() << std::endl;
		std::cout << "events_arr[i].get_vtx_2().y(): " << events_arr[i].get_vtx_2().y() << std::endl;
		std::cout << "events_arr[i].get_vtx_2().z(): " << events_arr[i].get_vtx_2().z() << std::endl;
		std::cout << "dz_arr_mj1[i]: " << dz_arr_mj1[i] << std::endl;
		std::cout << "dz_arr_mj2[i]: " << dz_arr_mj2[i] << std::endl;
		std::cout << "dzdiff_arr[i]: " << dzdiff_arr[i] << std::endl;
		std::cout << "chi2_dz: " << chi2_dz << std::endl <<std::endl;

		std::cout << "dxydiff_arr[i]: " << dxydiff_arr[i] << std::endl;
		std::cout << "chi2_dxy: " << chi2_dxy << std::endl;
		std::cout << "dxy1: " << dxy_arr_mj1[i] << std::endl;
		std::cout << "dxy2: " << dxy_arr_mj2[i] << std::endl;
		std::cout << "chi2_dxy1: " << chi2_dxy1 << std::endl;
		std::cout << "chi2_dxy2: " << chi2_dxy2 << std::endl;

		std::cout << "genMuonGroups[0][0]->vertex().x()-events_arr[i].get_vtx_1().x(): " << genMuonGroups[0][0]->vertex().x()-events_arr[i].get_vtx_1().x() << std::endl;
		std::cout << "genMuonGroups[0][0]->vertex().y()-events_arr[i].get_vtx_1().y(): " << genMuonGroups[0][0]->vertex().y()-events_arr[i].get_vtx_1().y() << std::endl;
		std::cout << "genMuonGroups[0][0]->vertex().z()-events_arr[i].get_vtx_1().z(): " << genMuonGroups[0][0]->vertex().z()-events_arr[i].get_vtx_1().z() << std::endl;
		std::cout << "genMuonGroups[1][0]->vertex().x()-events_arr[i].get_vtx_2().x(): " << genMuonGroups[1][0]->vertex().x()-events_arr[i].get_vtx_2().x() << std::endl;
		std::cout << "genMuonGroups[1][0]->vertex().y()-events_arr[i].get_vtx_2().y(): " << genMuonGroups[1][0]->vertex().y()-events_arr[i].get_vtx_2().y() << std::endl;
		std::cout << "genMuonGroups[1][0]->vertex().z()-events_arr[i].get_vtx_2().z(): " << genMuonGroups[1][0]->vertex().z()-events_arr[i].get_vtx_2().z() << std::endl;

		foundEvent = true;
	      }
	    }
	    else if (rec11 > 20 && rec22 < 20){

	      if (singleMatch){
		if ((fabs(events_arr[i].get_vtx_1().x() - genMuonGroups[1][0]->vertex().x())<0.1) && (fabs(events_arr[i].get_vtx_1().y() - genMuonGroups[1][0]->vertex().y())<0.1) && (fabs(events_arr[i].get_vtx_1().z() - genMuonGroups[1][0]->vertex().z())<0.1)){

		  std::cout << "************MUJET1 CLOSE TO GEN************" << std::endl;

		  std::cout << "events_arr[i].get_vtx_1().x(): " << events_arr[i].get_vtx_1().x() << std::endl;
		  std::cout << "events_arr[i].get_vtx_1().y(): " << events_arr[i].get_vtx_1().y() << std::endl;
		  std::cout << "events_arr[i].get_vtx_1().z(): " << events_arr[i].get_vtx_1().z() << std::endl;
		  std::cout << "events_arr[i].get_vtx_2().x(): " << events_arr[i].get_vtx_2().x() << std::endl;
		  std::cout << "events_arr[i].get_vtx_2().y(): " << events_arr[i].get_vtx_2().y() << std::endl;
		  std::cout << "events_arr[i].get_vtx_2().z(): " << events_arr[i].get_vtx_2().z() << std::endl;
		  std::cout << "dz_arr_mj1[i]: " << dz_arr_mj1[i] << std::endl;
		  std::cout << "dz_arr_mj2[i]: " << dz_arr_mj2[i] << std::endl;
		  std::cout << "dzdiff_arr[i]: " << dzdiff_arr[i] << std::endl;
		  std::cout << "chi2_dz: " << chi2_dz << std::endl <<std::endl;

		  std::cout << "dxydiff_arr[i]: " << dxydiff_arr[i] << std::endl;
		  std::cout << "chi2_dxy: " << chi2_dxy << std::endl;
		  std::cout << "dxy1: " << dxy_arr_mj1[i] << std::endl;
		  std::cout << "dxy2: " << dxy_arr_mj2[i] << std::endl;
		  std::cout << "chi2_dxy1: " << chi2_dxy1 << std::endl;
		  std::cout << "chi2_dxy2: " << chi2_dxy2 << std::endl;
		}

		if ((fabs(events_arr[i].get_vtx_2().x() - genMuonGroups[0][0]->vertex().x())<0.1) && (fabs(events_arr[i].get_vtx_2().y() - genMuonGroups[0][0]->vertex().y())<0.1) && (fabs(events_arr[i].get_vtx_2().z() - genMuonGroups[0][0]->vertex().z())<0.1)){

		  std::cout << "************MUJET2 CLOSE TO GEN************" << std::endl;

		  std::cout << "events_arr[i].get_vtx_1().x(): " << events_arr[i].get_vtx_1().x() << std::endl;
		  std::cout << "events_arr[i].get_vtx_1().y(): " << events_arr[i].get_vtx_1().y() << std::endl;
		  std::cout << "events_arr[i].get_vtx_1().z(): " << events_arr[i].get_vtx_1().z() << std::endl;
		  std::cout << "events_arr[i].get_vtx_2().x(): " << events_arr[i].get_vtx_2().x() << std::endl;
		  std::cout << "events_arr[i].get_vtx_2().y(): " << events_arr[i].get_vtx_2().y() << std::endl;
		  std::cout << "events_arr[i].get_vtx_2().z(): " << events_arr[i].get_vtx_2().z() << std::endl;
		  std::cout << "dz_arr_mj1[i]: " << dz_arr_mj1[i] << std::endl;
		  std::cout << "dz_arr_mj2[i]: " << dz_arr_mj2[i] << std::endl;
		  std::cout << "dzdiff_arr[i]: " << dzdiff_arr[i] << std::endl;
		  std::cout << "chi2_dz: " << chi2_dz << std::endl <<std::endl;

		  std::cout << "dxydiff_arr[i]: " << dxydiff_arr[i] << std::endl;
		  std::cout << "chi2_dxy: " << chi2_dxy << std::endl;
		  std::cout << "dxy1: " << dxy_arr_mj1[i] << std::endl;
		  std::cout << "dxy2: " << dxy_arr_mj2[i] << std::endl;
		  std::cout << "chi2_dxy1: " << chi2_dxy1 << std::endl;
		  std::cout << "chi2_dxy2: " << chi2_dxy2 << std::endl;

		}
	      }


	      if ((fabs(events_arr[i].get_vtx_2().x() - genMuonGroups[0][0]->vertex().x())<0.1) && (fabs(events_arr[i].get_vtx_2().y() - genMuonGroups[0][0]->vertex().y())<0.1) && (fabs(events_arr[i].get_vtx_2().z() - genMuonGroups[0][0]->vertex().z())<0.1) && (fabs(events_arr[i].get_vtx_1().x() - genMuonGroups[1][0]->vertex().x())<0.1) && (fabs(events_arr[i].get_vtx_1().y() - genMuonGroups[1][0]->vertex().y())<0.1) && (fabs(events_arr[i].get_vtx_1().z() - genMuonGroups[1][0]->vertex().z())<0.1) && goodGenDz && !foundEvent ){

		// if ((fabs(events_arr[i].get_vtx_2().x() - genMuonGroups[0][0]->vertex().x())<0.1) && (fabs(events_arr[i].get_vtx_2().y() - genMuonGroups[0][0]->vertex().y())<0.1) && (fabs(events_arr[i].get_vtx_2().z() - genMuonGroups[0][0]->vertex().z())<0.1) && (fabs(events_arr[i].get_vtx_1().x() - genMuonGroups[1][0]->vertex().x())<0.1) && (fabs(events_arr[i].get_vtx_1().y() - genMuonGroups[1][0]->vertex().y())<0.1) && (fabs(events_arr[i].get_vtx_1().z() - genMuonGroups[1][0]->vertex().z())<0.1) ){

		std::cout << "************BOTH JETS CLOSE TO GEN************" << std::endl;

		if (goodGenDz) std::cout << "gen dz < 0.3" << std::endl;
		if (!goodGenDz) std::cout << "gen dz > 0.3" << std::endl;

		std::cout << "events_arr[i].get_vtx_1().x(): " << events_arr[i].get_vtx_1().x() << std::endl;
		std::cout << "events_arr[i].get_vtx_1().y(): " << events_arr[i].get_vtx_1().y() << std::endl;
		std::cout << "events_arr[i].get_vtx_1().z(): " << events_arr[i].get_vtx_1().z() << std::endl;
		std::cout << "events_arr[i].get_vtx_2().x(): " << events_arr[i].get_vtx_2().x() << std::endl;
		std::cout << "events_arr[i].get_vtx_2().y(): " << events_arr[i].get_vtx_2().y() << std::endl;
		std::cout << "events_arr[i].get_vtx_2().z(): " << events_arr[i].get_vtx_2().z() << std::endl;
		std::cout << "dz_arr_mj1[i]: " << dz_arr_mj1[i] << std::endl;
		std::cout << "dz_arr_mj2[i]: " << dz_arr_mj2[i] << std::endl;
		std::cout << "dzdiff_arr[i]: " << dzdiff_arr[i] << std::endl;
		std::cout << "chi2_dz: " << chi2_dz << std::endl <<std::endl;

		std::cout << "dxydiff_arr[i]: " << dxydiff_arr[i] << std::endl;
		std::cout << "chi2_dxy: " << chi2_dxy << std::endl;
		std::cout << "dxy1: " << dxy_arr_mj1[i] << std::endl;
		std::cout << "dxy2: " << dxy_arr_mj2[i] << std::endl;
		std::cout << "chi2_dxy1: " << chi2_dxy1 << std::endl;
		std::cout << "chi2_dxy2: " << chi2_dxy2 << std::endl;


		std::cout << "genMuonGroups[1][0]->vertex().x()-events_arr[i].get_vtx_1().x(): " << genMuonGroups[1][0]->vertex().x()-events_arr[i].get_vtx_1().x() << std::endl;
		std::cout << "genMuonGroups[1][0]->vertex().y()-events_arr[i].get_vtx_1().y(): " << genMuonGroups[1][0]->vertex().y()-events_arr[i].get_vtx_1().y() << std::endl;
		std::cout << "genMuonGroups[1][0]->vertex().z()-events_arr[i].get_vtx_1().z(): " << genMuonGroups[1][0]->vertex().z()-events_arr[i].get_vtx_1().z() << std::endl;
		std::cout << "genMuonGroups[0][0]->vertex().x()-events_arr[i].get_vtx_2().x(): " << genMuonGroups[0][0]->vertex().x()-events_arr[i].get_vtx_2().x() << std::endl;
		std::cout << "genMuonGroups[0][0]->vertex().y()-events_arr[i].get_vtx_2().y(): " << genMuonGroups[0][0]->vertex().y()-events_arr[i].get_vtx_2().y() << std::endl;
		std::cout << "genMuonGroups[0][0]->vertex().z()-events_arr[i].get_vtx_2().z(): " << genMuonGroups[0][0]->vertex().z()-events_arr[i].get_vtx_2().z() << std::endl;

		foundEvent = true;
	      }




	    }
	      }
	    
	  
	    




	    if (verbose){
	      std::cout << " *****************" << std::endl;
	      std::cout << "dz_arr_mj1[i]: " << dz_arr_mj1[i] << std::endl;
	      std::cout << "dz_arr_mj2[i]: " << dz_arr_mj2[i] << std::endl;
	      std::cout << "dzdiff_arr[i]: " << dzdiff_arr[i] << std::endl;
	      std::cout << "((dzdiff_arr[i])/(0.0258375)): " << ((dzdiff_arr[i])/(0.0258375)) << std::endl;
	      std::cout << "chi2_dz: " << chi2_dz << std::endl;
	    }

	    verbose = false;
	    if (verbose){
	      std::cout << "dxy1_vs_chi2->Fill(" << dxy_arr_mj1[i] << "," << chi2_dxy1 << ");" << std::endl;
	      std::cout << "dxy2_vs_chi2->Fill(" << dxy_arr_mj2[i] << "," << chi2_dxy2 << ");" << std::endl;
	      std::cout << "dxy1_dxy2_chi2->Fill(" << dxy_arr_mj1[i] << "," << dxy_arr_mj2[i] << "," << chi2_dxy1+chi2_dxy2 << ");" << std::endl;
	      std::cout << "dxy1_plus_dxy2_chi2->Fill(" << dxy_arr_mj1[i]+dxy_arr_mj2[i] << "," << chi2_dxy1+chi2_dxy2 << ");" << std::endl;
	    }
	    verbose = false;


	    // std::cout << "dxy_vs_chi2->Fill(" << dxydiff_arr[i] << "," << chi2_dxy << ");" << std::endl;

	    // std::cout << "dz_dxy_chi2sum->SetBinContent(" << dzdiff_arr[i] << "," << dxydiff_arr[i] << "," << chi2_dxy+chi2_dz << ");" << std::endl;


	    verbose = false;
	    if (verbose){
	      std::cout << "dz_vs_chi2sum->Fill(" << dzdiff_arr[i] << "," << chi2_dz+chi2_dxy << ");" << std::endl;
	      std::cout << "dz_vs_chi2_dxy->Fill(" << dzdiff_arr[i] << "," << chi2_dxy << ");" << std::endl;
	      std::cout << "dz_vs_chi2_dz->Fill(" << dzdiff_arr[i] << "," << chi2_dz << ");" << std::endl;
	      std::cout << "dxy_vs_chi2->Fill(" << dxydiff_arr[i] << "," << chi2_dxy << ");" << std::endl;
	      std::cout << "mass_diff_vs_chi2_dz->Fill(" << events_arr[i].get_vtx_massmj1()-events_arr[i].get_vtx_massmj2() << "," << chi2_dz << ");" << std::endl;
	      std::cout << "mass_diff_vs_chi2sum->Fill(" << events_arr[i].get_vtx_massmj1()-events_arr[i].get_vtx_massmj2() << "," << chi2_dz+chi2_dxy << ");" << std::endl;

	    }
	    verbose = false;
	    // std::cout << "dz_vs_chi2->Fill(" << dzdiff_arr[i] << "," << chi2_dz1+chi2_dz2 << ");" << std::endl;

	    // chi2_dz = ((dzdiff_arr[i]-dzdiff_ave)/(sqrt(dzdiff_var)))*((dzdiff_arr[i]-dzdiff_ave)/(sqrt(dzdiff_var)));


	    // if (chi2_dz+chi2_dxy < threshold1 && dzdiff_arr[i] > -1 && dzdiff_arr[i] < 1){
	    // if (chi2_dz < threshold1){


	    if (chi2_dxy < threshold2){
	      threshold2 = chi2_dxy;
	      verbose = false;
	      if (verbose){
		std::cout << "******CHI2DXYMIN*****************" << std::endl;
		std::cout << "chi2_dxy: " << chi2_dxy << std::endl;
		std::cout << "((dxydiff_arr[i])/(0.00819594)): " << ((dxydiff_arr[i])/(0.00819594)) << std::endl;
		std::cout << "(dxydiff_arr[i]): " << (dxydiff_arr[i]) << std::endl << std::endl;

		std::cout << "chi2_dxy1: " << chi2_dxy1 << std::endl;
		std::cout << "((dxy_arr_mj1[i])/(0.00284177)): " << ((dxy_arr_mj1[i])/(0.00284177)) << std::endl;
		std::cout << "(dxy_arr_mj1[i]): " << (dxy_arr_mj1[i]) << std::endl << std::endl;

		std::cout << "chi2_dxy2: " << chi2_dxy2 << std::endl;
		std::cout << "((dxy_arr_mj2[i])/(0.00377043)): " << ((dxy_arr_mj2[i])/(0.00377043)) << std::endl;
		std::cout << "(dxy_arr_mj2[i]): " << (dxy_arr_mj2[i]) << std::endl << std::endl;

		std::cout << "dz_arr_mj1[i]: " << dz_arr_mj1[i] << std::endl;
		std::cout << "dz_arr_mj2[i]: " << dz_arr_mj2[i] << std::endl;
		std::cout << "dzdiff_arr[i]: " << dzdiff_arr[i] << std::endl;
		std::cout << "((dzdiff_arr[i])/(0.0258375)): " << ((dzdiff_arr[i])/(0.0258375)) << std::endl;
		std::cout << "chi2_dz: " << chi2_dz << std::endl;
		std::cout << "***********************" << std::endl;
	      }
	      verbose = false;
	    }

 	    // if (chi2_dz < threshold1){
	    if (chi2_dz+chi2_dxy < threshold1){

	      // threshold1 = chi2_dz;
	      threshold1 = chi2_dz+chi2_dxy;

	      final_chi2_dxy1 = chi2_dxy1;
	      final_chi2_dxy2 = chi2_dxy2;

	      final_chi2_dxy = chi2_dxy;
	      final_chi2_dz = chi2_dz;
	      final_chi2_dz_val = dzdiff_arr[i];
	      final_chi2_dxy_val = dxydiff_arr[i];
	      verbose = false;
	      if (verbose){
		std::cout << "******* MIN **********" << std::endl;
		std::cout << "******* MIN **********" << std::endl;
		std::cout << "******* MIN **********" << std::endl;
		std::cout << "******* MIN **********" << std::endl;
		std::cout << "dzdiff_arr[i]: " << dzdiff_arr[i] << std::endl;
		std::cout << "dxydiff_arr[i]: " << dxydiff_arr[i] << std::endl;
		std::cout << "chi2_dz: " << chi2_dz << std::endl;
		std::cout << "chi2_dxy: " << chi2_dxy << std::endl;
		std::cout << "xyz 1: " << events_arr[i].get_vtx_1().x() << ", " <<events_arr[i].get_vtx_1().y() << ", " << events_arr[i].get_vtx_1().z() << std::endl;
		std::cout << "xyz 2: " << events_arr[i].get_vtx_2().x() << ", " <<events_arr[i].get_vtx_2().y() << ", " << events_arr[i].get_vtx_2().z() << std::endl<<std::endl;
	      }
	      verbose = false;

	      massmj1 = events_arr[i].get_vtx_massmj1();
	      massmj2 = events_arr[i].get_vtx_massmj2();
	      dzmj1 =  events_arr[i].get_vtx_dzmj1();
	      dzmj2 =  events_arr[i].get_vtx_dzmj2();
	      dxymj1 =  events_arr[i].get_vtx_dxymj1();
	      dxymj2 =  events_arr[i].get_vtx_dxymj2();
	      
	      final_mj1_vtx = events_arr[i].get_vtx_mj1();
	      final_mj2_vtx = events_arr[i].get_vtx_mj2();

	      newx1 = events_arr[i].get_vtx_1().x();
	      newy1 = events_arr[i].get_vtx_1().y();
	      newz1 = events_arr[i].get_vtx_1().z();
	      newx2 = events_arr[i].get_vtx_2().x();
	      newy2 = events_arr[i].get_vtx_2().y();
	      newz2 = events_arr[i].get_vtx_2().z();


	      verbose = false;
	      if (rec11<20 && rec22 > 20){

		gen_match_x1diff = fabs(events_arr[i].get_vtx_1().x() - genMuonGroups[0][0]->vertex().x()); 
		if (verbose) std::cout << "gen_match_x1diff: " << gen_match_x1diff << ", events_arr[i].get_vtx_1().x(): " << events_arr[i].get_vtx_1().x() << ", genMuonGroups[0][0]->vertex().x(): " << genMuonGroups[0][0]->vertex().x() << std::endl;
		gen_match_y1diff = fabs(events_arr[i].get_vtx_1().y() - genMuonGroups[0][0]->vertex().y()); 
		if (verbose) std::cout << "gen_match_y1diff: " << gen_match_y1diff << ", events_arr[i].get_vtx_1().y(): " << events_arr[i].get_vtx_1().y() << ", genMuonGroups[0][0]->vertex().y(): " << genMuonGroups[0][0]->vertex().y() << std::endl;
		gen_match_z1diff = fabs(events_arr[i].get_vtx_1().z() - genMuonGroups[0][0]->vertex().z()); 
		if (verbose) std::cout << "gen_match_z1diff: " << gen_match_z1diff << ", events_arr[i].get_vtx_1().z(): " << events_arr[i].get_vtx_1().z() << ", genMuonGroups[0][0]->vertex().z(): " << genMuonGroups[0][0]->vertex().z() << std::endl;
		gen_match_x2diff = fabs(events_arr[i].get_vtx_2().x() - genMuonGroups[1][0]->vertex().x()); 
		if (verbose) std::cout << "gen_match_x2diff: " << gen_match_x2diff << ", events_arr[i].get_vtx_2().x(): " << events_arr[i].get_vtx_2().x() << ", genMuonGroups[1][0]->vertex().x(): " << genMuonGroups[1][0]->vertex().x() << std::endl;
		gen_match_y2diff = fabs(events_arr[i].get_vtx_2().y() - genMuonGroups[1][0]->vertex().y()); 
		if (verbose) std::cout << "gen_match_y2diff: " << gen_match_y2diff << ", events_arr[i].get_vtx_2().y(): " << events_arr[i].get_vtx_2().y() << ", genMuonGroups[1][0]->vertex().y(): " << genMuonGroups[1][0]->vertex().y() << std::endl;
		gen_match_z2diff = fabs(events_arr[i].get_vtx_2().z() - genMuonGroups[1][0]->vertex().z()); 
		if (verbose) std::cout << "gen_match_z2diff: " << gen_match_z2diff << ", events_arr[i].get_vtx_2().z(): " << events_arr[i].get_vtx_2().z() << ", genMuonGroups[1][0]->vertex().z(): " << genMuonGroups[1][0]->vertex().z() << std::endl;
	      }
	      else if (rec11 > 20 && rec22 < 20){
		gen_match_x1diff = fabs(events_arr[i].get_vtx_1().x() - genMuonGroups[1][0]->vertex().x()); 
		if (verbose) std::cout << "gen_match_x1diff: " << gen_match_x1diff << ", events_arr[i].get_vtx_1().x(): " << events_arr[i].get_vtx_1().x() << ", genMuonGroups[1][0]->vertex().x(): " << genMuonGroups[1][0]->vertex().x() << std::endl;
		gen_match_y1diff = fabs(events_arr[i].get_vtx_1().y() - genMuonGroups[1][0]->vertex().y()); 
		if (verbose) std::cout << "gen_match_y1diff: " << gen_match_y1diff << ", events_arr[i].get_vtx_1().y(): " << events_arr[i].get_vtx_1().y() << ", genMuonGroups[1][0]->vertex().y(): " << genMuonGroups[1][0]->vertex().y() << std::endl;
		gen_match_z1diff = fabs(events_arr[i].get_vtx_1().z() - genMuonGroups[1][0]->vertex().z()); 
		if (verbose) std::cout << "gen_match_z1diff: " << gen_match_z1diff << ", events_arr[i].get_vtx_1().z(): " << events_arr[i].get_vtx_1().z() << ", genMuonGroups[1][0]->vertex().z(): " << genMuonGroups[1][0]->vertex().z() << std::endl;
		gen_match_x2diff = fabs(events_arr[i].get_vtx_2().x() - genMuonGroups[0][0]->vertex().x()); 
		if (verbose) std::cout << "gen_match_x2diff: " << gen_match_x2diff << ", events_arr[i].get_vtx_2().x(): " << events_arr[i].get_vtx_2().x() << ", genMuonGroups[0][0]->vertex().x(): " << genMuonGroups[0][0]->vertex().x() << std::endl;
		gen_match_y2diff = fabs(events_arr[i].get_vtx_2().y() - genMuonGroups[0][0]->vertex().y()); 
		if (verbose) std::cout << "gen_match_y2diff: " << gen_match_y2diff << ", events_arr[i].get_vtx_2().y(): " << events_arr[i].get_vtx_2().y() << ", genMuonGroups[0][0]->vertex().y(): " << genMuonGroups[0][0]->vertex().y() << std::endl;
		gen_match_z2diff = fabs(events_arr[i].get_vtx_2().z() - genMuonGroups[0][0]->vertex().z()); 
		if (verbose) std::cout << "gen_match_z2diff: " << gen_match_z2diff << ", events_arr[i].get_vtx_2().z(): " << events_arr[i].get_vtx_2().z() << ", genMuonGroups[0][0]->vertex().z(): " << genMuonGroups[0][0]->vertex().z() << std::endl;
	      }


	      verbose = false;
	      if (verbose){
		std::cout << "events_arr[i].get_vtx_dzmj1(): " << events_arr[i].get_vtx_dzmj1() << std::endl;
		std::cout << "events_arr[i].get_vtx_dzmj2(): " << events_arr[i].get_vtx_dzmj2() << std::endl;
		std::cout << "diff: " << fabs(events_arr[i].get_vtx_dzmj1()-events_arr[i].get_vtx_dzmj2()) << std::endl;

		std::cout << "events_arr[i].get_vtx_dxymj1(): " << events_arr[i].get_vtx_dxymj1() << std::endl;
		std::cout << "events_arr[i].get_vtx_dxymj2(): " << events_arr[i].get_vtx_dxymj2() << std::endl;
		std::cout << "diff: " << fabs(events_arr[i].get_vtx_dxymj1()-events_arr[i].get_vtx_dxymj2()) << std::endl;

		std::cout << "events_arr[i].get_vtx_massmj1(): " << events_arr[i].get_vtx_massmj1() << std::endl;
		std::cout << "events_arr[i].get_vtx_massmj2(): " << events_arr[i].get_vtx_massmj2() << std::endl;
		std::cout << "diff: " << fabs(events_arr[i].get_vtx_massmj1()-events_arr[i].get_vtx_massmj2()) << std::endl;

		std::cout << "-----------------------------------------" << std::endl;
		std::cout << "events_arr[min_dice].get_vtx_dzmj1(): " << events_arr[min_dice].get_vtx_dzmj1() << std::endl;
		std::cout << "events_arr[min_dice].get_vtx_dzmj2(): " << events_arr[min_dice].get_vtx_dzmj2() << std::endl;
		std::cout << "diff: " << fabs(events_arr[min_dice].get_vtx_dzmj1()-events_arr[min_dice].get_vtx_dzmj2()) << std::endl;
		std::cout << "chi2 dz (min_dice): " << ((dzdiff_arr[min_dice]-dzdiff_ave)/(0.0258375))*((dzdiff_arr[min_dice]-dzdiff_ave)/(0.0258375)) << std::endl;
		std::cout << "dzdiff_ave: " << dzdiff_ave << std::endl;
		std::cout << "-----------------------------------------" << std::endl;
	      }


	      for ( unsigned int j = 0; j <= 1; j++ ) { 
		if (verbose) std::cout << " ---------------------------------------- " << std::endl;

		double isoTmp = 0.0;
		double n_dz;
		if ( j == 0 ){
		  diMuonTmp = final_mj1_vtx;
		  n_dz = dzmj1;
		  if (verbose) std::cout << "i = 0, mj1" << std::endl;
		}
		if ( j == 1 ){
		  diMuonTmp = final_mj2_vtx;
		  n_dz = dzmj2;
		  if (verbose) std::cout << "i = 0, mj2" << std::endl;
		}
	      
		for (reco::TrackCollection::const_iterator track = tracks->begin(); track != tracks->end(); ++track) {
		  bool trackIsMuon = false;
		  if ( sameTrack( &*track, &*((*muJets)[j].muon(0)->innerTrack()) ) || sameTrack( &*track, &*((*muJets)[j].muon(1)->innerTrack()) ) ) trackIsMuon = true;
		  if (!trackIsMuon) {
		    double dPhi = My_dPhi( diMuonTmp.phi(), track->phi() );
		    double dEta = diMuonTmp.eta() - track->eta();
		    double dR = sqrt( dPhi*dPhi + dEta*dEta ); 

		    if ( dR < 0.4 && track->pt() > 0.5 ) {
		      double dz = fabs( track->dz(beamSpot->position()) - n_dz);
		      if (verbose) std::cout << "u j: " << j << ", track->pt(): " << track->pt() << ", track->eta(): " << track->eta() << ", track->phi(): " << track->phi() << ", dPhi: " << dPhi << ", dEta: " << dEta << ", dR: " << dR << ", dz " << dz << std::endl;

		      if (verbose) std::cout << "u j: " << j << ", diMuonTmp.pt(): " << diMuonTmp.pt() << ", diMuonTmp.eta(): " << diMuonTmp.eta() << ", diMuonTmp.phi(): " << diMuonTmp.phi() << ", dPhi: " << dPhi << ", dEta: " << dEta << ", dR: " << dR << ", dz " << dz << ", n_dz: " << n_dz << std::endl;
		      
		      if ( dz < 0.1 ){
			isoTmp += track->pt();
			if (verbose) std::cout << "u isoTmp: " << isoTmp << std::endl;
		      }
		    }    
		  }
		}
		if ( j == 0 ) isomj1 = isoTmp;
		if ( j == 1 ) isomj2 = isoTmp;
	      }

	      if (verbose){
		std::cout << "events_arr[i].get_vtx_dzmj1(): " << events_arr[i].get_vtx_dzmj1() << std::endl;
		std::cout << "events_arr[i].get_vtx_dzmj2(): " << events_arr[i].get_vtx_dzmj2() << std::endl;
		std::cout << "diff: " << fabs(events_arr[i].get_vtx_dzmj1()-events_arr[i].get_vtx_dzmj2()) << std::endl;

		std::cout << "events_arr[i].get_vtx_dxymj1(): " << events_arr[i].get_vtx_dxymj1() << std::endl;
		std::cout << "events_arr[i].get_vtx_dxymj2(): " << events_arr[i].get_vtx_dxymj2() << std::endl;
		std::cout << "diff: " << fabs(events_arr[i].get_vtx_dxymj1()-events_arr[i].get_vtx_dxymj2()) << std::endl;

		std::cout << "events_arr[i].get_vtx_massmj1(): " << events_arr[i].get_vtx_massmj1() << std::endl;
		std::cout << "events_arr[i].get_vtx_massmj2(): " << events_arr[i].get_vtx_massmj2() << std::endl;
		std::cout << "diff: " << fabs(events_arr[i].get_vtx_massmj1()-events_arr[i].get_vtx_massmj2()) << std::endl;
	      }
	    

	    }
	  }

	verbose = false;
	if (verbose){
	  std::cout << "final_dz_dxy_chi2->Fill(" << final_chi2_dz_val << ","  << final_chi2_dxy_val << "," << final_chi2_dxy+final_chi2_dz <<  ");" << std::endl;

	  std::cout << "final_chi2_dz->Fill(" << final_chi2_dz  <<  ");" << std::endl;

	  std::cout << "final_dz_chi2sum->Fill(" << final_chi2_dz_val << ","  << chi2_dxy1+chi2_dxy2+chi2_dz <<  ");" << std::endl;
	}
	verbose = false;
	if (verbose){
	  std::cout << "************************" << std::endl;
	  std::cout << "final_chi2_dz: " << final_chi2_dz << std::endl;
	  std::cout << "final_chi2_dxy: " << final_chi2_dxy << std::endl;
	  std::cout << "final_chi2_dz1: " << final_chi2_dz1 << std::endl;
	  std::cout << "final_chi2_dz2: " << final_chi2_dz2 << std::endl;
	  std::cout << "final_chi2_dxy1: " << final_chi2_dxy1 << std::endl;
	  std::cout << "final_chi2_dxy2: " << final_chi2_dxy2 << std::endl;
	}
	verbose = false;

	if (verbose){
	  std::cout << std::endl << "dzdiff_ave: " << dzdiff_ave << std::endl;
	  std::cout << "dzdiff_var: " << dzdiff_var << std::endl;
	  std::cout << "sqrt(dzdiff_var): " << sqrt(dzdiff_var) << std::endl;
	  std::cout << "dxydiff_ave: " << dxydiff_ave << std::endl;
	  std::cout << "dxydiff_var: " << dxydiff_var << std::endl;
	  std::cout << "sqrt(dxydiff_var): " << sqrt(dxydiff_var) << std::endl;

	  std::cout << "dzdiff_abs_ave: " << dzdiff_abs_ave << std::endl;
	  std::cout << "dzdiff_abs_var: " << dzdiff_abs_var << std::endl;
	  std::cout << "sqrt(dzdiff_abs_var): " << sqrt(dzdiff_abs_var) << std::endl;
	  std::cout << "dxydiff_abs_ave: " << dxydiff_abs_ave << std::endl;
	  std::cout << "dxydiff_abs_var: " << dxydiff_abs_var << std::endl;
	  std::cout << "sqrt(dxydiff_abs_var): " << sqrt(dxydiff_abs_var) << std::endl;



	  std::cout << "dz_ave_mj1: " << dz_ave_mj1 << std::endl;
	  std::cout << "dz_ave_mj2: " << dz_ave_mj2 << std::endl;
	  std::cout << "dz_var_mj1: " << dz_var_mj1 << std::endl;
	  std::cout << "dz_var_mj2: " << dz_var_mj2 << std::endl;
	  std::cout << "sqrt(dz_var_mj1): " << sqrt(dz_var_mj1) << std::endl;
	  std::cout << "sqrt(dz_var_mj2): " << sqrt(dz_var_mj2) << std::endl << std::endl;

	  std::cout << "dxy_ave_mj1: " << dxy_ave_mj1 << std::endl;
	  std::cout << "dxy_ave_mj2: " << dxy_ave_mj2 << std::endl;
	  std::cout << "dxy_var_mj1: " << dxy_var_mj1 << std::endl;
	  std::cout << "dxy_var_mj2: " << dxy_var_mj2 << std::endl;
	  std::cout << "sqrt(dxy_var_mj1): " << sqrt(dxy_var_mj1) << std::endl;
	  std::cout << "sqrt(dxy_var_mj2): " << sqrt(dxy_var_mj2) << std::endl << std::endl;
	}

      }
      else std::cout << "muJet vertex not valid" << std::endl;
    }
  }

  verbose = false;

  if ( muJet1!= NULL && muJet2!= NULL ) {
    if (nMuJets == 2) {
      if (((*muJets)[0].vertexValid()) && ((*muJets)[1].vertexValid())){


	Global3DPoint muonjet1v = muJet1->vertexPoint();
	Global3DPoint muonjet2v = muJet2->vertexPoint();

	if (rec11<20 && rec22 > 20){


	  gen_match_x1diff_old = fabs(muonjet1v.x() - genMuonGroups[0][0]->vertex().x()); 
	  if (verbose) std::cout << "gen_match_x1diff_old: " << gen_match_x1diff_old << ", muonjet1v.x(): " << muonjet1v.x() << ", genMuonGroups[0][0]->vertex().x(): " << genMuonGroups[0][0]->vertex().x() << std::endl;
	  gen_match_y1diff_old = fabs(muonjet1v.y() - genMuonGroups[0][0]->vertex().y()); 
	  if (verbose) std::cout << "gen_match_y1diff_old: " << gen_match_y1diff_old << ", muonjet1v.y(): " << muonjet1v.y() << ", genMuonGroups[0][0]->vertex().y(): " << genMuonGroups[0][0]->vertex().y() << std::endl;
	  gen_match_z1diff_old = fabs(muonjet1v.z() - genMuonGroups[0][0]->vertex().z()); 
	  if (verbose) std::cout << "gen_match_z1diff_old: " << gen_match_z1diff_old << ", muonjet1v.z(): " << muonjet1v.z() << ", genMuonGroups[0][0]->vertex().z(): " << genMuonGroups[0][0]->vertex().z() << std::endl;
	  gen_match_x2diff_old = fabs(muonjet2v.x() - genMuonGroups[1][0]->vertex().x()); 
	  if (verbose) std::cout << "gen_match_x2diff_old: " << gen_match_x2diff_old << ", muonjet2v.x(): " << muonjet2v.x() << ", genMuonGroups[1][0]->vertex().x(): " << genMuonGroups[1][0]->vertex().x() << std::endl;
	  gen_match_y2diff_old = fabs(muonjet2v.y() - genMuonGroups[1][0]->vertex().y()); 
	  if (verbose) std::cout << "gen_match_y2diff_old: " << gen_match_y2diff_old << ", muonjet2v.y(): " << muonjet2v.y() << ", genMuonGroups[1][0]->vertex().y(): " << genMuonGroups[1][0]->vertex().y() << std::endl;
	  gen_match_z2diff_old = fabs(muonjet2v.z() - genMuonGroups[1][0]->vertex().z()); 
	  if (verbose) std::cout << "gen_match_z2diff_old: " << gen_match_z2diff_old << ", muonjet2v.z(): " << muonjet2v.z() << ", genMuonGroups[1][0]->vertex().z(): " << genMuonGroups[1][0]->vertex().z() << std::endl;
	}
	else if (rec11 > 20 && rec22 < 20){
	  gen_match_x1diff_old = fabs(muonjet1v.x() - genMuonGroups[1][0]->vertex().x()); 
	  if (verbose) std::cout << "gen_match_x1diff_old: " << gen_match_x1diff_old << ", muonjet1v.x(): " << muonjet1v.x() << ", genMuonGroups[1][0]->vertex().x(): " << genMuonGroups[1][0]->vertex().x() << std::endl;
	  gen_match_y1diff_old = fabs(muonjet1v.y() - genMuonGroups[1][0]->vertex().y()); 
	  if (verbose) std::cout << "gen_match_y1diff_old: " << gen_match_y1diff_old << ", muonjet1v.y(): " << muonjet1v.y() << ", genMuonGroups[1][0]->vertex().y(): " << genMuonGroups[1][0]->vertex().y() << std::endl;
	  gen_match_z1diff_old = fabs(muonjet1v.z() - genMuonGroups[1][0]->vertex().z()); 
	  if (verbose) std::cout << "gen_match_z1diff_old: " << gen_match_z1diff_old << ", muonjet1v.z(): " << muonjet1v.z() << ", genMuonGroups[1][0]->vertex().z(): " << genMuonGroups[1][0]->vertex().z() << std::endl;
	  gen_match_x2diff_old = fabs(muonjet2v.x() - genMuonGroups[0][0]->vertex().x()); 
	  if (verbose) std::cout << "gen_match_x2diff_old: " << gen_match_x2diff_old << ", muonjet2v.x(): " << muonjet2v.x() << ", genMuonGroups[0][0]->vertex().x(): " << genMuonGroups[0][0]->vertex().x() << std::endl;
	  gen_match_y2diff_old = fabs(muonjet2v.y() - genMuonGroups[0][0]->vertex().y()); 
	  if (verbose) std::cout << "gen_match_y2diff_old: " << gen_match_y2diff_old << ", muonjet2v.y(): " << muonjet2v.y() << ", genMuonGroups[0][0]->vertex().y(): " << genMuonGroups[0][0]->vertex().y() << std::endl;
	  gen_match_z2diff_old = fabs(muonjet2v.z() - genMuonGroups[0][0]->vertex().z()); 
	  if (verbose) std::cout << "gen_match_z2diff_old: " << gen_match_z2diff_old << ", muonjet2v.z(): " << muonjet2v.z() << ", genMuonGroups[0][0]->vertex().z(): " << genMuonGroups[0][0]->vertex().z() << std::endl;
	}
      }
    }
  }


  double u_isomj1 = isomj1;
  double u_isomj2 = isomj2;



  verbose = false;
  if (verbose){
    //   // if (min_dz_n > 0.5){
    // if (min_dz_n > 0.5 && min_dz_n < 0.7){
    if ( muJet1 != NULL && muJet2 != NULL ) {     
      if (((*muJets)[0].vertexValid()) && ((*muJets)[1].vertexValid())){
	// if (gen_match_z1diff_old < gen_match_z1diff && gen_match_x1diff_old < gen_match_x1diff && gen_match_y1diff_old < gen_match_y1diff){
	
	  std::cout << "************************************" << std::endl;

	  std::cout << "Event: " << event << std::endl;
	  std::cout << "fabs(dzmj1-dzmj2): " << temp_dz_diff << std::endl;	  
	  std::cout << "min_dz_n : " << min_dz_n << std::endl;
	  std::cout << "dzmj1: " << dzmj1 << ", dzmj2: " << dzmj2 << std::endl;
	  std::cout << "temp_dzmj1: " << temp_dzmj1 << ", temp_dzmj2: " << temp_dzmj2 << std::endl;
	  std::cout << "massmj1: " << massmj1 << ", massmj2: " << massmj2 << std::endl;

	  std::cout << "old xyz1: " << muJet1->vertexPoint().x() << ", " << muJet1->vertexPoint().y() << ", " << muJet1->vertexPoint().z() << std::endl;
	  std::cout << "new xyz1: " << newx1 << ", " << newy1 << ", " << newz1 << std::endl;
	  std::cout << "gen xyz1: " << genMuonGroups[0][0]->vertex().x() << ", " << genMuonGroups[0][0]->vertex().y() << ", " << genMuonGroups[0][0]->vertex().z() << std::endl ;
	  std::cout << "old xyz2: " << muJet2->vertexPoint().x() << ", " << muJet2->vertexPoint().y() << ", " << muJet2->vertexPoint().z() << std::endl;
	  std::cout << "new xyz2: " << newx2 << ", " << newy2 << ", " << newz2 << std::endl;
	  std::cout << "gen xyz2: " << genMuonGroups[1][0]->vertex().x() << ", " << genMuonGroups[1][0]->vertex().y() << ", " << genMuonGroups[1][0]->vertex().z() << std::endl;

	  std::cout <<" gen_match_x1diff = " << gen_match_x1diff << ", gen_match_x1diff_old: " << gen_match_x1diff_old << std::endl;
	  std::cout <<" gen_match_y1diff = " << gen_match_y1diff << ", gen_match_y1diff_old: " << gen_match_y1diff_old << std::endl;
	  std::cout <<" gen_match_z1diff = " << gen_match_z1diff << ", gen_match_z1diff_old: " << gen_match_z1diff_old << std::endl;
	  std::cout <<" gen_match_x2diff = " << gen_match_x2diff << ", gen_match_x2diff_old: " << gen_match_x2diff_old << std::endl;
	  std::cout <<" gen_match_y2diff = " << gen_match_y2diff << ", gen_match_y2diff_old: " << gen_match_y2diff_old << std::endl;
	  std::cout <<" gen_match_z2diff = " << gen_match_z2diff << ", gen_match_z2diff_old: " << gen_match_z2diff_old << std::endl;

	// }
      }
    }
  }
  verbose = false;

  edm::Handle<pat::TriggerEvent> triggerEvent;
  iEvent.getByLabel("patTriggerEvent", triggerEvent);
  
  //    std::cout<<"  trigger   "<<std::endl;
  bool isDiMuonHLTFired = false;
  if (    ( triggerEvent->path("HLT_Mu17_Mu8_v22") && triggerEvent->path("HLT_Mu17_Mu8_v22")->wasAccept() )
	  || ( triggerEvent->path("HLT_Mu17_Mu8_v21") && triggerEvent->path("HLT_Mu17_Mu8_v21")->wasAccept() )
	  || ( triggerEvent->path("HLT_Mu17_Mu8_v19") && triggerEvent->path("HLT_Mu17_Mu8_v19")->wasAccept() )
	  || ( triggerEvent->path("HLT_Mu17_Mu8_v18") && triggerEvent->path("HLT_Mu17_Mu8_v18")->wasAccept() )
	  || ( triggerEvent->path("HLT_Mu17_Mu8_v17") && triggerEvent->path("HLT_Mu17_Mu8_v17")->wasAccept() )
	  || ( triggerEvent->path("HLT_Mu17_Mu8_v16") && triggerEvent->path("HLT_Mu17_Mu8_v16")->wasAccept() ) ) {
    isDiMuonHLTFired = true;
  }  

  if(isDiMuonHLTFired) trigger=1;
  else trigger=0;


  
  // if ( muJet1 != NULL && muJet2 != NULL ) {

  //   if (((*muJets)[0].vertexValid()) && ((*muJets)[1].vertexValid())){

  //     const pat::MultiMuon *diMuonTmp = NULL;

  //     for ( unsigned int i = 1; i <= 2; i++ ) { 
  // 	double isoTmp = 0.0;
  // 	if ( i == 1 ){
  // 	  diMuonTmp = muJet1;
  // 	  if (verbose) std::cout << "i = 1, mj1" << std::endl;
  // 	}
  // 	if ( i == 2 ){
  // 	  diMuonTmp = muJet2;
  // 	  if (verbose) std::cout << "i = 2, mj2" << std::endl;
  // 	}

  // 	for (reco::TrackCollection::const_iterator track = tracks->begin(); track != tracks->end(); ++track) {
  // 	  bool trackIsMuon = false;
  // 	  if ( diMuonTmp->sameTrack( &*track, &*(diMuonTmp->muon(0)->innerTrack()) ) || diMuonTmp->sameTrack( &*track, &*(diMuonTmp->muon(1)->innerTrack()) ) ) trackIsMuon = true;
  // 	  if (!trackIsMuon) {
  // 	    double dPhi = My_dPhi( diMuonTmp->phi(), track->phi() );
  // 	    double dEta = diMuonTmp->eta() - track->eta();
  // 	    double dR = sqrt( dPhi*dPhi + dEta*dEta ); 
  // 	    if ( dR < 0.4 && track->pt() > 0.5 ) {
  // 	      double dz = fabs( track->dz(beamSpot->position()) - diMuonTmp->dz(beamSpot->position()) );
  // 	      if (verbose) std::cout << "n i: " << i << ", track->pt(): " << track->pt() << ", track->eta(): " << track->eta() << ", track->phi(): " << track->phi() << ", dPhi: " << dPhi << ", dEta: " << dEta << ", dR: " << dR << ", dz " << dz << std::endl;
	      
  // 	      double zpt =  pow(diMuonTmp->vertexMomentum().x()*diMuonTmp->vertexMomentum().x() + diMuonTmp->vertexMomentum().y()*diMuonTmp->vertexMomentum().y(),0.5);

  // 	      if (verbose) std::cout << "n i: " << i << ", diMuonTmp->pt(): " << zpt << ", diMuonTmp->eta(): " << diMuonTmp->eta() << ", diMuonTmp->phi(): " << diMuonTmp->phi() << ", dPhi: " << dPhi << ", dEta: " << dEta << ", dR: " << dR << ", dz " << dz << ", diMuonTmp->dz(beamSpot->position(): " << diMuonTmp->dz(beamSpot->position()) << std::endl;

  // 	      if ( dz < 0.1 ){
  // 		isoTmp += track->pt();
  // 		if (verbose) std::cout << "n isoTmp: " << isoTmp << std::endl;
  // 	      }
  // 	    }    
  // 	  }
  // 	}
  // 	if ( i == 1 ) isomj1 = isoTmp;
  // 	if ( i == 2 ) isomj2 = isoTmp;
  //     }
  //   }  
  // }

  verbose = false;
  if (verbose){
    if (u_isomj1 != isomj1 &&  u_isomj2 != isomj2){
      std::cout << "Event: " << event << ", lumi: " << lumi << ", run: " << run << std::endl;
      std::cout << "u_isomj1: " << u_isomj1 << ", isomj1: " << isomj1 << std::endl;
      std::cout << "u_isomj2: " << u_isomj2 << ", isomj2: " << isomj2 << std::endl;
    }
  }

  edm::Handle<reco::VertexCollection> primaryVertices;
  iEvent.getByLabel("offlinePrimaryVertices", primaryVertices);
  
  bool isVtxOk=false;
  for (reco::VertexCollection::const_iterator vertex = primaryVertices->begin();  vertex != primaryVertices->end();  ++vertex) {
    if (vertex->isValid() && !vertex->isFake() && vertex->tracksSize() > 3 && fabs(vertex->z()) < 24.) isVtxOk=true;
  }

  if(isVtxOk) isVtx=1;
  
  // if (dphi1 !=0 ) std::cout << "end dphi1: " << dphi1 << std::endl;
  // if (dphi2 !=0 ) std::cout << "end dphi2: " << dphi2 << std::endl;

  verbose = false;
  if (verbose){
    if ( muJet1!= NULL && muJet2!= NULL ) {
      if (nMuJets == 2) {
	if (((*muJets)[0].vertexValid()) && ((*muJets)[1].vertexValid())){
      
	  std::cout << "*************" << std::endl;
	  std::cout << "Event: " << event << std::endl;
	  std::cout << "dzmj1: " << dzmj1 << ", dzmj2: " << dzmj2 << std::endl;
	  std::cout << "fabs(dzmj1-dzmj2): " << fabs(dzmj1-dzmj2) << std::endl;
	  std::cout << "massmj1: " << massmj1 << ", massmj2: " << massmj2 << std::endl;
	  std::cout << "min_dz_n : " << min_dz_n << std::endl;
  
	}
      }
    }
  }



  AnaTree->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void 
Analysis::beginJob()
{
  std::cout << "BEGIN JOB" << std::endl;
  
  edm::Service<TFileService> tFileService;
  AnaTree = tFileService->make<TTree>("Events", "Events");


  // Event variables
  AnaTree->Branch("event",     &event,     "event/I");
  AnaTree->Branch("trigger",   &trigger,   "trigger/I");
  AnaTree->Branch("isVtx",     &isVtx,     "isVtx/I");
  AnaTree->Branch("genmuons",  &genmuons,  "genmuons/I");
  AnaTree->Branch("gengamD",   &gengamD,   "gengamD/I");
  AnaTree->Branch("genHiggs",  &genHiggs,  "genHiggs/I");
  AnaTree->Branch("recmuons",  &recmuons,  "recmuons/I");
  AnaTree->Branch("recmujets", &recmujets, "recmujets/I");

  // Generator Muons
  AnaTree->Branch("ptgengamD0",      &ptgengamD,      "ptgengamD[2]/F");
  AnaTree->Branch("ptgenHiggs",      &ptgenHiggs,     "ptgenHiggs[0]/F");

  AnaTree->Branch("ptgenMuons",     &ptgenMuons,    "ptgenMuons[genmuons]/F");
  AnaTree->Branch("etagenMuons",    &etagenMuons,   "etagenMuons[genmuons]/F");

  AnaTree->Branch("ptgenMuonsA0",    &ptgenMuonsA0,   "ptgenMuonsA0[2]/F");
  AnaTree->Branch("ptgenMuonsA1",    &ptgenMuonsA1,   "ptgenMuonsA1[2]/F");
  AnaTree->Branch("etagenMuonsA0",   &etagenMuonsA0,  "etagenMuonsA0[2]/F");
  AnaTree->Branch("etagenMuonsA1",   &etagenMuonsA1,  "etagenMuonsA1[2]/F");
  

  // Reconstructed Muons

  AnaTree->Branch("ptrecMuons",  &ptrecMuons,    "ptrecMuons[recmuons]/F");
  AnaTree->Branch("etarecMuons", &etarecMuons,   "etarecMuons[recmuons]/F");

  AnaTree->Branch("mj1muons",    &mj1muons,  "mj1muons/I");
  AnaTree->Branch("mj2muons",    &mj2muons,  "mj2muons/I");

  AnaTree->Branch("ptmj1muons",  &ptmj1muons, "ptmj1muons[mj1muons]/F");
  AnaTree->Branch("ptmj2muons",  &ptmj2muons, "ptmj2muons[mj2muons]/F");
  AnaTree->Branch("etamj1muons", &etamj1muons, "etamj1muons[mj1muons]/F");
  AnaTree->Branch("etamj2muons", &etamj2muons, "etamj2muons[mj2muons]/F");

  AnaTree->Branch("dzmj1",       &dzmj1,  "dzmj1/F");
  AnaTree->Branch("dzmj2",       &dzmj2,  "dzmj2/F");


  AnaTree->Branch("three_sig_1",       &three_sig_1,  "three_sig_1/I");
  AnaTree->Branch("three_sig_2",       &three_sig_2,  "three_sig_2/I");

  AnaTree->Branch("three_sig_x1",       &three_sig_x1,  "three_sig_x1/I");
  AnaTree->Branch("three_sig_y1",       &three_sig_y1,  "three_sig_y1/I");
  AnaTree->Branch("three_sig_z1",       &three_sig_z1,  "three_sig_z1/I");
  AnaTree->Branch("three_sig_x2",       &three_sig_x2,  "three_sig_x2/I");
  AnaTree->Branch("three_sig_y2",       &three_sig_y2,  "three_sig_y2/I");
  AnaTree->Branch("three_sig_z2",       &three_sig_z2,  "three_sig_z2/I");



  AnaTree->Branch("dxymj1",       &dxymj1,  "dxymj1/F");
  AnaTree->Branch("dxymj2",       &dxymj2,  "dxymj2/F");


  AnaTree->Branch("dzmj1_old",       &dzmj1_old,  "dzmj1_old/F");
  AnaTree->Branch("dzmj2_old",       &dzmj2_old,  "dzmj2_old/F");

  AnaTree->Branch("dphi1",       &dphi1,  "dphi1/F");
  AnaTree->Branch("dphi2",       &dphi2,  "dphi2/F");

  AnaTree->Branch("massmj1",     &massmj1, "massmj1/F");
  AnaTree->Branch("massmj2",     &massmj2, "massmj2/F");

  AnaTree->Branch("gen_match_dzmj1",       &gen_match_dzmj1,  "gen_match_dzmj1/F");
  AnaTree->Branch("gen_match_dzmj2",       &gen_match_dzmj2,  "gen_match_dzmj2/F");
  AnaTree->Branch("gen_match_massmj1",     &gen_match_massmj1, "gen_match_massmj1/F");
  AnaTree->Branch("gen_match_massmj2",     &gen_match_massmj2, "gen_match_massmj2/F");
  AnaTree->Branch("gen_match_chi2mj1",     &gen_match_chi2mj1, "gen_match_chi2mj1/F");
  AnaTree->Branch("gen_match_chi2mj2",     &gen_match_chi2mj2, "gen_match_chi2mj2/F");

  AnaTree->Branch("gen_match_x1diff",  &gen_match_x1diff, "gen_match_x1diff/F");
  AnaTree->Branch("gen_match_y1diff",  &gen_match_y1diff, "gen_match_y1diff/F");
  AnaTree->Branch("gen_match_z1diff",  &gen_match_z1diff, "gen_match_z1diff/F");
  AnaTree->Branch("gen_match_x2diff",  &gen_match_x2diff, "gen_match_x2diff/F");
  AnaTree->Branch("gen_match_y2diff",  &gen_match_y2diff, "gen_match_y2diff/F");
  AnaTree->Branch("gen_match_z2diff",  &gen_match_z2diff, "gen_match_z2diff/F");

  AnaTree->Branch("gen_match_x1diff_old",  &gen_match_x1diff_old, "gen_match_x1diff_old/F");
  AnaTree->Branch("gen_match_y1diff_old",  &gen_match_y1diff_old, "gen_match_y1diff_old/F");
  AnaTree->Branch("gen_match_z1diff_old",  &gen_match_z1diff_old, "gen_match_z1diff_old/F");
  AnaTree->Branch("gen_match_x2diff_old",  &gen_match_x2diff_old, "gen_match_x2diff_old/F");
  AnaTree->Branch("gen_match_y2diff_old",  &gen_match_y2diff_old, "gen_match_y2diff_old/F");
  AnaTree->Branch("gen_match_z2diff_old",  &gen_match_z2diff_old, "gen_match_z2diff_old/F");

  AnaTree->Branch("gen_mismatch_dzmj1",       &gen_mismatch_dzmj1,  "gen_mismatch_dzmj1/F");
  AnaTree->Branch("gen_mismatch_dzmj2",       &gen_mismatch_dzmj2,  "gen_mismatch_dzmj2/F");
  AnaTree->Branch("gen_mismatch_massmj1",     &gen_mismatch_massmj1, "gen_mismatch_massmj1/F");
  AnaTree->Branch("gen_mismatch_massmj2",     &gen_mismatch_massmj2, "gen_mismatch_massmj2/F");
  AnaTree->Branch("gen_mismatch_chi2mj1",     &gen_mismatch_chi2mj1, "gen_mismatch_chi2mj1/F");
  AnaTree->Branch("gen_mismatch_chi2mj2",     &gen_mismatch_chi2mj2, "gen_mismatch_chi2mj2/F");

  AnaTree->Branch("isomj1",     &isomj1, "isomj1/F");
  AnaTree->Branch("isomj2",     &isomj2, "isomj2/F");

}

// ------------ method called once each job just after ending the event loop  ------------
void 
Analysis::endJob() 
{
  
}

// ------------ method called when starting to processes a run  ------------
void 
Analysis::beginRun(edm::Run const & iRun, edm::EventSetup const & iSetup)
{
  
}

// ------------ method called when ending the processing of a run  ------------
void 
Analysis::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
Analysis::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
Analysis::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Analysis::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Analysis);
