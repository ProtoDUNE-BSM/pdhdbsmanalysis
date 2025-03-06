////////////////////////////////////////////////////////////////////////
// Class:       PandoraNeutrinoAnaProtoDUNE
// Plugin Type: analyzer (Unknown Unknown)
// File:        PandoraNeutrinoAnaProtoDUNE_module.cc
//
// Generated at Mon Apr  8 02:24:48 2024 by Ciaran Hasnip using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
//#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Additional Framework includes
#include "art_root_io/TFileService.h"

// ROOT includes
#include <TTree.h>
#include <TH1.h>

#include <string>

namespace ana {
  class PandoraNeutrinoAnaProtoDUNE;
}

//-----------------------------------------------
// Define analyser class
class ana::PandoraNeutrinoAnaProtoDUNE : public art::EDAnalyzer {
public:
  explicit PandoraNeutrinoAnaProtoDUNE(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PandoraNeutrinoAnaProtoDUNE(PandoraNeutrinoAnaProtoDUNE const&) = delete;
  PandoraNeutrinoAnaProtoDUNE(PandoraNeutrinoAnaProtoDUNE&&) = delete;
  PandoraNeutrinoAnaProtoDUNE& operator=(PandoraNeutrinoAnaProtoDUNE const&) = delete;
  PandoraNeutrinoAnaProtoDUNE& operator=(PandoraNeutrinoAnaProtoDUNE&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  void beginSubRun(art::SubRun const& subRun) override;
  void endSubRun(art::SubRun const& subRun) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  TTree *fSimulationNtuple;
  TTree *fRecoNtuple;
  TTree *fSubrunTree;

  // Declare member data here.
  unsigned int fEventID;
  unsigned int fRun;
  unsigned int fSubRun;

  art::InputTag fMCTruthLabel; ///< The name of the producer that tracked
                              ///< simulated particles through the detector
  art::InputTag fPFParticleLabel;
  double fSetPOT;
  
  // MCTruth information
  int fSimPDG;     ///< PDG ID of the particle being processed
  int fSimTrackID; ///< GEANT ID of the particle being processed
  int fCCNC; ///< Is neutrino interaction a CC or NC interaction
  unsigned int fTPCID; ///< TPC ID where neutrino interacts

  double fE;
  double fnuVertexX;
  double fnuVertexY;
  double fnuVertexZ;

  bool fInFV;

  TH1D *hMCNeutrinoEnergy;
  
  // Geomtry Information
  geo::GeometryCore const* fGeometryService; ///< pointer to Geometry provider
  std::vector<double> fFiducialBoundaries;

  // POT Summary data
  double fPOT;
  double fGoodPOT;

  // Pandora PFParticle information
  int fNPFParticles;
  int fNPrimaries;
  int fNPrimaryDaughters;
  TH1D *hTriggeredNeutrinoEnergy;
  TH1D *hRecoNeutrinoEnergy;
};


//-----------------------------------------------
// Analyser class constructor
ana::PandoraNeutrinoAnaProtoDUNE::PandoraNeutrinoAnaProtoDUNE(fhicl::ParameterSet const& p)
  : EDAnalyzer{p} 
  , fMCTruthLabel(p.get<std::string>("MCTruthLabel"))
  , fPFParticleLabel(p.get<std::string>("PFParticleLabel"))
  , fSetPOT(p.get<double>("SetPOT"))

  // More initializers here.
{
  
  // Get a pointer to the geometry service provider.
  fGeometryService = lar::providerFrom<geo::Geometry>();
  // TPC 1 is the first proper TPC - TPC 0 is for track stubs
  const geo::TPCGeo& tpc = fGeometryService->Cryostat().TPC(1);
  fFiducialBoundaries.push_back(0.); // central x
  fFiducialBoundaries.push_back(tpc.Width() - 0.05*tpc.Width()); // outer x
  fFiducialBoundaries.push_back(0.05*tpc.Height()); // bottom y
  fFiducialBoundaries.push_back(tpc.Height() - 0.05*tpc.Height()); // top y
  fFiducialBoundaries.push_back(0.05*(tpc.Length()*2.));
  fFiducialBoundaries.push_back((tpc.Length()*2.) - 0.05*(tpc.Length()*2));

  for (size_t i=0; i<fFiducialBoundaries.size(); i++) {
    std::cout << "\n bound = " << fFiducialBoundaries.at(i);
  }

  // Call appropriate consumes<>() for any products to be retrieved by this module.
  consumes<std::vector<simb::MCTruth>>(fMCTruthLabel);
}

//-----------------------------------------------
void ana::PandoraNeutrinoAnaProtoDUNE::analyze(art::Event const& e)
{
  // Implementation of required member function here.
  fEventID = e.id().event();
  fRun = e.run();
  fSubRun = e.subRun();

  fNPFParticles = 0;
  fNPrimaries = 0;
  fNPrimaryDaughters = 0;

  fInFV = false;

  // Define a "handle" to point to a vector of the objects.
  auto truthHandle = e.getValidHandle<std::vector<simb::MCTruth>>(fMCTruthLabel);

  for (auto const& truth : (*truthHandle)) {
    if (truth.NeutrinoSet()) {
      const auto &nu = truth.GetNeutrino();
      const auto &neutrino = nu.Nu();

      fSimPDG = neutrino.PdgCode();

      fCCNC = nu.CCNC();

      fE = neutrino.E();


      double fPrimaryStart[4];
      double fPrimaryVertex[4];

      const size_t numberTrajectoryPoints = neutrino.NumberTrajectoryPoints();
      const int last = numberTrajectoryPoints - 1;
      const TLorentzVector& positionStart = neutrino.Position(0);
      const TLorentzVector& positionEnd = neutrino.Position(last);
      // Set the vertex position - it should be the same value for each event	
      positionStart.GetXYZT(fPrimaryStart);
      positionEnd.GetXYZT(fPrimaryVertex);

      fnuVertexX = fPrimaryVertex[0];
      fnuVertexY = fPrimaryVertex[1];
      fnuVertexZ = fPrimaryVertex[2];

      if (std::fabs(fnuVertexX) < fFiducialBoundaries.at(1)) {
        if (fnuVertexY > fFiducialBoundaries.at(2) && 
            fnuVertexY < fFiducialBoundaries.at(3)) {
          if (fnuVertexZ > fFiducialBoundaries.at(4) && 
              fnuVertexZ < fFiducialBoundaries.at(5)) {
            fInFV = true;
          }
        }
      }

      geo::Point_t nuV_point(fnuVertexX, fnuVertexY, fnuVertexZ);
      fTPCID = fGeometryService->FindTPCAtPosition(nuV_point).TPC;
      if (fTPCID > 7 || fTPCID < 0) fTPCID = -999;
      
      // Store total event outputs in the TTree
      fSimulationNtuple->Fill();
    }
  }
  hMCNeutrinoEnergy->Fill(fE, fSetPOT/fPOT);

  art::Handle<std::vector<recob::PFParticle>> pfpHandle;
  //std::vector<art::Ptr<recob::PFParticle>> pfpVec;
  if (!e.getByLabel(fPFParticleLabel, pfpHandle)) { 
    return;
  }
  auto const &pfpVec = *pfpHandle;
  /*if (e.getByLabel(fPFParticleLabel, pfpHandle)) {
    art::fill_ptr_vector(pfpVec, pfpHandle);
  } else {
    return;
  }*/
  // Made it to here means that the event at least triggered so that Pandora was run
  hTriggeredNeutrinoEnergy->Fill(fE, fSetPOT/fPOT);

  // Check if there are any PFParticles
  if (pfpVec.empty()) {
    std::cout << "pfpVec.empty - should I be here?" << std::endl;
    return;
  }
  // Made it to here means that the event at least triggered so that Pandora was run
  //hTriggeredNeutrinoEnergy->Fill(fE, );

  // Initialise neutrino ID
  size_t neutrinoID(std::numeric_limits<size_t>::max());

  // Loop over PFParticles and look for a reconstructed neutrino
  //for (const art::Ptr<recob::PFParticle>& pfp : pfpVec) {
  for (const auto& pfp : pfpVec) {
    fNPFParticles++;

    //if (!pfp->IsPrimary()) continue;
    if (!pfp.IsPrimary()) continue;

    //neutrinoID = pfp->Self();
    //fNPrimaryDaughters = pfp->NumDaighters();
    neutrinoID = pfp.Self();
    fNPrimaryDaughters = pfp.NumDaughters();
    fNPrimaries++;
  }

  // Check that a reconstructed neutrino was found
  if (neutrinoID == std::numeric_limits<size_t>::max()) return;

  // Made to here so neutrino was reconstructed;
  hRecoNeutrinoEnergy->Fill(fE, fSetPOT/fPOT);

  fRecoNtuple->Fill();
}

//-----------------------------------------------
// Define outputs at start of the job
void ana::PandoraNeutrinoAnaProtoDUNE::beginJob() {
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;

  hMCNeutrinoEnergy = tfs->make<TH1D>("Total_MC_Nu_Energy", ";Energy (GeV);", 20, 0, 200);
  hTriggeredNeutrinoEnergy = tfs->make<TH1D>("Triggered_Nu_Energy", ";Energy (GeV);", 20, 0, 200);
  hRecoNeutrinoEnergy = tfs->make<TH1D>("Reco_Nu_Energy", ";Energy (GeV);", 20, 0, 200);

  // Get TFileService to create an output tree
  fSimulationNtuple = tfs->make<TTree>("GenieTruth", "GENIE Output Tree");

  // Add branches to TTree
  fSimulationNtuple->Branch("eventID", &fEventID);
  fSimulationNtuple->Branch("SubRun", &fSubRun, "SubRun/I");
  fSimulationNtuple->Branch("Run", &fRun, "Run/I");
  fSimulationNtuple->Branch("PDG", &fSimPDG, "PDG/I");
  fSimulationNtuple->Branch("CCNC", &fCCNC, "CCNC/I");
  fSimulationNtuple->Branch("TPCID", &fTPCID);

  fSimulationNtuple->Branch("E", &fE, "E/D");
  fSimulationNtuple->Branch("nuVertexX", &fnuVertexX, "nuVertexX/D");
  fSimulationNtuple->Branch("nuVertexY", &fnuVertexY, "nuVertexY/D");
  fSimulationNtuple->Branch("nuVertexZ", &fnuVertexZ, "nuVertexZ/D");
  fSimulationNtuple->Branch("InFV", &fInFV, "InFV/B");
  
  fRecoNtuple = tfs->make<TTree>("Reco", "Pandora Reco Output Tree");
  fRecoNtuple->Branch("eventID", &fEventID);
  fRecoNtuple->Branch("SubRun", &fSubRun, "SubRun/I");
  fRecoNtuple->Branch("Run", &fRun, "Run/I");
  fRecoNtuple->Branch("nPFParticles", &fNPFParticles);
  fRecoNtuple->Branch("nPrimaries", &fNPrimaries);
  fRecoNtuple->Branch("nPrimaryDaughters", &fNPrimaryDaughters);

  fSubrunTree = tfs->make<TTree>("SubRunTree", "SubRun-level information");
  fSubrunTree->Branch("POT", &fPOT, "POT/D");
  fSubrunTree->Branch("GoodPOT", &fGoodPOT, "GoodPOT/D");

}

//-----------------------------------------------
void ana::PandoraNeutrinoAnaProtoDUNE::beginSubRun(art::SubRun const& subRun) {
    const auto potSummaryHandle = subRun.getValidHandle<sumdata::POTSummary>("generator");
    const auto &potSummary = *potSummaryHandle;
    fPOT = potSummary.totpot;
    fGoodPOT = potSummary.totgoodpot;
     
    std::cout << "POTSummary content: totpot = " << potSummary.totpot 
      << ", totgoodpot = " << potSummary.totgoodpot << std::endl;
      
    // Fill the TTree with the current subrun's POT information
    fSubrunTree->Fill();
}

//-----------------------------------------------
void ana::PandoraNeutrinoAnaProtoDUNE::endSubRun(art::SubRun const& subRun) {}

//-----------------------------------------------
void ana::PandoraNeutrinoAnaProtoDUNE::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(ana::PandoraNeutrinoAnaProtoDUNE)
