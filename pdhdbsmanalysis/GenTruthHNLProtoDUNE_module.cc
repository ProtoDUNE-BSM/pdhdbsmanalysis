////////////////////////////////////////////////////////////////////////
// Class:       GenTruthHNLProtoDUNE
// Plugin Type: analyzer (Unknown Unknown)
// File:        GenTruthHNLProtoDUNE_module.cc
// Author:      Hamza Amar
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "larcoreobj/SummaryData/POTSummary.h"

#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larsim/Utils/TruthMatchUtils.h"

// additional Framework includes
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

// ROOT includes
#include "TTree.h"

namespace ana {
  class GenTruthHNLProtoDUNE;
}

/**
 * @brief GenTruthHNLProtoDUNE class is an EDAnalyzer that performs a validation of the true HNL decay generation
 * This module is designed to provide the true HNL decay information.
 * The module produces two trees: fTreeHNLTrue and fTreePOT.
 * The fTreeHNLTrue tree contains the true HNL decay information, including the energy, vertex position, momentum and direction.
 * The fTreePOT tree contains the POT summary information.
 */
class ana::GenTruthHNLProtoDUNE : public art::EDAnalyzer {
public:
  explicit GenTruthHNLProtoDUNE(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  GenTruthHNLProtoDUNE(GenTruthHNLProtoDUNE const&) = delete;
  GenTruthHNLProtoDUNE(GenTruthHNLProtoDUNE&&) = delete;
  GenTruthHNLProtoDUNE& operator=(GenTruthHNLProtoDUNE const&) = delete;
  GenTruthHNLProtoDUNE& operator=(GenTruthHNLProtoDUNE&&) = delete;

  // Custom functions
  void ProcessTruthInformation(const std::vector<simb::MCTruth>& truthCollection, TTree* tree);
  void CreateBranches(TTree* tree, unsigned int& eventID, int& run, int& subRun, int& pdg, int& mother,
                      double& T0,
                      double& vertexPosX, double& vertexPosY, double& vertexPosZ,
                      double& energy, double& kineticEnergy,
                      double& px, double& py, double& pz, double& p,
                      double& dirCosX, double& dirCosY, double& dirCosZ);

  // Required functions.
  void analyze(art::Event const& e) override;

  void beginSubRun(art::SubRun const& subRun) override;
  void endSubRun(art::SubRun const& subRun) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:
  // Declare member data here
  std::string fMCTruthLabel;
  std::string fCosmicLabel;
  std::string fAr39Label, fAr42Label, fK42Label, fKr85Label;

  // POT Summary data
  const bool fRollUpUnsavedIDs = true;
  double fPOT = 0;
  double fTotalPOT = 0;
  double fGoodPOT = 0;

  // General event information
  int fRun;
  int fSubRun;
  unsigned int fEventID;

  // Trees
  TTree* fTreeHNLTrue; // Signal
  TTree* fTreeCosmicsTrue; // Cosmics
  TTree *fTreeAr39True, *fTreeAr42True, *fTreeK42True, *fTreeKr85True; // Radiologicals
  TTree* fTreePOT;

  // Variables for fTreeHNLTrue
  double fT0; // Time of the event in ns
  double fEnergy, fKineticEnergy; // HNL product energy in GeV
  int fPDG; // HNL product PDG Code
  int fmother; // To check if the HNL decay product is a primary particle, i.e. mother = -1
  double fVertexPosX, fVertexPosY, fVertexPosZ; // HNL product vertex position in cm
  double fPX, fPY, fPZ; // HNL product module components at the vertex in GeV
  double fP; // Modulus of HNL product momentum at the vertex in GeV
  double fDirCosX, fDirCosY, fDirCosZ; // HNL product direction cosines  

};

/**
 * @brief Construct the GenTruthHNLProtoDUNE module.
 * @param p fhicl::ParameterSet.
 * The constructor of the GenTruthHNLProtoDUNE module initializes the module by getting the configuration parameters from the fhicl file.
 */
ana::GenTruthHNLProtoDUNE::GenTruthHNLProtoDUNE(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}, 
  // More initializers here.
  fMCTruthLabel(p.get<std::string>("MCTruthLabel")),
  fCosmicLabel(p.get<std::string>("CosmicLabel")),
  fAr39Label(p.get<std::string>("Ar39Label")),
  fAr42Label(p.get<std::string>("Ar42Label")),
  fK42Label(p.get<std::string>("K42Label")),
  fKr85Label(p.get<std::string>("Kr85Label"))
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

/**
 * @brief Analyze the event.
 * @param e art::Event.
 * The analyze function is the main function of the GenTruthHNLProtoDUNE module.
 * The function is called for each event.
 * The function gets the HNL decay truth information.
 * The function then fills the fTreeHNLTrue and fTreePOT.
 */
void ana::GenTruthHNLProtoDUNE::analyze(art::Event const& e)
{
  // Set all general event information
  fRun     = e.run();
  fSubRun  = e.subRun();
  fEventID = e.id().event();

  // Get the HNL decay truth information
  auto truthHandle = e.getValidHandle<std::vector<simb::MCTruth>>(fMCTruthLabel);
  ana::GenTruthHNLProtoDUNE::ProcessTruthInformation(*truthHandle, fTreeHNLTrue);

  // Get the cosmics truth information
  auto cosmicsHandle = e.getValidHandle<std::vector<simb::MCTruth>>(fCosmicLabel);
  ana::GenTruthHNLProtoDUNE::ProcessTruthInformation(*cosmicsHandle, fTreeCosmicsTrue);

  // Get the radiologicals truth information
  // Ar39
  auto Ar39Handle = e.getValidHandle<std::vector<simb::MCTruth>>(fAr39Label);
  ana::GenTruthHNLProtoDUNE::ProcessTruthInformation(*Ar39Handle, fTreeAr39True);
  // Ar42
  auto Ar42Handle = e.getValidHandle<std::vector<simb::MCTruth>>(fAr42Label);
  ana::GenTruthHNLProtoDUNE::ProcessTruthInformation(*Ar42Handle, fTreeAr42True);
  // K42
  auto K42Handle = e.getValidHandle<std::vector<simb::MCTruth>>(fK42Label);
  ana::GenTruthHNLProtoDUNE::ProcessTruthInformation(*K42Handle, fTreeK42True);
  // Kr85
  auto Kr85Handle = e.getValidHandle<std::vector<simb::MCTruth>>(fKr85Label);
  ana::GenTruthHNLProtoDUNE::ProcessTruthInformation(*Kr85Handle, fTreeKr85True);

}

void ana::GenTruthHNLProtoDUNE::ProcessTruthInformation(const std::vector<simb::MCTruth>& truthCollection, TTree* tree) {
  for (auto const& truth : truthCollection) {
    for (int i = 0; i < truth.NParticles(); ++i) {
      const auto& particle = truth.GetParticle(i);
      fmother = particle.Mother();
      fPDG = particle.PdgCode();
      fT0 = particle.T();
      fEnergy = particle.E();
      fKineticEnergy = fEnergy - particle.Mass();
      fPX = particle.Px();
      fPY = particle.Py();
      fPZ = particle.Pz();
      fP = particle.P();
      fDirCosX = fPX / fP;
      fDirCosY = fPY / fP;
      fDirCosZ = fPZ / fP;
      fVertexPosX = particle.Vx();
      fVertexPosY = particle.Vy();
      fVertexPosZ = particle.Vz();
      // Fill the TTree
      tree->Fill();
    }
  }
}

void ana::GenTruthHNLProtoDUNE::CreateBranches(TTree* tree, unsigned int& eventID, int& run, int& subRun, int& pdg, int& mother,
                    double& T0,
                    double& vertexPosX, double& vertexPosY, double& vertexPosZ,
                    double& energy, double& kineticEnergy,
                    double& px, double& py, double& pz, double& p,
                    double& dirCosX, double& dirCosY, double& dirCosZ) {
  // General event information
  tree->Branch("Event", &eventID, "Event/I");
  tree->Branch("Run", &run, "Run/I");
  tree->Branch("SubRun", &subRun, "SubRun/I");
  // Particle true information
  tree->Branch("PDG", &pdg, "PDG/I");
  tree->Branch("Mother", &mother, "Mother/I");
  tree->Branch("T0", &T0, "T0/D"); // ns
  tree->Branch("VertexPositionX", &vertexPosX, "VertexPosX/D"); // cm
  tree->Branch("VertexPositionY", &vertexPosY, "VertexPosY/D"); // cm
  tree->Branch("VertexPositionZ", &vertexPosZ, "VertexPosZ/D"); // cm
  tree->Branch("Energy", &energy, "Energy/D"); // GeV
  tree->Branch("KineticEnergy", &kineticEnergy, "KineticEnergy/D"); // GeV
  tree->Branch("MomentumX", &px, "MomentumX/D"); // GeV
  tree->Branch("MomentumY", &py, "MomentumY/D"); // GeV
  tree->Branch("MomentumZ", &pz, "MomentumZ/D"); // GeV
  tree->Branch("Momentum", &p, "Momentum/D"); // GeV
  tree->Branch("DirectionX", &dirCosX, "DirCosX/D"); // Director cosines
  tree->Branch("DirectionY", &dirCosY, "DirCosY/D");
  tree->Branch("DirectionZ", &dirCosZ, "DirCosZ/D");
}

/**
 * @brief Begin the job.
 * The beginJob function initializes the output TTree for the true HNL decay information, 
 * and the POT information.
 */
void ana::GenTruthHNLProtoDUNE::beginJob() {
  // Make our handle to the TFileService
  art::ServiceHandle<art::TFileService> tfs;

  // Create TTree and branches for TTreeTrue
  fTreeHNLTrue = tfs->make<TTree>("HNLTree", "Tree containing HNL decay event information");
  ana::GenTruthHNLProtoDUNE::CreateBranches(fTreeHNLTrue, fEventID, fRun, fSubRun, fPDG, fmother,
                 fT0,
                 fVertexPosX, fVertexPosY, fVertexPosZ,
                 fEnergy, fKineticEnergy,
                 fPX, fPY, fPZ, fP,
                 fDirCosX, fDirCosY, fDirCosZ);

  // Create TTree and branches for TTreeCosmicsTrue
  fTreeCosmicsTrue = tfs->make<TTree>("CosmicsTree", "Tree containing cosmics event information");
  ana::GenTruthHNLProtoDUNE::CreateBranches(fTreeCosmicsTrue, fEventID, fRun, fSubRun, fPDG, fmother,
                 fT0,
                 fVertexPosX, fVertexPosY, fVertexPosZ,
                 fEnergy, fKineticEnergy,
                 fPX, fPY, fPZ, fP,
                 fDirCosX, fDirCosY, fDirCosZ);

  // Create TTree and branches for TTreeAr39True
  fTreeAr39True = tfs->make<TTree>("Ar39Tree", "Tree containing Ar39 event information");
  ana::GenTruthHNLProtoDUNE::CreateBranches(fTreeAr39True, fEventID, fRun, fSubRun, fPDG, fmother,
                 fT0,
                 fVertexPosX, fVertexPosY, fVertexPosZ,
                 fEnergy, fKineticEnergy,
                 fPX, fPY, fPZ, fP,
                 fDirCosX, fDirCosY, fDirCosZ);

  // Create TTree and branches for TTreeAr42True
  fTreeAr42True = tfs->make<TTree>("Ar42Tree", "Tree containing Ar42 event information");
  ana::GenTruthHNLProtoDUNE::CreateBranches(fTreeAr42True, fEventID, fRun, fSubRun, fPDG, fmother,
                 fT0,
                 fVertexPosX, fVertexPosY, fVertexPosZ,
                 fEnergy, fKineticEnergy,
                 fPX, fPY, fPZ, fP,
                 fDirCosX, fDirCosY, fDirCosZ);
  
  // Create TTree and branches for TTreeK42True
  fTreeK42True = tfs->make<TTree>("K42Tree", "Tree containing K42 event information");
  ana::GenTruthHNLProtoDUNE::CreateBranches(fTreeK42True, fEventID, fRun, fSubRun, fPDG, fmother,
                 fT0,
                 fVertexPosX, fVertexPosY, fVertexPosZ,
                 fEnergy, fKineticEnergy,
                 fPX, fPY, fPZ, fP,
                 fDirCosX, fDirCosY, fDirCosZ);
  
  // Create TTree and branches for TTreeKr85True
  fTreeKr85True = tfs->make<TTree>("Kr85Tree", "Tree containing Kr85 event information");
  ana::GenTruthHNLProtoDUNE::CreateBranches(fTreeKr85True, fEventID, fRun, fSubRun, fPDG, fmother,
                 fT0,
                 fVertexPosX, fVertexPosY, fVertexPosZ,
                 fEnergy, fKineticEnergy,
                 fPX, fPY, fPZ, fP,
                 fDirCosX, fDirCosY, fDirCosZ);

  // Create TTree and branches for TTreePOT
  fTreePOT = tfs->make<TTree>("POTTree", "Tree containing POT information");
  // PoT information
  fTreePOT -> Branch("POT", &fPOT, "POT/D");
  fTreePOT -> Branch("GoodPOT", &fGoodPOT, "GoodPOT/D");

}

/**
 * @brief Begin the subrun.
 * @param subRun art::SubRun.
 * The beginSubRun function is called at the beginning of each subrun.
 * The function gets the POT information from the POTSummary and fills the fTreePOT tree.
 */
void ana::GenTruthHNLProtoDUNE::beginSubRun(art::SubRun const& subRun)
{
  const auto potSummaryHandle = subRun.getValidHandle<sumdata::POTSummary>("generator");
  const auto &potSummary = *potSummaryHandle;
  // Get the POT information
  fPOT = potSummary.totpot;
  fGoodPOT = potSummary.totgoodpot;
  // Fill the TTree fTreePOT
  fTreePOT->Fill();
}

/**
 * @brief End the subrun.
 * @param subRun art::SubRun.
 * The endSubRun function is called at the end of each subrun.
 */
void ana::GenTruthHNLProtoDUNE::endSubRun(art::SubRun const& subRun) {}

/**
 * @brief End the job.
 * The endJob function is called at the end of the job.
 */
void ana::GenTruthHNLProtoDUNE::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(ana::GenTruthHNLProtoDUNE)
