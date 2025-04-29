////////////////////////////////////////////////////////////////////////
// Class:       NeutrinoValidation
// Plugin Type: analyzer (Unknown Unknown)
// File:        NeutrinoValidation_module.cc
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

#include "detdataformats/trigger/TriggerObjectOverlay.hpp"
#include "detdataformats/trigger/TriggerPrimitive.hpp"
#include "detdataformats/trigger/TriggerActivityData.hpp"
#include "detdataformats/trigger/TriggerCandidateData.hpp"

#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "dunereco/AnaUtils/DUNEAnaPFParticleUtils.h"
#include "dunereco/AnaUtils/DUNEAnaEventUtils.h"
#include "dunereco/AnaUtils/DUNEAnaShowerUtils.h"
#include "dunereco/AnaUtils/DUNEAnaTrackUtils.h"
#include "dunereco/AnaUtils/DUNEAnaSliceUtils.h"
#include "dunereco/AnaUtils/DUNEAnaUtilsBase.h"
#include "dunereco/FDSensOpt/NeutrinoAngularRecoAlg/NeutrinoAngularRecoAlg.h"
#include "dunereco/FDSensOpt/NeutrinoEnergyRecoAlg/NeutrinoEnergyRecoAlg.h"
#include "larsim/Utils/TruthMatchUtils.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

// additional Framework includes
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

// ROOT includes
#include "TTree.h"

#include <fstream>

namespace ana {
  class NeutrinoValidation;
}

/**
 * @brief NeutrinoValidation class is an EDAnalyzer that performs a validation of the neutrino reconstruction
 * This module is designed to validate the neutrino reconstruction by comparing the reconstructed neutrino energy, position and direction with other neutrino reconstruction algorithms or the true neutrino information.
 * The module uses the NeutrinoAngularRecoAlg and NeutrinoEnergyRecoAlg to reconstruct the neutrino direction and energy, respectively.
 * The module also provides a detailed information of the neutrino daughters, including the track and shower information.
 * The module produces tree trees: fTreeTrue, fTreePOT and fTreeReco.
 * The fTreeTrue tree contains the true neutrino information, including the neutrino energy, vertex position, momentum and direction.
 * The fTreePOT tree contains the POT summary information.
 * The fTreeReco tree contains the reconstructed neutrino information, including the neutrino energy, vertex position, direction and the daughter information.
 */
class ana::NeutrinoValidation : public art::EDAnalyzer {
public:
  explicit NeutrinoValidation(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  NeutrinoValidation(NeutrinoValidation const&) = delete;
  NeutrinoValidation(NeutrinoValidation&&) = delete;
  NeutrinoValidation& operator=(NeutrinoValidation const&) = delete;
  NeutrinoValidation& operator=(NeutrinoValidation&&) = delete;

  double GetSliceCaloEnergy(const art::Ptr<recob::Slice>& slicePtr, art::Event const& e);
  unsigned int GetNuDaughterInfo(
    const art::Ptr<recob::PFParticle>& pfparticlePtr, art::Event const& e,
    std::vector<unsigned int>& nuDaughterIsTrack, std::vector<unsigned int>& nuDaughterIsShower,
    std::vector<double>& TrackDirectionX, std::vector<double>& TrackDirectionY, std::vector<double>& TrackDirectionZ,
    std::vector<double>& ShowerDirectionX, std::vector<double>& ShowerDirectionY, std::vector<double>& ShowerDirectionZ, 
    std::vector<double>& KineticEnergyTrack, std::vector<double>& ShowerEnergy, std::vector<int>& DaughterPDG
  );

  // Required functions.
  void analyze(art::Event const& e) override;

  void beginSubRun(art::SubRun const& subRun) override;
  void endSubRun(art::SubRun const& subRun) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here
  std::string fSliceLabel;
  std::string fPFParticleLabel;
  std::string fVertexLabel;
  std::string fTrackLabel;
  std::string fShowerLabel;
  std::string fHitLabel;
  std::string fMCParticleLabel;
  std::string fMCTruthLabel;
  std::string fCalorimetryLabel;
  std::string fWireLabel;
  std::string fTALabel;
  bool fSliceCaloEnergy;

  // Trigger activity flag 
  int fTA;

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
  TTree* fTreeTrue;
  TTree* fTreePOT;
  TTree* fTreeReco;

  // Variables for fTreeTrue
  double fnuEnergy; // Neutrino energy in GeV
  int fnuPDG; // Neutrino PDG Code
  int fmother; // To check if the neutrino is a primary particle, i.e. mother = -1
  double fnuVertexPosX, fnuVertexPosY, fnuVertexPosZ; // Neutrino vertex position in cm
  double fnuPX, fnuPY, fnuPZ; // Neutrino module components at the vertex in GeV
  double fnuP; // Modulus of neutrino momentum at the vertex in GeV
  double fnuDirCosX, fnuDirCosY, fnuDirCosZ; // Neutrino direction cosines

  dune::NeutrinoAngularRecoAlg fNeutrinoRecoAngle;
  dune::NeutrinoEnergyRecoAlg fNeutrinoRecoEnergy;

  // Variables for fTreeReco
  // Jagged arrays for thorough primary neutrino event information
  // If only one primary neutrino per event, e.g. within the most energetic slice, 
  // then the vectors will have only one element (could be a scalar)
  std::vector<int> fRecoPDG, fRecoDaughterPDG; // PDG code
  std::vector<long unsigned int> fnuHits; // Number of hits
  std::vector<double> fRecoEnergy; // Energy in GeV
  std::vector<double> fRecoVertexPosX, fRecoVertexPosY, fRecoVertexPosZ; // Vertex position in cm
  std::vector<double> fRecoDirCosX, fRecoDirCosY, fRecoDirCosZ; // Direction cosines
  double fSliceEnergy; // Slice energy in GeV

  // Jagged arrays for thorough neutrino daughters event information
  std::vector<unsigned int> fDaughterIsTrack, fDaughterIsShower;
  std::vector<double> fTrackDirectionX, fTrackDirectionY, fTrackDirectionZ; // Track direction cosines
  std::vector<double> fShowerDirectionX, fShowerDirectionY, fShowerDirectionZ; // Shower direction cosines
  std::vector<double> fKineticEnergyTrack, fShowerEnergy; // Energy in MeV

  // Pandora PFParticle information
  // Number of PFParticles, neutrinos, daughters
  unsigned int fNoPFParticles = 0;
  unsigned int fNoNeutrinos = 0;
  unsigned int fNoDaughters = 0; 
  // Metadata
  std::vector<int> fnuSliceKey, fSliceIndex, fnuID, fnuScore;
};

/**
 * @brief Construct the NeutrinoValidation module.
 * @param p fhicl::ParameterSet.
 * The constructor of the NeutrinoValidation module initializes the module by getting the configuration parameters from the fhicl file.
 * The module also initializes the NeutrinoAngularRecoAlg and NeutrinoEnergyRecoAlg algorithms.
 */
ana::NeutrinoValidation::NeutrinoValidation(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}, 
  // More initializers here.
  fSliceLabel(p.get<std::string>("SliceLabel")),
  fPFParticleLabel(p.get<std::string>("PFParticleLabel")),
  fVertexLabel(p.get<std::string>("VertexLabel")),
  fTrackLabel(p.get<std::string>("TrackLabel")),
  fShowerLabel(p.get<std::string>("ShowerLabel")),
  fHitLabel(p.get<std::string>("HitLabel")),
  fMCParticleLabel(p.get<std::string>("MCParticleLabel")),
  fMCTruthLabel(p.get<std::string>("MCTruthLabel")),
  fCalorimetryLabel(p.get<std::string>("CalorimetryLabel")),
  fWireLabel(p.get<std::string>("WireLabel")),
  fTALabel(p.get<std::string>("TALabel")),
  fSliceCaloEnergy(p.get<bool>("SliceCaloEnergy")),
  fNeutrinoRecoAngle(p, fTrackLabel, fShowerLabel, fHitLabel,
      fWireLabel, fTrackLabel, fShowerLabel, fHitLabel),
  fNeutrinoRecoEnergy(p, fTrackLabel, fShowerLabel, fHitLabel,
      fWireLabel, fTrackLabel, fShowerLabel, fHitLabel)
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

/**
 * @brief Analyze the event.
 * @param e art::Event.
 * The analyze function is the main function of the NeutrinoValidation module.
 * The function is called for each event and performs the validation of the neutrino reconstruction.
 * The function gets the trigger activity information flag, the neutrino truth information and the slices and PFParticles from Pandora reconstructed data.
 * The function then fills the fTreeTrue, fTreePOT and fTreeReco trees.
 */
void ana::NeutrinoValidation::analyze(art::Event const& e)
{
  // Set all general event information
  fRun     = e.run();
  fSubRun  = e.subRun();
  fEventID = e.id().event();

  // Get the trigger information
  art::Handle<std::vector<dunedaq::trgdataformats::TriggerActivityData>> taHandle;
  fTA = (!e.getByLabel(fTALabel, taHandle) || taHandle->empty()) ? 0 : 1;

  // Get the neutrino truth information
  auto truthHandle = e.getValidHandle<std::vector<simb::MCTruth>>(fMCTruthLabel);
  for (auto const& truth : (*truthHandle)) {
    if (truth.NeutrinoSet()) {
      const auto &nu = truth.GetNeutrino();
      const auto &neutrino = nu.Nu();
      fnuEnergy = neutrino.E();
      fnuPDG = neutrino.PdgCode();
      fmother = neutrino.Mother();
      fnuVertexPosX = neutrino.Vx();
      fnuVertexPosY = neutrino.Vy();
      fnuVertexPosZ = neutrino.Vz();
      fnuPX = neutrino.Px();
      fnuPY = neutrino.Py();
      fnuPZ = neutrino.Pz();
      fnuP = neutrino.P();
      fnuDirCosX = fnuPX/fnuP;
      fnuDirCosY = fnuPY/fnuP;
      fnuDirCosZ = fnuPZ/fnuP;    
      // Fill the TTree fTreeTrue
      fTreeTrue->Fill();
    }
  }

  // Handle slices and PFParticles
  art::Handle<std::vector<recob::Slice>> sliceHandle = e.getHandle<std::vector<recob::Slice>>(fSliceLabel);
  if (sliceHandle.isValid()) {
    std::vector<art::Ptr<recob::Slice>> slicePtrVector;
    art::fill_ptr_vector(slicePtrVector, sliceHandle);
    // Get the most energetic slice
    auto max_energy_slice_it = slicePtrVector.begin();
    if (fSliceCaloEnergy) {
      max_energy_slice_it = std::max_element(slicePtrVector.begin(), slicePtrVector.end(), [&](const auto& slicePtr1, const auto& slicePtr2) {
        return GetSliceCaloEnergy(slicePtr1, e) < GetSliceCaloEnergy(slicePtr2, e);
      });
    } else {
      max_energy_slice_it = std::max_element(slicePtrVector.begin(), slicePtrVector.end(), [&](const auto& slicePtr1, const auto& slicePtr2) {
        return fNeutrinoRecoEnergy.CalculateNeutrinoEnergy(e, slicePtr1, true).fNuLorentzVector.E() < fNeutrinoRecoEnergy.CalculateNeutrinoEnergy(e, slicePtr2, true).fNuLorentzVector.E();
      });
    }
    art::Ptr<recob::Slice> most_energetic_slice = *max_energy_slice_it;
    if (fSliceCaloEnergy) {
      // std::cerr << "Slice energy Calo: " << GetSliceCaloEnergy(most_energetic_slice, e) / 1000 << std::endl;
      fSliceEnergy = GetSliceCaloEnergy(most_energetic_slice, e) / 1000; // Convert to GeV
    } else {
      // std::cerr << "Slice energy Neutrino: " << fNeutrinoRecoEnergy.CalculateNeutrinoEnergy(e, most_energetic_slice, true).fNuLorentzVector.E() << std::endl;
      fSliceEnergy = fNeutrinoRecoEnergy.CalculateNeutrinoEnergy(e, most_energetic_slice, true).fNuLorentzVector.E();
    }
    // Check if the slice is a neutrino
    art::FindManyP<recob::PFParticle> slicePFPAssoc(sliceHandle, e, fPFParticleLabel);
    auto pfparticlePtrVector = slicePFPAssoc.at(most_energetic_slice.key());
    auto pfparticleHandle = e.getHandle<std::vector<recob::PFParticle>>(fPFParticleLabel);
    fNoPFParticles = pfparticlePtrVector.size();
    for (const auto& pfparticlePtr : pfparticlePtrVector) {
      if (pfparticlePtr->IsPrimary() and dune_ana::DUNEAnaPFParticleUtils::IsNeutrino(pfparticlePtr)) {
        ++fNoNeutrinos;
        fRecoPDG.push_back(pfparticlePtr->PdgCode());
        fnuSliceKey.push_back(most_energetic_slice.key());
        fnuID.push_back(pfparticlePtr->Self());
        // Get metadata
        art::Ptr<larpandoraobj::PFParticleMetadata> pandoraMetaData = dune_ana::DUNEAnaPFParticleUtils::GetMetadata(pfparticlePtr, e, fPFParticleLabel);
        std::map<std::string, float> fPFPPropertiesMap = pandoraMetaData->GetPropertiesMap();
        fnuScore.push_back(fPFPPropertiesMap["NuScore"]);
        fSliceIndex.push_back(fPFPPropertiesMap["SliceIndex"]);

        // Get number of hits
        auto hits = dune_ana::DUNEAnaSliceUtils::GetHits(most_energetic_slice, e, fSliceLabel);
        fnuHits.push_back(hits.size());

        // Get vertex position
        auto nu_vertex = dune_ana::DUNEAnaPFParticleUtils::GetVertex(pfparticlePtr, e, fVertexLabel);
        fRecoVertexPosX.push_back(nu_vertex->position().X());
        fRecoVertexPosY.push_back(nu_vertex->position().Y());
        fRecoVertexPosZ.push_back(nu_vertex->position().Z());

        // Get neutrino direction from interaction (hits within the slice, 3 planes)
        dune::Point_t v_point;
        v_point.SetCoordinates(fRecoVertexPosX.back(), fRecoVertexPosY.back(), fRecoVertexPosZ.back());
        auto nu_angle = fNeutrinoRecoAngle.CalculateNeutrinoAngle(e, most_energetic_slice, v_point);
        fRecoDirCosX.push_back(nu_angle.fRecoDirection.X());
        fRecoDirCosY.push_back(nu_angle.fRecoDirection.Y());
        fRecoDirCosZ.push_back(nu_angle.fRecoDirection.Z());
        // Get neutrino energy from interaction (hits within the slice at collection plane)
        auto energy_output = fNeutrinoRecoEnergy.CalculateNeutrinoEnergy(e, most_energetic_slice, true);
        fRecoEnergy.push_back(energy_output.fNuLorentzVector.E());

        // Get daughter/child information
        fNoDaughters = GetNuDaughterInfo(
          pfparticlePtr, e, fDaughterIsTrack, fDaughterIsShower, 
          fTrackDirectionX, fTrackDirectionY, fTrackDirectionZ, 
          fShowerDirectionX, fShowerDirectionY, fShowerDirectionZ, 
          fKineticEnergyTrack, fShowerEnergy, fRecoDaughterPDG
        );

        // Fill the TTree fTreeReco
        fTreeReco->Fill();
        // Reset counters
        fNoPFParticles = 0;
        fNoNeutrinos = 0;
        fNoDaughters = 0;
        // Clear the vectors
        fRecoPDG.clear();
        fnuHits.clear();
        fRecoEnergy.clear();
        fRecoVertexPosX.clear();
        fRecoVertexPosY.clear();
        fRecoVertexPosZ.clear();
        fRecoDirCosX.clear();
        fRecoDirCosY.clear();
        fRecoDirCosZ.clear();
        fDaughterIsTrack.clear();
        fDaughterIsShower.clear();
        fTrackDirectionX.clear();
        fTrackDirectionY.clear();
        fTrackDirectionZ.clear();
        fShowerDirectionX.clear();
        fShowerDirectionY.clear();
        fShowerDirectionZ.clear();
        fKineticEnergyTrack.clear();
        fShowerEnergy.clear();
        fRecoDaughterPDG.clear();
        fnuSliceKey.clear();
        fSliceIndex.clear();
        fnuID.clear();
        fnuScore.clear();

        // Only a primary neutrino is considered per slice (the most energetic one in this case)
        break;
      } else if (pfparticlePtr == pfparticlePtrVector.back()) {
        // Fill the TTree fTreeReco with -999 values to create a useful mask for validation
        fRecoPDG.push_back(-999);
        fnuHits.push_back(0);
        fRecoEnergy.push_back(-999);
        fRecoVertexPosX.push_back(-999);
        fRecoVertexPosY.push_back(-999);
        fRecoVertexPosZ.push_back(-999);
        fRecoDirCosX.push_back(-999);
        fRecoDirCosY.push_back(-999);
        fRecoDirCosZ.push_back(-999);
        fDaughterIsTrack.push_back(-999);
        fDaughterIsShower.push_back(-999);
        fTrackDirectionX.push_back(-999);
        fTrackDirectionY.push_back(-999);
        fTrackDirectionZ.push_back(-999);
        fShowerDirectionX.push_back(-999);
        fShowerDirectionY.push_back(-999);
        fShowerDirectionZ.push_back(-999);
        fKineticEnergyTrack.push_back(-999);
        fShowerEnergy.push_back(-999);
        fRecoDaughterPDG.push_back(-999);
        fnuSliceKey.push_back(-999);
        fSliceIndex.push_back(-999);
        fnuID.push_back(-999);
        fnuScore.push_back(-999);

        // Fill the TTree fTreeReco
        fTreeReco->Fill();
        // Reset counters
        fNoPFParticles = 0;
        fNoNeutrinos = 0;
        fNoDaughters = 0;
        // Clear the vectors
        fRecoPDG.clear();
        fnuHits.clear();
        fRecoEnergy.clear();
        fRecoVertexPosX.clear();
        fRecoVertexPosY.clear();
        fRecoVertexPosZ.clear();
        fRecoDirCosX.clear();
        fRecoDirCosY.clear();
        fRecoDirCosZ.clear();
        fDaughterIsTrack.clear();
        fDaughterIsShower.clear();
        fTrackDirectionX.clear();
        fTrackDirectionY.clear();
        fTrackDirectionZ.clear();
        fShowerDirectionX.clear();
        fShowerDirectionY.clear();
        fShowerDirectionZ.clear();
        fKineticEnergyTrack.clear();
        fShowerEnergy.clear();
        fRecoDaughterPDG.clear();
        fnuSliceKey.clear();
        fSliceIndex.clear();
        fnuID.clear();
        fnuScore.clear();    
      }
    }
  }
  else {
    std::cerr << "Slice handle is not valid" << std::endl;
  }
  std::cerr << "End of event" << std::endl;
}

/**
 * @brief Get the total energy of the slice.
 * @param slicePtr art::Ptr<recob::Slice>.
 * @param e art::Event.
 * @return double, total energy of the slice.
 */
double ana::NeutrinoValidation::GetSliceCaloEnergy(const art::Ptr<recob::Slice>& slicePtr, art::Event const& e)
{
  art::ValidHandle<std::vector<recob::Slice>> sliceHandle = e.getValidHandle<std::vector<recob::Slice>>(fSliceLabel);
  art::FindManyP<recob::PFParticle> slicePFPAssoc(sliceHandle, e, fPFParticleLabel);
  std::vector<art::Ptr<recob::PFParticle>> pfparticlePtrVector = slicePFPAssoc.at(slicePtr.key());
  double total_energy = 0;

  if (!pfparticlePtrVector.empty()) {
    for (const art::Ptr<recob::PFParticle>& pfparticlePtr : pfparticlePtrVector) {
      if (dune_ana::DUNEAnaPFParticleUtils::IsTrack(pfparticlePtr, e, fPFParticleLabel, fTrackLabel)) {
        auto track = dune_ana::DUNEAnaPFParticleUtils::GetTrack(pfparticlePtr, e, fPFParticleLabel, fTrackLabel);
        auto track_calo = dune_ana::DUNEAnaTrackUtils::GetCalorimetry(track, e, fTrackLabel, fCalorimetryLabel);
        if (track_calo) {
          total_energy += track_calo->KineticEnergy();
        }
      } else if (dune_ana::DUNEAnaPFParticleUtils::IsShower(pfparticlePtr, e, fPFParticleLabel, fShowerLabel)) {
        auto shower = dune_ana::DUNEAnaPFParticleUtils::GetShower(pfparticlePtr, e, fPFParticleLabel, fShowerLabel);
        auto showerHandle = e.getHandle<std::vector<recob::Shower>>(fShowerLabel);
        art::FindManyP<anab::Calorimetry> showerCaloAssoc(showerHandle, e, "pandoraShowercalonosce"); // Same output if "pandoraShowercalonosce" -> "pandoraShowercalo"
        auto shower_calo = showerCaloAssoc.at(shower.key());
        if (!shower_calo.empty()) {
          total_energy += shower_calo.at(2)->KineticEnergy();
        }
      }
    }
  }

  return total_energy;
}

/**
 * @brief Get the daughter information of the neutrino.
 * @param pfparticlePtr art::Ptr<recob::PFParticle>.
 * @param e art::Event.
 * @param nuDaughterIsTrack std::vector<unsigned int>&.
 * @param nuDaughterIsShower std::vector<unsigned int>&.
 * @param TrackDirectionX std::vector<double>&.
 * @param TrackDirectionY std::vector<double>&.
 * @param TrackDirectionZ std::vector<double>&.
 * @param ShowerDirectionX std::vector<double>&.
 * @param ShowerDirectionY std::vector<double>&.
 * @param ShowerDirectionZ std::vector<double>&.
 * @param KineticEnergyTrack std::vector<double>&.
 * @param ShowerEnergy std::vector<double>&.
 * @param DaughterPDG std::vector<int>&.
 * @return unsigned int, number of daughters.
 */
unsigned int ana::NeutrinoValidation::GetNuDaughterInfo(
  const art::Ptr<recob::PFParticle>& pfparticlePtr, art::Event const& e,
  std::vector<unsigned int>& nuDaughterIsTrack, std::vector<unsigned int>& nuDaughterIsShower,
  std::vector<double>& TrackDirectionX, std::vector<double>& TrackDirectionY, std::vector<double>& TrackDirectionZ,
  std::vector<double>& ShowerDirectionX, std::vector<double>& ShowerDirectionY, std::vector<double>& ShowerDirectionZ, 
  std::vector<double>& KineticEnergyTrack, std::vector<double>& ShowerEnergy, std::vector<int>& DaughterPDG 
) {

  std::vector<art::Ptr<recob::PFParticle>> daughtersPFP = dune_ana::DUNEAnaPFParticleUtils::GetChildParticles(pfparticlePtr, e, fPFParticleLabel);
  for (const auto& daughterPFP : daughtersPFP) {
    if (dune_ana::DUNEAnaPFParticleUtils::IsTrack(daughterPFP, e, fPFParticleLabel, fTrackLabel)) {
      nuDaughterIsTrack.push_back(1);
      nuDaughterIsShower.push_back(0);
      auto track = dune_ana::DUNEAnaPFParticleUtils::GetTrack(daughterPFP, e, fPFParticleLabel, fTrackLabel);
      auto track_calo = dune_ana::DUNEAnaTrackUtils::GetCalorimetry(track, e, fTrackLabel, fCalorimetryLabel); // Same output for pandoracalo & pandoracalonosce
      TrackDirectionX.push_back(track->VertexDirection().X());
      TrackDirectionY.push_back(track->VertexDirection().Y());
      TrackDirectionZ.push_back(track->VertexDirection().Z());
      KineticEnergyTrack.push_back(track_calo->KineticEnergy());
      DaughterPDG.push_back(daughterPFP->PdgCode());
    } else if (dune_ana::DUNEAnaPFParticleUtils::IsShower(daughterPFP, e, fPFParticleLabel, fShowerLabel)) {
      nuDaughterIsShower.push_back(1);
      nuDaughterIsTrack.push_back(0);
      auto shower = dune_ana::DUNEAnaPFParticleUtils::GetShower(daughterPFP, e, fPFParticleLabel, fShowerLabel);
      auto showerHandle = e.getHandle<std::vector<recob::Shower>>(fShowerLabel);
      art::FindManyP<anab::Calorimetry> showerCaloAssoc(showerHandle, e, "pandoraShowercalo"); // Same output if "pandoraShowercalonosce" -> "pandoraShowercalo"
      auto shower_calo = showerCaloAssoc.at(shower.key());
      ShowerDirectionX.push_back(shower->Direction().X());
      ShowerDirectionY.push_back(shower->Direction().Y());
      ShowerDirectionZ.push_back(shower->Direction().Z());
      ShowerEnergy.push_back(shower_calo.at(2)->KineticEnergy());
      DaughterPDG.push_back(daughterPFP->PdgCode());
    //   if(!shower->Energy().empty()) {
    //   shower->Energy().at(0);
    //   shower->Energy().at(2);
    //   shower->Energy().at(shower->best_plane()); 
    } else {
        nuDaughterIsTrack.push_back(0);
        nuDaughterIsShower.push_back(0);
    }
  }

  return daughtersPFP.size();
}

/**
 * @brief Begin the job.
 * The beginJob function initializes the output TTree for the true neutrino information, 
 * the POT information and the reconstructed neutrino information.
 */
void ana::NeutrinoValidation::beginJob()
{
  // Make our handle to the TFileService
  art::ServiceHandle<art::TFileService> tfs;

  // Create TTree and branches for TTreeTrue
  fTreeTrue = tfs->make<TTree>("NeutrinoTree", "Tree containing neutrino event information");
  // General event information
  fTreeTrue -> Branch( "Event" , &fEventID, "Event/I"  );
  fTreeTrue -> Branch( "Run"   , &fRun    , "Run/I"    );
  fTreeTrue -> Branch( "SubRun", &fSubRun , "SubRun/I" );
  // Trigger info
  fTreeTrue -> Branch("TA", &fTA, "TA/I");
  // Neutrino true information
  fTreeTrue -> Branch("PDG", &fnuPDG, "PDG/I");
  fTreeTrue -> Branch("Mother", &fmother, "Mother/I");
  fTreeTrue -> Branch("VertexPositionX", &fnuVertexPosX, "VertexPosX/D"); // cm
  fTreeTrue -> Branch("VertexPositionY", &fnuVertexPosY, "VertexPosY/D"); // cm
  fTreeTrue -> Branch("VertexPositionZ", &fnuVertexPosZ, "VertexPosZ/D"); // cm
  fTreeTrue -> Branch("Energy", &fnuEnergy, "Energy/D"); // GeV
  fTreeTrue -> Branch("MomentumX", &fnuPX, "MomentumX/D"); // GeV
  fTreeTrue -> Branch("MomentumY", &fnuPY, "MomentumY/D"); // GeV
  fTreeTrue -> Branch("MomentumZ", &fnuPZ, "MomentumZ/D"); // GeV
  fTreeTrue -> Branch("Momentum", &fnuP, "Momentum/D"); // GeV
  fTreeTrue -> Branch("DirectionX", &fnuDirCosX, "DirCosX/D"); // Director cosines
  fTreeTrue -> Branch("DirectionY", &fnuDirCosY, "DirCosY/D");
  fTreeTrue -> Branch("DirectionZ", &fnuDirCosZ, "DirCosZ/D");

  // Create TTree and branches for TTreePOT
  fTreePOT = tfs->make<TTree>("POTTree", "Tree containing POT information");
  // PoT information
  fTreePOT -> Branch("POT", &fPOT, "POT/D");
  fTreePOT -> Branch("GoodPOT", &fGoodPOT, "GoodPOT/D");

  // Create TTree and branches for TTreeReco
  fTreeReco = tfs->make<TTree>("NeutrinoRecoTree", "Tree containing selected reconstructed neutrino event information");
  // General event information
  fTreeReco -> Branch( "Event" , &fEventID, "Event/I"  );
  fTreeReco -> Branch( "Run"   , &fRun    , "Run/I"    );
  fTreeReco -> Branch( "SubRun", &fSubRun , "SubRun/I" );
  // Trigger info
  fTreeReco -> Branch("TA", &fTA, "TA/I");
  // Neutrino reco event information
  fTreeReco -> Branch("No. Neutrinos", &fNoNeutrinos, "No. Neutrinos/I");
  fTreeReco -> Branch("No. PFParticles", &fNoPFParticles, "No. PFParticles/I");
  fTreeReco -> Branch("No. Daughters", &fNoDaughters, "No. Daughters/I");
  fTreeReco -> Branch("SliceEnergy", &fSliceEnergy);
  fTreeReco -> Branch("PDG", &fRecoPDG);
  fTreeReco -> Branch("DaughterPDG", &fRecoDaughterPDG);
  fTreeReco -> Branch("Hits", &fnuHits);
  fTreeReco -> Branch("Energy", &fRecoEnergy); 
  fTreeReco -> Branch("VertexPositionX", &fRecoVertexPosX);
  fTreeReco -> Branch("VertexPositionY", &fRecoVertexPosY);
  fTreeReco -> Branch("VertexPositionZ", &fRecoVertexPosZ);
  fTreeReco -> Branch("DirectionX", &fRecoDirCosX);
  fTreeReco -> Branch("DirectionY", &fRecoDirCosY);
  fTreeReco -> Branch("DirectionZ", &fRecoDirCosZ);
  // Neutrino daughter reco event information
  fTreeReco -> Branch("DaughterIsTrack", &fDaughterIsTrack);
  fTreeReco -> Branch("DaughterIsShower", &fDaughterIsShower);
  fTreeReco -> Branch("TrackDirectionX", &fTrackDirectionX);
  fTreeReco -> Branch("TrackDirectionY", &fTrackDirectionY);
  fTreeReco -> Branch("TrackDirectionZ", &fTrackDirectionZ);
  fTreeReco -> Branch("ShowerDirectionX", &fShowerDirectionX);
  fTreeReco -> Branch("ShowerDirectionY", &fShowerDirectionY);
  fTreeReco -> Branch("ShowerDirectionZ", &fShowerDirectionZ);
  fTreeReco -> Branch("KineticEnergyTrack", &fKineticEnergyTrack);
  fTreeReco -> Branch("ShowerEnergy", &fShowerEnergy);
  // Metadata information
  fTreeReco -> Branch("NuSliceKey", &fnuSliceKey);
  fTreeReco -> Branch("NuSliceIndex", &fSliceIndex);
  fTreeReco -> Branch("NuID", &fnuID);
  fTreeReco -> Branch("NuScore", &fnuScore);
}

/**
 * @brief Begin the subrun.
 * @param subRun art::SubRun.
 * The beginSubRun function is called at the beginning of each subrun.
 * The function gets the POT information from the POTSummary and fills the fTreePOT tree.
 */
void ana::NeutrinoValidation::beginSubRun(art::SubRun const& subRun)
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
void ana::NeutrinoValidation::endSubRun(art::SubRun const& subRun) {}

/**
 * @brief End the job.
 * The endJob function is called at the end of the job.
 */
void ana::NeutrinoValidation::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(ana::NeutrinoValidation)
