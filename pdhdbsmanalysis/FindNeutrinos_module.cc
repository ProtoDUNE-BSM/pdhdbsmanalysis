////////////////////////////////////////////////////////////////////////
// Class:       FindNeutrinos
// Plugin Type: analyzer (Unknown Unknown)
// File:        FindNeutrinos_module.cc
//
// Generated at Fri Feb 14 05:35:05 2025 by Dario Pullia using cetskelgen
// from cetlib version 3.18.02.
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


// #include "messagefacility/MessageLogger/MessageLogger.h"
// #include "art/Framework/Services/Registry/ServiceHandle.h"

#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "dunereco/AnaUtils/DUNEAnaPFParticleUtils.h"
#include "dunereco/AnaUtils/DUNEAnaEventUtils.h"
#include "dunereco/AnaUtils/DUNEAnaShowerUtils.h"
#include "dunereco/AnaUtils/DUNEAnaTrackUtils.h"
#include "dunereco/AnaUtils/DUNEAnaUtilsBase.h"
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




namespace NeutrinoAna {
  class FindNeutrinos;
}


class NeutrinoAna::FindNeutrinos : public art::EDAnalyzer {
public:
  explicit FindNeutrinos(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  FindNeutrinos(FindNeutrinos const&) = delete;
  FindNeutrinos(FindNeutrinos&&) = delete;
  FindNeutrinos& operator=(FindNeutrinos const&) = delete;
  FindNeutrinos& operator=(FindNeutrinos&&) = delete;

  double GetTotalEnergy(const art::Ptr<recob::Slice>& slicePtr, art::Event const& e);

  std::vector<double> GetDaugtherInfoDFS(const art::Ptr<recob::PFParticle>& pfparticlePtr, art::Event const& e);
  std::vector<TruthMatchUtils::G4ID> GetTrueInfoChainDFS(const art::Ptr<recob::PFParticle>& pfparticlePtr, art::Event const& e);
  TruthMatchUtils::G4ID GetTrueInfo(const art::Ptr<recob::PFParticle>& pfparticlePtr, art::Event const& e);


  // Required functions.
  void analyze(art::Event const& e) override;
  void beginSubRun(art::SubRun const& subRun) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Create output TTree
  TTree *fTree_reco;
  TTree *fTree_truth;
  
  // Tree variables reco
  unsigned int fEventID_reco;
  double fVx_reco;
  double fVy_reco;
  double fVz_reco;
  int fPdgCode_reco;
  double fEnergy_reco;
  double fDirectionX_reco;
  double fDirectionY_reco;
  double fDirectionZ_reco;
  int fNHits_reco;
  int fNPFPs_reco;
  int fTrueOriginID_reco;
  int fEventNumber_reco;
  
  // Tree variables true
  unsigned int fEventID_true;
  double fVx_true;
  double fVy_true;
  double fVz_true;
  double fE_true;
  double fPx_true;
  double fPy_true;
  double fPz_true;
  int fPdgCode_true;
  int fMother_true;
  double fPOT_true;
  double fGoodPOT_true;
  double fTotalPOT_true;
  int fTA_true;
  int fEventNumber_true;



  // Declare member data here.
  std::string fTrackLabel;
  std::string fShowerLabel;
  std::string fVertexLabel;
  std::string fSliceLabel;
  std::string fPFParticleLabel;
  std::string fHitsModuleLabel;
  std::string fMCParticleLabel;
  std::string fMCTruthLabel;
  std::string fCalorimetryLabel;
  std::string fTALabel;
  std::string fOutputFolder;
  std::string fOutputFileName;
  bool fGetTruth;

  int fTA;
  

  // 
  const bool fRollUpUnsavedIDs = true;
  double fPOT = 0;
  double fTotalPOT = 0;
  double fGoodPOT = 0;
  //
  int fEventNumber = 0;

};


NeutrinoAna::FindNeutrinos::FindNeutrinos(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}, 
  // More initializers here.
  fTrackLabel(p.get<std::string>("TrackLabel")),
  fShowerLabel(p.get<std::string>("ShowerLabel")),
  fVertexLabel(p.get<std::string>("VertexLabel")),
  fSliceLabel(p.get<std::string>("SliceLabel")),
  fPFParticleLabel(p.get<std::string>("PFParticleLabel")),
  fHitsModuleLabel(p.get<std::string>("HitsModuleLabel")),
  fMCParticleLabel(p.get<std::string>("MCParticleLabel")),
  fMCTruthLabel(p.get<std::string>("MCTruthLabel")),
  fCalorimetryLabel(p.get<std::string>("CalorimetryLabel")),
  fTALabel(p.get<std::string>("TALabel")),
  fOutputFolder(p.get<std::string>("OutputFolder")),
  fOutputFileName(p.get<std::string>("OutputFileName")),
  fGetTruth(p.get<bool>("GetTruth"))
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void NeutrinoAna::FindNeutrinos::analyze(art::Event const& e)
{
  // Implementation of required member function here.
  // The plan here is to go over the events and find if a neutrino is present.
  // I will iterate over the slices to get the most energetic one.
  // Then I will get the PFParticles associated with the slice. 
  // If the PFParticle is a neutrino, I will get the vertex and the track associated with it. 
  // The vertex should be well inside the detector, and I want a number of daughters of the neutrino.
  
  // Initialize the variables to 0
  fEventID_reco = 0;
  fVx_reco = 0;
  fVy_reco = 0;
  fVz_reco = 0;
  fPdgCode_reco = 0;
  fEnergy_reco = 0;
  fDirectionX_reco = 0;
  fDirectionY_reco = 0;
  fDirectionZ_reco = 0;
  fNHits_reco = 0;
  fNPFPs_reco = 0;
  fTrueOriginID_reco = 0;
  fEventNumber_reco = 0;  

  fEventID_true = 0;
  fVx_true = 0;
  fVy_true = 0;
  fVz_true = 0;
  fE_true = 0;
  fPx_true = 0;
  fPy_true = 0;
  fPz_true = 0;
  fPdgCode_true = 0;
  fMother_true = 0;
  fPOT_true = 0;
  fGoodPOT_true = 0;
  fTotalPOT_true = 0;
  fTA_true = 0;
  fEventNumber_true = 0;
  // open the output files
  std::ofstream outfile_truth;
  outfile_truth.open(fOutputFolder + "/" + fOutputFileName+"_truth.txt", std::ios_base::app);
  // outfile_truth<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<std::endl;

  std::ofstream outfile_recos;
  outfile_recos.open(fOutputFolder + "/" + fOutputFileName+"_pfparticle.txt", std::ios_base::app);
  fEventNumber++;
  outfile_recos<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<fEventNumber<<std::endl;  
  

  // get the trigger information
  art::Handle<std::vector<dunedaq::trgdataformats::TriggerActivityData>> taHandle;
  if (!taHandle.isValid()) {
    fTA = -1;
  } 

 if (!e.getByLabel(fTALabel, taHandle)) {
    fTA = 0;
  } else {
    if (taHandle->size() == 0) {
      fTA = 0;
    } else {
      fTA = 1;
    }
  } 
  double fE = 0;

  if (fGetTruth) {
    auto truthHandle = e.getValidHandle<std::vector<simb::MCTruth>>(fMCTruthLabel);
    for (auto const& truth : (*truthHandle)) {
      if (truth.NeutrinoSet()) {
        const auto &nu = truth.GetNeutrino();
        const auto &neutrino = nu.Nu();
        fE = neutrino.E();
        outfile_truth << neutrino.Vx() << " " << neutrino.Vy() << " " << neutrino.Vz() << " " << fE << " " << neutrino.Px() << " " << neutrino.Py() << " " << neutrino.Pz() << " " << neutrino.PdgCode() << " " << neutrino.Mother() << " " << fPOT << " " << fGoodPOT << " " << fTotalPOT << " " << fTA << std::endl;
        // fill the variables
        fEventID_true = e.id().event();
        fVx_true = neutrino.Vx();
        fVy_true = neutrino.Vy();
        fVz_true = neutrino.Vz();
        fE_true = fE;
        fPx_true = neutrino.Px();
        fPy_true = neutrino.Py();
        fPz_true = neutrino.Pz();
        fPdgCode_true = neutrino.PdgCode();
        fMother_true = neutrino.Mother();
        fPOT_true = fPOT;
        fGoodPOT_true = fGoodPOT;
        fTotalPOT_true = fTotalPOT;
        fTA_true = fTA;
        fEventNumber_true = fEventNumber;
      }
    }
  }
  fTree_truth->Fill();


  // -------------------------------------------------------------------
  
  art::Handle<std::vector<recob::Slice>> sliceHandle = e.getHandle<std::vector<recob::Slice>>(fSliceLabel);
  if (sliceHandle.isValid()){
    std::vector<art::Ptr<recob::Slice>> slicePtrVector;
    art::fill_ptr_vector(slicePtrVector, sliceHandle);
    std::cout << "Number of Slices: " << slicePtrVector.size() << std::endl;

    // get the most energetic slice
    double max_energy = 0;
    double energy = 0;
    art::Ptr<recob::Slice> most_energetic_slice;
    for (const art::Ptr<recob::Slice>& slicePtr : slicePtrVector) {
      energy = GetTotalEnergy(slicePtr, e);
      if (energy > max_energy) {
        max_energy = energy;
        most_energetic_slice = slicePtr;
      }
    }
    std::cout << "Most energetic slice: " << max_energy << std::endl;
    // check if the slice is a neutrino
    art::FindManyP<recob::PFParticle> slicePFPAssoc(sliceHandle, e, fSliceLabel);
    std::vector<art::Ptr<recob::PFParticle>> pfparticlePtrVector = slicePFPAssoc.at(most_energetic_slice.key());
    art::Handle<std::vector<recob::PFParticle>> pfparticleHandle = e.getHandle<std::vector<recob::PFParticle>>(fPFParticleLabel);

    for (const art::Ptr<recob::PFParticle>& pfparticlePtr : pfparticlePtrVector) {
      if (pfparticlePtr->IsPrimary() and (pfparticlePtr->PdgCode() == 12 or pfparticlePtr->PdgCode() == 14 or pfparticlePtr->PdgCode() == 16 or pfparticlePtr->PdgCode() == -12 or pfparticlePtr->PdgCode() == -14 or pfparticlePtr->PdgCode() == -16)) {
        std::cout << "Neutrino found" << std::endl;
      }


      art::FindManyP<recob::Vertex> pfVertexAssoc(pfparticleHandle, e, fVertexLabel);
      std::vector<art::Ptr<recob::Vertex>> vertices = pfVertexAssoc.at(pfparticlePtr.key());
      
      std::vector<double> info = GetDaugtherInfoDFS(pfparticlePtr, e);
      
      int trueOriginID = -999;
      if (fGetTruth){
        trueOriginID = GetTrueInfo(pfparticlePtr, e);
      }

      double vertex_x, vertex_y, vertex_z;
      if (vertices.size() > 0) {
        vertex_x = vertices.at(0)->position().X();
        vertex_y = vertices.at(0)->position().Y();
        vertex_z = vertices.at(0)->position().Z();
      }
      else {
        vertex_x = -999;
        vertex_y = -999;
        vertex_z = -999;
      }

      int pdg = pfparticlePtr->PdgCode();

      outfile_recos << vertex_x << " " << vertex_y << " " << vertex_z <<  " " << pdg << " " << info[0] << " " << info[1] << " " << info[2] << " " << info[3] << " " << info[4] << " " << info[5] << " " << trueOriginID << std::endl;
      // fill the variables
      fEventID_reco = e.id().event();
      fVx_reco = vertex_x;
      fVy_reco = vertex_y;
      fVz_reco = vertex_z;
      fPdgCode_reco = pdg;
      fEnergy_reco = info[2];
      fDirectionX_reco = info[3];
      fDirectionY_reco = info[4];
      fDirectionZ_reco = info[5];
      fNHits_reco = info[0];
      fNPFPs_reco = info[1]; 
      fTrueOriginID_reco = trueOriginID;
      fEventNumber_reco = fEventNumber;
      fTree_reco->Fill();

    }
  }
  else {
    std::cout << "Slice handle is not valid" << std::endl;
  }

  std::cout << "End of event" << std::endl;


  outfile_truth.close();
  outfile_recos.close();

}

double NeutrinoAna::FindNeutrinos::GetTotalEnergy(const art::Ptr<recob::Slice>& slicePtr, art::Event const& e){

  art::ValidHandle<std::vector<recob::Slice>> sliceHandle = e.getValidHandle<std::vector<recob::Slice>>(fSliceLabel);

  art::FindManyP<recob::PFParticle> slicePFPAssoc(sliceHandle, e, fSliceLabel);
  std::vector<art::Ptr<recob::PFParticle>> pfparticlePtrVector = slicePFPAssoc.at(slicePtr.key());


  double total_energy = 0;

  for (const art::Ptr<recob::PFParticle>& pfparticlePtr : pfparticlePtrVector) {
    if (dune_ana::DUNEAnaPFParticleUtils::IsTrack(pfparticlePtr, e, fPFParticleLabel, fTrackLabel)) {
      // get the track
      art::Ptr<recob::Track> this_track =  dune_ana::DUNEAnaPFParticleUtils::GetTrack(pfparticlePtr, e, fPFParticleLabel, fTrackLabel);
      art::Ptr<anab::Calorimetry> this_calo = dune_ana::DUNEAnaTrackUtils::GetCalorimetry(this_track, e, fTrackLabel, fCalorimetryLabel);
      total_energy += this_calo->KineticEnergy();

    }
    else if (dune_ana::DUNEAnaPFParticleUtils::IsShower(pfparticlePtr, e, fPFParticleLabel, fShowerLabel)) {
      // get the shower
      art::Ptr<recob::Shower> this_shower =  dune_ana::DUNEAnaPFParticleUtils::GetShower(pfparticlePtr, e, fPFParticleLabel, fShowerLabel);
      // do not use the utility, use the association
      art::Handle<std::vector<recob::Shower>> showerHandle = e.getHandle<std::vector<recob::Shower>>(fShowerLabel);
      art::FindManyP<anab::Calorimetry> showerCaloAssoc(showerHandle, e, "pandoraShowercalonosce");
      std::vector<art::Ptr<anab::Calorimetry>> calos = showerCaloAssoc.at(this_shower.key());
      total_energy += calos.at(0)->KineticEnergy();
    }
  }

  return total_energy;
}

int NeutrinoAna::FindNeutrinos::GetTrueInfo(const art::Ptr<recob::PFParticle>& pfparticlePtr, art::Event const& e) {
  TruthMatchUtils::G4ID trueID = 0;
  std::vector<TruthMatchUtils::G4ID> trueOriginIDs_vector;
  // get the daughter pfparticles
  std::cout << "Getting true info" << std::endl;
  trueOriginIDs_vector = GetTrueInfoChainDFS(pfparticlePtr, e);
  std::cout << "Got true info" << std::endl;

  for (TruthMatchUtils::G4ID trueOriginID : trueOriginIDs_vector) {
    if (trueOriginID == 1) {
      trueID++;
    }
  }
  

  return trueID;
}

std::vector<TruthMatchUtils::G4ID> NeutrinoAna::FindNeutrinos::GetTrueInfoChainDFS(const art::Ptr<recob::PFParticle>& pfparticlePtr, art::Event const& e) {
  std::vector<TruthMatchUtils::G4ID> trueIDs_vector;
  // get the daughter pfparticles
  std::vector<art::Ptr<recob::PFParticle>> daughters = dune_ana::DUNEAnaPFParticleUtils::GetChildParticles(pfparticlePtr, e, fPFParticleLabel);

  if (daughters.size() == 0) {
    // get the hits
    art::Handle<std::vector<recob::PFParticle>> pfparticleHandle = e.getHandle<std::vector<recob::PFParticle>>(fPFParticleLabel);
    std::vector<art::Ptr<recob::Hit>> hits = dune_ana::DUNEAnaPFParticleUtils::GetHits(pfparticlePtr, e, fPFParticleLabel);
    // get the true info
    const auto clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);
    TruthMatchUtils::G4ID hitID(TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData, hits, fRollUpUnsavedIDs));


    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    const simb::MCParticle* mcparticle = pi_serv->TrackIdToParticle_P(hitID);
    // check if the pointer is valid
    if (mcparticle == nullptr) {
      trueIDs_vector.push_back(0);
      
    } else{
      int trackID = mcparticle->TrackId();
      simb::MCTruth mcTruth = pi_serv->TrackIdToMCTruth(trackID);
      trueIDs_vector.push_back(mcTruth.Origin());  
    }

  }
  else {
    for (const art::Ptr<recob::PFParticle>& daughter : daughters) {
      std::vector<TruthMatchUtils::G4ID> daughter_trueIDs = GetTrueInfoChainDFS(daughter, e);
      trueIDs_vector.insert(trueIDs_vector.end(), daughter_trueIDs.begin(), daughter_trueIDs.end());
    }
  }
  return trueIDs_vector;
}

std::vector<double> NeutrinoAna::FindNeutrinos::GetDaugtherInfoDFS(const art::Ptr<recob::PFParticle>& pfparticlePtr, art::Event const& e) {
  // i want the total number of hits, the total number of pfp, the total energy
  std::vector<double> info;
  info.push_back(0); // Number of hits
  info.push_back(0); // Number of pfparticles
  info.push_back(0); // Energy
  info.push_back(0); // Direction X
  info.push_back(0); // Direction Y
  info.push_back(0); // Direction Z



  // get the daughter pfparticles
  std::vector<art::Ptr<recob::PFParticle>> daughters = dune_ana::DUNEAnaPFParticleUtils::GetChildParticles(pfparticlePtr, e, fPFParticleLabel);

  if (daughters.size() == 0) {
    // get the hits
    art::Handle<std::vector<recob::PFParticle>> pfparticleHandle = e.getHandle<std::vector<recob::PFParticle>>(fPFParticleLabel);
    std::vector<art::Ptr<recob::Hit>> hits = dune_ana::DUNEAnaPFParticleUtils::GetHits(pfparticlePtr, e, fPFParticleLabel);
    info[0] = hits.size();

    // get the number of pfparticles
    info[1] = 1;

    
    // get the energy
    if (dune_ana::DUNEAnaPFParticleUtils::IsTrack(pfparticlePtr, e, fPFParticleLabel, fTrackLabel)) {
      // get the track
      art::Ptr<recob::Track> this_track =  dune_ana::DUNEAnaPFParticleUtils::GetTrack(pfparticlePtr, e, fPFParticleLabel, fTrackLabel);
      art::Ptr<anab::Calorimetry> this_calo = dune_ana::DUNEAnaTrackUtils::GetCalorimetry(this_track, e, fTrackLabel, fCalorimetryLabel);
      info[2] = this_calo->KineticEnergy();

      // get the direction
      info[3] = this_track->VertexDirection().X();
      info[4] = this_track->VertexDirection().Y();
      info[5] = this_track->VertexDirection().Z();

    }
    else if (dune_ana::DUNEAnaPFParticleUtils::IsShower(pfparticlePtr, e, fPFParticleLabel, fShowerLabel)) {
      // get the shower
      art::Ptr<recob::Shower> this_shower =  dune_ana::DUNEAnaPFParticleUtils::GetShower(pfparticlePtr, e, fPFParticleLabel, fShowerLabel);
      // do not use the utility, use the association
      art::Handle<std::vector<recob::Shower>> showerHandle = e.getHandle<std::vector<recob::Shower>>(fShowerLabel);
      art::FindManyP<anab::Calorimetry> showerCaloAssoc(showerHandle, e, "pandoraShowercalonosce");
      std::vector<art::Ptr<anab::Calorimetry>> calos = showerCaloAssoc.at(this_shower.key());
      info[2] = calos.at(0)->KineticEnergy();
      
      // get the direction
      info[3] = this_shower->Direction().X();
      info[4] = this_shower->Direction().Y();
      info[5] = this_shower->Direction().Z();
    }

  }
  else {
    for (const art::Ptr<recob::PFParticle>& daughter : daughters) {
      std::vector<double> daughter_info = GetDaugtherInfoDFS(daughter, e);
      info[0] += daughter_info[0];
      info[1] += daughter_info[1];
      info[2] += daughter_info[2];
      info[3] += daughter_info[3];  
      info[4] += daughter_info[4];
      info[5] += daughter_info[5];
    }
  }
  return info;
}

void NeutrinoAna::FindNeutrinos::beginSubRun(art::SubRun const& subRun) {
  
  if (fGetTruth){

    const auto potSummaryHandle = subRun.getValidHandle<sumdata::POTSummary>("generator");
    const auto &potSummary = *potSummaryHandle;
    fPOT = potSummary.totpot;
    fGoodPOT = potSummary.totgoodpot;
    fTotalPOT += fPOT;
  
  }

}

void NeutrinoAna::FindNeutrinos::beginJob()
{
  // Implementation of optional member function here.
  // Get the TFileService to create the output TTree for us
  art::ServiceHandle<art::TFileService> tfs;
  // fTree = tfs->make<TTree>("tree", "Output TTree");
  fTree_reco = tfs->make<TTree>("tree_reco", "Output TTree reco");
  fTree_truth = tfs->make<TTree>("tree_truth", "Output TTree truth");

  // Add branches to TTree
  fTree_reco->Branch("eventID", &fEventID_reco);
  fTree_reco->Branch("vx", &fVx_reco);
  fTree_reco->Branch("vy", &fVy_reco);
  fTree_reco->Branch("vz", &fVz_reco);
  fTree_reco->Branch("pdgCode", &fPdgCode_reco);
  fTree_reco->Branch("energy", &fEnergy_reco);
  fTree_reco->Branch("directionX", &fDirectionX_reco);
  fTree_reco->Branch("directionY", &fDirectionY_reco);
  fTree_reco->Branch("directionZ", &fDirectionZ_reco);
  fTree_reco->Branch("nHits", &fNHits_reco);
  fTree_reco->Branch("nPFPs", &fNPFPs_reco);
  fTree_reco->Branch("trueOriginID", &fTrueOriginID_reco);

  fTree_truth->Branch("eventID", &fEventID_true);
  fTree_truth->Branch("vx", &fVx_true);
  fTree_truth->Branch("vy", &fVy_true);
  fTree_truth->Branch("vz", &fVz_true);
  fTree_truth->Branch("E", &fE_true);
  fTree_truth->Branch("Px", &fPx_true);
  fTree_truth->Branch("Py", &fPy_true);
  fTree_truth->Branch("Pz", &fPz_true);
  fTree_truth->Branch("pdgCode", &fPdgCode_true);
  fTree_truth->Branch("mother", &fMother_true);
  fTree_truth->Branch("POT", &fPOT_true);
  fTree_truth->Branch("goodPOT", &fGoodPOT_true);
  fTree_truth->Branch("totalPOT", &fTotalPOT_true);
  fTree_truth->Branch("TA", &fTA_true);


}

void NeutrinoAna::FindNeutrinos::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(NeutrinoAna::FindNeutrinos)
