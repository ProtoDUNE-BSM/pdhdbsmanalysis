////////////////////////////////////////////////////////////////////////
// Class:       AnalyzeEvents
// Plugin Type: analyzer (Unknown Unknown)
// File:        AnalyzeEvents_module.cc
//
// Generated at Fri Nov 29 04:45:27 2024 by Dario Pullia using cetskelgen
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




namespace test {
  class AnalyzeEvents;
}

class test::AnalyzeEvents : public art::EDAnalyzer {
public:
  explicit AnalyzeEvents(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  AnalyzeEvents(AnalyzeEvents const&) = delete;
  AnalyzeEvents(AnalyzeEvents&&) = delete;
  AnalyzeEvents& operator=(AnalyzeEvents const&) = delete;
  AnalyzeEvents& operator=(AnalyzeEvents&&) = delete;

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

  // define input labels
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


  // PDG vector
  std::vector<int> pdg_vector;
  std::map<int, int> pdg_count;

  // 
  const bool fRollUpUnsavedIDs = true;

  // PoT variables
  double fPOT = 0;
  double fTotalPOT = 0;
  double fGoodPOT = 0;



};

test::AnalyzeEvents::AnalyzeEvents(fhicl::ParameterSet const& p)
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

void test::AnalyzeEvents::analyze(art::Event const& e)
{
  // Implementation of required member function here.
  // Set the event ID
  // unsigned int EventID = e.id().event();
  // std::cout << "Event ID: " << EventID << std::endl;
  // auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);

  // int fSimPDG = 0;
  // int fCCNC = 0;
  double fE = 0;
  std::ofstream outfile_truth;
  outfile_truth.open(fOutputFolder + "/" + fOutputFileName+"_truth.txt", std::ios_base::app);

  std::ofstream outfile_recos;
  outfile_recos.open(fOutputFolder + "/" + fOutputFileName+"_pfparticle.txt", std::ios_base::app);
  outfile_recos<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<std::endl;

  int fTA = 0;

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

  if (fGetTruth) {
    auto truthHandle = e.getValidHandle<std::vector<simb::MCTruth>>(fMCTruthLabel);
    for (auto const& truth : (*truthHandle)) {
      if (truth.NeutrinoSet()) {
        const auto &nu = truth.GetNeutrino();
        const auto &neutrino = nu.Nu();

        // fSimPDG = neutrino.PdgCode();
        // fCCNC = nu.CCNC();
        fE = neutrino.E();
        
        // std::cout << "Neutrino truth Origin " << truth.Origin() << std::endl;
        // std::cout << "PDG: " << fSimPDG << " CCNC: " << fCCNC << " E: " << fE << " Vertex (XYZ): " << neutrino.Vx() << " " << neutrino.Vy() << " " << neutrino.Vz() << " Mother: " << neutrino.Mother() << " Process: " << neutrino.Process() << std::endl;

        // print to file
        // std::cout<<"Output file: "<<fOutputFolder + "/" + fOutputFileName+"_truth.txt"<<std::endl;
        outfile_truth << neutrino.Vx() << " " << neutrino.Vy() << " " << neutrino.Vz() << " " << fE << " " << neutrino.Px() << " " << neutrino.Py() << " " << neutrino.Pz() << " " << neutrino.PdgCode() << " " << neutrino.Mother() << " " << fPOT << " " << fGoodPOT << " " << fTotalPOT << " " << fTA << std::endl;
      }
    }
  }


  // get the vector of pointers to pfparticles
  art::Handle<std::vector<recob::PFParticle>> pfparticleHandle = e.getHandle<std::vector<recob::PFParticle>>(fPFParticleLabel);
  if (pfparticleHandle.isValid()){
    std::vector<art::Ptr<recob::PFParticle>> pfparticlePtrVector;
    art::fill_ptr_vector(pfparticlePtrVector, pfparticleHandle);
    // std::cout << "Number of PFParticles: " << pfparticlePtrVector.size() << std::endl;


    for (const art::Ptr<recob::PFParticle>& pfparticlePtr : pfparticlePtrVector) {
      if (pfparticlePtr->IsPrimary() and (pfparticlePtr->PdgCode() == 12 or pfparticlePtr->PdgCode() == 14 or pfparticlePtr->PdgCode() == 16 or pfparticlePtr->PdgCode() == -12 or pfparticlePtr->PdgCode() == -14 or pfparticlePtr->PdgCode() == -16)) {
        // std::cout << "Neutrino found" << std::endl;
        // std::cout << "PDG: " << pfparticlePtr->PdgCode() << " Parent: " << pfparticlePtr->Parent() << " Self: " << pfparticlePtr->Self() << std::endl;
        // get the vertex
        art::FindManyP<recob::Vertex> pfVertexAssoc(pfparticleHandle, e, fVertexLabel);
        std::vector<art::Ptr<recob::Vertex>> vertices = pfVertexAssoc.at(pfparticlePtr.key());
        // std::cout << "Vertex position: (" << vertices.at(0)->position().X() << ", " << vertices.at(0)->position().Y() << ", " << vertices.at(0)->position().Z() << ")" << std::endl;

        // get the track
        std::vector<double> info = GetDaugtherInfoDFS(pfparticlePtr, e);
        

        int trueOriginID = -999;
        if (fGetTruth){
          trueOriginID = GetTrueInfo(pfparticlePtr, e);
        }
        outfile_recos << vertices.at(0)->position().X() << " " << vertices.at(0)->position().Y() << " " << vertices.at(0)->position().Z() <<  " " << pfparticlePtr->PdgCode() << " " << info[0] << " " << info[1] << " " << info[2] << " " << info[3] << " " << info[4] << " " << info[5] << " " << trueOriginID << std::endl;
      }
    }

  }
  else {
    // std::cout << "PFParticle handle is not valid" << std::endl;
  }
  outfile_truth.close();
  outfile_recos.close();
}

// parameters for the subrun functions



void test::AnalyzeEvents::beginSubRun(art::SubRun const& subRun) {
  
  const auto potSummaryHandle = subRun.getValidHandle<sumdata::POTSummary>("generator");
  const auto &potSummary = *potSummaryHandle;
  fPOT = potSummary.totpot;
  fGoodPOT = potSummary.totgoodpot;
  fTotalPOT += fPOT;

}

int test::AnalyzeEvents::GetTrueInfo(const art::Ptr<recob::PFParticle>& pfparticlePtr, art::Event const& e) {
  TruthMatchUtils::G4ID trueID = 0;
  std::vector<TruthMatchUtils::G4ID> trueOriginIDs_vector;
  // get the daughter pfparticles
  trueOriginIDs_vector = GetTrueInfoChainDFS(pfparticlePtr, e);

  for (TruthMatchUtils::G4ID trueOriginID : trueOriginIDs_vector) {
    if (trueOriginID == 1) {
      trueID++;
    }
  }
  

  return trueID;
}

std::vector<TruthMatchUtils::G4ID> test::AnalyzeEvents::GetTrueInfoChainDFS(const art::Ptr<recob::PFParticle>& pfparticlePtr, art::Event const& e) {
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

std::vector<double> test::AnalyzeEvents::GetDaugtherInfoDFS(const art::Ptr<recob::PFParticle>& pfparticlePtr, art::Event const& e) {
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



void test::AnalyzeEvents::beginJob()
{
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;

}

void test::AnalyzeEvents::endJob()
{
  // Implementation of optional member function here.

}

DEFINE_ART_MODULE(test::AnalyzeEvents)
