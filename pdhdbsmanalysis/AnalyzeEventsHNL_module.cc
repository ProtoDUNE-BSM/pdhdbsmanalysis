////////////////////////////////////////////////////////////////////////
// Class:       AnalyzeEventsHNL
// Plugin Type: analyzer (Unknown Unknown)
// File:        AnalyzeEventsNeutrino_module.cc
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
#include "dunereco/FDSensOpt/NeutrinoAngularRecoAlg/NeutrinoAngularRecoAlg.h"
#include "dunereco/FDSensOpt/NeutrinoEnergyRecoAlg/NeutrinoEnergyRecoAlg.h"

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




namespace hnlAna {
  class AnalyzeEventsHNL;
}


class hnlAna::AnalyzeEventsHNL : public art::EDAnalyzer {
public:
  explicit AnalyzeEventsHNL(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  AnalyzeEventsHNL(AnalyzeEventsHNL const&) = delete;
  AnalyzeEventsHNL(AnalyzeEventsHNL&&) = delete;
  AnalyzeEventsHNL& operator=(AnalyzeEventsHNL const&) = delete;
  AnalyzeEventsHNL& operator=(AnalyzeEventsHNL&&) = delete;

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
  TTree *fTree_noCuts;
  TTree *fTree_fiducialPass;
  TTree *fTree_dirZPass;
  TTree *fTree_isTrackPass;
  TTree *fTree_isShowerPass;

  // dune reco functions
  dune::NeutrinoAngularRecoAlg fNeutrinoRecoAngle;
  dune::NeutrinoEnergyRecoAlg fNeutrinoRecoEnergy;

  
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
  int fPassCut;
  int fPassCut_after_fiducial;
  int fPassCut_after_directioncut;
  int fPassCut_after_trackcut;

  // Tree variables true
  unsigned int fEventID_true;
  double fVx_true;
  double fVy_true;
  double fVz_true;
  double fE_true;
  double fPx_true;
  double fPy_true;
  double fPz_true;
  double fP_true;
  int fPdgCode_true;
  int fMother_true;
  double fPOT_true;
  double fGoodPOT_true;
  double fTotalPOT_true;
  int fTA_true;
  int fEventNumber_true;
  double fDirCosX_true;
  double fDirCosY_true;
  double fDirCosZ_true;
  double fKineticEnergy_true; 



  unsigned int fEventID_noCuts;
  std::vector<double> fVx_noCuts;
  std::vector<double> fVy_noCuts;
  std::vector<double> fVz_noCuts;
  std::vector<double> fEnergy_noCuts;
  std::vector<double> fDirectionX_noCuts;
  std::vector<double> fDirectionY_noCuts;
  std::vector<double> fDirectionZ_noCuts;
  std::vector<int> fTrueOriginID_noCuts;

  unsigned int fEventID_fiducialPass;
  std::vector<double> fVx_fiducialPass; 
  std::vector<double> fVy_fiducialPass;
  std::vector<double> fVz_fiducialPass;
  std::vector<double> fEnergy_fiducialPass;
  std::vector<double> fDirectionX_fiducialPass;
  std::vector<double> fDirectionY_fiducialPass;
  std::vector<double> fDirectionZ_fiducialPass;
  std::vector<int> fTrueOriginID_fiducialPass;

  unsigned int fEventID_dirZPass;
  std::vector<double> fVx_dirZPass;
  std::vector<double> fVy_dirZPass;
  std::vector<double> fVz_dirZPass;
  std::vector<double> fEnergy_dirZPass;
  std::vector<double> fDirectionX_dirZPass;
  std::vector<double> fDirectionY_dirZPass;
  std::vector<double> fDirectionZ_dirZPass;
  std::vector<int> fTrueOriginID_dirZPass;

  unsigned int fEventID_isTrackPass;
  std::vector<double> fVx_isTrackPass;
  std::vector<double> fVy_isTrackPass;
  std::vector<double> fVz_isTrackPass;
  std::vector<double> fEnergy_isTrackPass;
  std::vector<double> fDirectionX_isTrackPass;
  std::vector<double> fDirectionY_isTrackPass;
  std::vector<double> fDirectionZ_isTrackPass;
  std::vector<int> fTrueOriginID_isTrackPass;

  unsigned int fEventID_isShowerPass;
  std::vector<double> fVx_isShowerPass;
  std::vector<double> fVy_isShowerPass;
  std::vector<double> fVz_isShowerPass;
  std::vector<double> fEnergy_isShowerPass;
  std::vector<double> fDirectionX_isShowerPass;
  std::vector<double> fDirectionY_isShowerPass;
  std::vector<double> fDirectionZ_isShowerPass;
  std::vector<int> fTrueOriginID_isShowerPass;
  std::vector<int> fIsTrack;
  std::vector<int> fIsShower;

  




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
  int total_events = 0;
  int pass_energy = 0;
  int pass_track_multiplicity = 0;
  int pass_all_cuts = 0;

  

  // 
  const bool fRollUpUnsavedIDs = true;
  double fPOT = 0;
  double fTotalPOT = 0;
  double fGoodPOT = 0;
  //
  int fEventNumber = 0;

};


hnlAna::AnalyzeEventsHNL::AnalyzeEventsHNL(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}, 
  // More initializers here.
  fNeutrinoRecoAngle(p, "pandoraTrack", "pandoraShower", "pandora", "wclsdatahd", "pandoraTrack", "pandoraShower", "pandora"),
  fNeutrinoRecoEnergy(p, "pandoraTrack", "pandoraShower", "pandora", "wclsdatahd", "pandoraTrack", "pandoraShower", "pandora"),
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

void hnlAna::AnalyzeEventsHNL::analyze(art::Event const& e)
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
  fPassCut = 0;
  fPassCut_after_fiducial = 0;
  fPassCut_after_directioncut = 0;
  fPassCut_after_trackcut = 0;

  fEventID_true = 0;
  fVx_true = 0;
  fVy_true = 0;
  fVz_true = 0;
  fE_true = 0;
  fPx_true = 0;
  fPy_true = 0;
  fPz_true = 0;
  fP_true = 0;
  fDirCosX_true = 0;
  fDirCosY_true = 0;
  fDirCosZ_true = 0;
  fPdgCode_true = 0;
  fMother_true = 0;
  fPOT_true = 0;
  fGoodPOT_true = 0;
  fTotalPOT_true = 0;
  fTA_true = 0;
  fEventNumber_true = 0;
  fKineticEnergy_true = 0; 



  fEventNumber++;
  

  // get the trigger information

  art::Handle<std::vector<dunedaq::trgdataformats::TriggerActivityData>> taHandle;
  if (!e.getByLabel(fTALabel, taHandle)) {
      // std::cout << "Warning: No TriggerActivityData found for label " << fTALabel << std::endl;
      fTA_true = 0;
  } else if (taHandle->empty()) {
      // std::cout << "Warning: TriggerActivityData retrieved but is empty." << std::endl;
      fTA_true = 0;
  } else {
      fTA_true = 1;
  }

  //   art::Handle<std::vector<dunedaq::trgdataformats::TriggerActivityData>> taHandle;
  //   if (!taHandle.isValid()) {
  //     std::cerr << "Failed to retrieve TriggerActivityData from event. "
  //                  << "This may be expected if no trigger primitives were found." 
  //                  << std::endl;
  //     fTA = -1;
  //   } 

  //  if (!e.getByLabel(fTALabel, taHandle)) {
  //     fTA = 0;
  //   } else {
  //     if (taHandle->size() == 0) {
  //       fTA = 0;
  //     } else {
  //       fTA = 1;
  //     }
  //   } 
    // double fE = 0;

  if (fGetTruth) {
    auto truthHandle = e.getValidHandle<std::vector<simb::MCTruth>>(fMCTruthLabel);
    for (auto const& truth : (*truthHandle)) {
      TLorentzVector tempVector;
      TLorentzVector totalVector;
      tempVector.SetXYZM(0, 0, 0, 0);
      totalVector.SetXYZM(0, 0, 0, 0);
      for (int i = 0; i < truth.NParticles(); i++) {
        // std::cout << "Number of particles: " << truth.NParticles() << std::endl;
        // std::cout << "Particle " << i << std::endl;
        const auto& particle = truth.GetParticle(i);
        fE_true = particle.E();
        fEventID_true = e.id().event();

        // Production vertex (where HNL is created) : Mother = -1 (HNL is not identified?)
        // fVx_prod_true = particle.Vx();
        // fVy_prod_true = particle.Vy();
        // fVz_prod_true = particle.Vz();
        // decay vertex (where HNL decays)
        fVx_true = particle.EndX();
        fVy_true = particle.EndY();
        fVz_true = particle.EndZ();
        fE_true = particle.E();
        fPx_true = particle.Px();
        fPy_true = particle.Py();
        fPz_true = particle.Pz();
        fPdgCode_true = particle.PdgCode();
        fMother_true = particle.Mother();
        // Kinetic Energy = E - mass
        double mass = particle.Mass();
        fKineticEnergy_true = fE_true - mass;

        fP_true = particle.P();
        fDirCosX_true = fPx_true / fP_true;
        fDirCosY_true = fPy_true / fP_true;
        fDirCosZ_true = fPz_true / fP_true;

        fTree_truth->Fill();
        tempVector.SetXYZM(particle.Px(), particle.Py(), particle.Pz(), particle.Mass());
        totalVector += tempVector;

        // std::cout << "Event number: " << fEventID_true << std::endl;
        // std::cout << "True vertex position: (" << fVx_true << ", " << fVy_true << ", " << fVz_true << ")" << std::endl;
        // std::cout << "True origin ID: " << fTrueOriginID_reco << std::endl;
        // std::cout << "True PDG code: " << fPdgCode_true << std::endl;
        // std::cout << "True energy: " << fE_true << std::endl;
        // std::cout << "True momentum: (" << fPx_true << ", " << fPy_true << ", " << fPz_true << ")" << std::endl;
        // std::cout << "True mother: " << fMother_true << std::endl;
        // std::cout << "Kinetic Energy: " << fKineticEnergy_true << " for pdg code: " << fPdgCode_true << " and mass: " << mass << "and Energy: " << fE_true << std::endl;

      }// end of loop over mc particles
    }// loop of truthHandle
  }// end of if fGetTruth




  // -------------------------------------------------------------------
  
  art::Handle<std::vector<recob::Slice>> sliceHandle = e.getHandle<std::vector<recob::Slice>>(fSliceLabel);

  if (sliceHandle.isValid()){
    std::vector<art::Ptr<recob::Slice>> slicePtrVector;
    art::fill_ptr_vector(slicePtrVector, sliceHandle);

    // get the most energetic slice
    double max_energy = 0;
    double energy = 0;
    art::Ptr<recob::Slice> most_energetic_slice;
    for (const art::Ptr<recob::Slice>& slicePtr : slicePtrVector) {
      // energy = GetTotalEnergy(slicePtr, e);
      // use the slice energy
      dune::EnergyRecoOutput this_energy_output = fNeutrinoRecoEnergy.CalculateNeutrinoEnergy(e, slicePtr, true);
      energy = this_energy_output.fNuLorentzVector.E();

      if (energy > max_energy) {
        max_energy = energy;
        most_energetic_slice = slicePtr;
      }
    }
    

    art::FindManyP<recob::PFParticle> slicePFPAssoc(sliceHandle, e, fSliceLabel);
    std::vector<art::Ptr<recob::PFParticle>> pfparticlePtrVector = slicePFPAssoc.at(most_energetic_slice.key());
    art::Handle<std::vector<recob::PFParticle>> pfparticleHandle = e.getHandle<std::vector<recob::PFParticle>>(fPFParticleLabel);
    


    std::vector<double> selectedvertexX;
    std::vector<double> selectedvertexY;
    std::vector<double> selectedvertexZ;
    std::vector<double> selectedEnergy;
    std::vector<double> selectedDirectionX;
    std::vector<double> selectedDirectionY;
    std::vector<double> selectedDirectionZ;
    std::vector<int> selectedTrueOriginID;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for (const art::Ptr<recob::PFParticle>& pfparticlePtr : pfparticlePtrVector) {

      art::FindManyP<recob::Vertex> pfVertexAssoc(pfparticleHandle, e, fVertexLabel);
      std::vector<art::Ptr<recob::Vertex>> vertices = pfVertexAssoc.at(pfparticlePtr.key());
      
      std::vector<double> info = GetDaugtherInfoDFS(pfparticlePtr, e);
      
      // True Origin
      int trueOriginID = -999;
      if (fGetTruth){
        trueOriginID = GetTrueInfo(pfparticlePtr, e);
      }

      // Vertex
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
      selectedTrueOriginID.push_back(trueOriginID);
      selectedvertexX.push_back(vertex_x);
      selectedvertexY.push_back(vertex_y);
      selectedvertexZ.push_back(vertex_z);
      dune::Point_t default_v_point;
      default_v_point.SetCoordinates(vertex_x, vertex_y, vertex_z);
      dune::AngularRecoOutput nu_angle = fNeutrinoRecoAngle.CalculateNeutrinoAngle(e, most_energetic_slice, default_v_point);
      selectedDirectionX.push_back(nu_angle.fRecoDirection.x()/(sqrt(nu_angle.fRecoDirection.x()*nu_angle.fRecoDirection.x() + nu_angle.fRecoDirection.y()*nu_angle.fRecoDirection.y() + nu_angle.fRecoDirection.z()*nu_angle.fRecoDirection.z())));
      selectedDirectionY.push_back(nu_angle.fRecoDirection.y()/(sqrt(nu_angle.fRecoDirection.x()*nu_angle.fRecoDirection.x() + nu_angle.fRecoDirection.y()*nu_angle.fRecoDirection.y() + nu_angle.fRecoDirection.z()*nu_angle.fRecoDirection.z())));
      selectedDirectionZ.push_back(nu_angle.fRecoDirection.z()/(sqrt(nu_angle.fRecoDirection.x()*nu_angle.fRecoDirection.x() + nu_angle.fRecoDirection.y()*nu_angle.fRecoDirection.y() + nu_angle.fRecoDirection.z()*nu_angle.fRecoDirection.z())));
      selectedEnergy.push_back(fNeutrinoRecoEnergy.CalculateNeutrinoEnergy(e, most_energetic_slice, true).fNuLorentzVector.E());
    }// end of loop over pfparticles
    fEventID_noCuts = e.id().event();
    fTrueOriginID_noCuts = selectedTrueOriginID;
    fVx_noCuts = selectedvertexX;
    fVy_noCuts = selectedvertexY;
    fVz_noCuts = selectedvertexZ;
    fEnergy_noCuts = selectedEnergy;
    fDirectionX_noCuts = selectedDirectionX;
    fDirectionY_noCuts = selectedDirectionY;
    fDirectionZ_noCuts = selectedDirectionZ;
    fTree_noCuts->Fill();
    selectedvertexX.clear();
    selectedvertexY.clear();
    selectedvertexZ.clear();
    selectedEnergy.clear();
    selectedDirectionX.clear();
    selectedDirectionY.clear();
    selectedDirectionZ.clear();
    selectedTrueOriginID.clear();
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for (const art::Ptr<recob::PFParticle>& pfparticlePtr : pfparticlePtrVector) {

      art::FindManyP<recob::Vertex> pfVertexAssoc(pfparticleHandle, e, fVertexLabel);
      std::vector<art::Ptr<recob::Vertex>> vertices = pfVertexAssoc.at(pfparticlePtr.key());
      
      std::vector<double> info = GetDaugtherInfoDFS(pfparticlePtr, e);
      
      // True Origin
      int trueOriginID = -999;
      if (fGetTruth){
        trueOriginID = GetTrueInfo(pfparticlePtr, e);
      }

      // Vertex
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


      if (!((vertex_y <= 550) && (vertex_z >= 20))) {
        continue;
      }


      selectedTrueOriginID.push_back(trueOriginID);
      selectedvertexX.push_back(vertex_x);
      selectedvertexY.push_back(vertex_y);
      selectedvertexZ.push_back(vertex_z);
      dune::Point_t default_v_point;
      default_v_point.SetCoordinates(vertex_x, vertex_y, vertex_z);
      dune::AngularRecoOutput nu_angle = fNeutrinoRecoAngle.CalculateNeutrinoAngle(e, most_energetic_slice, default_v_point);
      selectedDirectionX.push_back(nu_angle.fRecoDirection.x()/(sqrt(nu_angle.fRecoDirection.x()*nu_angle.fRecoDirection.x() + nu_angle.fRecoDirection.y()*nu_angle.fRecoDirection.y() + nu_angle.fRecoDirection.z()*nu_angle.fRecoDirection.z())));
      selectedDirectionY.push_back(nu_angle.fRecoDirection.y()/(sqrt(nu_angle.fRecoDirection.x()*nu_angle.fRecoDirection.x() + nu_angle.fRecoDirection.y()*nu_angle.fRecoDirection.y() + nu_angle.fRecoDirection.z()*nu_angle.fRecoDirection.z())));
      selectedDirectionZ.push_back(nu_angle.fRecoDirection.z()/(sqrt(nu_angle.fRecoDirection.x()*nu_angle.fRecoDirection.x() + nu_angle.fRecoDirection.y()*nu_angle.fRecoDirection.y() + nu_angle.fRecoDirection.z()*nu_angle.fRecoDirection.z())));
      selectedEnergy.push_back(fNeutrinoRecoEnergy.CalculateNeutrinoEnergy(e, most_energetic_slice, true).fNuLorentzVector.E());
    }// end of loop over pfparticles
    fEventID_fiducialPass = e.id().event();
    fTrueOriginID_fiducialPass = selectedTrueOriginID;
    fVx_fiducialPass = selectedvertexX;
    fVy_fiducialPass = selectedvertexY;
    fVz_fiducialPass = selectedvertexZ;
    fEnergy_fiducialPass = selectedEnergy;
    fDirectionX_fiducialPass = selectedDirectionX;
    fDirectionY_fiducialPass = selectedDirectionY;
    fDirectionZ_fiducialPass = selectedDirectionZ;
    fTree_fiducialPass->Fill();
    selectedvertexX.clear();
    selectedvertexY.clear();
    selectedvertexZ.clear();
    selectedEnergy.clear();
    selectedDirectionX.clear();
    selectedDirectionY.clear();
    selectedDirectionZ.clear();
    selectedTrueOriginID.clear();
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for (const art::Ptr<recob::PFParticle>& pfparticlePtr : pfparticlePtrVector) {

      art::FindManyP<recob::Vertex> pfVertexAssoc(pfparticleHandle, e, fVertexLabel);
      std::vector<art::Ptr<recob::Vertex>> vertices = pfVertexAssoc.at(pfparticlePtr.key());
      
      std::vector<double> info = GetDaugtherInfoDFS(pfparticlePtr, e);
      
      // True Origin
      int trueOriginID = -999;
      if (fGetTruth){
        trueOriginID = GetTrueInfo(pfparticlePtr, e);
      }

      // Vertex
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


      if (!((vertex_y <= 550) && (vertex_z >= 20))) {
        continue;
      }


      dune::Point_t default_v_point;
      default_v_point.SetCoordinates(vertex_x, vertex_y, vertex_z);
      dune::AngularRecoOutput nu_angle = fNeutrinoRecoAngle.CalculateNeutrinoAngle(e, most_energetic_slice, default_v_point);


      double tempDirectionZ = nu_angle.fRecoDirection.z()/(sqrt(nu_angle.fRecoDirection.x()*nu_angle.fRecoDirection.x() + nu_angle.fRecoDirection.y()*nu_angle.fRecoDirection.y() + nu_angle.fRecoDirection.z()*nu_angle.fRecoDirection.z()));
      if (tempDirectionZ < 0.7) {
        continue;
      }


      selectedTrueOriginID.push_back(trueOriginID);
      selectedvertexX.push_back(vertex_x);
      selectedvertexY.push_back(vertex_y);
      selectedvertexZ.push_back(vertex_z);
      selectedDirectionX.push_back(nu_angle.fRecoDirection.x()/(sqrt(nu_angle.fRecoDirection.x()*nu_angle.fRecoDirection.x() + nu_angle.fRecoDirection.y()*nu_angle.fRecoDirection.y() + nu_angle.fRecoDirection.z()*nu_angle.fRecoDirection.z())));
      selectedDirectionY.push_back(nu_angle.fRecoDirection.y()/(sqrt(nu_angle.fRecoDirection.x()*nu_angle.fRecoDirection.x() + nu_angle.fRecoDirection.y()*nu_angle.fRecoDirection.y() + nu_angle.fRecoDirection.z()*nu_angle.fRecoDirection.z())));
      selectedDirectionZ.push_back(nu_angle.fRecoDirection.z()/(sqrt(nu_angle.fRecoDirection.x()*nu_angle.fRecoDirection.x() + nu_angle.fRecoDirection.y()*nu_angle.fRecoDirection.y() + nu_angle.fRecoDirection.z()*nu_angle.fRecoDirection.z())));
      selectedEnergy.push_back(fNeutrinoRecoEnergy.CalculateNeutrinoEnergy(e, most_energetic_slice, true).fNuLorentzVector.E());
    }// end of loop over pfparticles
    fEventID_dirZPass = e.id().event();
    fTrueOriginID_dirZPass = selectedTrueOriginID;
    fVx_dirZPass = selectedvertexX;
    fVy_dirZPass = selectedvertexY;
    fVz_dirZPass = selectedvertexZ;
    fEnergy_dirZPass = selectedEnergy;
    fDirectionX_dirZPass = selectedDirectionX;
    fDirectionY_dirZPass = selectedDirectionY;
    fDirectionZ_dirZPass = selectedDirectionZ;
    fTree_dirZPass->Fill();
    selectedvertexX.clear();
    selectedvertexY.clear();
    selectedvertexZ.clear();
    selectedEnergy.clear();
    selectedDirectionX.clear();
    selectedDirectionY.clear();
    selectedDirectionZ.clear();
    selectedTrueOriginID.clear();
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for (const art::Ptr<recob::PFParticle>& pfparticlePtr : pfparticlePtrVector) {

      art::FindManyP<recob::Vertex> pfVertexAssoc(pfparticleHandle, e, fVertexLabel);
      std::vector<art::Ptr<recob::Vertex>> vertices = pfVertexAssoc.at(pfparticlePtr.key());
      
      std::vector<double> info = GetDaugtherInfoDFS(pfparticlePtr, e);
      
      // True Origin
      int trueOriginID = -999;
      if (fGetTruth){
        trueOriginID = GetTrueInfo(pfparticlePtr, e);
      }

      // Vertex
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


      if (!((vertex_y <= 550) && (vertex_z >= 20))) {
        continue;
      }


      dune::Point_t default_v_point;
      default_v_point.SetCoordinates(vertex_x, vertex_y, vertex_z);
      dune::AngularRecoOutput nu_angle = fNeutrinoRecoAngle.CalculateNeutrinoAngle(e, most_energetic_slice, default_v_point);


      // double tempDirectionZ = nu_angle.fRecoDirection.z()/(sqrt(nu_angle.fRecoDirection.x()*nu_angle.fRecoDirection.x() + nu_angle.fRecoDirection.y()*nu_angle.fRecoDirection.y() + nu_angle.fRecoDirection.z()*nu_angle.fRecoDirection.z()));
      // if (tempDirectionZ < 0.7) {
      //   continue;
      // }

      std::vector<art::Ptr<recob::PFParticle>> daughters = dune_ana::DUNEAnaPFParticleUtils::GetChildParticles(pfparticlePtr, e, fPFParticleLabel);
      if(daughters.size() != 0){
        continue;
      }

      if(!dune_ana::DUNEAnaPFParticleUtils::IsTrack(pfparticlePtr, e, fPFParticleLabel, fTrackLabel)){
        continue;
      }


      selectedTrueOriginID.push_back(trueOriginID);
      selectedvertexX.push_back(vertex_x);
      selectedvertexY.push_back(vertex_y);
      selectedvertexZ.push_back(vertex_z);
      selectedDirectionX.push_back(nu_angle.fRecoDirection.x()/(sqrt(nu_angle.fRecoDirection.x()*nu_angle.fRecoDirection.x() + nu_angle.fRecoDirection.y()*nu_angle.fRecoDirection.y() + nu_angle.fRecoDirection.z()*nu_angle.fRecoDirection.z())));
      selectedDirectionY.push_back(nu_angle.fRecoDirection.y()/(sqrt(nu_angle.fRecoDirection.x()*nu_angle.fRecoDirection.x() + nu_angle.fRecoDirection.y()*nu_angle.fRecoDirection.y() + nu_angle.fRecoDirection.z()*nu_angle.fRecoDirection.z())));
      selectedDirectionZ.push_back(nu_angle.fRecoDirection.z()/(sqrt(nu_angle.fRecoDirection.x()*nu_angle.fRecoDirection.x() + nu_angle.fRecoDirection.y()*nu_angle.fRecoDirection.y() + nu_angle.fRecoDirection.z()*nu_angle.fRecoDirection.z())));
      selectedEnergy.push_back(fNeutrinoRecoEnergy.CalculateNeutrinoEnergy(e, most_energetic_slice, true).fNuLorentzVector.E());
    }// end of loop over pfparticles
    fEventID_isTrackPass = e.id().event();
    fTrueOriginID_isTrackPass = selectedTrueOriginID;
    fVx_isTrackPass = selectedvertexX;
    fVy_isTrackPass = selectedvertexY;
    fVz_isTrackPass = selectedvertexZ;
    fEnergy_isTrackPass = selectedEnergy;
    fDirectionX_isTrackPass = selectedDirectionX;
    fDirectionY_isTrackPass = selectedDirectionY;
    fDirectionZ_isTrackPass = selectedDirectionZ;
    fTree_isTrackPass->Fill();
    selectedvertexX.clear();
    selectedvertexY.clear();
    selectedvertexZ.clear();
    selectedEnergy.clear();
    selectedDirectionX.clear();
    selectedDirectionY.clear();
    selectedDirectionZ.clear();
    selectedTrueOriginID.clear();


    std::vector<int> selectedIsTrack;
    std::vector<int> selectedIsShower;

    for (const art::Ptr<recob::PFParticle>& pfparticlePtr : pfparticlePtrVector) {

      art::FindManyP<recob::Vertex> pfVertexAssoc(pfparticleHandle, e, fVertexLabel);
      std::vector<art::Ptr<recob::Vertex>> vertices = pfVertexAssoc.at(pfparticlePtr.key());
      
      std::vector<double> info = GetDaugtherInfoDFS(pfparticlePtr, e);
      
      // True Origin
      int trueOriginID = -999;
      if (fGetTruth){
        trueOriginID = GetTrueInfo(pfparticlePtr, e);
      }

      // Vertex
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


      if (!((vertex_y <= 550) && (vertex_z >= 20))) {
        continue;
      }


      dune::Point_t default_v_point;
      default_v_point.SetCoordinates(vertex_x, vertex_y, vertex_z);
      dune::AngularRecoOutput nu_angle = fNeutrinoRecoAngle.CalculateNeutrinoAngle(e, most_energetic_slice, default_v_point);


      // double tempDirectionZ = nu_angle.fRecoDirection.z()/(sqrt(nu_angle.fRecoDirection.x()*nu_angle.fRecoDirection.x() + nu_angle.fRecoDirection.y()*nu_angle.fRecoDirection.y() + nu_angle.fRecoDirection.z()*nu_angle.fRecoDirection.z()));
      // if (tempDirectionZ < 0.7) {
      //   continue;
      // }

      std::vector<art::Ptr<recob::PFParticle>> daughters = dune_ana::DUNEAnaPFParticleUtils::GetChildParticles(pfparticlePtr, e, fPFParticleLabel);
      if(daughters.size() != 0){
        continue;
      }

      if(dune_ana::DUNEAnaPFParticleUtils::IsTrack(pfparticlePtr, e, fPFParticleLabel, fTrackLabel)){
        selectedIsTrack.push_back(1);
      }

      if(dune_ana::DUNEAnaPFParticleUtils::IsShower(pfparticlePtr, e, fPFParticleLabel, fShowerLabel)){
        selectedIsShower.push_back(1);
      }
      
      selectedTrueOriginID.push_back(trueOriginID);
      selectedvertexX.push_back(vertex_x);
      selectedvertexY.push_back(vertex_y);
      selectedvertexZ.push_back(vertex_z);
      selectedDirectionX.push_back(nu_angle.fRecoDirection.x()/(sqrt(nu_angle.fRecoDirection.x()*nu_angle.fRecoDirection.x() + nu_angle.fRecoDirection.y()*nu_angle.fRecoDirection.y() + nu_angle.fRecoDirection.z()*nu_angle.fRecoDirection.z())));
      selectedDirectionY.push_back(nu_angle.fRecoDirection.y()/(sqrt(nu_angle.fRecoDirection.x()*nu_angle.fRecoDirection.x() + nu_angle.fRecoDirection.y()*nu_angle.fRecoDirection.y() + nu_angle.fRecoDirection.z()*nu_angle.fRecoDirection.z())));
      selectedDirectionZ.push_back(nu_angle.fRecoDirection.z()/(sqrt(nu_angle.fRecoDirection.x()*nu_angle.fRecoDirection.x() + nu_angle.fRecoDirection.y()*nu_angle.fRecoDirection.y() + nu_angle.fRecoDirection.z()*nu_angle.fRecoDirection.z())));
      selectedEnergy.push_back(fNeutrinoRecoEnergy.CalculateNeutrinoEnergy(e, most_energetic_slice, true).fNuLorentzVector.E());
    }// end of loop over pfparticles
    fEventID_isShowerPass = e.id().event();
    fTrueOriginID_isShowerPass = selectedTrueOriginID;
    fVx_isShowerPass = selectedvertexX;
    fVy_isShowerPass = selectedvertexY;
    fVz_isShowerPass = selectedvertexZ;
    fEnergy_isShowerPass = selectedEnergy;
    fDirectionX_isShowerPass = selectedDirectionX;
    fDirectionY_isShowerPass = selectedDirectionY;
    fDirectionZ_isShowerPass = selectedDirectionZ;
    fIsTrack = selectedIsTrack;
    fIsShower = selectedIsShower;

    if((fIsTrack.size()==fIsShower.size())){
      fTree_isShowerPass->Fill();
    }
    selectedvertexX.clear();
    selectedvertexY.clear();
    selectedvertexZ.clear();
    selectedEnergy.clear();
    selectedDirectionX.clear();
    selectedDirectionY.clear();
    selectedDirectionZ.clear();
    selectedTrueOriginID.clear();
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // uint32_t timeHigh_ns = e.time().timeHigh();
    // uint32_t timeLow_ns = e.time().timeLow();
    // fTime_aggregate = timeHigh_ns*1e9 + timeLow_ns;

  }// end of if sliceHandle is valid
}// end of analyze function
//------------------------------------------------------------------------------------------------------------------




//--------------------------------------------------------------------------------------------------------------------------
double hnlAna::AnalyzeEventsHNL::GetTotalEnergy(const art::Ptr<recob::Slice>& slicePtr, art::Event const& e){
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





int hnlAna::AnalyzeEventsHNL::GetTrueInfo(const art::Ptr<recob::PFParticle>& pfparticlePtr, art::Event const& e) {
  // needs to be rethinked. What's the point of having the total number of elements with the same origin?
  TruthMatchUtils::G4ID trueID = 0;
  std::vector<TruthMatchUtils::G4ID> trueOriginIDs_vector;
  // get the daughter pfparticles
  // std::cout << "Getting true info" << std::endl;
  trueOriginIDs_vector = GetTrueInfoChainDFS(pfparticlePtr, e);
  // std::cout << "Got true info" << std::endl;

  for (TruthMatchUtils::G4ID trueOriginID : trueOriginIDs_vector) {
    if (trueOriginID == 1) {
      trueID++;
    }
  }
  return trueID;
}

std::vector<TruthMatchUtils::G4ID> hnlAna::AnalyzeEventsHNL::GetTrueInfoChainDFS(const art::Ptr<recob::PFParticle>& pfparticlePtr, art::Event const& e) {
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

std::vector<double> hnlAna::AnalyzeEventsHNL::GetDaugtherInfoDFS(const art::Ptr<recob::PFParticle>& pfparticlePtr, art::Event const& e) {
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
      // get dedx which is vector float and print the sum with units
      // std::vector<float> dEdx = this_calo->dEdx();
      // for (size_t i = 0; i < dEdx.size(); i++) {
      //   std::cout << "dEdx[" << i << "] = " << dEdx[i] << " MeV/cm" << std::endl;
      // }

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

void hnlAna::AnalyzeEventsHNL::beginSubRun(art::SubRun const& subRun) {
  
  if (fGetTruth){

    const auto potSummaryHandle = subRun.getValidHandle<sumdata::POTSummary>("generator");
    const auto &potSummary = *potSummaryHandle;
    fPOT = potSummary.totpot;
    fGoodPOT = potSummary.totgoodpot;
    fTotalPOT += fPOT;
  
  }

}

void hnlAna::AnalyzeEventsHNL::beginJob()
{
  // Implementation of optional member function here.
  // Get the TFileService to create the output TTree for us
  art::ServiceHandle<art::TFileService> tfs;
  fTree_reco = tfs->make<TTree>("tree_reco", "Output TTree reco");
  fTree_truth = tfs->make<TTree>("tree_truth", "Output TTree truth");
  fTree_noCuts = tfs->make<TTree>("tree_noCuts", "Output TTree noCuts");
  fTree_fiducialPass = tfs->make<TTree>("tree_fiducialPass", "Output TTree fiducialPass");
  fTree_dirZPass = tfs->make<TTree>("tree_dirZPass", "Output TTree dirZPass");
  fTree_isTrackPass = tfs->make<TTree>("tree_isTrackPass", "Output TTree isTrackPass");
  fTree_isShowerPass = tfs->make<TTree>("tree_isShowerPass", "Output TTree isShowerPass");


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
  fTree_reco->Branch("eventNumber", &fEventNumber_reco);
  fTree_reco->Branch("passCut", &fPassCut);
  fTree_reco->Branch("passCut_after_fiducial", &fPassCut_after_fiducial);
  fTree_reco->Branch("passCut_after_directioncut", &fPassCut_after_directioncut);
  fTree_reco->Branch("passCut_after_trackcut", &fPassCut_after_trackcut);

  fTree_truth->Branch("eventID", &fEventID_true);
  fTree_truth->Branch("vx", &fVx_true);
  fTree_truth->Branch("vy", &fVy_true);
  fTree_truth->Branch("vz", &fVz_true);
  fTree_truth->Branch("E", &fE_true);
  fTree_truth->Branch("Px", &fPx_true);
  fTree_truth->Branch("Py", &fPy_true);
  fTree_truth->Branch("Pz", &fPz_true);
  fTree_truth->Branch("P", &fP_true);
  fTree_truth->Branch("directionX", &fDirCosX_true);
  fTree_truth->Branch("directionY", &fDirCosY_true);
  fTree_truth->Branch("directionZ", &fDirCosZ_true);
  fTree_truth->Branch("pdgCode", &fPdgCode_true);
  fTree_truth->Branch("mother", &fMother_true);
  fTree_truth->Branch("POT", &fPOT_true);
  fTree_truth->Branch("goodPOT", &fGoodPOT_true);
  fTree_truth->Branch("totalPOT", &fTotalPOT_true);
  fTree_truth->Branch("TA", &fTA_true);
  fTree_truth->Branch("eventNumber", &fEventNumber_true);


  fTree_noCuts->Branch("eventID", &fEventID_noCuts);
  fTree_noCuts->Branch("vx", &fVx_noCuts);
  fTree_noCuts->Branch("vy", &fVy_noCuts);
  fTree_noCuts->Branch("vz", &fVz_noCuts);
  fTree_noCuts->Branch("energy", &fEnergy_noCuts);
  fTree_noCuts->Branch("directionX", &fDirectionX_noCuts);
  fTree_noCuts->Branch("directionY", &fDirectionY_noCuts);
  fTree_noCuts->Branch("directionZ", &fDirectionZ_noCuts);
  fTree_noCuts->Branch("trueOriginID", &fTrueOriginID_noCuts);

  fTree_fiducialPass->Branch("eventID", &fEventID_fiducialPass);
  fTree_fiducialPass->Branch("vx", &fVx_fiducialPass);
  fTree_fiducialPass->Branch("vy", &fVy_fiducialPass);
  fTree_fiducialPass->Branch("vz", &fVz_fiducialPass);
  fTree_fiducialPass->Branch("energy", &fEnergy_fiducialPass);
  fTree_fiducialPass->Branch("directionX", &fDirectionX_fiducialPass);
  fTree_fiducialPass->Branch("directionY", &fDirectionY_fiducialPass);
  fTree_fiducialPass->Branch("directionZ", &fDirectionZ_fiducialPass);
  fTree_fiducialPass->Branch("trueOriginID", &fTrueOriginID_fiducialPass);

  fTree_dirZPass->Branch("eventID", &fEventID_dirZPass);
  fTree_dirZPass->Branch("vx", &fVx_dirZPass);
  fTree_dirZPass->Branch("vy", &fVy_dirZPass);
  fTree_dirZPass->Branch("vz", &fVz_dirZPass);
  fTree_dirZPass->Branch("energy", &fEnergy_dirZPass);
  fTree_dirZPass->Branch("directionX", &fDirectionX_dirZPass);
  fTree_dirZPass->Branch("directionY", &fDirectionY_dirZPass);
  fTree_dirZPass->Branch("directionZ", &fDirectionZ_dirZPass);
  fTree_dirZPass->Branch("trueOriginID", &fTrueOriginID_dirZPass);


  fTree_isTrackPass->Branch("eventID", &fEventID_isTrackPass);
  fTree_isTrackPass->Branch("vx", &fVx_isTrackPass);
  fTree_isTrackPass->Branch("vy", &fVy_isTrackPass);
  fTree_isTrackPass->Branch("vz", &fVz_isTrackPass);
  fTree_isTrackPass->Branch("energy", &fEnergy_isTrackPass);
  fTree_isTrackPass->Branch("directionX", &fDirectionX_isTrackPass);
  fTree_isTrackPass->Branch("directionY", &fDirectionY_isTrackPass);
  fTree_isTrackPass->Branch("directionZ", &fDirectionZ_isTrackPass);
  fTree_isTrackPass->Branch("trueOriginID", &fTrueOriginID_isTrackPass);

  fTree_isShowerPass->Branch("eventID", &fEventID_isShowerPass);
  fTree_isShowerPass->Branch("vx", &fVx_isShowerPass);
  fTree_isShowerPass->Branch("vy", &fVy_isShowerPass);
  fTree_isShowerPass->Branch("vz", &fVz_isShowerPass);
  fTree_isShowerPass->Branch("energy", &fEnergy_isShowerPass);
  fTree_isShowerPass->Branch("directionX", &fDirectionX_isShowerPass);
  fTree_isShowerPass->Branch("directionY", &fDirectionY_isShowerPass);
  fTree_isShowerPass->Branch("directionZ", &fDirectionZ_isShowerPass);
  fTree_isShowerPass->Branch("trueOriginID", &fTrueOriginID_isShowerPass);
  fTree_isShowerPass->Branch("isTrack", &fIsTrack);
  fTree_isShowerPass->Branch("isShower", &fIsShower);







}

void hnlAna::AnalyzeEventsHNL::endJob()
{
  // // Implementation of optional member function here.
  // std::cout << "==== Cut Flow Summary ====" << std::endl;
  // std::cout << "Total events: " << total_events << std::endl;
  // std::cout << "After fiducial cut: " << pass_fiducial << " (" << (100.0 * pass_fiducial / total_events) << "%)" << std::endl;
  // std::cout << "After direction cut: " << pass_direction << " (" << (100.0 * pass_direction / total_events) << "%)" << std::endl;
  // // std::cout << "After energy cut: " << pass_energy << " (" << (100.0 * pass_energy / total_events) << "%)" << std::endl;
  // std::cout << "After track multiplicity cut: " << pass_track_multiplicity << " (" << (100.0 * pass_track_multiplicity / total_events) << "%)" << std::endl;
  // std::cout << "After all cuts: " << pass_all_cuts << " (" << (100.0 * pass_all_cuts / total_events) << "%)" << std::endl;
  // std::cout << "=========================" << std::endl;
}

DEFINE_ART_MODULE(hnlAna::AnalyzeEventsHNL)