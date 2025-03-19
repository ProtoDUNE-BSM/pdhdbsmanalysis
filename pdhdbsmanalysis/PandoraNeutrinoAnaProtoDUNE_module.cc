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
// #include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Slice.h"
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
#include "canvas/Persistency/Common/FindManyP.h"

#include "detdataformats/trigger/TriggerObjectOverlay.hpp"
#include "detdataformats/trigger/TriggerPrimitive.hpp"
#include "detdataformats/trigger/TriggerActivityData.hpp"
#include "detdataformats/trigger/TriggerCandidateData.hpp"

#include "dunereco/FDSensOpt/NeutrinoAngularRecoAlg/NeutrinoAngularRecoAlg.h"
#include "dunereco/FDSensOpt/NeutrinoEnergyRecoAlg/NeutrinoEnergyRecoAlg.h"
#include "dunereco/AnaUtils/DUNEAnaEventUtils.h"
#include "dunereco/AnaUtils/DUNEAnaPFParticleUtils.h"
#include "dunereco/AnaUtils/DUNEAnaTrackUtils.h"
#include "dunereco/AnaUtils/DUNEAnaShowerUtils.h"

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
  int fDecaynum;
 
  unsigned int fTA;

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

  // Geomtry Information
  geo::GeometryCore const* fGeometryService; ///< pointer to Geometry provider
  std::vector<double> fFiducialBoundaries;

  // POT Summary data
  double fPOT;
  double fGoodPOT;

  // Pandora PFParticle information
  std::vector<int> fnuSliceKey;
  std::vector<int> fnuID;
  std::vector<int> fRecoPDG;
  std::vector<int> fnuScore;
  std::vector<int> fSliceIndex;
  unsigned int fNPFParticles;
  unsigned int fNPrimaryChildren;
  unsigned int fNNeutrinos;

  std::vector<double> fRecoVertexX;
  std::vector<double> fRecoVertexY;
  std::vector<double> fRecoVertexZ;
 
  dune::NeutrinoAngularRecoAlg fNeutrinoRecoAngle;
  dune::NeutrinoEnergyRecoAlg fNeutrinoRecoEnergy;

  std::vector<double> fNuAngleRecoX;
  std::vector<double> fNuAngleRecoY;
  std::vector<double> fNuAngleRecoZ;

  std::vector<double> fNuEnergy;

  std::vector<std::vector<unsigned int>> vChildIsTrack;
  std::vector<std::vector<unsigned int>> vChildIsShower;
  std::vector<std::vector<double>> vTrackDirectionX;
  std::vector<std::vector<double>> vTrackDirectionY;
  std::vector<std::vector<double>> vTrackDirectionZ;
  std::vector<std::vector<double>> vShowerDirectionX;
  std::vector<std::vector<double>> vShowerDirectionY;
  std::vector<std::vector<double>> vShowerDirectionZ;
  std::vector<std::vector<double>> vKineticEnergyTrack;
  std::vector<std::vector<double>> vShowerEnergy;

};


//-----------------------------------------------
// Analyser class constructor
ana::PandoraNeutrinoAnaProtoDUNE::PandoraNeutrinoAnaProtoDUNE(fhicl::ParameterSet const& p)
  : EDAnalyzer{p} 
  , fMCTruthLabel(p.get<std::string>("MCTruthLabel"))
  , fPFParticleLabel(p.get<std::string>("PFParticleLabel"))
  , fSetPOT(p.get<double>("SetPOT"))
  , fDecaynum(p.get<int>("decay"))
  , fNeutrinoRecoAngle(p, "pandoraTrack", "pandoraShower", "pandora",
      "wclsdatahd", "pandoraTrack", "pandoraShower", "pandora")
  , fNeutrinoRecoEnergy(p, "pandoraTrack", "pandoraShower", "pandora",
      "wclsdatahd", "pandoraTrack", "pandoraShower", "pandora")
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

  fInFV = false;
  
  art::Handle<std::vector<dunedaq::trgdataformats::TriggerActivityData>> taHandle;
  if (!e.getByLabel("tamakerTPC", taHandle)) {
      fTA = 0;
  } else {
    if (taHandle->size() == 0) {
      fTA = 0;
    } else {
      std::cout << ">>> Found " << taHandle->size() << " TAs in Event " << fEventID << std::endl;
      fTA = 1;
    }
  }

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

  art::ValidHandle<std::vector<recob::Slice>> sliceHandle = e.getValidHandle<std::vector<recob::Slice>>("pandora");
  std::vector<art::Ptr<recob::Slice>> sliceVector;

  if (sliceHandle.isValid()) art::fill_ptr_vector(sliceVector, sliceHandle);

  art::FindManyP<recob::PFParticle> slicePFPAssoc(sliceHandle, e, "pandora");

  fnuSliceKey.clear();
  fnuID.clear();
  fRecoPDG.clear();
  fNPFParticles = 0;
  fNNeutrinos = 0;
  fNPrimaryChildren = 0;
  fRecoVertexX.clear();
  fRecoVertexY.clear();
  fRecoVertexZ.clear();
  fNuAngleRecoX.clear();
  fNuAngleRecoY.clear();
  fNuAngleRecoZ.clear();
  fNuEnergy.clear();
  vChildIsTrack.clear();
  vChildIsShower.clear();
  vTrackDirectionX.clear();
  vTrackDirectionY.clear();
  vTrackDirectionZ.clear();
  vShowerDirectionX.clear();
  vShowerDirectionY.clear();
  vShowerDirectionZ.clear();

  for (const art::Ptr<recob::Slice> &slice : sliceVector) {
    std::vector<art::Ptr<recob::PFParticle>> slicePFPs(slicePFPAssoc.at(slice.key()));
    for (const art::Ptr<recob::PFParticle> &slicePFP : slicePFPs) {
      bool isNeutrino = dune_ana::DUNEAnaPFParticleUtils::IsNeutrino(slicePFP);
      bool isPrimary = slicePFP->IsPrimary();

      if (!(isNeutrino && isPrimary)) continue;

      fNNeutrinos++;
      fnuSliceKey.push_back(slice.key());
      fnuID.push_back(slicePFP->Self());
      fNPFParticles = slicePFPs.size();
      fNPrimaryChildren = slicePFP->NumDaughters();
      fRecoPDG.push_back(slicePFP->PdgCode());

      art::Ptr<larpandoraobj::PFParticleMetadata> pandoraMetaData = dune_ana::DUNEAnaPFParticleUtils::GetMetadata(slicePFP, e, "pandora");
      std::map<std::string, float> fPFPPropertiesMap = pandoraMetaData->GetPropertiesMap();
      fnuScore.push_back(fPFPPropertiesMap["NuScore"]);
      fSliceIndex.push_back(fPFPPropertiesMap["SliceIndex"]);

      art::Ptr<recob::Vertex> nu_vertex = dune_ana::DUNEAnaPFParticleUtils::GetVertex(slicePFP, e, "pandora");
      auto v_point = nu_vertex->position();
      dune::Point_t default_v_point;
      default_v_point.SetCoordinates(v_point.x(), v_point.y(), v_point.z());

      fRecoVertexX.push_back(default_v_point.x());
      fRecoVertexY.push_back(default_v_point.y());
      fRecoVertexZ.push_back(default_v_point.z());

      // Get direction of overall shower from atmospheric reco code
      dune::AngularRecoOutput nu_angle = fNeutrinoRecoAngle.CalculateNeutrinoAngle(e, slice, default_v_point);
      fNuAngleRecoX.push_back(nu_angle.fRecoDirection.x());
      fNuAngleRecoY.push_back(nu_angle.fRecoDirection.y());
      fNuAngleRecoZ.push_back(nu_angle.fRecoDirection.z());

      dune::EnergyRecoOutput energy_output = fNeutrinoRecoEnergy.CalculateNeutrinoEnergy(e, slice, true);

      fNuEnergy.push_back(energy_output.fNuLorentzVector.E());

      std::vector<art::Ptr<recob::PFParticle>> vNuChildren = dune_ana::DUNEAnaPFParticleUtils::GetChildParticles(slicePFP, e, "pandora");
  
      std::vector<unsigned int> vNuChildIsTrack;
      std::vector<unsigned int> vNuChildIsShower;
      std::vector<double> vNuTrackDirectionX;
      std::vector<double> vNuTrackDirectionY;
      std::vector<double> vNuTrackDirectionZ;
      std::vector<double> vNuShowerDirectionX;
      std::vector<double> vNuShowerDirectionY;
      std::vector<double> vNuShowerDirectionZ;
      std::vector<double> vNuKineticEnergyTrack;
      std::vector<double> vNuShowerEnergy;

      for (const art::Ptr<recob::PFParticle> &childPFP : vNuChildren) {
        if (dune_ana::DUNEAnaPFParticleUtils::IsTrack(childPFP, e, "pandora", "pandoraTrack")) {
          vNuChildIsTrack.push_back(1);
          art::Ptr<recob::Track> track = dune_ana::DUNEAnaPFParticleUtils::GetTrack(childPFP, e, "pandora", "pandoraTrack");
          vNuTrackDirectionX.push_back(track->StartDirection().x());
          vNuTrackDirectionY.push_back(track->StartDirection().y());
          vNuTrackDirectionZ.push_back(track->StartDirection().z());

          art::Ptr<anab::Calorimetry> track_cal = dune_ana::DUNEAnaTrackUtils::GetCalorimetry(track, e, "pandoraTrack", "pandoracalonosce");
          vNuKineticEnergyTrack.push_back(track_cal->KineticEnergy());
        } else {
          vNuChildIsTrack.push_back(0);
        }
        if (dune_ana::DUNEAnaPFParticleUtils::IsShower(childPFP, e, "pandora", "pandoraShower")) {
          vNuChildIsShower.push_back(1);
          art::Ptr<recob::Shower> shower = dune_ana::DUNEAnaPFParticleUtils::GetShower(childPFP, e, "pandora", "pandoraShower");
          vNuShowerDirectionX.push_back(shower->Direction().x());
          vNuShowerDirectionY.push_back(shower->Direction().y());
          vNuShowerDirectionZ.push_back(shower->Direction().z());
          if (!shower->Energy().empty()) {
            vNuShowerEnergy.push_back(shower->Energy().at(0));
          } else {
            vNuShowerEnergy.push_back(0);
          }
        } else {
          vNuChildIsShower.push_back(0);
        }

      }
      // Fill child information into neutrino vectors
      vChildIsTrack.push_back(vNuChildIsTrack);
      vChildIsShower.push_back(vNuChildIsShower);
      vTrackDirectionX.push_back(vNuTrackDirectionX);
      vTrackDirectionY.push_back(vNuTrackDirectionY);
      vTrackDirectionZ.push_back(vNuTrackDirectionZ);
      vShowerDirectionX.push_back(vNuShowerDirectionX);
      vShowerDirectionY.push_back(vNuShowerDirectionY);
      vShowerDirectionZ.push_back(vNuShowerDirectionZ);
      vKineticEnergyTrack.push_back(vNuKineticEnergyTrack);
      vShowerEnergy.push_back(vNuShowerEnergy);

      // Try connecting to MCTruth
      //std::vector<art::Ptr<recob::Hit>> recoHits = dune_ana::DUNEAnaPFParticleUtils::GetHits(slicePFP, e, "pandora");
    }

  }

  fRecoNtuple->Fill();

}

//-----------------------------------------------
// Define outputs at start of the job
void ana::PandoraNeutrinoAnaProtoDUNE::beginJob() {
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;

  // Get TFileService to create an output tree
  fSimulationNtuple = tfs->make<TTree>("GenieTruth", "GENIE Output Tree");

  // Add branches to TTree
  fSimulationNtuple->Branch("eventID", &fEventID);
  fSimulationNtuple->Branch("Decaynum", &fDecaynum);
  fSimulationNtuple->Branch("TA", &fTA);
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
  fRecoNtuple->Branch("Decaynum", &fDecaynum);
  fRecoNtuple->Branch("TA", &fTA);
  fRecoNtuple->Branch("SubRun", &fSubRun, "SubRun/I");
  fRecoNtuple->Branch("Run", &fRun, "Run/I");
  fRecoNtuple->Branch("nPFParticles", &fNPFParticles);
  fRecoNtuple->Branch("nNeutrinos", &fNNeutrinos);
  fRecoNtuple->Branch("nPrimaryChildren", &fNPrimaryChildren);
  fRecoNtuple->Branch("nuSliceKey", &fnuSliceKey);
  fRecoNtuple->Branch("SliceIndex", &fSliceIndex);
  fRecoNtuple->Branch("nuID", &fnuID);
  fRecoNtuple->Branch("nuScore", &fnuScore);
  fRecoNtuple->Branch("RecoPDG", &fRecoPDG);
  fRecoNtuple->Branch("RecoVertexX", &fRecoVertexX);
  fRecoNtuple->Branch("RecoVertexY", &fRecoVertexY);
  fRecoNtuple->Branch("RecoVertexZ", &fRecoVertexZ);
  fRecoNtuple->Branch("NuAngleRecoX", &fNuAngleRecoX);
  fRecoNtuple->Branch("NuAngleRecoY", &fNuAngleRecoY);
  fRecoNtuple->Branch("NuAngleRecoZ", &fNuAngleRecoZ);
  fRecoNtuple->Branch("NuEnergy", &fNuEnergy);
  fRecoNtuple->Branch("vChildIsTrack", &vChildIsTrack);
  fRecoNtuple->Branch("vChildIsShower", &vChildIsShower);
  fRecoNtuple->Branch("vTrackDirectionX", &vTrackDirectionX);
  fRecoNtuple->Branch("vTrackDirectionY", &vTrackDirectionY);
  fRecoNtuple->Branch("vTrackDirectionZ", &vTrackDirectionZ);
  fRecoNtuple->Branch("vShowerDirectionX", &vShowerDirectionX);
  fRecoNtuple->Branch("vShowerDirectionY", &vShowerDirectionY);
  fRecoNtuple->Branch("vShowerDirectionZ", &vShowerDirectionZ);
  fRecoNtuple->Branch("KineticEnergyTrack", &vKineticEnergyTrack);
  fRecoNtuple->Branch("ShowerEnergy", &vShowerEnergy);

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
