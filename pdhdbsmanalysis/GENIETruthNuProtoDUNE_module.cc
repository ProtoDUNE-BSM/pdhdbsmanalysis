////////////////////////////////////////////////////////////////////////
// Class:       GENIETruthNuProtoDUNE
// Plugin Type: analyzer (Unknown Unknown)
// File:        GENIETruthNuProtoDUNE_module.cc
//
// Generated at Mon Apr  8 02:24:48 2024 by Ciaran Hasnip using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "detdataformats/trigger/TriggerObjectOverlay.hpp"
#include "detdataformats/trigger/TriggerPrimitive.hpp"
#include "detdataformats/trigger/TriggerActivityData.hpp"
#include "detdataformats/trigger/TriggerCandidateData.hpp"

// Additional Framework includes
#include "art_root_io/TFileService.h"

// ROOT includes
#include <TTree.h>
#include <TH1.h>

#include <string>

namespace ana {
  class GENIETruthNuProtoDUNE;
}

using timestamp_t = dunedaq::trgdataformats::timestamp_t;
using channel_t = dunedaq::trgdataformats::channel_t;
using triggerprimitive_t = dunedaq::trgdataformats::TriggerPrimitive;

// Define analyser class
class ana::GENIETruthNuProtoDUNE : public art::EDAnalyzer {
public:
  explicit GENIETruthNuProtoDUNE(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  GENIETruthNuProtoDUNE(GENIETruthNuProtoDUNE const&) = delete;
  GENIETruthNuProtoDUNE(GENIETruthNuProtoDUNE&&) = delete;
  GENIETruthNuProtoDUNE& operator=(GENIETruthNuProtoDUNE const&) = delete;
  GENIETruthNuProtoDUNE& operator=(GENIETruthNuProtoDUNE&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  void beginSubRun(art::SubRun const& subRun) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:
    
  std::pair<channel_t, channel_t> pCollectionAPA1IDs;
  std::pair<channel_t, channel_t> pCollectionAPA2IDs;
  std::pair<channel_t, channel_t> pCollectionAPA3IDs;
  std::pair<channel_t, channel_t> pCollectionAPA4IDs;

  TTree *fSimulationNtuple;
  TTree *fSubrunTree;

  // Declare member data here.
  unsigned int fEventID;
  unsigned int fRun;
  unsigned int fSubRun;

  art::InputTag fMCTruthLabel; ///< The name of the producer that tracked
                                          ///< simulated particles through the detector
  
  TH1D *hMCNeutrinoEnergy;

  double fSetPOT;

  int fSimPDG;     ///< PDG ID of the particle being processed
  int fSimTrackID; ///< GEANT ID of the particle being processed
  int fCCNC; ///< Is neutrino interaction a CC or NC interaction
  unsigned int fTPCID; ///< TPC ID where neutrino interacts

  double fE;
  double fnuStartX;
  double fnuStartY;
  double fnuStartZ;
  double fnuVertexX;
  double fnuVertexY;
  double fnuVertexZ;

  bool fInFV;
  
  geo::GeometryCore const* fGeometryService; ///< pointer to Geometry provider
  std::vector<double> fFiducialBoundaries;
  double fZedge;

  double fDecaynum;
  double fPOT;
  double fGoodPOT;
  double fTotalPOT;

  bool fTA;
  int fnTAs;
  int fAPA_id;
  std::vector<int> fAPA_ids;
  std::vector<double> fTPTAADCIntSum;

};


// Analyser class constructor
ana::GENIETruthNuProtoDUNE::GENIETruthNuProtoDUNE(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
  , fMCTruthLabel(p.get<std::string>("MCTruthLabel"))
  , fSetPOT(p.get<double>("SetPOT"))
  , fDecaynum(p.get<int>("decay"))
  // More initializers here.
{
  
  fAPA_id = 0;
  pCollectionAPA1IDs = std::make_pair(2080, 2559);
  pCollectionAPA2IDs = std::make_pair(7200, 7680);  
  pCollectionAPA3IDs = std::make_pair(4160, 4639);
  pCollectionAPA4IDs = std::make_pair(9280, 9759); 
 
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

  fTotalPOT = 0;
  fZedge = tpc.Length()*2.;

  // Call appropriate consumes<>() for any products to be retrieved by this module.
  consumes<std::vector<simb::MCTruth>>(fMCTruthLabel);
}

void ana::GENIETruthNuProtoDUNE::analyze(art::Event const& e)
{
  // Implementation of required member function here.
  fEventID = e.id().event();
  fRun = e.run();
  fSubRun = e.subRun();

  fInFV = false;

  art::Handle<std::vector<dunedaq::trgdataformats::TriggerActivityData>> taHandle;
  if (!e.getByLabel("tamakerTPC", taHandle)) {
      fTA = false;
  } else {
    if (taHandle->size() == 0) {
      fTA = false;
    } else {
      std::cout << ">>> Found " << taHandle->size() << " TAs in Event " << fEventID << std::endl;
      fTA = true;
    }
  }

  fnTAs = taHandle->size();

  const art::FindManyP<triggerprimitive_t> findTPsInTAs(taHandle, e, "tamakerTPC");
  if ( ! findTPsInTAs.isValid() ) {
    std::cout << " [WARNING] TPs not found in TA." << std::endl;
  }                                                                                                
  
  for (size_t ta = 0; ta < taHandle->size(); ta++) {
    std::cout << "START TA " << ta << " out of " << taHandle->size() << std::endl;
    auto fTPs = findTPsInTAs.at(ta);

    std::cout << "Found " << fTPs.size() << " TPs in TA " << ta << std::endl;

    timestamp_t first_tick = taHandle->at(ta).time_start;
    timestamp_t last_tick = taHandle->at(ta).time_end;
  
    timestamp_t TAWindow = last_tick - first_tick;
    if (TAWindow < 20e3) TAWindow = 20e3;

    std::cout << ">>> TAWindow = " << TAWindow << std::endl;
    double ADCIntSum = std::accumulate(fTPs.begin(), fTPs.end(), 0,
        [](double sum, const art::Ptr<triggerprimitive_t> &tp) { return sum + tp->adc_integral; });

    fTPTAADCIntSum.push_back(ADCIntSum);

    // Now sort in channel number
    std::sort(fTPs.begin(), fTPs.end(),
        [] (const art::Ptr<triggerprimitive_t> &lh, const art::Ptr<triggerprimitive_t> &rh) -> bool { return lh->channel < rh->channel; });

    channel_t current_chan = fTPs.at(0)->channel;
    
    std::cout << "First tick = " << first_tick << ", last tick = " << last_tick << std::endl;
    std::cout << "First channel = " << current_chan << std::endl;

    std::string title("");

    if (current_chan >= pCollectionAPA1IDs.first) {
      if (current_chan <= pCollectionAPA1IDs.second) {
        // APA 1, TPC 1
        fAPA_id = 1;
      }
    }
    if (current_chan >= pCollectionAPA3IDs.first) {
      if (current_chan <= pCollectionAPA3IDs.second) {
        // APA 3, TPC 2
        fAPA_id = 3;
      }
    }
    if (current_chan >= pCollectionAPA2IDs.first) {
      if (current_chan <= pCollectionAPA2IDs.second) {
        // APA 2, TPC 5
        fAPA_id = 2;
      }
    }
    if (current_chan >= pCollectionAPA4IDs.first) {
      if (current_chan <= pCollectionAPA4IDs.second) {
        // APA 4, TPC 6
        fAPA_id = 4;
      }
    }

    std::cout << "APA ID = " << fAPA_id << std::endl;
    fAPA_ids.push_back(fAPA_id);
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

      std::cout << "E = " << fE << std::endl;
      if (fE < 5. && fnTAs > 0) {
        std::cout << ">>> TRIGGERED LOW ENERGY E = " << fE << std::endl;
      }

      double fPrimaryStart[4];
      double fPrimaryVertex[4];

      const size_t numberTrajectoryPoints = neutrino.NumberTrajectoryPoints();
      const int last = numberTrajectoryPoints - 1;
      const TLorentzVector& positionStart = neutrino.Position(0);
      const TLorentzVector& positionEnd = neutrino.Position(last);
      // Set the vertex position - it should be the same value for each event	
      positionStart.GetXYZT(fPrimaryStart);
      positionEnd.GetXYZT(fPrimaryVertex);

      fnuStartX = fPrimaryStart[0];
      fnuStartY = fPrimaryStart[1];
      fnuStartZ = fPrimaryStart[2];
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
  hMCNeutrinoEnergy->Fill(fE);
}

// Define outputs at start of the job
void ana::GENIETruthNuProtoDUNE::beginJob() {
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;

  hMCNeutrinoEnergy = tfs->make<TH1D>("Total_MC_Nu_Energy", ";Energy (GeV);", 20, 0, 200);

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
  fSimulationNtuple->Branch("Decaynum", &fDecaynum, "Decaynum/D");
  fSimulationNtuple->Branch("POT", &fPOT, "POT/D");
  fSimulationNtuple->Branch("nuStartX", &fnuStartX, "nuStartX/D");
  fSimulationNtuple->Branch("nuStartY", &fnuStartY, "nuStartY/D");
  fSimulationNtuple->Branch("nuStartZ", &fnuStartZ, "nuStartZ/D");
  fSimulationNtuple->Branch("nuVertexX", &fnuVertexX, "nuVertexX/D");
  fSimulationNtuple->Branch("nuVertexY", &fnuVertexY, "nuVertexY/D");
  fSimulationNtuple->Branch("nuVertexZ", &fnuVertexZ, "nuVertexZ/D");
  fSimulationNtuple->Branch("InFV", &fInFV, "InFV/B");
  fSimulationNtuple->Branch("TA", &fTA, "TA/B");
  fSimulationNtuple->Branch("nTAs", &fnTAs, "nTAs/I");
  fSimulationNtuple->Branch("APA_ids", &fAPA_ids);
  fSimulationNtuple->Branch("fTPTAADCIntSum", &fTPTAADCIntSum);

  fSubrunTree = tfs->make<TTree>("SubRunTree", "SubRun-level information");
  fSubrunTree->Branch("POT", &fPOT, "POT/D");
  fSubrunTree->Branch("GoodPOT", &fGoodPOT, "GoodPOT/D");
}

void ana::GENIETruthNuProtoDUNE::beginSubRun(art::SubRun const& subRun) {
  
  const auto potSummaryHandle = subRun.getValidHandle<sumdata::POTSummary>("generator");
  const auto &potSummary = *potSummaryHandle;
  fPOT = potSummary.totpot;
  fGoodPOT = potSummary.totgoodpot;
 
  fTotalPOT += fPOT;
  std::cout << "POTSummary content: totpot = " << potSummary.totpot 
    << ", totgoodpot = " << potSummary.totgoodpot << std::endl;
      
  // Fill the TTree with the current subrun's POT information
  fSubrunTree->Fill();
}

void ana::GENIETruthNuProtoDUNE::endJob()
{
  // Implementation of optional member function here.
  std::cout << "Total POT = " << fTotalPOT << std::endl;
  hMCNeutrinoEnergy->Scale(fSetPOT / fTotalPOT);
}

DEFINE_ART_MODULE(ana::GENIETruthNuProtoDUNE)
