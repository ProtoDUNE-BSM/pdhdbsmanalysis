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
#include "larcore/Geometry/WireReadout.h"
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
  TH1D *hMCNumuEnergy;
  TH1D *hMCNueEnergy;
  TH1D *hMCNumubarEnergy;
  TH1D *hMCNuebarEnergy;

  double fSetPOT;

  int fSimPDG;     ///< PDG ID of the particle being processed
  int fSimTrackID; ///< GEANT ID of the particle being processed
  int fCCNC; ///< Is neutrino interaction a CC or NC interaction
  int fTarget;
  unsigned int fTPCID; ///< TPC ID where neutrino interacts

  double fE;
  double fnuVertexX;
  double fnuVertexY;
  double fnuVertexZ;
  double fnuPx;
  double fnuPy;
  double fnuPz;

  bool fInFV;
  
  geo::GeometryCore const* fGeometryService; ///< pointer to Geometry provider
  //geo::WireReadoutGeom const& wireReadout;
  std::vector<double> fFiducialBoundaries;
  double fZedge;

  double fPOT;
  double fGoodPOT;
  double fTotalPOT;

  bool fTA;
  bool fTAcollection;
  int fnTAs;
  int fROP;
  int fAPA_id;
  std::vector<int> fROPs;
  std::vector<int> fAPA_ids;
  std::vector<double> fTPTAADCIntSum;

};


// Analyser class constructor
ana::GENIETruthNuProtoDUNE::GENIETruthNuProtoDUNE(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
  , fMCTruthLabel(p.get<std::string>("MCTruthLabel"))
  , fSetPOT(p.get<double>("SetPOT"))
  // More initializers here.
{
  
  fROP = 0;
  fAPA_id = 0;
  fROPs.clear();
  fAPA_ids.clear();
  pCollectionAPA1IDs = std::make_pair(2080, 2559);
  pCollectionAPA2IDs = std::make_pair(7200, 7680);  
  pCollectionAPA3IDs = std::make_pair(4160, 4639);
  pCollectionAPA4IDs = std::make_pair(9280, 9759); 
 
  // Get a pointer to the geometry service provider.
  fGeometryService = lar::providerFrom<geo::Geometry>();
  std::string info = fGeometryService->Info();
  //std::cout << info;
  // TPC 1 is the first proper TPC - TPC 0 is for track stubs
  const geo::TPCGeo& tpc = fGeometryService->Cryostat().TPC(1);
  //geo::WireReadoutGeom const& wireReadout = art::ServiceHandle<geo::WireReadout>()->Get();
  //wireReadout = art::ServiceHandle<geo::WireReadout>()->Get();
  //fFiducialBoundaries.push_back(0.); // central x
  //fFiducialBoundaries.push_back(tpc.Width() - 0.05*tpc.Width()); // outer x
  //fFiducialBoundaries.push_back(0.05*tpc.Height()); // bottom y
  //fFiducialBoundaries.push_back(tpc.Height() - 0.05*tpc.Height()); // top y
  //fFiducialBoundaries.push_back(0.05*(tpc.Length()*2.));
  //fFiducialBoundaries.push_back((tpc.Length()*2.) - 0.05*(tpc.Length()*2));
  
  fFiducialBoundaries.push_back(0.); // central x
  fFiducialBoundaries.push_back(tpc.Width()); // outer x
  fFiducialBoundaries.push_back(0.); // bottom y
  fFiducialBoundaries.push_back(tpc.Height()); // top y
  fFiducialBoundaries.push_back(0.);
  fFiducialBoundaries.push_back(tpc.Length()*2.);

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
  
  geo::WireReadoutGeom const& wireReadout = art::ServiceHandle<geo::WireReadout>()->Get();
  
  fROPs.clear();
  fAPA_ids.clear();
  
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

  fTAcollection = false;

  for (size_t ta = 0; ta < taHandle->size(); ta++) {
    const art::FindManyP<triggerprimitive_t> findTPsInTAs(taHandle, e, "tamakerTPC");
    if ( ! findTPsInTAs.isValid() ) {
      std::cout << " [WARNING] TPs not found in TA." << std::endl;
    }                                                                                                
    auto fTPs = findTPsInTAs.at(ta);

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

    auto rop = wireReadout.ChannelToROP(current_chan);
    auto tpc = rop.parentID().TPCset;
    fROP = rop.ROP;
    fAPA_id = tpc;
   
    // Often only interested if there is a collection TA
    if (fROP == 2 || fROP == 3) fTAcollection = true;
    else if (fAPA_id == 0 && fROP == 1) fTAcollection = true;
    //if (faPA_id == 0 && fROP == 3)

    fROPs.push_back(fROP);
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

      fTarget = nu.Target();

      fE = neutrino.E();

      double fPrimaryVertex[4];

      const TLorentzVector& positionStart = neutrino.Position(0);
      // Set the vertex position - it should be the same value for each event	
      positionStart.GetXYZT(fPrimaryVertex);

      fnuVertexX = fPrimaryVertex[0];
      fnuVertexY = fPrimaryVertex[1];
      fnuVertexZ = fPrimaryVertex[2];

      fnuPx = neutrino.Px();
      fnuPy = neutrino.Py();
      fnuPz = neutrino.Pz();

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
      if (fTPCID > 7 || fTPCID < 0) fTPCID = -1;
      
      // Store total event outputs in the TTree
      fSimulationNtuple->Fill();
    }
  }
  hMCNeutrinoEnergy->Fill(fE);
 
  switch(fSimPDG) {
    case 14:
      hMCNumuEnergy->Fill(fE);
      break;
    case 12:
      hMCNueEnergy->Fill(fE);
      break;
    case -14:
      hMCNumubarEnergy->Fill(fE);
      break;
    case -12:
      hMCNuebarEnergy->Fill(fE);
      break;
    default:
      std::cout << "Warning - no pdg recognised!" << std::endl;
  }

}

// Define outputs at start of the job
void ana::GENIETruthNuProtoDUNE::beginJob() {
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;

  hMCNeutrinoEnergy = tfs->make<TH1D>("Total_MC_Nu_Energy", ";Energy (GeV);", 20, 0, 200);
  hMCNumuEnergy = tfs->make<TH1D>("Total_MC_Numu_Energy", ";Energy (GeV);", 20, 0, 200);
  hMCNueEnergy = tfs->make<TH1D>("Total_MC_Nue_Energy", ";Energy (GeV);", 20, 0, 200);
  hMCNumubarEnergy = tfs->make<TH1D>("Total_MC_Numubar_Energy", ";Energy (GeV);", 20, 0, 200);
  hMCNuebarEnergy = tfs->make<TH1D>("Total_MC_Nuebar_Energy", ";Energy (GeV);", 20, 0, 200);

  // Get TFileService to create an output tree
  fSimulationNtuple = tfs->make<TTree>("GenieTruth", "GENIE Output Tree");

  // Add branches to TTree
  fSimulationNtuple->Branch("eventID", &fEventID);
  fSimulationNtuple->Branch("SubRun", &fSubRun, "SubRun/I");
  fSimulationNtuple->Branch("Run", &fRun, "Run/I");
  fSimulationNtuple->Branch("PDG", &fSimPDG, "PDG/I");
  fSimulationNtuple->Branch("CCNC", &fCCNC, "CCNC/I");
  fSimulationNtuple->Branch("Target", &fTarget, "Target/I");
  fSimulationNtuple->Branch("TPCID", &fTPCID);

  fSimulationNtuple->Branch("E", &fE, "E/D");
  fSimulationNtuple->Branch("POT", &fPOT, "POT/D");
  fSimulationNtuple->Branch("nuVertexX", &fnuVertexX, "nuVertexX/D");
  fSimulationNtuple->Branch("nuVertexY", &fnuVertexY, "nuVertexY/D");
  fSimulationNtuple->Branch("nuVertexZ", &fnuVertexZ, "nuVertexZ/D");
  fSimulationNtuple->Branch("nuPx", &fnuPx, "nuPx/D");
  fSimulationNtuple->Branch("nuPy", &fnuPy, "nuPy/D");
  fSimulationNtuple->Branch("nuPz", &fnuPz, "nuPz/D");
  fSimulationNtuple->Branch("InFV", &fInFV, "InFV/B");
  fSimulationNtuple->Branch("TA", &fTA, "TA/B");
  fSimulationNtuple->Branch("TAcollection", &fTAcollection, "TAcollection/B");
  fSimulationNtuple->Branch("nTAs", &fnTAs, "nTAs/I");
  fSimulationNtuple->Branch("ROPs", &fROPs);
  fSimulationNtuple->Branch("APA_ids", &fAPA_ids);
  fSimulationNtuple->Branch("fTPTAADCIntSum", &fTPTAADCIntSum);

  fSubrunTree = tfs->make<TTree>("SubRunTree", "SubRun-level information");
  fSubrunTree->Branch("POT", &fPOT, "POT/D");
  fSubrunTree->Branch("GoodPOT", &fGoodPOT, "GoodPOT/D");
}

void ana::GENIETruthNuProtoDUNE::beginSubRun(art::SubRun const& subRun) {
  
  const auto potSummaryHandle = subRun.getValidHandle<sumdata::POTSummary>(fMCTruthLabel);
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
  hMCNumuEnergy->Scale(fSetPOT / fTotalPOT);
  hMCNueEnergy->Scale(fSetPOT / fTotalPOT);
  hMCNumubarEnergy->Scale(fSetPOT / fTotalPOT);
  hMCNuebarEnergy->Scale(fSetPOT / fTotalPOT);

  double total_events = hMCNeutrinoEnergy->Integral();
  double numu_events = hMCNumuEnergy->Integral();
  double nue_events = hMCNueEnergy->Integral();
  double numubar_events = hMCNumubarEnergy->Integral();
  double nuebar_events = hMCNuebarEnergy->Integral();

  std::string total_title = "Total #nu: " + std::to_string(total_events) + " in " + std::to_string(fSetPOT) + " POT";
  std::string numu_title = "#nu_{#mu}: " + std::to_string(numu_events) + " in " + std::to_string(fSetPOT) + " POT";
  std::string nue_title = "#nu_{e}: " + std::to_string(nue_events) + " in " + std::to_string(fSetPOT) + " POT";
  std::string numubar_title = "#bar{#nu}_{#mu}: " + std::to_string(numubar_events) + " in " + std::to_string(fSetPOT) + " POT";
  std::string nuebar_title = "#bar{#nu}_{e}: " + std::to_string(nuebar_events) + " in " + std::to_string(fSetPOT) + " POT";

  hMCNeutrinoEnergy->SetTitle(total_title.c_str());
  hMCNumuEnergy->SetTitle(numu_title.c_str());
  hMCNueEnergy->SetTitle(nue_title.c_str());
  hMCNumubarEnergy->SetTitle(numubar_title.c_str());
  hMCNuebarEnergy->SetTitle(nuebar_title.c_str());

  std::cout << total_title << std::endl;
  std::cout << numu_title << std::endl;
  std::cout << nue_title << std::endl;
  std::cout << numubar_title << std::endl;
  std::cout << nuebar_title << std::endl;

}

DEFINE_ART_MODULE(ana::GENIETruthNuProtoDUNE)
