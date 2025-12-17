/**
 * @file GeneralProtoDUNETriggerTrainingDataMaker_module.cc
 *
 * @brief This module is an analyzer that displays information about trigger primitives, trigger activities, and trigger candidates.
 *
 * It reads trigger primitive, trigger activity, and trigger candidate data from the input event and fills three separate TTree objects with the data.
 * The module also provides general event information such as run number, subrun number, and event ID.
 *
 * The module takes three input tags to specify the collections of trigger primitive, trigger activity, and trigger candidate data.
 * It also supports an optional verbosity level to control the amount of output printed to the console.
 *
 * The filled TTree objects can be used for further analysis or visualization of the trigger data.
 */
////////////////////////////////////////////////////////////////////////
// Class:       GeneralProtoDUNETriggerTrainingDataMaker
// Plugin Type: analyzer (Unknown Unknown)
// File:        GeneralProtoDUNETriggerTrainingDataMaker_module.cc
//
// Generated at Mon Apr 29 11:24:28 2024 by Hamza Amar Es-sghir using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

#include "detdataformats/trigger/TriggerCandidateData.hpp"
#include "detdataformats/trigger/TriggerActivityData.hpp"
#include "detdataformats/trigger/TriggerPrimitive.hpp"
#include "detdataformats/DetID.hpp"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/RDTimeStamp.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataalg/DetectorInfo/DetectorProperties.h"

#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcore/Geometry/WireReadout.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/make_tool.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/Assns.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Additional framework includes
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
//#include "dunetrigger/TriggerSim/Verbosity.hh"

// ROOT includes
#include <TH1I.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TTree.h>
#include <TFile.h>
#include <TGraph.h>
#include <TRandom.h>
#include <TVector3.h>
#include <fcntl.h>

#include <memory>
#include <algorithm>
#include <iostream>

namespace duneana {
  class GeneralProtoDUNETriggerTrainingDataMaker;
  enum Verbosity{
    // these would all implicitly have these values
    // but it's best to explicitly define things
    // for anyone in the future looking
    kQuiet = 0,
    kInfo = 1,
    kDebug = 2,
    kVerbose = 3
  };
}


class duneana::GeneralProtoDUNETriggerTrainingDataMaker : public art::EDAnalyzer {
public:
  explicit GeneralProtoDUNETriggerTrainingDataMaker(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  GeneralProtoDUNETriggerTrainingDataMaker(GeneralProtoDUNETriggerTrainingDataMaker const&) = delete;
  GeneralProtoDUNETriggerTrainingDataMaker(GeneralProtoDUNETriggerTrainingDataMaker&&) = delete;
  GeneralProtoDUNETriggerTrainingDataMaker& operator=(GeneralProtoDUNETriggerTrainingDataMaker const&) = delete;
  GeneralProtoDUNETriggerTrainingDataMaker& operator=(GeneralProtoDUNETriggerTrainingDataMaker&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;

private:

  // Create output Tree
  TTree *fTPWindowTree;
  TTree *fTPNuWindowTree;
  TTree *fMCTruthTree;

  // General event information
  int fRun;
  int fSubRun;
  unsigned int fEventID;
  int fEventIterator;
  
  using timestamp_t = uint64_t;
  using channel_t = uint32_t;
  using version_t = uint16_t;
  using detid_t = uint16_t;
  
  ////////////////////////
  // Truth variables //
  ///////////////////////
  double fE;
  int fAPA; // >>> APA number of interaction
  int fTA; // >>> Was there a trigger?
  double fnuVertexX;
  double fnuVertexY;
  double fnuVertexZ;
  double fDriftDistance;
  double fVertexTime;
  timestamp_t fVertexTick;
  timestamp_t fVertexTick_new_evd;
  timestamp_t fVertexTick_new_weight;
  
  ////////////////////////
  // fTPTree variables //
  ///////////////////////
  channel_t fChannelID;
  timestamp_t fStart_time;
  timestamp_t fTime_over_threshold;
  timestamp_t fTime_peak;
  uint32_t fADC_integral;
  uint16_t fADC_peak;
  int fROP_ID;
  detid_t fDetID;
  int fType;
  int fAlgorithm;
  
  //////////////////////////
  // Generic time windows for cosmics
  /////////////////////////
  std::vector<std::vector<int>> fWindow_apacrp;
  std::vector<std::vector<int>> fWindow_planeid;
  std::vector<std::vector<timestamp_t>> fWindow_timepeak;
  std::vector<std::vector<channel_t>> fWindow_channelid;
  std::vector<std::vector<uint32_t>> fWindow_adcintegral;
  std::vector<std::vector<uint64_t>> fWindow_tot;
  std::vector<std::vector<uint64_t>> fWindow_adcpeak;

  //////////////////////////
  // Special time windows for signal events
  /////////////////////////
  std::vector<std::vector<int>> fNuWindow_apacrp;
  std::vector<std::vector<int>> fNuWindow_planeid;
  std::vector<std::vector<timestamp_t>> fNuWindow_timepeak;
  std::vector<std::vector<channel_t>> fNuWindow_channelid;
  std::vector<std::vector<uint32_t>> fNuWindow_adcintegral;
  std::vector<std::vector<uint64_t>> fNuWindow_tot;
  std::vector<std::vector<uint64_t>> fNuWindow_adcpeak;
  
  art::InputTag tp_tag_;
  art::InputTag ta_tag_;
  art::InputTag tc_tag_;
  int verbosity_;
  art::InputTag fMCTruthLabel; ///< The name of the producer that tracked

  // Trigger information
  std::vector<dunedaq::trgdataformats::TriggerActivityData> fTriggerActivity;
  std::vector<dunedaq::trgdataformats::TriggerPrimitive> fTriggerPrimitive;
  std::vector<dunedaq::trgdataformats::TriggerCandidateData> fTriggerCandidate;
  
};


duneana::GeneralProtoDUNETriggerTrainingDataMaker::GeneralProtoDUNETriggerTrainingDataMaker(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  , tp_tag_(p.get<art::InputTag>("tp_tag"))
  , ta_tag_(p.get<art::InputTag>("ta_tag"))
  , tc_tag_(p.get<art::InputTag>("tc_tag"))
  , verbosity_(p.get<int>("verbosity",0))
  , fMCTruthLabel(p.get<std::string>("MCTruthLabel"))
{
  consumes<std::vector<dunedaq::trgdataformats::TriggerPrimitive>>(tp_tag_);
  consumes<art::Assns<dunedaq::trgdataformats::TriggerPrimitive, dunedaq::trgdataformats::TriggerActivityData> >(ta_tag_);
  consumes<std::vector<dunedaq::trgdataformats::TriggerCandidateData>>(tc_tag_);
}

void duneana::GeneralProtoDUNETriggerTrainingDataMaker::analyze(art::Event const& e)
{

  fWindow_apacrp.clear();
  fWindow_planeid.clear();
  fWindow_timepeak.clear();
  fWindow_channelid.clear();
  fWindow_adcintegral.clear();
  fWindow_tot.clear();
  fWindow_adcpeak.clear();

  fWindow_apacrp.assign(9, std::vector<int>());
  fWindow_planeid.assign(9, std::vector<int>());
  fWindow_timepeak.assign(9, std::vector<timestamp_t>());
  fWindow_channelid.assign(9, std::vector<channel_t>());
  fWindow_adcintegral.assign(9, std::vector<uint32_t>());
  fWindow_tot.assign(9, std::vector<uint64_t>());
  fWindow_adcpeak.assign(9, std::vector<uint64_t>());

  // Clear Neutrino time windows
  fNuWindow_apacrp.clear();
  fNuWindow_planeid.clear();
  fNuWindow_timepeak.clear();
  fNuWindow_channelid.clear();
  fNuWindow_adcintegral.clear();
  fNuWindow_tot.clear();
  fNuWindow_adcpeak.clear();
  
  fNuWindow_apacrp.assign(1, std::vector<int>());
  fNuWindow_planeid.assign(1, std::vector<int>());
  fNuWindow_timepeak.assign(1, std::vector<timestamp_t>());
  fNuWindow_channelid.assign(1, std::vector<channel_t>());
  fNuWindow_adcintegral.assign(1, std::vector<uint32_t>());
  fNuWindow_tot.assign(1, std::vector<uint64_t>());
  fNuWindow_adcpeak.assign(1, std::vector<uint64_t>());

  // Set all general event information
  fRun    = e.run();
  fSubRun = e.subRun();
  fEventID = e.id().event();
    
  // Load the geometry service
  art::ServiceHandle<geo::Geometry> geom;
  geo::WireReadoutGeom const& wireReadout = art::ServiceHandle<geo::WireReadout>()->Get();
  // space charge and neutrino vertex correction.
  // Method taken from 
  // https://code-doc.larsoft.org/docs/latest/html/classevd_1_1SimulationDrawer.html#a8be0b3e65554ff1b8414c82baca11249
  const spacecharge::SpaceCharge* sce = lar::providerFrom<spacecharge::SpaceChargeService>();
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e);

  // Get truth information about neutrino if there is one
  art::Handle<std::vector<simb::MCTruth>> truthHandle;
  bool neutrinoMC(false);
  bool plane_exists(false);
  if (e.getByLabel(fMCTruthLabel, truthHandle)) {
    neutrinoMC = true;
    for (auto const& truth : (*truthHandle)) {
      if (truth.NeutrinoSet()) {
        const auto &nu = truth.GetNeutrino();
        const auto &neutrino = nu.Nu();

        fE = neutrino.E();
      
        double fPrimaryVertex[4];

        const TLorentzVector& positionStart = neutrino.Position(0);
        // Set the vertex position - it should be the same value for each event	
        positionStart.GetXYZT(fPrimaryVertex);
        fnuVertexX = fPrimaryVertex[0];
        fnuVertexY = fPrimaryVertex[1];
        fnuVertexZ = fPrimaryVertex[2];
      }
    }

    auto nuV_point = geo::Point_t(fnuVertexX, fnuVertexY, fnuVertexZ);
    
    geo::Point_t sceOffset{0, 0, 0};
    if (sce->EnableCorrSCE()) sceOffset = sce->GetPosOffsets(nuV_point);

    geo::Point_t const corr_nuV_point{fnuVertexX - sceOffset.X(), fnuVertexY + sceOffset.Y(), fnuVertexZ + sceOffset.Z()};

    geo::GeometryCore const* fGeometryService = lar::providerFrom<geo::Geometry>();
    auto plane = wireReadout.Plane(fGeometryService->FindTPCAtPosition(corr_nuV_point), geo::View_t::kW);
    plane_exists = wireReadout.HasPlane(plane.ID());
    fAPA = static_cast<int>(wireReadout.TPCtoTPCset(plane.ID().parentID()).deepestIndex()) + 1; // TPCset + 1 = APA/CRP
    std::cout << "APA/CRP number = " << fAPA << std::endl;

    if (plane_exists) {
      double time = detProp.ConvertXToTicks(corr_nuV_point.X(), plane.ID());
      fVertexTick_new_evd = static_cast<timestamp_t>(time);
      fVertexTick_new_weight = static_cast<timestamp_t>(time * 31.25);

      fDriftDistance = plane.DistanceFromPlane(nuV_point);
      fVertexTime = fDriftDistance / 0.16; //>> drift speed  = 0.16cm/us
      fVertexTick = static_cast<timestamp_t>((fVertexTime * 1000.) / 16.); //>> 16ns per tick
    } else { // only calculate drift distance for real TPCs
      fDriftDistance = 0.;
      fVertexTime = 0.;
      fVertexTick = 0;
      fVertexTick_new_evd = 0;
      fVertexTick_new_weight = 0;
    }
  }

  // Define neutrino time window if it exists
  timestamp_t nuWindowStart = 0;
  timestamp_t nuWindowEnd = 0;
  bool doNuWindow = false;

  std::cout << "VertexTick = " << fVertexTick << ", Plane = " << fAPA << std::endl;
  if (plane_exists && neutrinoMC) {
    timestamp_t vertex_tick_copy = fVertexTick_new_weight;
    int int_nuWindowStart = static_cast<int>(vertex_tick_copy) - 8000;
    if (int_nuWindowStart < 0) int_nuWindowStart = 0;
    nuWindowStart = static_cast<timestamp_t>(int_nuWindowStart);

    nuWindowEnd = nuWindowStart + 20000;
    //if (nuWindowEnd > 200000) nuWindowEnd = 200000;
    doNuWindow = true;
    std::cout << "doNuWindow TRUE: APA = " << fAPA << std::endl;
    std::cout << "start window = " << nuWindowStart << ", end = " << nuWindowEnd << std::endl;
  }

  if (neutrinoMC && !doNuWindow) {
    return;
  }
  // Done with MC Tree now fill once for each neutrino
  fMCTruthTree -> Fill();

  // Take TPs from event
  auto tp_handle = e.getValidHandle< std::vector<dunedaq::trgdataformats::TriggerPrimitive> >(tp_tag_);  
  fTriggerPrimitive = *tp_handle;

  // Take TAs from event
  auto ta_handle = e.getValidHandle< std::vector<dunedaq::trgdataformats::TriggerActivityData> >(ta_tag_);
  fTriggerActivity = *ta_handle;

  // Take TPs from TP-TA association from event
  auto tpfromtpta_handle = e.getValidHandle< art::Assns<dunedaq::trgdataformats::TriggerPrimitive, dunedaq::trgdataformats::TriggerActivityData> >(ta_tag_);
  auto fTPfromTPTA = *tpfromtpta_handle;

  // Take TAs from TP-TA association from event
  auto tafromtpta_handle = e.getValidHandle< art::Assns<dunedaq::trgdataformats::TriggerActivityData, dunedaq::trgdataformats::TriggerPrimitive> >(ta_tag_);
  auto fTAfromTPTA = *tafromtpta_handle;

  if(verbosity_ >= Verbosity::kInfo)
  {
    std::cout << "Found " << fTriggerPrimitive.size() << " TPs" << std::endl;
    std::cout << "Found " << fTriggerActivity.size() << " TAs" << std::endl;
  }
 
  fTA = (int)fTriggerActivity.size();

  // Fill TP tree
  for(long unsigned int i=0; i < fTriggerPrimitive.size(); i++)
  {
    fChannelID = fTriggerPrimitive[i].channel;
    fStart_time = fTriggerPrimitive[i].time_start;
    fTime_over_threshold = fTriggerPrimitive[i].time_over_threshold;
    fTime_peak = fTriggerPrimitive[i].time_peak;
    fADC_integral = fTriggerPrimitive[i].adc_integral;
    fADC_peak = fTriggerPrimitive[i].adc_peak;
    fDetID = fTriggerPrimitive[i].detid;
    fType = static_cast<int>(fTriggerPrimitive[i].type);
    fAlgorithm = static_cast<int>(fTriggerPrimitive[i].algorithm);
    
    int apa = 0;
    auto rop = wireReadout.ChannelToROP(fChannelID);
    auto tpc = rop.parentID().TPCset;
    apa = tpc + 1;

    //if (rop.ROP == 0 || rop.ROP == 1) continue; // Only look at collection plane for now
    int fROP = static_cast<int>(rop.ROP);
    if (fROP == 3) fROP = 2; // rop of 3 is just collection plane in NP04

    // Determine the time window index (0 to 9).
    int windowIndex = fTime_peak / 20000;
    //if (windowIndex >= 10) windowIndex = 9;  // safeguard for times at the upper edge
    if (windowIndex >= 9) continue;  // safeguard for times at the upper edge and only accept 9 windows per readout
        
    fWindow_apacrp[windowIndex].push_back(apa); 
    fWindow_planeid[windowIndex].push_back(fROP); 
    fWindow_timepeak[windowIndex].push_back(fTime_peak); 
    fWindow_channelid[windowIndex].push_back(fChannelID); 
    fWindow_adcintegral[windowIndex].push_back(fADC_integral); 
    fWindow_tot[windowIndex].push_back(fTime_over_threshold); 
    fWindow_adcpeak[windowIndex].push_back(fADC_peak);

    // Fill neutrino window if defined
    if (doNuWindow && fTime_peak >= nuWindowStart && fTime_peak <= nuWindowEnd && apa == fAPA) {
      fNuWindow_apacrp[0].push_back(fAPA);
      fNuWindow_planeid[0].push_back(fROP);
      fNuWindow_timepeak[0].push_back(fTime_peak);
      fNuWindow_channelid[0].push_back(fChannelID);
      fNuWindow_adcintegral[0].push_back(fADC_integral);
      fNuWindow_tot[0].push_back(fTime_over_threshold); 
      fNuWindow_adcpeak[0].push_back(fADC_peak);
    }

  }

  fTPWindowTree -> Fill();

  if (doNuWindow && neutrinoMC) {
    fTPNuWindowTree -> Fill(); 
  }
}

void duneana::GeneralProtoDUNETriggerTrainingDataMaker::beginJob()
{
  // Make our handle to the TFileService
  art::ServiceHandle<art::TFileService> tfs;
  // The TTrees
  fTPWindowTree = tfs->make<TTree>("TPWindowTree", "time windows for cosmic images");
  
  fTPNuWindowTree = tfs->make<TTree>("TPNuWindowTree", "time windows for neutrino images");
  fMCTruthTree = tfs->make<TTree>("GenieTruth", "GENIE Output Tree");

  // Add branches to TTree
  fMCTruthTree -> Branch( "eventID", &fEventID);
  fMCTruthTree -> Branch( "SubRun", &fSubRun, "SubRun/I");
  fMCTruthTree -> Branch( "Run", &fRun, "Run/I");
  fMCTruthTree -> Branch( "APA_CRP", &fAPA);
  fMCTruthTree -> Branch( "TA" , &fTA, "TA/I");
  fMCTruthTree -> Branch( "DriftDistance", &fDriftDistance);
  fMCTruthTree -> Branch( "VertexTime", &fVertexTime);
  fMCTruthTree -> Branch( "VertexTick", &fVertexTick);
  fMCTruthTree -> Branch( "VertexTick_new_evd", &fVertexTick_new_evd);
  fMCTruthTree -> Branch( "VertexTick_new_weight", &fVertexTick_new_weight);
  fMCTruthTree -> Branch( "E", &fE, "E/D");
  fMCTruthTree -> Branch( "nuVertexX", &fnuVertexX, "nuVertexX/D");
  fMCTruthTree -> Branch( "nuVertexY", &fnuVertexY, "nuVertexY/D");
  fMCTruthTree -> Branch( "nuVertexZ", &fnuVertexZ, "nuVertexZ/D");
  
  ////////////////////////////////////////
  // fTriggerPrimitive tree information //
  ////////////////////////////////////////
  fTPWindowTree -> Branch( "Window_apacrp", &fWindow_apacrp);
  fTPWindowTree -> Branch( "Window_planeid", &fWindow_planeid);
  fTPWindowTree -> Branch( "Window_timepeak", &fWindow_timepeak);
  fTPWindowTree -> Branch( "Window_channelid", &fWindow_channelid);
  fTPWindowTree -> Branch( "Window_adcintegral", &fWindow_adcintegral);
  fTPWindowTree -> Branch( "Window_tot", &fWindow_tot);
  fTPWindowTree -> Branch( "Window_adcpeak", &fWindow_adcpeak);
  
  fTPNuWindowTree -> Branch( "APA" , &fAPA, "APA/I"  );
  fTPNuWindowTree -> Branch( "TA" , &fTA, "TA/I"  );
  fTPNuWindowTree -> Branch( "Window_apacrp", &fNuWindow_apacrp);
  fTPNuWindowTree -> Branch( "Window_planeid", &fNuWindow_planeid);
  fTPNuWindowTree -> Branch( "Window_timepeak", &fNuWindow_timepeak);
  fTPNuWindowTree -> Branch( "Window_channelid", &fNuWindow_channelid);
  fTPNuWindowTree -> Branch( "Window_adcintegral", &fNuWindow_adcintegral);
  fTPNuWindowTree -> Branch( "Window_tot", &fNuWindow_tot);
  fTPNuWindowTree -> Branch( "Window_adcpeak", &fNuWindow_adcpeak);
  
}

DEFINE_ART_MODULE(duneana::GeneralProtoDUNETriggerTrainingDataMaker)
