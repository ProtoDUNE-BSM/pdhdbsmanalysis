/**
 * @file SmallTriggerTPCInfoDisplay_module.cc
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
// Class:       SmallTriggerTPCInfoDisplay
// Plugin Type: analyzer (Unknown Unknown)
// File:        SmallTriggerTPCInfoDisplay_module.cc
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
  class SmallTriggerTPCInfoDisplay;
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


class duneana::SmallTriggerTPCInfoDisplay : public art::EDAnalyzer {
public:
  explicit SmallTriggerTPCInfoDisplay(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SmallTriggerTPCInfoDisplay(SmallTriggerTPCInfoDisplay const&) = delete;
  SmallTriggerTPCInfoDisplay(SmallTriggerTPCInfoDisplay&&) = delete;
  SmallTriggerTPCInfoDisplay& operator=(SmallTriggerTPCInfoDisplay const&) = delete;
  SmallTriggerTPCInfoDisplay& operator=(SmallTriggerTPCInfoDisplay&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;

private:

  // Create output Tree
  TTree *fTPTree;
  TTree *fTPWindowAPA1Tree;
  TTree *fTPWindowAPA2Tree;
  TTree *fTPWindowAPA3Tree;
  TTree *fTPWindowAPA4Tree;
  TTree *fTPNuWindowTree;
  TTree *fTAWindowTree;
  TTree *fTPfromTPTATree;
  TTree *fTATree;
  TTree *fTCTree;
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
  
  double fE;
  int fTPCID; ///< TPC ID where neutrino interacts
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
  std::vector<std::vector<timestamp_t>> fAPA1Window_timepeak;
  std::vector<std::vector<channel_t>> fAPA1Window_channelid;
  std::vector<std::vector<uint32_t>> fAPA1Window_adcintegral;
  std::vector<std::vector<uint64_t>> fAPA1Window_tot;
  std::vector<std::vector<uint64_t>> fAPA1Window_adcpeak;
  
  std::vector<std::vector<timestamp_t>> fAPA2Window_timepeak;
  std::vector<std::vector<channel_t>> fAPA2Window_channelid;
  std::vector<std::vector<uint32_t>> fAPA2Window_adcintegral;
  std::vector<std::vector<uint64_t>> fAPA2Window_tot;
  std::vector<std::vector<uint64_t>> fAPA2Window_adcpeak;

  std::vector<std::vector<timestamp_t>> fAPA3Window_timepeak;
  std::vector<std::vector<channel_t>> fAPA3Window_channelid;
  std::vector<std::vector<uint32_t>> fAPA3Window_adcintegral;
  std::vector<std::vector<uint64_t>> fAPA3Window_tot;
  std::vector<std::vector<uint64_t>> fAPA3Window_adcpeak;
  
  std::vector<std::vector<timestamp_t>> fAPA4Window_timepeak;
  std::vector<std::vector<channel_t>> fAPA4Window_channelid;
  std::vector<std::vector<uint32_t>> fAPA4Window_adcintegral;
  std::vector<std::vector<uint64_t>> fAPA4Window_tot;
  std::vector<std::vector<uint64_t>> fAPA4Window_adcpeak;
  
  //////////////////////////
  // Special time windows for signal events
  /////////////////////////
  std::vector<std::vector<timestamp_t>> fNuWindow_timepeak;
  std::vector<std::vector<channel_t>> fNuWindow_channelid;
  std::vector<std::vector<uint32_t>> fNuWindow_adcintegral;
  std::vector<std::vector<uint64_t>> fNuWindow_tot;
  std::vector<std::vector<uint64_t>> fNuWindow_adcpeak;
  
  //////////////////////////
  // Special time windows for TAs
  /////////////////////////
  std::vector<int> fapaTA;
  std::vector<std::vector<timestamp_t>> fTAWindow_timepeak;
  std::vector<std::vector<channel_t>> fTAWindow_channelid;
  std::vector<std::vector<uint32_t>> fTAWindow_adcintegral;
  std::vector<std::vector<uint64_t>> fTAWindow_tot;
  std::vector<std::vector<uint64_t>> fTAWindow_adcpeak;
  
  ////////////////////////
  // fTATree variables //
  ///////////////////////
  channel_t fChannel_start_TA;
  channel_t fChannel_end_TA;
  channel_t fChannel_peak_TA;
  int fTAROP_ID;
  timestamp_t fTime_start_TA;
  timestamp_t fTime_end_TA;
  timestamp_t fTime_peak_TA;
  timestamp_t fTime_activity;
  uint32_t fADC_integral_TA;
  uint16_t fADC_peak_TA;
  int fAlgorithm_TA;

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


duneana::SmallTriggerTPCInfoDisplay::SmallTriggerTPCInfoDisplay(fhicl::ParameterSet const& p)
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

void duneana::SmallTriggerTPCInfoDisplay::analyze(art::Event const& e)
{

  fAPA1Window_timepeak.clear();
  fAPA2Window_timepeak.clear();
  fAPA3Window_timepeak.clear();
  fAPA4Window_timepeak.clear();
  
  fAPA1Window_channelid.clear();
  fAPA2Window_channelid.clear();
  fAPA3Window_channelid.clear();
  fAPA4Window_channelid.clear();
  
  fAPA1Window_adcintegral.clear();
  fAPA2Window_adcintegral.clear();
  fAPA3Window_adcintegral.clear();
  fAPA4Window_adcintegral.clear();
  
  fAPA1Window_tot.clear();
  fAPA2Window_tot.clear();
  fAPA3Window_tot.clear();
  fAPA4Window_tot.clear();
  
  fAPA1Window_adcpeak.clear();
  fAPA2Window_adcpeak.clear();
  fAPA3Window_adcpeak.clear();
  fAPA4Window_adcpeak.clear();

  fAPA1Window_timepeak.assign(9, std::vector<timestamp_t>());
  fAPA2Window_timepeak.assign(9, std::vector<timestamp_t>());
  fAPA3Window_timepeak.assign(9, std::vector<timestamp_t>());
  fAPA4Window_timepeak.assign(9, std::vector<timestamp_t>());

  fAPA1Window_channelid.assign(9, std::vector<channel_t>());
  fAPA2Window_channelid.assign(9, std::vector<channel_t>());
  fAPA3Window_channelid.assign(9, std::vector<channel_t>());
  fAPA4Window_channelid.assign(9, std::vector<channel_t>());

  fAPA1Window_adcintegral.assign(9, std::vector<uint32_t>());
  fAPA2Window_adcintegral.assign(9, std::vector<uint32_t>());
  fAPA3Window_adcintegral.assign(9, std::vector<uint32_t>());
  fAPA4Window_adcintegral.assign(9, std::vector<uint32_t>());

  fAPA1Window_tot.assign(9, std::vector<uint64_t>());
  fAPA2Window_tot.assign(9, std::vector<uint64_t>());
  fAPA3Window_tot.assign(9, std::vector<uint64_t>());
  fAPA4Window_tot.assign(9, std::vector<uint64_t>());
  
  fAPA1Window_adcpeak.assign(9, std::vector<uint64_t>());
  fAPA2Window_adcpeak.assign(9, std::vector<uint64_t>());
  fAPA3Window_adcpeak.assign(9, std::vector<uint64_t>());
  fAPA4Window_adcpeak.assign(9, std::vector<uint64_t>());

  // Clear Neutrino time windows
  fNuWindow_timepeak.clear();
  fNuWindow_channelid.clear();
  fNuWindow_adcintegral.clear();
  fNuWindow_tot.clear();
  fNuWindow_adcpeak.clear();
  
  fNuWindow_timepeak.assign(1, std::vector<timestamp_t>());
  fNuWindow_channelid.assign(1, std::vector<channel_t>());
  fNuWindow_adcintegral.assign(1, std::vector<uint32_t>());
  fNuWindow_tot.assign(1, std::vector<uint64_t>());
  fNuWindow_adcpeak.assign(1, std::vector<uint64_t>());

  // assign TA windows later
  fTAWindow_timepeak.clear();
  fTAWindow_channelid.clear();
  fTAWindow_adcintegral.clear();
  fTAWindow_tot.clear();
  fTAWindow_adcpeak.clear();
  fapaTA.clear();


  // Set all general event information
  fRun    = e.run();
  fSubRun = e.subRun();
  fEventID = e.id().event();
  fEventIterator++;
    
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
  //bool inFV = false;
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

    //bool inFV = false;
    /*if (fnuVertexX < 350 && fnuVertexX > -350 
        && fnuVertexY > 0 && fnuVertexY < 607 
        && fnuVertexZ > 0 && fnuVertexZ < 460) {
      inFV = true;
    }*/

    auto nuV_point = geo::Point_t(fnuVertexX, fnuVertexY, fnuVertexZ);
    
    geo::Point_t sceOffset{0, 0, 0};
    if (sce->EnableCorrSCE()) sceOffset = sce->GetPosOffsets(nuV_point);

    geo::Point_t const corr_nuV_point{fnuVertexX - sceOffset.X(), fnuVertexY + sceOffset.Y(), fnuVertexZ + sceOffset.Z()};

    geo::GeometryCore const* fGeometryService = lar::providerFrom<geo::Geometry>();
    fTPCID = static_cast<int>(fGeometryService->FindTPCAtPosition(nuV_point).TPC);
    std::cout << "TPCID = " << fTPCID << std::endl;
    if (fTPCID > 7 || fTPCID < 0) fTPCID = -1;
 
    switch (fTPCID) {
      case 1:
        fAPA = 1;
        break;
      case 2:
        fAPA = 3;
        break;
      case 5:
        fAPA = 2;
        break;
      case 6:
        fAPA = 4;
        break;
      default:
        fAPA = 0;
        break;
    }

    //if (fAPA == 1 || fAPA == 2 || fAPA == 3 || fAPA == 4) {
      auto plane = wireReadout.Plane(fGeometryService->FindTPCAtPosition(corr_nuV_point), geo::View_t::kW);

      double time = detProp.ConvertXToTicks(corr_nuV_point.X(), plane.ID());
      fVertexTick_new_evd = static_cast<timestamp_t>(time);
      fVertexTick_new_weight = static_cast<timestamp_t>(time * 31.25);

      fDriftDistance = plane.DistanceFromPlane(nuV_point);
      fVertexTime = fDriftDistance / 0.16; //>> drift speed  = 0.16cm/us
      fVertexTick = static_cast<timestamp_t>((fVertexTime * 1000.) / 16.); //>> 16ns per tick
    /*} else { // only calculate drift distance for real TPCs
      fDriftDistance = 0.;
      fVertexTime = 0.;
      fVertexTick = 0;
      fVertexTick_new_evd = 0;
      fVertexTick_new_weight = 0;
    }*/
  }

  // Define neutrino time window if it exists
  timestamp_t nuWindowStart = 0;
  timestamp_t nuWindowEnd = 0;
  bool doNuWindow = false;

  std::cout << "VertexTick = " << fVertexTick << ", APA = " << fAPA << std::endl;
  //if (fAPA > 0 && fVertexTick > 0 && neutrinoMC && inFV) {
  if (fVertexTick > 0 && neutrinoMC) {
    //timestamp_t vertex_tick_copy = fVertexTick;
    timestamp_t vertex_tick_copy = fVertexTick_new_weight;
    //int int_nuWindowStart = static_cast<int>(vertex_tick_copy) - 10000;
    int int_nuWindowStart = static_cast<int>(vertex_tick_copy) - 8000;
    if (int_nuWindowStart < 0) int_nuWindowStart = 0;
    nuWindowStart = static_cast<timestamp_t>(int_nuWindowStart);

    nuWindowEnd = nuWindowStart + 20000;
    if (nuWindowEnd > 200000) nuWindowEnd = 200000;
    doNuWindow = true;
    std::cout << "doNuWindow TRUE: APA = " << fAPA << std::endl;
    std::cout << "start window = " << nuWindowStart << ", end = " << nuWindowEnd << std::endl;
  }

  if (neutrinoMC && !doNuWindow) {
    fEventIterator--;
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

  // Take TCs from event
  //auto tc_handle = e.getValidHandle< std::vector<dunedaq::trgdataformats::TriggerCandidateData> >(tc_tag_);
  //fTriggerCandidate = *tc_handle;
  
  if(verbosity_ >= Verbosity::kInfo)
  {
    //std::cout << "Found " << rawdigit_vec.size() << " raw::RawDigits" << std::endl;
    std::cout << "Found " << fTriggerPrimitive.size() << " TPs" << std::endl;
    std::cout << "Found " << fTriggerActivity.size() << " TAs" << std::endl;
    //std::cout << "Found " << fTriggerCandidate.size() << " TCs" << std::endl;
  }
 
  fTA = (int)fTriggerActivity.size();
  fTAWindow_timepeak.assign(fTA, std::vector<timestamp_t>());
  fTAWindow_channelid.assign(fTA, std::vector<channel_t>());
  fTAWindow_adcintegral.assign(fTA, std::vector<uint32_t>());
  fTAWindow_tot.assign(fTA, std::vector<uint64_t>());
  fTAWindow_adcpeak.assign(fTA, std::vector<uint64_t>());

  // Fill TA tree
  for(long unsigned int i=0; i < fTriggerActivity.size(); i++)
  {
    fChannel_start_TA = fTriggerActivity[i].channel_start;
    fChannel_end_TA = fTriggerActivity[i].channel_end;
    fChannel_peak_TA = fTriggerActivity[i].channel_peak;
    fTime_start_TA = fTriggerActivity[i].time_start;
    fTime_end_TA = fTriggerActivity[i].time_end;
    fTime_peak_TA = fTriggerActivity[i].time_peak;
    fTime_activity = fTriggerActivity[i].time_activity;
    fADC_integral_TA = fTriggerActivity[i].adc_integral;
    fADC_peak_TA = fTriggerActivity[i].adc_peak;
    fAlgorithm_TA = static_cast<int>(fTriggerActivity[i].algorithm);
    
    auto rop = wireReadout.ChannelToROP(fChannelID);
    auto tpc = rop.parentID().TPCset;
    //auto first_channel_rop = wireReadout.FirstChannelInROP(rop);
    //auto n_channels_rop = wireReadout.Nchannels(rop);
    //fROP_ID = rop.ROP;
    //auto tpcid = wireReadout.ROPtoTPCs(rop);
    //int apaTA = rop.TPCset;
    int apaTA = tpc;
    //std::cout << ">>> TA in APA " << apaTA << "\n";
    //int rop_index = rop.deepestIndex();
    //std::cout << "rop: " << rop_index << ", 1st ch: " << first_channel_rop << ", n chan: " << n_channels_rop << "\n";

    /*
    int apaTA = 0;
    for (const auto &t : tpcid) {
      std::cout << "TA in TPC " << t.TPC << "\n";
      if (t.TPC == 0 || t.TPC == 1) {
        apaTA = 1;
      } else if (t.TPC == 2 || t.TPC == 3) {
        apaTA = 3;
      } else if (t.TPC == 4 || t.TPC == 5) {
        apaTA = 2;
      } else if (t.TPC == 6 || t.TPC == 7) {
        apaTA = 4;
      }
    }
    */
/*
    int apaTA = 0;
    //if (fTriggerActivity[i].channel_start >= 2080 && fTriggerActivity[i].channel_end <= 2559) {
    if (fTriggerActivity[i].channel_end <= 2560) {
      apaTA = 1;
    } else if (fTriggerActivity[i].channel_start >= 7200 && fTriggerActivity[i].channel_end <= 7680) {
      apaTA = 2;
    } else if (fTriggerActivity[i].channel_start >= 4160 && fTriggerActivity[i].channel_end <= 4640) {
      apaTA = 3;
    } else if (fTriggerActivity[i].channel_start >= 9280 && fTriggerActivity[i].channel_end <= 9760) {
      apaTA = 4;
    } else {
      // do nothing
    }
    //std::cout << "TA in APA " << apaTA << "; test APA = " << test_apaTA << "\n";
   */ 
    fapaTA.push_back(apaTA);
    
    // Fill tree
    fTATree -> Fill();
  }

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
    
    // Get ROP ID (ReadOut Plane ID)
    //auto rop = wireReadout.ChannelToROP(fChannelID);
    //fROP_ID = rop.ROP;
    //auto tpcid = wireReadout.ROPtoTPCs(rop);
    //int apa = rop.TPCset;
    /*int apa = 0;
    for (const auto &t : tpcid) {
      //std::cout << "TA in TPC " << t.TPC << "\n";
      if (t.TPC == 0 || t.TPC == 1) {
        apa = 1;
      } else if (t.TPC == 2 || t.TPC == 3) {
        apa = 3;
      } else if (t.TPC == 4 || t.TPC == 5) {
        apa = 2;
      } else if (t.TPC == 6 || t.TPC == 7) {
        apa = 4;
      }
    }*/
     
    // Fill tree
    //fTPTree -> Fill();
    int apa = 0;
    auto rop = wireReadout.ChannelToROP(fChannelID);
    auto tpc = rop.parentID().TPCset;
    //auto first_channel_rop = wireReadout.FirstChannelInROP(rop);
    //auto n_channels_rop = wireReadout.Nchannels(rop);
    apa = tpc;
    //std::cout << "TP in APA " << apa << "\n";
    //int rop_index = rop.deepestIndex();
    //std::cout << "rop: " << rop_index << ", 1st ch: " << first_channel_rop << ", n chan: " << n_channels_rop << "\n";
    //if (fChannelID >= 2080 && fChannelID <= 2559) {
    /*if (fChannelID <= 2560) {
      apa = 1;
    } else if (fChannelID >= 7200 && fChannelID <= 7680) {
      apa = 2;
    } else if (fChannelID >= 4160 && fChannelID <= 4640) {
      apa = 3;
    } else if (fChannelID >= 9280 && fChannelID <= 9760) {
      apa = 4;
    } else {
      // do nothing
    }*/
    

    // Determine the time window index (0 to 9).
    int windowIndex = fTime_peak / 20000;
    //if (windowIndex >= 10) windowIndex = 9;  // safeguard for times at the upper edge
    if (windowIndex >= 9) continue;  // safeguard for times at the upper edge and only accept 9 windows per readout

    switch(apa) {
      case 1:
        fAPA1Window_timepeak[windowIndex].push_back(fTime_peak); 
        fAPA1Window_channelid[windowIndex].push_back(fChannelID); 
        fAPA1Window_adcintegral[windowIndex].push_back(fADC_integral); 
        fAPA1Window_tot[windowIndex].push_back(fTime_over_threshold); 
        fAPA1Window_adcpeak[windowIndex].push_back(fADC_peak);
        break;
      case 2:
        fAPA2Window_timepeak[windowIndex].push_back(fTime_peak); 
        fAPA2Window_channelid[windowIndex].push_back(fChannelID); 
        fAPA2Window_adcintegral[windowIndex].push_back(fADC_integral); 
        fAPA2Window_tot[windowIndex].push_back(fTime_over_threshold); 
        fAPA2Window_adcpeak[windowIndex].push_back(fADC_peak);
        break;
      case 3:
        fAPA3Window_timepeak[windowIndex].push_back(fTime_peak); 
        fAPA3Window_channelid[windowIndex].push_back(fChannelID); 
        fAPA3Window_adcintegral[windowIndex].push_back(fADC_integral); 
        fAPA3Window_tot[windowIndex].push_back(fTime_over_threshold); 
        fAPA3Window_adcpeak[windowIndex].push_back(fADC_peak);
        break;
      case 4:
        fAPA4Window_timepeak[windowIndex].push_back(fTime_peak); 
        fAPA4Window_channelid[windowIndex].push_back(fChannelID); 
        fAPA4Window_adcintegral[windowIndex].push_back(fADC_integral); 
        fAPA4Window_tot[windowIndex].push_back(fTime_over_threshold); 
        fAPA4Window_adcpeak[windowIndex].push_back(fADC_peak);
        break;
      default: 
        break;
    }

    // Fill neutrino window if defined
    if (doNuWindow && fTime_peak >= nuWindowStart && fTime_peak <= nuWindowEnd && apa == fAPA) {
      std::cout << "Neutrino in APA " << apa << "\n";
      fNuWindow_timepeak[0].push_back(fTime_peak);
      fNuWindow_channelid[0].push_back(fChannelID);
      fNuWindow_adcintegral[0].push_back(fADC_integral);
      fNuWindow_tot[0].push_back(fTime_over_threshold); 
      fNuWindow_adcpeak[0].push_back(fADC_peak);
    }

    if (fTA > 0) {
      for(long unsigned int ta=0; ta < fTriggerActivity.size(); ta++) {
        if (doNuWindow && fTime_peak >= fTriggerActivity[ta].time_start && fTime_peak <= fTriggerActivity[ta].time_end && apa == fapaTA[ta]) {
          fTAWindow_timepeak[ta].push_back(fTime_peak);
          fTAWindow_channelid[ta].push_back(fChannelID);
          fTAWindow_adcintegral[ta].push_back(fADC_integral);
          fTAWindow_tot[ta].push_back(fTime_over_threshold); 
          fTAWindow_adcpeak[ta].push_back(fADC_peak);
        }
      }
    }
  }

  fTPWindowAPA1Tree -> Fill();
  fTPWindowAPA2Tree -> Fill();
  fTPWindowAPA3Tree -> Fill();
  fTPWindowAPA4Tree -> Fill();

  if (doNuWindow && neutrinoMC) {
    fTPNuWindowTree -> Fill(); 
  } else if (neutrinoMC) {
    fEventIterator--;
  }

  if (fTA > 0) {
    fTAWindowTree -> Fill(); 
  }
  
  // Fill TA tree
  for(long unsigned int i=0; i < fTriggerActivity.size(); i++)
  {
    auto rop = wireReadout.ChannelToROP(fTriggerActivity[i].channel_start);
    fTAROP_ID = rop.ROP;
    fChannel_start_TA = fTriggerActivity[i].channel_start;
    fChannel_end_TA = fTriggerActivity[i].channel_end;
    fChannel_peak_TA = fTriggerActivity[i].channel_peak;
    fTime_start_TA = fTriggerActivity[i].time_start;
    fTime_end_TA = fTriggerActivity[i].time_end;
    fTime_peak_TA = fTriggerActivity[i].time_peak;
    fTime_activity = fTriggerActivity[i].time_activity;
    fADC_integral_TA = fTriggerActivity[i].adc_integral;
    fADC_peak_TA = fTriggerActivity[i].adc_peak;
    fAlgorithm_TA = static_cast<int>(fTriggerActivity[i].algorithm);

    // Fill tree
    fTATree -> Fill();
  }

}

void duneana::SmallTriggerTPCInfoDisplay::beginJob()
{
  fEventIterator = 0;
  // Make our handle to the TFileService
  art::ServiceHandle<art::TFileService> tfs;
  // The TTrees
  fTPWindowAPA1Tree = tfs->make<TTree>("TPWindowAPA1Tree", "time windows in apa 1");
  fTPWindowAPA2Tree = tfs->make<TTree>("TPWindowAPA2Tree", "time windows in apa 2");
  fTPWindowAPA3Tree = tfs->make<TTree>("TPWindowAPA3Tree", "time windows in apa 3");
  fTPWindowAPA4Tree = tfs->make<TTree>("TPWindowAPA4Tree", "time windows in apa 4");
  
  fTPNuWindowTree = tfs->make<TTree>("TPNuWindowTree", "nu time windows");
  fTAWindowTree = tfs->make<TTree>("TAWindowTree", "TA time windows in");
  fTATree = tfs->make<TTree>("TATree", "DAQ trigger activity maker tree");
  fMCTruthTree = tfs->make<TTree>("GenieTruth", "GENIE Output Tree");

  // Add branches to TTree
  fMCTruthTree -> Branch("eventID", &fEventID);
  fMCTruthTree -> Branch( "EventIterator" , &fEventIterator, "EventIterator/I"  );
  fMCTruthTree -> Branch("SubRun", &fSubRun, "SubRun/I");
  fMCTruthTree -> Branch("Run", &fRun, "Run/I");
  fMCTruthTree -> Branch("TPCID", &fTPCID);
  fMCTruthTree -> Branch("APA", &fAPA);
  fMCTruthTree -> Branch("DriftDistance", &fDriftDistance);
  fMCTruthTree -> Branch("VertexTime", &fVertexTime);
  fMCTruthTree -> Branch("VertexTick", &fVertexTick);
  fMCTruthTree -> Branch("VertexTick_new_evd", &fVertexTick_new_evd);
  fMCTruthTree -> Branch("VertexTick_new_weight", &fVertexTick_new_weight);
  fMCTruthTree -> Branch("E", &fE, "E/D");
  fMCTruthTree -> Branch("nuVertexX", &fnuVertexX, "nuVertexX/D");
  fMCTruthTree -> Branch("nuVertexY", &fnuVertexY, "nuVertexY/D");
  fMCTruthTree -> Branch("nuVertexZ", &fnuVertexZ, "nuVertexZ/D");
  
  ////////////////////////////////////////
  // fTriggerPrimitive tree information //
  ////////////////////////////////////////
  fTPWindowAPA1Tree -> Branch( "EventIterator" , &fEventIterator, "EventIterator/I"  );
  fTPWindowAPA1Tree -> Branch( "APA1Window_timepeak", &fAPA1Window_timepeak);
  fTPWindowAPA1Tree -> Branch( "APA1Window_channelid", &fAPA1Window_channelid);
  fTPWindowAPA1Tree -> Branch( "APA1Window_adcintegral", &fAPA1Window_adcintegral);
  fTPWindowAPA1Tree -> Branch( "APA1Window_tot", &fAPA1Window_tot);
  fTPWindowAPA1Tree -> Branch( "APA1Window_adcpeak", &fAPA1Window_adcpeak);

  fTPWindowAPA2Tree -> Branch( "EventIterator" , &fEventIterator, "EventIterator/I"  );
  fTPWindowAPA2Tree -> Branch( "APA2Window_timepeak", &fAPA2Window_timepeak);
  fTPWindowAPA2Tree -> Branch( "APA2Window_channelid", &fAPA2Window_channelid);
  fTPWindowAPA2Tree -> Branch( "APA2Window_adcintegral", &fAPA2Window_adcintegral);
  fTPWindowAPA2Tree -> Branch( "APA2Window_tot", &fAPA2Window_tot);
  fTPWindowAPA2Tree -> Branch( "APA2Window_adcpeak", &fAPA2Window_adcpeak);
  
  fTPWindowAPA3Tree -> Branch( "EventIterator" , &fEventIterator, "EventIterator/I"  );
  fTPWindowAPA3Tree -> Branch( "APA3Window_timepeak", &fAPA3Window_timepeak);
  fTPWindowAPA3Tree -> Branch( "APA3Window_channelid", &fAPA3Window_channelid);
  fTPWindowAPA3Tree -> Branch( "APA3Window_adcintegral", &fAPA3Window_adcintegral);
  fTPWindowAPA3Tree -> Branch( "APA3Window_tot", &fAPA3Window_tot);
  fTPWindowAPA3Tree -> Branch( "APA3Window_adcpeak", &fAPA3Window_adcpeak);
  
  fTPWindowAPA4Tree -> Branch( "EventIterator" , &fEventIterator, "EventIterator/I"  );
  fTPWindowAPA4Tree -> Branch( "APA4Window_timepeak", &fAPA4Window_timepeak);
  fTPWindowAPA4Tree -> Branch( "APA4Window_channelid", &fAPA4Window_channelid);
  fTPWindowAPA4Tree -> Branch( "APA4Window_adcintegral", &fAPA4Window_adcintegral);
  fTPWindowAPA4Tree -> Branch( "APA4Window_tot", &fAPA4Window_tot);
  fTPWindowAPA4Tree -> Branch( "APA4Window_adcpeak", &fAPA4Window_adcpeak);
  
  fTPNuWindowTree -> Branch( "EventIterator" , &fEventIterator, "EventIterator/I"  );
  fTPNuWindowTree -> Branch( "APA" , &fAPA, "APA/I"  );
  fTPNuWindowTree -> Branch( "TA" , &fTA, "TA/I"  );
  fTPNuWindowTree -> Branch( "NuWindow_timepeak", &fNuWindow_timepeak);
  fTPNuWindowTree -> Branch( "NuWindow_channelid", &fNuWindow_channelid);
  fTPNuWindowTree -> Branch( "NuWindow_adcintegral", &fNuWindow_adcintegral);
  fTPNuWindowTree -> Branch( "NuWindow_tot", &fNuWindow_tot);
  fTPNuWindowTree -> Branch( "NuWindow_adcpeak", &fNuWindow_adcpeak);
  
  fTAWindowTree -> Branch( "EventIterator" , &fEventIterator, "EventIterator/I"  );
  fTAWindowTree -> Branch( "APATA" , &fapaTA );
  fTAWindowTree -> Branch( "TA" , &fTA, "TA/I"  );
  fTAWindowTree -> Branch( "TAWindow_timepeak", &fTAWindow_timepeak);
  fTAWindowTree -> Branch( "TAWindow_channelid", &fTAWindow_channelid);
  fTAWindowTree -> Branch( "TAWindow_adcintegral", &fTAWindow_adcintegral);
  fTAWindowTree -> Branch( "TAWindow_tot", &fTAWindow_tot);
  fTAWindowTree -> Branch( "TAWindow_adcpeak", &fTAWindow_adcpeak);

  ////////////////////////////////////////
  // fTriggerActivity tree information //
  ///////////////////////////////////////
  // General event information
  fTATree -> Branch( "Event" , &fEventID, "Event/I"  );
  fTATree -> Branch( "EventIterator" , &fEventIterator, "EventIterator/I"  );
  fTATree -> Branch( "Run"   , &fRun    , "Run/I"    );
  fTATree -> Branch( "SubRun", &fSubRun , "SubRun/I" );
  // Trigger activity information
  fTATree -> Branch( "Channel_start" , &fChannel_start_TA);
  fTATree -> Branch( "Channel_end" , &fChannel_end_TA);
  fTATree -> Branch( "Channel_peak"  , &fChannel_peak_TA);
  fTATree -> Branch( "ROP"  , &fTAROP_ID);
  fTATree -> Branch( "Time_start" , &fTime_start_TA);
  fTATree -> Branch( "Time_end" , &fTime_end_TA);
  fTATree -> Branch( "Time_peak" , &fTime_peak_TA);
  fTATree -> Branch( "Time_activity" , &fTime_activity);
  fTATree -> Branch( "ADC_integral" , &fADC_integral_TA);
  fTATree -> Branch( "ADC_peak" , &fADC_peak_TA);
  fTATree -> Branch( "DetID" , &fDetID);
  fTATree -> Branch( "Type" , &fType);
  fTATree -> Branch( "Algorithm" , &fAlgorithm_TA);

}

DEFINE_ART_MODULE(duneana::SmallTriggerTPCInfoDisplay)
