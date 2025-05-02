////////////////////////////////////////////////////////////////////////
// Class:       AnalyzeReco
// Plugin Type: analyzer (Unknown Unknown)
// File:        AnalyzeReco_module.cc
//
// Generated at Sun Jan 19 10:46:24 2025 by Ciaran Hasnip using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

#include "larcoreobj/SummaryData/POTSummary.h"

#include "art_root_io/TFileService.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "detdataformats/trigger/Types.hpp"

// ROOT includes
#include <TTree.h>

using timestamp_t = dunedaq::trgdataformats::timestamp_t;

namespace ana {
  class AnalyzeReco;
}


class ana::AnalyzeReco : public art::EDAnalyzer {
public:
  explicit AnalyzeReco(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  AnalyzeReco(AnalyzeReco const&) = delete;
  AnalyzeReco(AnalyzeReco&&) = delete;
  AnalyzeReco& operator=(AnalyzeReco const&) = delete;
  AnalyzeReco& operator=(AnalyzeReco&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginSubRun(art::SubRun const& subRun) override;
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.

  TTree *fAnaTree;

  std::string fSPSBeamData;
  std::ifstream fInputData;

  unsigned int fEventID;
  unsigned int fRun;
  unsigned int fSubRun;

  uint64_t fPoT_threshold;
  uint64_t fFirstTimeStamp;
  uint64_t fLastTimeStamp;
  double fPoT;
  double fTotalPoT;

  uint64_t sum_PoT;
  timestamp_t run_livetime;
  timestamp_t run_spill_livetime;
  std::vector<std::pair<timestamp_t, uint64_t>> vSpillClockPoT; // Vector to store the spill clock times and PoT values
};


ana::AnalyzeReco::AnalyzeReco(fhicl::ParameterSet const& p)
  : EDAnalyzer{p} ,
    fSPSBeamData(p.get<std::string>("sps_beamdata")) ,
    fPoT_threshold(p.get<uint64_t>("PoT_threshold"))
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void ana::AnalyzeReco::analyze(art::Event const& e)
{
  // Implementation of required member function here.
  fEventID = e.id().event();
  fRun = e.run();
  fSubRun = e.subRun();

  uint64_t timeHigh_ns = e.time().timeHigh() * 1e9;
  uint64_t timeLow_ns = e.time().timeLow();
  uint64_t fEventTimeStamp = (timeHigh_ns + timeLow_ns) * 1e-6;
  if (fEventTimeStamp < fFirstTimeStamp) fFirstTimeStamp = fEventTimeStamp;
  if (fEventTimeStamp > fLastTimeStamp) fLastTimeStamp = fEventTimeStamp;

  // Fill trr
  fAnaTree->Fill();
}

void ana::AnalyzeReco::beginSubRun(art::SubRun const& subRun) {
 
  fPoT = 0;
  // Only neutrino MC has or needs this
  /*if (!evt.isRealData()) {
    const auto potSummaryHandle = subRun.getValidHandle<sumdata::POTSummary>("generator");
    const auto &potSummary = *potSummaryHandle;
    fPoT = potSummary.totgoodpot;
 
    fTotalPoT += fPoT;
    std::cout << "POTSummary content: totpot = " << potSummary.totpot 
      << ", totgoodpot = " << potSummary.totgoodpot << std::endl;
      
    // Fill the TTree with the current subrun's POT information
    //fSubrunTree->Fill();
  }*/
}

void ana::AnalyzeReco::beginJob()
{

  fTotalPoT = 0;

  art::ServiceHandle<art::TFileService> tfs;
  fAnaTree = tfs->make<TTree>("ana", "Analysis Tree");

  fAnaTree->Branch("eventID", &fEventID); 


  // Implementation of optional member function here.
  std::cout << "SPS beam data file: " << fSPSBeamData << "\n";
  fInputData.open(fSPSBeamData);

  if (!fInputData.good()) {
    throw std::runtime_error("Input csv file " + fSPSBeamData + " cannot be read.");
  }

  std::string line;
  // Read and discard the header line
  if (std::getline(fInputData, line)) {
    // Optionally, you can store the header if needed
    // std::vector<std::string> header;
    // std::stringstream headerStream(line);
    // std::string cell;
    // while (std::getline(headerStream, cell, ',')) {
    //     header.push_back(cell);
    // }
  }

  sum_PoT = 0; // running sum of PoT in SPS data file

  while (std::getline(fInputData, line)) {
    std::vector<std::string> data;
    std::stringstream lineStream(line);
    std::string cell;
    while (std::getline(lineStream, cell, ',')) {
      data.push_back(cell);
    }
        
    try {
      timestamp_t clock = static_cast<timestamp_t>(std::stod(data[0])*1e3); // Convert the string in ms to a double and then to a timestamp_t
      uint64_t PoT = static_cast<uint64_t>(std::stoull(data[1])); // Convert the string to an unsigned long long and then to a uint64_t
      if (PoT >= fPoT_threshold) {
        vSpillClockPoT.push_back(std::make_pair(clock, PoT)); // Store the spill clock time and PoT value
        sum_PoT += PoT;
      }
    } catch (const std::invalid_argument& e) {
            // Handle the case where the string is not a valid double
    } catch (const std::out_of_range& e) {
            // Handle the case where the double is out of range
    }
        
  }

  

  // Get total time between first spill and last spill to get length of run
  run_livetime = vSpillClockPoT.end()->first - vSpillClockPoT.begin()->first;
  // Get the amount of time the spill has been on
  run_spill_livetime = static_cast<timestamp_t>(vSpillClockPoT.size() * 4785); 

  std::cout << "In " << fSPSBeamData << " there are " << vSpillClockPoT.size() << " SPS beam spills.\n\n";
  std::cout << "In " << fSPSBeamData << " total PoT = " << sum_PoT << ", run livetime = " << run_livetime << 
    " and amount of time with spill ON = " << run_spill_livetime << "\n\n";

  // Set these defaults to be changed in analyze function
  fFirstTimeStamp = vSpillClockPoT.end()->first;
  fLastTimeStamp = vSpillClockPoT.begin()->first;
}

void ana::AnalyzeReco::endJob()
{
  // Implementation of optional member function here.

  uint64_t sum_event_PoT(0);
  int number_of_spills(0);
  for (size_t spill = 0; spill < vSpillClockPoT.size(); ++spill) {
    if (vSpillClockPoT[spill].first > fFirstTimeStamp &&
        vSpillClockPoT[spill].first < fLastTimeStamp) {
      sum_event_PoT += vSpillClockPoT[spill].second;
      number_of_spills++;
    }
  }
  timestamp_t run_event_livetime = fLastTimeStamp - fFirstTimeStamp;
  timestamp_t run_event_spill_livetime = static_cast<timestamp_t>(number_of_spills * 4785);
  
  std::cout << "In " << fSPSBeamData << " there are " << vSpillClockPoT.size() << " SPS beam spills.\n\n";
  std::cout << "In " << fSPSBeamData << " total PoT = " << sum_PoT << ", run livetime = " << run_livetime << 
    " and amount of time with spill ON = " << run_spill_livetime << "\n\n";
 
  std::cout << "When looking at first/last event total PoT = " << sum_event_PoT << ", run livetime = " << run_event_livetime <<  
    " and time with spill on = " << run_event_spill_livetime << "\n\n";
}

DEFINE_ART_MODULE(ana::AnalyzeReco)
