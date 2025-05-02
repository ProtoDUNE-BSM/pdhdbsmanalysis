////////////////////////////////////////////////////////////////////////
// Class:       PDHDRawDataImage
// Plugin Type: analyzer (Unknown Unknown)
// File:        PDHDRawDataImage_module.cc
//
// Generated at Thu Feb  6 10:26:27 2025 by Ciaran Hasnip using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

// ART includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
// LArSoft includes
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/WireReadout.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/WireReadoutStandardGeom.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/RawData/RawDigit.h"
// DUNE includes
#include "detdataformats/trigger/TriggerObjectOverlay.hpp"
#include "detdataformats/trigger/TriggerPrimitive.hpp"
#include "detdataformats/trigger/TriggerActivityData.hpp"
// Additional Framework includes
#include "art_root_io/TFileService.h"
// ROOT includes
#include "TTree.h"

namespace ana {
  class PDHDRawDataImage;
}

constexpr int kMaxNumberCh = 10240;
constexpr int kMaxTicks= 8000;

class ana::PDHDRawDataImage : public art::EDAnalyzer {
public:
  explicit PDHDRawDataImage(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PDHDRawDataImage(PDHDRawDataImage const&) = delete;
  PDHDRawDataImage(PDHDRawDataImage&&) = delete;
  PDHDRawDataImage& operator=(PDHDRawDataImage const&) = delete;
  PDHDRawDataImage& operator=(PDHDRawDataImage&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

  void reset();

private:

  // Declare member data here.
  art::InputTag fRawLabel;
  art::InputTag fTALabel;

  geo::GeometryCore const* fGeometryService; ///< pointer to Geometry provider

  TTree *fTree;

  int fRun;
  int fSubrun;
  int fEvent;

  int fTA;
  int fW_plane[kMaxNumberCh];
  int fW_ch[kMaxNumberCh];
  float fW_signal[kMaxNumberCh][kMaxTicks];

};


ana::PDHDRawDataImage::PDHDRawDataImage(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
  , fRawLabel(p.get< art::InputTag >("RawLabel"))
  , fTALabel(p.get< art::InputTag >("TALabel"))
  // More initializers here.
{
  fGeometryService = lar::providerFrom<geo::Geometry>();
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void ana::PDHDRawDataImage::analyze(art::Event const& e)
{
  // Implementation of required member function here.
  fRun = e.run();
  fSubrun = e.subRun();
  fEvent = e.id().event();
  
  auto const& wireReadoutGeom = art::ServiceHandle<geo::WireReadout const>()->Get();

  auto const& rawdigits = *(e.getValidHandle<std::vector<raw::RawDigit>>(fRawLabel));

  art::Handle<std::vector<dunedaq::trgdataformats::TriggerActivityData>> taHandle;
  if (!e.getByLabel(fTALabel, taHandle)) {
      fTA = 0;
  } else {
    if (taHandle->size() == 0) {
      fTA = 0;
    } else {
      std::cout << ">>> Found " << taHandle->size() << " TAs in Event " << fEvent << std::endl;
      fTA = 1;
    }
  }

  // Don't bother continuing if the event does not trigger
  if (fTA != 1) return;
  
  unsigned int idx = 0;
  for (raw::RawDigit const& rawdigit: rawdigits) {
    std::cout << "Compression: " << rawdigit.Compression() << std::endl;
    if (idx >= kMaxNumberCh) {
      mf::LogWarning("RawDigitAna") << "Number of channels in rawdigit exceeds kMaxNumberCh. Only saving the first " << kMaxNumberCh << " channels.";
      break;
    }
    unsigned int jdx = 0;
    std::cout << "Number of ticks = " << rawdigit.NADC() << std::endl;
    for ( int adc = 0; adc < int(rawdigit.NADC()); ++adc){
      if (jdx >= kMaxTicks) {
        mf::LogWarning("RawDigitAna") << "Number of ticks in rawdigit exceeds kMaxTicks. Only saving the first " << kMaxTicks << " ticks.";
        break;
      }
      fW_signal[idx][jdx]= rawdigit.ADC(adc);
      ++ jdx;
    }
    fW_plane[idx] = wireReadoutGeom.ChannelToROP(rawdigit.Channel()).ROP;
    fW_ch[idx] = rawdigit.Channel();
    ++ idx;
  }
  fTree->Fill();
}

void ana::PDHDRawDataImage::beginJob()
{
  // Implementation of optional member function here.
  reset();
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("rawdigitTree","Tree with rawdigit info");
  fTree->Branch("event", &fEvent);
  fTree->Branch("subrun", &fSubrun);
  fTree->Branch("run", &fRun);
  fTree->Branch("TA", &fTA);
  fTree->Branch("w_ch", &fW_ch, Form("w_ch[%d]/I", kMaxNumberCh));
  fTree->Branch("w_plane", &fW_plane, Form("w_plane[%d]/I", kMaxNumberCh));
  fTree->Branch("w_signal", &fW_signal, Form("w_signal[%d][%d]/F", kMaxNumberCh, kMaxTicks));
}

void ana::PDHDRawDataImage::endJob()
{
  // Implementation of optional member function here.
}

void ana::PDHDRawDataImage::reset()
{

  for ( unsigned int i=0; i<kMaxNumberCh; ++i){
    fW_plane[i] = -1;   
    fW_ch[i] = -1;   
    for ( unsigned int j=0; j<kMaxTicks; ++j){
      fW_signal[i][j] = 0;
    }
  }
}

DEFINE_ART_MODULE(ana::PDHDRawDataImage)
