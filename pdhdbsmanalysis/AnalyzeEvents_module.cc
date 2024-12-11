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
#include "messagefacility/MessageLogger/MessageLogger.h"

// additional Framework includes
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"


#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"

// ROOT includes
#include "TTree.h"

#include <fstream>




namespace test {
  class AnalyzeEvents;
}
// pdhdkeepupstage2 | pandora................ | ....................... | std::vector<recob::Cluster>................................................................................ | .1131
// pdhdkeepupstage2 | pandora................ | ....................... | std::vector<recob::PFParticle>............................................................................. | ..578
// pdhdkeepupstage2 | pandoraShower.......... | ....................... | std::vector<recob::Shower>................................................................................. | ..260

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

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.
  TTree* fTrackTree;
  TTree* fShowerTree;
  int n_failed_tracks = 0;
  int n_failed_showers = 0;
  int n_success_tracks = 0;
  int n_success_showers = 0;
  // Tree variables
  std::vector<unsigned int> fTrackEventID;
  std::vector<double> fTrackVertexX;
  std::vector<double> fTrackVertexY;
  std::vector<double> fTrackVertexZ;
  std::vector<double> fTrackDirectionX;
  std::vector<double> fTrackDirectionY;
  std::vector<double> fTrackDirectionZ;
  std::vector<double> fTrackLength;

  std::vector<unsigned int> fShowerEventID;
  std::vector<double> fShowerDirectionX;
  std::vector<double> fShowerDirectionY;
  std::vector<double> fShowerDirectionZ;
  std::vector<double> fShowerStartX;
  std::vector<double> fShowerStartY;
  std::vector<double> fShowerStartZ;
  std::vector<double> fShowerLength;



  std::vector<recob::Track> fTrackVector;
  std::vector<recob::Shower> fShowerVector;
  // define input labels
  std::string fTrackLabel;
  std::string fShowerLabel;
};

test::AnalyzeEvents::AnalyzeEvents(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  // More initializers here.
  fTrackLabel(p.get<std::string>("TrackLabel")),
  fShowerLabel(p.get<std::string>("ShowerLabel"))
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void test::AnalyzeEvents::analyze(art::Event const& e)
{
  // Implementation of required member function here.
  // Set the event ID
  unsigned int EventID = e.id().event();
  // if (EventID != 1331) {
  //   return;
  // }
  std::cout << "Event ID: " << EventID << std::endl;

  art::Handle<std::vector<recob::Track>> trackHandle = e.getHandle<std::vector<recob::Track>>(fTrackLabel);
  if (trackHandle.isValid()) {
    std::cout << "Number of tracks: " << trackHandle->size() << std::endl;
    // Loop over tracks
    for (size_t trackIndex = 0; trackIndex < trackHandle->size(); trackIndex++) {
      std::cout << "Track vertex position: " << trackHandle->at(trackIndex).Vertex().X() << ", " << trackHandle->at(trackIndex).Vertex().Y() << ", " << trackHandle->at(trackIndex).Vertex().Z() << std::endl;
      std::cout << "Track direction: " << trackHandle->at(trackIndex).StartDirection().X() << ", " << trackHandle->at(trackIndex).StartDirection().Y() << ", " << trackHandle->at(trackIndex).StartDirection().Z() << std::endl;
      std::cout << "Track length: " << trackHandle->at(trackIndex).Length() << std::endl;
      fTrackVertexX.push_back(trackHandle->at(trackIndex).Vertex().X());
      fTrackVertexY.push_back(trackHandle->at(trackIndex).Vertex().Y());
      fTrackVertexZ.push_back(trackHandle->at(trackIndex).Vertex().Z());
      fTrackDirectionX.push_back(trackHandle->at(trackIndex).StartDirection().X());
      fTrackDirectionY.push_back(trackHandle->at(trackIndex).StartDirection().Y());
      fTrackDirectionZ.push_back(trackHandle->at(trackIndex).StartDirection().Z());
      fTrackLength.push_back(trackHandle->at(trackIndex).Length());
      fTrackEventID.push_back(EventID);
      std::cout << "Event ID: " << fTrackEventID.at(fTrackEventID.size() - 1) << std::endl;
      std::cout << "Event ID size: " << fTrackEventID.size() << std::endl;
    }
    n_success_tracks++;
  }
    else {
    std::cout << "Track handle is not valid" << std::endl;
    n_failed_tracks++;
  }

  art::Handle<std::vector<recob::Shower>> showerHandle = e.getHandle<std::vector<recob::Shower>>(fShowerLabel);
  if (showerHandle.isValid()) {
    std::cout << "Number of showers: " << showerHandle->size() << std::endl;
    // Loop over showers
    for (size_t showerIndex = 0; showerIndex < showerHandle->size(); showerIndex++) {
      std::cout << "Fillin shower tree" << std::endl;
      std::cout << "Shower direction: " << showerHandle->at(showerIndex).Direction().X() << ", " << showerHandle->at(showerIndex).Direction().Y() << ", " << showerHandle->at(showerIndex).Direction().Z() << std::endl;
      std::cout << "Shower start position: " << showerHandle->at(showerIndex).ShowerStart().X() << ", " << showerHandle->at(showerIndex).ShowerStart().Y() << ", " << showerHandle->at(showerIndex).ShowerStart().Z() << std::endl;
      std::cout << "Shower length: " << showerHandle->at(showerIndex).Length() << std::endl;
      fShowerDirectionX.push_back(showerHandle->at(showerIndex).Direction().X());
      fShowerDirectionY.push_back(showerHandle->at(showerIndex).Direction().Y());
      fShowerDirectionZ.push_back(showerHandle->at(showerIndex).Direction().Z());
      fShowerStartX.push_back(showerHandle->at(showerIndex).ShowerStart().X());
      fShowerStartY.push_back(showerHandle->at(showerIndex).ShowerStart().Y());
      fShowerStartZ.push_back(showerHandle->at(showerIndex).ShowerStart().Z());
      fShowerLength.push_back(showerHandle->at(showerIndex).Length());
      fShowerEventID.push_back(EventID);
      std::cout << "Event ID: " << fShowerEventID.at(fShowerEventID.size() - 1) << std::endl;
      std::cout << "Event ID size: " << fShowerEventID.size() << std::endl;
    }
    n_success_showers++;
  }
    else {
    std::cout << "Shower handle is not valid" << std::endl;
    n_failed_showers++;
  }


  // Fill tree
  fTrackTree->Fill();
  fShowerTree->Fill();


  // Clear vectors
  fTrackVertexX.clear();
  fTrackVertexY.clear();
  fTrackVertexZ.clear();
  fTrackDirectionX.clear();
  fTrackDirectionY.clear();
  fTrackDirectionZ.clear();
  fTrackLength.clear();
  fTrackEventID.clear();

  fShowerDirectionX.clear();
  fShowerDirectionY.clear();
  fShowerDirectionZ.clear();
  fShowerStartX.clear();
  fShowerStartY.clear();
  fShowerStartZ.clear();
  fShowerLength.clear();
  fShowerEventID.clear();

  std::cout << "Number of failed tracks: " << n_failed_tracks << std::endl;
  std::cout << "Number of failed showers: " << n_failed_showers << std::endl;
  std::cout << "Number of successful tracks: " << n_success_tracks << std::endl;
  std::cout << "Number of successful showers: " << n_success_showers << std::endl;
  // write these numbers to a file


}

void test::AnalyzeEvents::beginJob()
{
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;

  // Create tree for tracks
  fTrackTree = tfs->make<TTree>("trackTree", "Track Output TTree");
  fTrackTree->Branch("TrackEventID", &fTrackEventID);
  fTrackTree->Branch("TrackVertexX", &fTrackVertexX);
  fTrackTree->Branch("TrackVertexY", &fTrackVertexY);
  fTrackTree->Branch("TrackVertexZ", &fTrackVertexZ);
  fTrackTree->Branch("TrackDirectionX", &fTrackDirectionX);
  fTrackTree->Branch("TrackDirectionY", &fTrackDirectionY);
  fTrackTree->Branch("TrackDirectionZ", &fTrackDirectionZ);
  fTrackTree->Branch("TrackLength", &fTrackLength);

  // Create tree for showers
  fShowerTree = tfs->make<TTree>("showerTree", "Shower Output TTree");
  fShowerTree->Branch("ShowerEventID", &fShowerEventID);
  fShowerTree->Branch("ShowerDirectionX", &fShowerDirectionX);
  fShowerTree->Branch("ShowerDirectionY", &fShowerDirectionY);
  fShowerTree->Branch("ShowerDirectionZ", &fShowerDirectionZ);
  fShowerTree->Branch("ShowerStartX", &fShowerStartX);
  fShowerTree->Branch("ShowerStartY", &fShowerStartY);
  fShowerTree->Branch("ShowerStartZ", &fShowerStartZ);
  fShowerTree->Branch("ShowerLength", &fShowerLength);
  

}

void test::AnalyzeEvents::endJob()
{
  // Implementation of optional member function here.
  std::ofstream myfile;
  myfile.open("output.txt", std::ios_base::app);
  myfile << n_failed_showers << " " << n_failed_tracks << " " << n_success_showers << " " << n_success_tracks << std::endl;
  myfile.close();
}

DEFINE_ART_MODULE(test::AnalyzeEvents)
