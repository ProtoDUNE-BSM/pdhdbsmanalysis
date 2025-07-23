// FindNeutrinos_module.cc
//
// Refactored for clarity and performance:
//  - Modularized ADC pedestal subtraction
//  - Cached Handles & associations
//  - Guard clauses & early exits
//  - Range-based loops & std algos
//  - MessageLogger in place of std::cout
////////////////////////////////////////////////////////////////////////

// FindNeutrinos_module.cc — Corrected Include Section
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


#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"


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
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

// ROOT includes
#include "TTree.h"

#include <vector>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <cmath>


namespace NeutrinoAna {
class FindNeutrinos;
}  // namespace NeutrinoAna

class NeutrinoAna::FindNeutrinos : public art::EDAnalyzer {
public:
    explicit FindNeutrinos(fhicl::ParameterSet const& cfg);
    void analyze(art::Event const& event) override;
    void beginSubRun(art::SubRun const& subRun) override;
    void beginJob() override;
    void endJob() override {}

    // Disable copy / move
    FindNeutrinos(FindNeutrinos const&) = delete;
    FindNeutrinos(FindNeutrinos&&) = delete;
    FindNeutrinos& operator=(FindNeutrinos const&) = delete;
    FindNeutrinos& operator=(FindNeutrinos&&) = delete;

private:
    // --------------------------------------------------------------------------
    //  Utility helpers
    // --------------------------------------------------------------------------
    std::vector<std::vector<short>>
    subtractPedestal(std::vector<art::Ptr<raw::RawDigit>> const& rawDigits) const;

    bool detectGroundShake(std::vector<art::Ptr<raw::RawDigit>> const& rawDigits) const;

    int sumADCInRange(std::vector<art::Ptr<raw::RawDigit>> const& rawDigits,
                      int startTick,
                      int endTick) const;

    std::vector<double> computeDaughterSummary(
        art::Ptr<recob::PFParticle> const& pfParticle,
        art::Event const& event) const;
    std::vector<double> getInformation(
        art::Ptr<recob::PFParticle> const& pfParticle,
        art::Event const& event) const;

    int computeTrueOriginIdentifier(
        art::Ptr<recob::PFParticle> const& pfParticle,
        art::Event const& event) const;
    double computeZFromChannel(int channel) const;
    double computeEnergyDepositedInRange(
        std::vector<art::Ptr<recob::Hit>> const& hitPtrs,
        art::Event const& event,
        double startZ,
        double endZ) const;

    // --------------------------------------------------------------------------
    //  Reconstruction algorithms
    // --------------------------------------------------------------------------
    dune::NeutrinoAngularRecoAlg fNeutrinoAngularAlg;
    dune::NeutrinoEnergyRecoAlg  fNeutrinoEnergyAlg;

    // --------------------------------------------------------------------------
    //  Configuration labels
    // --------------------------------------------------------------------------
    std::string fTrackLabel;
    std::string fShowerLabel;
    std::string fSliceLabel;
    std::string fPFParticleLabel;
    std::string fVertexLabel;
    std::string fHitLabel;
    std::string fCalorimetryLabel;
    std::string fTriggerActivityLabel;
    std::string fMCTruthLabel;
    bool        fEnableTruth;

    // --------------------------------------------------------------------------
    //  TTree pointers
    // --------------------------------------------------------------------------
    TTree* fRecoTree      {nullptr};
    TTree* fTruthTree     {nullptr};
    TTree* fAggregateTree {nullptr};

    // --------------------------------------------------------------------------
    //  Reco-tree branches
    // --------------------------------------------------------------------------
    unsigned int fRecoEventID                 {0};
    double       fRecoVertexX                 {0.};
    double       fRecoVertexY                 {0.};
    double       fRecoVertexZ                 {0.};
    int          fRecoNeutrinoPdgCode         {0};
    double       fRecoNeutrinoEnergy          {0.};
    double       fRecoDirectionX              {0.};
    double       fRecoDirectionY              {0.};
    double       fRecoDirectionZ              {0.};
    int          fRecoNumberOfHits            {0};
    int          fRecoNumberOfPFParticles     {0};
    int          fRecoTrueOriginID            {0};
    int          fRecoEventSequenceNumber     {0};
    int          fRecoPassSelectionCriterion  {0};

    // --------------------------------------------------------------------------
    //  Truth-tree branches
    // --------------------------------------------------------------------------
    unsigned int fTruthEventID                {0};
    double       fTruthVertexX                {0.};
    double       fTruthVertexY                {0.};
    double       fTruthVertexZ                {0.};
    double       fTruthNeutrinoEnergy         {0.};
    double       fTruthNeutrinoMomentumX      {0.};
    double       fTruthNeutrinoMomentumY      {0.};
    double       fTruthNeutrinoMomentumZ      {0.};
    int          fTruthNeutrinoPdgCode        {0};
    int          fTruthNeutrinoMotherPdgCode  {0};
    double       fTruthSubrunTotalPOT         {0.};
    double       fTruthSubrunGoodPOT          {0.};
    double       fAccumulatedPOT              {0.};
    int          fTruthTriggerActivityFlag    {0};
    int          fTruthEventSequenceNumber    {0};

    // --------------------------------------------------------------------------
    //  Aggregate-tree branches
    // --------------------------------------------------------------------------
    unsigned int fAggregateEventID                {0};
    double       fAggregateVertexX                {0.};
    double       fAggregateVertexY                {0.};
    double       fAggregateVertexZ                {0.};
    double       fAggregateReconstructedEnergy    {0.};
    double       fAggregateDirectionX             {0.};
    double       fAggregateDirectionY             {0.};
    double       fAggregateDirectionZ             {0.};
    double       fAggregateDirectionX2            {0.};
    double       fAggregateDirectionY2            {0.};
    double       fAggregateDirectionZ2            {0.};
    int          fAggregateNumberOfHits           {0};
    int          fAggregateNumberOfPFParticles    {0};
    int          fAggregateTrueOriginID           {0};
    int          fAggregateEventSequenceNumber    {0};
    int          fAggregatePassSelectionCriterion {0};
    int          fAggregateSpillStatusFlag        {0};
    double       fAggregateEventTimestamp         {0.};
    int          fAggregateTotalNumberOfHits      {0};
    int          fAggregateTriggerCandidateCount  {0};
    int          fAggregateGroundShakeCount       {0};
    int          fAggregateSumOfLastADCTicks      {0};
    int          fAggregateSumOfTriggeredADCTicks {0};
    int          fAggregatePassSecondSelectionCriterion {0};
    double       fAggregateEnergyDepositedInFirst10cm {0.};
    double       fAggregateEnergyDepositedInSecond10cm {0.};
    double       fAggregateEnergyDepositedInThird10cm {0.};
    double       fAggregateEnergyDepositedInFourth10cm {0.};
    double       fAggregateEnergyDepositedInFifth10cm {0.};
    double       fAggregateEnergyDepositedInSixth10cm {0.};
    double       fAggregateEnergyDepositedInSeventh10cm {0.};   
    double       fAggregateEnergyDepositedInEighth10cm {0.};
    double       fAggregateEnergyDepositedInNinth10cm {0.};
    double       fAggregateEnergyDepositedInTenth10cm {0.};
    double       fAggregateEnergyDepositedInEleventh10cm {0.};
    double       fAggregateEnergyDepositedInTwelfth10cm {0.};
    double       fAggregateEnergyDepositedInThirteenth10cm {0.};
    double       fAggregateEnergyDepositedInFourteenth10cm {0.};
    double       fAggregateEnergyDepositedInFifteenth10cm {0.};
    double       fAggregateEnergyDepositedInFirst10cmBefore {0.};
    double       fAggregateEnergyDepositedInSecond10cmBefore {0.};

    // --------------------------------------------------------------------------
    //  Internal run / sub-run bookkeeping
    // --------------------------------------------------------------------------
    double fCurrentSubrunTotalPOT   {0.};
    double fCurrentSubrunGoodPOT    {0.};
    double fTotalAccumulatedPOT     {0.};
    int    fTriggerActivityPresent  {0};
    int    fGlobalEventCounter      {0};
};

// ============================================================================
//  Constructor
// ============================================================================
NeutrinoAna::FindNeutrinos::FindNeutrinos(fhicl::ParameterSet const& cfg)
  : EDAnalyzer{cfg}
  , fNeutrinoAngularAlg{cfg,
        "pandoraTrack", "pandoraShower", "pandora",
        "wclsdatahd", "pandoraTrack", "pandoraShower", "pandora"}
  , fNeutrinoEnergyAlg{cfg,
        "pandoraTrack", "pandoraShower", "pandora",
        "wclsdatahd", "pandoraTrack", "pandoraShower", "pandora"}
  , fTrackLabel           {cfg.get<std::string>("TrackLabel")}
  , fShowerLabel          {cfg.get<std::string>("ShowerLabel")}
  , fSliceLabel           {cfg.get<std::string>("SliceLabel")}
  , fPFParticleLabel      {cfg.get<std::string>("PFParticleLabel")}
  , fVertexLabel          {cfg.get<std::string>("VertexLabel")}
  , fHitLabel             {cfg.get<std::string>("HitsModuleLabel")}
  , fCalorimetryLabel     {cfg.get<std::string>("CalorimetryLabel")}
  , fTriggerActivityLabel {cfg.get<std::string>("TALabel")}
  , fMCTruthLabel         {cfg.get<std::string>("MCTruthLabel")}
  , fEnableTruth          {cfg.get<bool>("GetTruth")}
{}

// ============================================================================
//  beginJob – create TTree branches
// ============================================================================
void NeutrinoAna::FindNeutrinos::beginJob()
{
    art::ServiceHandle<art::TFileService> tfs;

    fRecoTree = tfs->make<TTree>("tree_reco", "Reconstructed neutrino candidates");
    fRecoTree->Branch("eventID",                &fRecoEventID);
    fRecoTree->Branch("vertexX",                &fRecoVertexX);
    fRecoTree->Branch("vertexY",                &fRecoVertexY);
    fRecoTree->Branch("vertexZ",                &fRecoVertexZ);
    fRecoTree->Branch("pdgCode",                &fRecoNeutrinoPdgCode);
    fRecoTree->Branch("energy",                 &fRecoNeutrinoEnergy);
    fRecoTree->Branch("directionX",             &fRecoDirectionX);
    fRecoTree->Branch("directionY",             &fRecoDirectionY);
    fRecoTree->Branch("directionZ",             &fRecoDirectionZ);
    fRecoTree->Branch("numberOfHits",           &fRecoNumberOfHits);
    fRecoTree->Branch("numberOfPFParticles",    &fRecoNumberOfPFParticles);
    fRecoTree->Branch("trueOriginID",           &fRecoTrueOriginID);
    fRecoTree->Branch("eventSequenceNumber",    &fRecoEventSequenceNumber);
    fRecoTree->Branch("passSelectionCriterion", &fRecoPassSelectionCriterion);

    fTruthTree = tfs->make<TTree>("tree_truth", "True neutrino information");
    fTruthTree->Branch("eventID",             &fTruthEventID);
    fTruthTree->Branch("vertexX",             &fTruthVertexX);
    fTruthTree->Branch("vertexY",             &fTruthVertexY);
    fTruthTree->Branch("vertexZ",             &fTruthVertexZ);
    fTruthTree->Branch("energy",              &fTruthNeutrinoEnergy);
    fTruthTree->Branch("momentumX",           &fTruthNeutrinoMomentumX);
    fTruthTree->Branch("momentumY",           &fTruthNeutrinoMomentumY);
    fTruthTree->Branch("momentumZ",           &fTruthNeutrinoMomentumZ);
    fTruthTree->Branch("pdgCode",             &fTruthNeutrinoPdgCode);
    fTruthTree->Branch("motherPdgCode",       &fTruthNeutrinoMotherPdgCode);
    fTruthTree->Branch("subrunTotalPOT",      &fTruthSubrunTotalPOT);
    fTruthTree->Branch("subrunGoodPOT",       &fTruthSubrunGoodPOT);
    fTruthTree->Branch("accumulatedPOT",      &fAccumulatedPOT);
    fTruthTree->Branch("triggerActivityFlag", &fTruthTriggerActivityFlag);
    fTruthTree->Branch("eventSequenceNumber", &fTruthEventSequenceNumber);

    fAggregateTree = tfs->make<TTree>("tree_aggregate", "One entry per event");
    fAggregateTree->Branch("eventID",                 &fAggregateEventID);
    fAggregateTree->Branch("vertexX",                 &fAggregateVertexX);
    fAggregateTree->Branch("vertexY",                 &fAggregateVertexY);
    fAggregateTree->Branch("vertexZ",                 &fAggregateVertexZ);
    fAggregateTree->Branch("reconstructedEnergy",     &fAggregateReconstructedEnergy);
    fAggregateTree->Branch("directionX",              &fAggregateDirectionX);
    fAggregateTree->Branch("directionY",              &fAggregateDirectionY);
    fAggregateTree->Branch("directionZ",              &fAggregateDirectionZ);
    fAggregateTree->Branch("directionX2",              &fAggregateDirectionX2);
    fAggregateTree->Branch("directionY2",              &fAggregateDirectionY2);
    fAggregateTree->Branch("directionZ2",              &fAggregateDirectionZ2);
    fAggregateTree->Branch("numberOfHits",            &fAggregateNumberOfHits);
    fAggregateTree->Branch("numberOfPFParticles",     &fAggregateNumberOfPFParticles);
    fAggregateTree->Branch("trueOriginID",            &fAggregateTrueOriginID);
    fAggregateTree->Branch("eventSequenceNumber",     &fAggregateEventSequenceNumber);
    fAggregateTree->Branch("passSelectionCriterion",  &fAggregatePassSelectionCriterion);
    fAggregateTree->Branch("spillStatusFlag",         &fAggregateSpillStatusFlag);
    fAggregateTree->Branch("eventTimestamp",          &fAggregateEventTimestamp);
    fAggregateTree->Branch("totalNumberOfHits",       &fAggregateTotalNumberOfHits);
    fAggregateTree->Branch("triggerCandidateCount",   &fAggregateTriggerCandidateCount);
    fAggregateTree->Branch("groundShakeCount",        &fAggregateGroundShakeCount);
    fAggregateTree->Branch("sumOfLastADCTicks",       &fAggregateSumOfLastADCTicks);
    fAggregateTree->Branch("sumOfTriggeredADCTicks",  &fAggregateSumOfTriggeredADCTicks);
    fAggregateTree->Branch("passSecondSelectionCriterion", &fAggregatePassSecondSelectionCriterion);
    fAggregateTree->Branch("energyDepositedInFirst10cm", &fAggregateEnergyDepositedInFirst10cm);
    fAggregateTree->Branch("energyDepositedInSecond10cm", &fAggregateEnergyDepositedInSecond10cm);
    fAggregateTree->Branch("energyDepositedInThird10cm", &fAggregateEnergyDepositedInThird10cm);
    fAggregateTree->Branch("energyDepositedInFourth10cm", &fAggregateEnergyDepositedInFourth10cm);
    fAggregateTree->Branch("energyDepositedInFifth10cm", &fAggregateEnergyDepositedInFifth10cm);
    fAggregateTree->Branch("energyDepositedInSixth10cm", &fAggregateEnergyDepositedInSixth10cm);
    fAggregateTree->Branch("energyDepositedInSeventh10cm", &fAggregateEnergyDepositedInSeventh10cm);
    fAggregateTree->Branch("energyDepositedInEighth10cm", &fAggregateEnergyDepositedInEighth10cm);
    fAggregateTree->Branch("energyDepositedInNinth10cm", &fAggregateEnergyDepositedInNinth10cm);
    fAggregateTree->Branch("energyDepositedInTenth10cm", &fAggregateEnergyDepositedInTenth10cm);
    fAggregateTree->Branch("energyDepositedInEleventh10cm", &fAggregateEnergyDepositedInEleventh10cm);
    fAggregateTree->Branch("energyDepositedInTwelfth10cm", &fAggregateEnergyDepositedInTwelfth10cm);
    fAggregateTree->Branch("energyDepositedInThirteenth10cm", &fAggregateEnergyDepositedInThirteenth10cm);
    fAggregateTree->Branch("energyDepositedInFourteenth10cm", &fAggregateEnergyDepositedInFourteenth10cm);
    fAggregateTree->Branch("energyDepositedInFifteenth10cm", &fAggregateEnergyDepositedInFifteenth10cm);
    fAggregateTree->Branch("energyDepositedInFirst10cmBefore", &fAggregateEnergyDepositedInFirst10cmBefore);
    fAggregateTree->Branch("energyDepositedInSecond10cmBefore", &fAggregateEnergyDepositedInSecond10cmBefore);

}

// ============================================================================
//  beginSubRun – cache POT information
// ============================================================================
void NeutrinoAna::FindNeutrinos::beginSubRun(art::SubRun const& subRun)
{
    if (!fEnableTruth) return;

    auto const& potSummaryHandle =
        subRun.getValidHandle<sumdata::POTSummary>("generator");

    fCurrentSubrunTotalPOT   = potSummaryHandle->totpot;
    fCurrentSubrunGoodPOT    = potSummaryHandle->totgoodpot;
    fTotalAccumulatedPOT    += fCurrentSubrunTotalPOT;
}

// ============================================================================
//  Pedestal-subtraction helper
// ============================================================================
std::vector<std::vector<short>>
NeutrinoAna::FindNeutrinos::subtractPedestal(
    std::vector<art::Ptr<raw::RawDigit>> const& rawDigits) const
{
    std::vector<std::vector<short>> pedSubtracted;
    pedSubtracted.reserve(rawDigits.size());

    for (auto const& rdPtr : rawDigits) {
        std::vector<short> samples = rdPtr->ADCs();
        if (samples.empty()) continue;

        // Build simple frequency map → mode
        std::unordered_map<short, int> frequency;
        frequency.reserve(512);
        for (short s : samples)
            ++frequency[s];

        short mode = std::max_element(
            frequency.begin(), frequency.end(),
            [](auto const& a, auto const& b) {
                return a.second < b.second;
            })->first;

        for (short& s : samples)
            s -= mode;

        pedSubtracted.push_back(std::move(samples));
    }

    return pedSubtracted;
}

// ============================================================================
//  Ground-shake detection
// ============================================================================
bool NeutrinoAna::FindNeutrinos::detectGroundShake(
    std::vector<art::Ptr<raw::RawDigit>> const& rawDigits) const
{
    auto waveforms = subtractPedestal(rawDigits);
    if (waveforms.empty()) return false;

    const int nChannels = waveforms.size();
    const int nSamples  = waveforms[0].size();
    const int threshold = 10;

    for (int tick = 0; tick < nSamples; ++tick) {
        int channelsAbove = 0;
        for (auto const& wf : waveforms) {
            if (wf[tick] > threshold)
                ++channelsAbove;
        }
        if (channelsAbove > 0.85 * nChannels)
            return true;
    }

    return false;
}

// ============================================================================
//  ADC integration helper
// ============================================================================
int NeutrinoAna::FindNeutrinos::sumADCInRange(
    std::vector<art::Ptr<raw::RawDigit>> const& rawDigits,
    int startTick,
    int endTick) const
{
    auto waveforms = subtractPedestal(rawDigits);
    int sum = 0;

    for (auto const& wf : waveforms) {
        const int size = wf.size();
        for (int tick = startTick; tick < endTick && tick < size; ++tick) {
            sum += wf[tick];
        }
    }

    return sum;
}

// ============================================================================
//  Daughter summary (total hits, total PFP count)
// ============================================================================
std::vector<double> NeutrinoAna::FindNeutrinos::computeDaughterSummary(
    art::Ptr<recob::PFParticle> const& pfParticle,
    art::Event const&                 event) const
{
    std::vector<double> summary(2, 0.0);  // hits, PFPs
    std::vector<art::Ptr<recob::PFParticle>> daughters = dune_ana::DUNEAnaPFParticleUtils::GetChildParticles(pfParticle, event, fPFParticleLabel);

    if (daughters.empty()) {
        // Hits
        auto hits = dune_ana::DUNEAnaPFParticleUtils::GetHits(
            pfParticle, event, fPFParticleLabel);
        summary[0] = static_cast<double>(hits.size());
        summary[1] = 1.0;  // this single PFParticle
    }
    else {
        // Hits and PFPs
        for (auto const& daughter : daughters) {
            auto daughter_info = computeDaughterSummary(daughter, event);
            summary[0] += daughter_info[0];
            summary[1] += daughter_info[1];
        }
    }

    return summary;
}

// ============================================================================
//  GetInformations (energy, direction)
// ============================================================================
std::vector<double> NeutrinoAna::FindNeutrinos::getInformation(
    art::Ptr<recob::PFParticle> const& pfParticle,
    art::Event const&                 event) const
{
    std::vector<double> summary(4, 0.0);  // E, dirX, dirY, dirZ
    // Track or shower?
    if (dune_ana::DUNEAnaPFParticleUtils::IsTrack(
            pfParticle, event, fPFParticleLabel, fTrackLabel)) {
        auto trackPtr = dune_ana::DUNEAnaPFParticleUtils::GetTrack(
            pfParticle, event, fPFParticleLabel, fTrackLabel);
        auto caloPtr = dune_ana::DUNEAnaTrackUtils::GetCalorimetry(
            trackPtr, event, fTrackLabel, fCalorimetryLabel);

        summary[0] = caloPtr->KineticEnergy();
        summary[1] = trackPtr->VertexDirection().X();
        summary[2] = trackPtr->VertexDirection().Y();
        summary[3] = trackPtr->VertexDirection().Z();

    } else if (dune_ana::DUNEAnaPFParticleUtils::IsShower(
                   pfParticle, event, fPFParticleLabel, fShowerLabel)) {
        auto showerPtr = dune_ana::DUNEAnaPFParticleUtils::GetShower(
            pfParticle, event, fPFParticleLabel, fShowerLabel);
        art::Handle<std::vector<recob::Shower>> showerHandle;
        event.getByLabel(fShowerLabel, showerHandle);

        art::FindManyP<anab::Calorimetry> caloAssn(
            showerHandle, event, "pandoraShowercalonosce");
        auto caloVec = caloAssn.at(showerPtr.key());

        summary[0] = caloVec.front()->KineticEnergy();
        summary[1] = showerPtr->Direction().X();
        summary[2] = showerPtr->Direction().Y();
        summary[3] = showerPtr->Direction().Z();
    }

    return summary;
}

// ============================================================================
//  True-origin identifier helper
// ============================================================================
int NeutrinoAna::FindNeutrinos::computeTrueOriginIdentifier(
    art::Ptr<recob::PFParticle> const& pfParticle,
    art::Event const&                 event) const
{
    auto const& clockData =
        art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(event);

    TruthMatchUtils::G4ID g4ID =
        TruthMatchUtils::TrueParticleIDFromTotalRecoHits(
            clockData,
            dune_ana::DUNEAnaPFParticleUtils::GetHits(
                pfParticle, event, fPFParticleLabel),
            /* rollUpUnsavedIDs */ true);

    art::ServiceHandle<cheat::ParticleInventoryService> pis;
    auto const* mcPtr = pis->TrackIdToParticle_P(g4ID);
    if (!mcPtr) return 0;

    return pis->TrackIdToMCTruth(mcPtr->TrackId()).Origin();
}

// ============================================================================
//  Get the energy deposited in a specific range
// ============================================================================
double NeutrinoAna::FindNeutrinos::computeEnergyDepositedInRange(
    std::vector<art::Ptr<recob::Hit>> const& hitPtrs,
    art::Event const&             event,
    double startZ, double endZ) const
{
    double energyDeposited = 0.0;

    // Loop over the hits and sum the energy deposited in the specified range
    for (auto const& hit : hitPtrs) {
        if (!hit) continue;
        // Include only collection hits
        if (hit->View() != 2) continue;

        // Check if the hit is within the first 10 cm
        auto channelID = hit->Channel();
        double z_pos = computeZFromChannel(channelID);
        if (z_pos < endZ && z_pos > startZ) {
            energyDeposited += hit->Integral();
        }
    }

    return energyDeposited;
}


// ============================================================================
//  Compute Z from channel
// ============================================================================
double NeutrinoAna::FindNeutrinos::computeZFromChannel(int channel) const
{
    // Constants (in cm)
    constexpr double apa_length_in_cm = 230.0;
    constexpr double wire_pitch_in_cm_collection = 0.479;
    // constexpr double wire_pitch_in_cm_induction_diagonal = 0.4669;
    // constexpr double apa_angle_deg = 90.0 - 35.7;
    // constexpr double offset_between_apa_in_cm = 2.4;
    // constexpr double apa_height_in_cm = 598.4;
    // constexpr double time_tick_in_cm = 0.0805;
    // constexpr double apa_width_in_cm = 4.7;
    // constexpr double backtracker_error_margin = 4.0;

    // // Calculate induction wire pitch (projected)
    // double apa_angle_rad = apa_angle_deg * M_PI / 180.0;
    // double wire_pitch_in_cm_induction = wire_pitch_in_cm_induction_diagonal / std::sin(apa_angle_rad);
    // double apa_angular_coeff = std::tan(apa_angle_rad);

    // first, get the APA number
    int apa_number = channel / 2560;
    // second, subtract the inductions
    int channel_number = (channel % 2560) - 1600;
    // third, consider both sides
    channel_number = channel_number%480;

    double offset = (apa_number < 2) ? 0.0 : apa_length_in_cm;

    double z_pos = channel_number * wire_pitch_in_cm_collection + offset;

    return z_pos;
}

// ============================================================================
//  Main Analysis
// ============================================================================
void NeutrinoAna::FindNeutrinos::analyze(art::Event const& event)
{
    ++fGlobalEventCounter;

    // ------------------------------------------------------------------------
    //  Cache handles once
    // ------------------------------------------------------------------------
    auto hitHandle =
        event.getValidHandle<std::vector<recob::Hit>>(fHitLabel);
    auto sliceHandle =
        event.getValidHandle<std::vector<recob::Slice>>(fSliceLabel);
    auto pfpHandle =
        event.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleLabel);

    std::vector<art::Ptr<recob::Slice>> slicePtrs;
    art::fill_ptr_vector(slicePtrs, sliceHandle);

    std::vector<art::Ptr<recob::Hit>> allHitPtrs;
    art::fill_ptr_vector(allHitPtrs, hitHandle);

    // ------------------------------------------------------------------------
    //  TriggerActivity flag
    // ------------------------------------------------------------------------
    art::Handle<std::vector<dunedaq::trgdataformats::TriggerActivityData>> taHandle;
    event.getByLabel(fTriggerActivityLabel, taHandle);
    fTriggerActivityPresent =
        (taHandle.isValid() && !taHandle->empty()) ? 1 : 0;
    std::cout << "TriggerActivity present: "
              << (fTriggerActivityPresent ? "Yes" : "No") << std::endl;

    // ------------------------------------------------------------------------
    //  Truth branch fill (if enabled)
    // ------------------------------------------------------------------------
    if (fEnableTruth) {
        fTruthEventID              = event.id().event();
        fTruthEventSequenceNumber  = fGlobalEventCounter;
        fTruthTriggerActivityFlag  = fTriggerActivityPresent;
        fTruthSubrunTotalPOT       = fCurrentSubrunTotalPOT;
        fTruthSubrunGoodPOT        = fCurrentSubrunGoodPOT;
        fAccumulatedPOT            = fTotalAccumulatedPOT;

        art::Handle<std::vector<simb::MCTruth>> truthHandle;
        if (event.getByLabel(fMCTruthLabel, truthHandle)) {
            for (auto const& mcTruth : *truthHandle) {
                if (!mcTruth.NeutrinoSet()) continue;
                auto const& nu = mcTruth.GetNeutrino().Nu();

                fTruthVertexX           = nu.Vx();
                fTruthVertexY           = nu.Vy();
                fTruthVertexZ           = nu.Vz();
                fTruthNeutrinoEnergy    = nu.E();
                fTruthNeutrinoMomentumX = nu.Px();
                fTruthNeutrinoMomentumY = nu.Py();
                fTruthNeutrinoMomentumZ = nu.Pz();
                fTruthNeutrinoPdgCode   = nu.PdgCode();
                fTruthNeutrinoMotherPdgCode = nu.Mother();
            }
        }
        if (fTruthTree) fTruthTree->Fill();
    }

    // ------------------------------------------------------------------------
    //  Find most-energetic slice
    // ------------------------------------------------------------------------
    art::FindManyP<recob::PFParticle> sliceToPFP(
        sliceHandle, event, fSliceLabel);

    double maxRecoE  = 0.;
    art::Ptr<recob::Slice> bestSlice;
    for (auto const& slicePtr : slicePtrs) {
        auto eOut  = fNeutrinoEnergyAlg.CalculateNeutrinoEnergy(event, slicePtr, true);
        double thisE = eOut.fNuLorentzVector.E();
        if (thisE > maxRecoE) {
            maxRecoE  = thisE;
            bestSlice = slicePtr;
        }
    }
    if (!bestSlice) {
        mf::LogWarning("FindNeutrinos")
            << "No slice found in event " << event.id().event();
        return;
    }

    // ------------------------------------------------------------------------
    //  PFParticles in the best slice
    // ------------------------------------------------------------------------
    auto pfpVec = sliceToPFP.at(bestSlice.key());

    // Reset aggregate counters
    fAggregateEventID                = event.id().event();
    fAggregateEventSequenceNumber    = fGlobalEventCounter;
    fAggregateVertexX                = -1.;
    fAggregateVertexY                = -1.;
    fAggregateVertexZ                = -1.;
    fAggregateReconstructedEnergy    = -1.;
    fAggregateDirectionX             = -1.;
    fAggregateDirectionY             = -1.;
    fAggregateDirectionZ             = -1.;
    fAggregateDirectionX2            = -1.;
    fAggregateDirectionY2            = -1.;
    fAggregateDirectionZ2            = -1.;
    fAggregateNumberOfHits           = -1;
    fAggregateNumberOfPFParticles    = -1;
    fAggregateTrueOriginID           = -1;
    fAggregatePassSelectionCriterion = -1;
    fAggregateSpillStatusFlag        = -1;
    fAggregateTriggerCandidateCount  = -1;
    fAggregateGroundShakeCount       = -1;
    fAggregateSumOfLastADCTicks      = -1;
    fAggregateSumOfTriggeredADCTicks = -1;
    fAggregatePassSecondSelectionCriterion = -1;
    fAggregateEnergyDepositedInFirst10cm = -1;
    fAggregateEnergyDepositedInSecond10cm = -1;
    fAggregateEnergyDepositedInThird10cm = -1;
    fAggregateEnergyDepositedInFourth10cm = -1;
    fAggregateEnergyDepositedInFifth10cm = -1;
    fAggregateEnergyDepositedInSixth10cm = -1;
    fAggregateEnergyDepositedInSeventh10cm = -1;
    fAggregateEnergyDepositedInEighth10cm = -1;
    fAggregateEnergyDepositedInNinth10cm = -1;
    fAggregateEnergyDepositedInTenth10cm = -1;
    fAggregateEnergyDepositedInEleventh10cm = -1;
    fAggregateEnergyDepositedInTwelfth10cm = -1;
    fAggregateEnergyDepositedInThirteenth10cm = -1;
    fAggregateEnergyDepositedInFourteenth10cm = -1;
    fAggregateEnergyDepositedInFifteenth10cm = -1;
    fAggregateEnergyDepositedInFirst10cmBefore = -1;
    fAggregateEnergyDepositedInSecond10cmBefore = -1;

    // ------------------------------------------------------------------------
    //  Loop PFParticles
    // ------------------------------------------------------------------------
    art::FindManyP<recob::Vertex> pfpToVertex(
        pfpHandle, event, fVertexLabel);

    for (auto const& pfp : pfpVec) {
        auto vertexPtrs = pfpToVertex.at(pfp.key());
        double vtxX = -999., vtxY = -999., vtxZ = -999.;
        if (!vertexPtrs.empty()) {
            vtxX = vertexPtrs.front()->position().X();
            vtxY = vertexPtrs.front()->position().Y();
            vtxZ = vertexPtrs.front()->position().Z();
        }

        auto summary       = getInformation(pfp, event);
        auto daughterSummary = computeDaughterSummary(pfp, event);
        int  trueOriginID  = fEnableTruth
                             ? computeTrueOriginIdentifier(pfp, event)
                             : -1;

        // Individual PFParticle → reco tree
        fRecoEventID                = event.id().event();
        fRecoVertexX                = vtxX;
        fRecoVertexY                = vtxY;
        fRecoVertexZ                = vtxZ;
        fRecoNeutrinoPdgCode        = pfp->PdgCode();
        fRecoNeutrinoEnergy         = summary[0];
        fRecoDirectionX             = summary[1];
        fRecoDirectionY             = summary[2];
        fRecoDirectionZ             = summary[3];
        fRecoNumberOfHits           = static_cast<int>(daughterSummary[0]);
        fRecoNumberOfPFParticles    = static_cast<int>(daughterSummary[1]);
        fRecoTrueOriginID           = trueOriginID;
        fRecoEventSequenceNumber    = fGlobalEventCounter;
        fRecoPassSelectionCriterion = 0;  // set below for primary nu

        // If this PFParticle is a primary neutrino candidate
        if (pfp->IsPrimary() &&
            (std::abs(pfp->PdgCode()) == 12 ||
             std::abs(pfp->PdgCode()) == 14 ||
             std::abs(pfp->PdgCode()) == 16)) {
            if (vtxY <= 550. && vtxZ >= 20. && daughterSummary[1] >= 4.0){
                fRecoPassSelectionCriterion = 1;
                fAggregatePassSelectionCriterion = 1;
            }
            else {
                fRecoPassSelectionCriterion = 0;
                fAggregatePassSelectionCriterion = 0;
            }

            fAggregateVertexX             = vtxX;
            fAggregateVertexY             = vtxY;
            fAggregateVertexZ             = vtxZ;
            fAggregateNumberOfPFParticles = static_cast<int>(daughterSummary[1]);
            fAggregateTrueOriginID        = trueOriginID;
            fAggregateNumberOfHits        = static_cast<int>(daughterSummary[0]);
            std::cout << "Found primary neutrino candidate in event" << std::endl;
        }
        fRecoTree->Fill();
    }

    // ------------------------------------------------------------------------
    //  Refine aggregate energy and direction with dedicated algs
    // ------------------------------------------------------------------------
    dune::Point_t seedVertex;
    seedVertex.SetCoordinates(
        fAggregateVertexX, fAggregateVertexY, fAggregateVertexZ);

    auto angOut = fNeutrinoAngularAlg.CalculateNeutrinoAngle(
        event, bestSlice, seedVertex);
    fAggregateDirectionX = angOut.fRecoDirection.x();
    fAggregateDirectionY = angOut.fRecoDirection.y();
    fAggregateDirectionZ = angOut.fRecoDirection.z();

    auto eOut = fNeutrinoEnergyAlg.CalculateNeutrinoEnergy(
        event, bestSlice, true);
    fAggregateReconstructedEnergy = eOut.fNuLorentzVector.E();

    // ------------------------------------------------------------------------
    //  Slice hit times (first & last)
    // ------------------------------------------------------------------------
    art::FindManyP<recob::Hit> sliceToHit(
        sliceHandle, event, fSliceLabel);
    auto sliceHits = sliceToHit.at(bestSlice.key());

    if (!sliceHits.empty()) {
        auto cmp = [](auto const& a, auto const& b) {
            return a->PeakTime() < b->PeakTime();
        };
        auto [minIt, maxIt] = std::minmax_element(
            sliceHits.begin(), sliceHits.end(), cmp);

        MF_LOG_DEBUG("FindNeutrinos")
            << "Slice hit-time window: "
            << (*minIt)->PeakTime() << " → " << (*maxIt)->PeakTime();
    }

    // // ------------------------------------------------------------------------
    // //  Get the max and min x,y,z coordinates of all the spacepoints
    // // ------------------------------------------------------------------------
    // double minX = 1e6, minY = 1e6, minZ = 1e6;
    // double maxX = -1e6, maxY = -1e6, maxZ = -1e6;

    // // Get PFParticles associated with slices
    // art::FindManyP<recob::PFParticle> sliceToPFParticle(
    //     sliceHandle, event, "pandora");
    // auto pfParticles = sliceToPFParticle.at(bestSlice.key());

    // // Get PFParticle handle for the SpacePoint association
    // art::Handle<std::vector<recob::PFParticle>> pfParticleHandle;
    // event.getByLabel("pandora", pfParticleHandle);

    // // Get SpacePoints associated with PFParticles
    // art::FindManyP<recob::SpacePoint> pfParticleToSpacePoint(
    //     pfParticleHandle, event, "pandora");

    // std::vector<art::Ptr<recob::SpacePoint>> allSpacePoints;

    // // Loop through PFParticles in the slice to collect all space points
    // if (!pfParticles.empty()) {
    //     for (const auto& pfp : pfParticles) {
    //         auto spacePoints = pfParticleToSpacePoint.at(pfp.key());
            
    //         // Add all space points from this PFParticle to our collection
    //         for (const auto& sp : spacePoints) {
    //             allSpacePoints.push_back(sp);
    //         }
    //     }
    // }

    // // Now calculate the bounding box from all collected space points
    // if (!allSpacePoints.empty()) {
    //     for (auto const& sp : allSpacePoints) {
    //         auto xyz = sp->XYZ();
    //         if (xyz[0] < minX) minX = xyz[0];
    //         if (xyz[1] < minY) minY = xyz[1];
    //         if (xyz[2] < minZ) minZ = xyz[2];
    //         if (xyz[0] > maxX) maxX = xyz[0];
    //         if (xyz[1] > maxY) maxY = xyz[1];
    //         if (xyz[2] > maxZ) maxZ = xyz[2];
    //     }
        
    //     std::cout << "Space points in event " << event.id().event() << ": "
    //             << allSpacePoints.size() << std::endl;
    //     std::cout << "X range: [" << minX << ", " << maxX << "]" << std::endl;
    //     std::cout << "Y range: [" << minY << ", " << maxY << "]" << std::endl;
    //     std::cout << "Z range: [" << minZ << ", " << maxZ << "]" << std::endl;
        
    // } else {
    //     mf::LogWarning("FindNeutrinos")
    //         << "No space points found in event " << event.id().event();
    //     return;
    // }

    // ------------------------------------------------------------------------
    //  Loop over the tracks not in the slice to see if there is a track parallel to the
    //  neutrino direction
    // ------------------------------------------------------------------------
    art::Handle<std::vector<recob::Track>> trackHandle;
    if (!event.getByLabel(fTrackLabel, trackHandle)) {
        mf::LogWarning("FindNeutrinos")
            << "No tracks found in event " << event.id().event();
        return;
    }
    std::vector<art::Ptr<recob::Track>> trackPtrs;
    art::fill_ptr_vector(trackPtrs, trackHandle);
    
    // get tracks ids in the slice
    std::set<int> sliceTracksIDs;
    for (auto const& PFParticle : pfpVec) {
        if (PFParticle->IsPrimary()) continue;  // skip primary neutrinos
        // check if PFParticle is a track
        if (!dune_ana::DUNEAnaPFParticleUtils::IsTrack(
                PFParticle, event, fPFParticleLabel, fTrackLabel)) continue;

        auto trackPtr = dune_ana::DUNEAnaPFParticleUtils::GetTrack(
            PFParticle, event, fPFParticleLabel, fTrackLabel);
        if (trackPtr.isNonnull()) {
            sliceTracksIDs.insert(trackPtr->ID());
        }
    }

    fAggregatePassSecondSelectionCriterion = 1;  // assume true until proven otherwise
    double maxTrackDot = 0.7;
    for (auto const& trackPtr : trackPtrs) {
        // Skip tracks that are part of the slice
        if (sliceTracksIDs.find(trackPtr->ID()) != sliceTracksIDs.end()) continue;
    
        if (trackPtr->Length() < 10.) continue;  // skip short tracks
        // Skip tracks if vertex z is bigger than neutrino vertex z
        if (trackPtr->Vertex().Z() > fAggregateVertexZ) continue;

        // Calculate the dot product with the neutrino direction
        // Normalize the track direction vector
        double tx = trackPtr->VertexDirection().X();
        // double ty = trackPtr->VertexDirection().Y();
        double tz = trackPtr->VertexDirection().Z();
        // double tnorm = std::sqrt(tx * tx + ty * ty + tz * tz);
        // if (tnorm == 0) continue; // skip invalid direction

        // double dot = (tx / tnorm) * fAggregateDirectionX +
        //          (ty / tnorm) * fAggregateDirectionY +
        //          (tz / tnorm) * fAggregateDirectionZ;
        
        double tnorm = std::sqrt(tx * tx + tz * tz);
        if (tnorm == 0) continue; // skip invalid direction

        double dot = (tx / tnorm) * fAggregateDirectionX +
                 (tz / tnorm) * fAggregateDirectionZ;
        
        if (std::abs(dot) > maxTrackDot) {
            // check how far the track is from the neutrino aggregate vertex
            double trackVtxX = trackPtr->Vertex().X();
            double trackVtxY = trackPtr->Vertex().Y();
            double trackVtxZ = trackPtr->Vertex().Z();
            // first, calculate the vector from the neutrino vertex to the track vertex
            double dx = trackVtxX - fAggregateVertexX;
            double dy = trackVtxY - fAggregateVertexY;
            double dz = trackVtxZ - fAggregateVertexZ;
            // then, calculate the vector product with the neutrino direction
            double vector_product_x = dy * fAggregateDirectionZ - dz * fAggregateDirectionY;
            double vector_product_y = dz * fAggregateDirectionX - dx * fAggregateDirectionZ;
            double vector_product_z = dx * fAggregateDirectionY - dy * fAggregateDirectionX;
            double distance = std::sqrt(vector_product_x * vector_product_x +
                                        vector_product_y * vector_product_y +
                                        vector_product_z * vector_product_z);
            distance /= std::sqrt(fAggregateDirectionX * fAggregateDirectionX +
                                fAggregateDirectionY * fAggregateDirectionY +
                                fAggregateDirectionZ * fAggregateDirectionZ);

            if (distance < 20.) { 
                fAggregatePassSecondSelectionCriterion = 0;
                std::cout << "Found a track parallel to the neutrino direction with dot product "
                          << dot << " and distance " << distance
                          << " from the neutrino vertex." << std::endl;
            }
        }
    }

    // ------------------------------------------------------------------------
    //  Compute the energy deposited in the first 10 cm of the track
    // ------------------------------------------------------------------------

    // get the hits associated with the best slice
    art::FindManyP<recob::Hit> hitAssn(sliceHandle, event, fSliceLabel);
    auto hitPtrs = hitAssn.at(bestSlice.key());

    fAggregateEnergyDepositedInFirst10cm = computeEnergyDepositedInRange(hitPtrs, event, fAggregateVertexZ, fAggregateVertexZ + 10.0);
    fAggregateEnergyDepositedInSecond10cm = computeEnergyDepositedInRange(hitPtrs, event, fAggregateVertexZ + 10.0, fAggregateVertexZ + 20.0);
    fAggregateEnergyDepositedInThird10cm = computeEnergyDepositedInRange(hitPtrs, event, fAggregateVertexZ + 20.0, fAggregateVertexZ + 30.0);
    fAggregateEnergyDepositedInFourth10cm = computeEnergyDepositedInRange(hitPtrs, event, fAggregateVertexZ + 30.0, fAggregateVertexZ + 40.0);
    fAggregateEnergyDepositedInFifth10cm = computeEnergyDepositedInRange(hitPtrs, event, fAggregateVertexZ + 40.0, fAggregateVertexZ + 50.0);
    fAggregateEnergyDepositedInSixth10cm = computeEnergyDepositedInRange(hitPtrs, event, fAggregateVertexZ + 50.0, fAggregateVertexZ + 60.0);
    fAggregateEnergyDepositedInSeventh10cm = computeEnergyDepositedInRange(hitPtrs, event, fAggregateVertexZ + 60.0, fAggregateVertexZ + 70.0);
    fAggregateEnergyDepositedInEighth10cm = computeEnergyDepositedInRange(hitPtrs, event, fAggregateVertexZ + 70.0, fAggregateVertexZ + 80.0);
    fAggregateEnergyDepositedInNinth10cm = computeEnergyDepositedInRange(hitPtrs, event, fAggregateVertexZ + 80.0, fAggregateVertexZ + 90.0);
    fAggregateEnergyDepositedInTenth10cm = computeEnergyDepositedInRange(hitPtrs, event, fAggregateVertexZ + 90.0, fAggregateVertexZ + 100.0);
    fAggregateEnergyDepositedInEleventh10cm = computeEnergyDepositedInRange(hitPtrs, event, fAggregateVertexZ + 100.0, fAggregateVertexZ + 110.0);
    fAggregateEnergyDepositedInTwelfth10cm = computeEnergyDepositedInRange(hitPtrs, event, fAggregateVertexZ + 110.0, fAggregateVertexZ + 120.0);
    fAggregateEnergyDepositedInThirteenth10cm = computeEnergyDepositedInRange(hitPtrs, event, fAggregateVertexZ + 120.0, fAggregateVertexZ + 130.0);
    fAggregateEnergyDepositedInFourteenth10cm = computeEnergyDepositedInRange(hitPtrs, event, fAggregateVertexZ + 130.0, fAggregateVertexZ + 140.0);
    fAggregateEnergyDepositedInFifteenth10cm = computeEnergyDepositedInRange(hitPtrs, event, fAggregateVertexZ + 140.0, fAggregateVertexZ + 150.0);
    // ------------------------------------------------------------------------
    //  Average time peak time of all hits in the slice
    // ------------------------------------------------------------------------
    double totalTime = 0.0;
    int hitCount = 0;
    int hitsInTPC0 = 0, hitsInTPC1 = 0;
    for (auto const& hitPtr : hitPtrs) {
        if (!hitPtr) continue;
        // Only consider collection view hits
        if (hitPtr->View() != 2) continue;
        totalTime += hitPtr->PeakTime();
        ++hitCount;
        if (hitPtr->Channel() / 2560 == 0 || hitPtr->Channel() / 2560 == 2) {
            ++hitsInTPC0;
        } else {
            ++hitsInTPC1;
        }
    }
    double averageTime = (hitCount > 0) ? totalTime / hitCount : -1.0;

    // get the hits in the same TPC as the best slice and close in time to the average time
    std::vector<art::Ptr<recob::Hit>> HitsCloseToSlice;
    int target_tpc = (hitsInTPC0 > hitsInTPC1) ? 0 : 1; // TPC with more hits in the slice
    for (auto const& hitPtr : allHitPtrs) {
        if (!hitPtr) continue;
        // Only consider collection view hits
        if (hitPtr->View() != 2) continue;
        // Check if the hit is close in time to the average time
        if (std::abs(hitPtr->PeakTime() - averageTime) > 32*100) continue;
        // Skip hits that are in the best slice
        if (std::find(hitPtrs.begin(), hitPtrs.end(), hitPtr) != hitPtrs.end()) continue;
        // Check if the hit is in the same TPC as the best slice
        int hit_tpc = (hitPtr->Channel() / 2560 == 0 || hitPtr->Channel() / 2560 == 2) ? 0 : 1;
        if (hit_tpc == target_tpc) {
            HitsCloseToSlice.push_back(hitPtr);
        }
    }
    fAggregateEnergyDepositedInFirst10cmBefore = computeEnergyDepositedInRange(HitsCloseToSlice, event, fAggregateVertexZ - 10.0, fAggregateVertexZ);
    fAggregateEnergyDepositedInSecond10cmBefore = computeEnergyDepositedInRange(HitsCloseToSlice, event, fAggregateVertexZ- 20.0, fAggregateVertexZ - 10.0);


    // ------------------------------------------------------------------------
    //  Trigger candidate count (non-ground-shake)
    // ------------------------------------------------------------------------
    art::Handle<std::vector<dunedaq::trgdataformats::TriggerCandidateData>> triggerCandHandle;
    if (event.getByLabel("triggerrawdecoder:daq", triggerCandHandle) && triggerCandHandle.isValid()) {
        for (auto const& tc : *triggerCandHandle) {
            if (tc.type != dunedaq::trgdataformats::TriggerCandidateData::Type::kADCSimpleWindow)
                ++fAggregateTriggerCandidateCount;
        }
    }

    // ------------------------------------------------------------------------
    //  Raw digits: split by APA & check for ground shakes
    // ------------------------------------------------------------------------
    std::vector<art::Ptr<raw::RawDigit>> rawDigitPtrs;
    art::Handle<std::vector<raw::RawDigit>> rawHandle;
    if (event.getByLabel("tpcrawdecoder:daq", rawHandle)) {
        art::fill_ptr_vector(rawDigitPtrs, rawHandle);

        std::vector<art::Ptr<raw::RawDigit>> apa[4];
        for (auto const& rd : rawDigitPtrs) {
            if (rd->Channel() % 2560 < 1600) continue;  // collection only
            int idx = rd->Channel() / 2560;
            if (idx >= 0 && idx < 4)
                apa[idx].push_back(rd);
        }

        for (int a = 0; a < 4; ++a) {
            if (detectGroundShake(apa[a]))
                ++fAggregateGroundShakeCount;

            int nSamples = apa[a].empty() ? 0 : apa[a][0]->Samples();
            fAggregateSumOfLastADCTicks +=
                sumADCInRange(apa[a], nSamples - 1500, nSamples);

            int startTick = (nSamples > 9000) ? 4100 : 100;
            fAggregateSumOfTriggeredADCTicks +=
                sumADCInRange(apa[a], startTick, startTick + 3000);
        }
    }

    // ------------------------------------------------------------------------
    //  Spill flag
    // ------------------------------------------------------------------------
    art::Handle<std::vector<bool>> spillHandle;
    if (event.getByLabel("spillflag:sps", spillHandle) && !spillHandle->empty())
        fAggregateSpillStatusFlag = (*spillHandle)[0];

    // ------------------------------------------------------------------------
    //  Event timestamp (timeHigh * 1 ns + timeLow ns)
    // ------------------------------------------------------------------------
    auto evtTime = event.time();
    fAggregateEventTimestamp =
        static_cast<double>(evtTime.timeHigh()) * 1e9 + evtTime.timeLow();

    // ------------------------------------------------------------------------
    //  Final direction normalisation
    // ------------------------------------------------------------------------
    double norm = std::sqrt(
        fAggregateDirectionX * fAggregateDirectionX +
        fAggregateDirectionY * fAggregateDirectionY +
        fAggregateDirectionZ * fAggregateDirectionZ);
    if (norm > 0.) {
        fAggregateDirectionX /= norm;
        fAggregateDirectionY /= norm;
        fAggregateDirectionZ /= norm;
    }
    // ------------------------------------------------------------------------
    //  Get the second direction
    // ------------------------------------------------------------------------
    fAggregateDirectionX2 = 0;
    fAggregateDirectionY2 = 0;
    fAggregateDirectionZ2 = 0;
    for (auto const& pfp : pfpVec) {
        auto summary = getInformation(pfp, event);
        fAggregateDirectionX2 += summary[1];
        fAggregateDirectionY2 += summary[2];
        fAggregateDirectionZ2 += summary[3];
    }
    // Normalize the second direction
    double norm2 = std::sqrt(
        fAggregateDirectionX2 * fAggregateDirectionX2 +
        fAggregateDirectionY2 * fAggregateDirectionY2 +
        fAggregateDirectionZ2 * fAggregateDirectionZ2);
    if (norm2 > 0.) {
        fAggregateDirectionX2 /= norm2;
        fAggregateDirectionY2 /= norm2;
        fAggregateDirectionZ2 /= norm2;
    }

    // Total number of hits in event
    fAggregateTotalNumberOfHits = static_cast<int>(allHitPtrs.size());

    // ------------------------------------------------------------------------
    //  Fill aggregate tree
    // ------------------------------------------------------------------------
    fAggregateTree->Fill();

        // Print aggregate neutrino candidate info if present
        if (fAggregateVertexX!=0 and fAggregateVertexX!=-1) // selected reco
        {
            std::cout << "Aggregate Neutrino Candidate Info:" << std::endl;
            std::cout << "Vertex: (" << fAggregateVertexX << ", " << fAggregateVertexY << ", " << fAggregateVertexZ << ")" << std::endl;
            std::cout << "Direction: (" << fAggregateDirectionX << ", " << fAggregateDirectionY << ", " << fAggregateDirectionZ << ")" << std::endl;
            std::cout << "Direction2: (" << fAggregateDirectionX2 << ", " << fAggregateDirectionY2 << ", " << fAggregateDirectionZ2 << ")" << std::endl;
            std::cout << "Energy: " << fAggregateReconstructedEnergy << " GeV" << std::endl;
            std::cout << "TrueOriginID: " << fAggregateTrueOriginID << std::endl;
            std::cout << "Number of PFParticles: " << fAggregateNumberOfPFParticles << std::endl;
            std::cout << "Energy deposited in the first 10 cm: " << fAggregateEnergyDepositedInFirst10cm << " ADC" << std::endl;
            std::cout << "Energy deposited in the second 10 cm: " << fAggregateEnergyDepositedInSecond10cm << " ADC" << std::endl;
            std::cout << "Energy deposited in the third 10 cm: " << fAggregateEnergyDepositedInThird10cm << " ADC" << std::endl;
            std::cout << "Energy deposited in the fourth 10 cm: " << fAggregateEnergyDepositedInFourth10cm << " ADC" << std::endl;
            std::cout << "Energy deposited in the fifth 10 cm: " << fAggregateEnergyDepositedInFifth10cm << " ADC" << std::endl;
            std::cout << "Energy deposited in the sixth 10 cm: " << fAggregateEnergyDepositedInSixth10cm << " ADC" << std::endl;
            std::cout << "Energy deposited in the seventh 10 cm: " << fAggregateEnergyDepositedInSeventh10cm << " ADC" << std::endl;
            std::cout << "Energy deposited in the eighth 10 cm: " << fAggregateEnergyDepositedInEighth10cm << " ADC" << std::endl;
            std::cout << "Energy deposited in the ninth 10 cm: " << fAggregateEnergyDepositedInNinth10cm << " ADC" << std::endl;
            std::cout << "Energy deposited in the tenth 10 cm: " << fAggregateEnergyDepositedInTenth10cm << " ADC" << std::endl;
            std::cout << "Energy deposited in the eleventh 10 cm: " << fAggregateEnergyDepositedInEleventh10cm << " ADC" << std::endl;
            std::cout << "Energy deposited in the twelfth 10 cm: " << fAggregateEnergyDepositedInTwelfth10cm << " ADC" << std::endl;
            std::cout << "Energy deposited in the thirteenth 10 cm: " << fAggregateEnergyDepositedInThirteenth10cm << " ADC" << std::endl;
            std::cout << "Energy deposited in the fourteenth 10 cm: " << fAggregateEnergyDepositedInFourteenth10cm << " ADC" << std::endl;
            std::cout << "Energy deposited in the fifteenth 10 cm: " << fAggregateEnergyDepositedInFifteenth10cm << " ADC" << std::endl;
            std::cout << "Energy deposited in the first 10 cm before: " << fAggregateEnergyDepositedInFirst10cmBefore << " ADC" << std::endl;
            std::cout << "Energy deposited in the second 10 cm before: " << fAggregateEnergyDepositedInSecond10cmBefore << " ADC" << std::endl;
            std::cout << "Spill Status Flag: " << fAggregateSpillStatusFlag << std::endl;
            std::cout << "PassSelectionCriterion: " << fAggregatePassSelectionCriterion << std::endl;
            std::cout << "PassSecondSelectionCriterion: " << fAggregatePassSecondSelectionCriterion << std::endl;
            // print the direction of all the daughters
            std::cout << "Daughter Directions: " << std::endl;
            for (auto const& pfp : pfpVec) {
                auto summary = getInformation(pfp, event);
                std::cout << "  Direction: (" << summary[1] << ", " << summary[2] << ", " << summary[3] << ")" << std::endl;
            }

        }

}

// ============================================================================
//  Module registration macro
// ============================================================================
DEFINE_ART_MODULE(NeutrinoAna::FindNeutrinos)
