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
        // print to be removed
        std::cout << "Daughter summary: Hits = " << summary[0]<< ", PFPs = " << summary[1] << std::endl;
        std::cout << "Daughter PFParticle PDG code: " << pfParticle->PdgCode() << std::endl;
        std::cout << "Daughter PFParticle isPrimary: " << pfParticle->IsPrimary() << std::endl;
        
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
//  Main analyze method
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
    fAggregateNumberOfHits           = -1;
    fAggregateNumberOfPFParticles    = -1;
    fAggregateTrueOriginID           = -1;
    fAggregatePassSelectionCriterion = -1;
    fAggregateSpillStatusFlag        = -1;
    fAggregateTriggerCandidateCount  = -1;
    fAggregateGroundShakeCount       = -1;
    fAggregateSumOfLastADCTicks      = -1;
    fAggregateSumOfTriggeredADCTicks = -1;

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
        
        // To be deleted in the future
        if (fRecoNumberOfPFParticles>40){
            std::cout << "Event " << event.id().event() << " has " << fRecoNumberOfPFParticles << " PFParticles." << std::endl;
            std::cout << "Vertex: (" << vtxX << ", " << vtxY << ", " << vtxZ << ")" << std::endl;
            std::cout << "Direction: (" << fRecoDirectionX << ", " << fRecoDirectionY << ", " << fRecoDirectionZ << ")" << std::endl;   
            std::cout << "Energy: " << fRecoNeutrinoEnergy << " GeV" << std::endl;
            std::cout << "PdgCode: " << fRecoNeutrinoPdgCode << std::endl;

        }

        // If this PFParticle is a primary neutrino candidate
        if (pfp->IsPrimary() &&
            (std::abs(pfp->PdgCode()) == 12 ||
             std::abs(pfp->PdgCode()) == 14 ||
             std::abs(pfp->PdgCode()) == 16)) {
            if (vtxY <= 550. && vtxZ >= 20. && summary[1] >= 4.0)
                fRecoPassSelectionCriterion = 1;

            fAggregateVertexX             = vtxX;
            fAggregateVertexY             = vtxY;
            fAggregateVertexZ             = vtxZ;
            fAggregateNumberOfPFParticles = static_cast<int>(daughterSummary[1]);
            fAggregateTrueOriginID        = trueOriginID;
            fAggregateNumberOfHits        = static_cast<int>(daughterSummary[0]);
            // print to be removed
            std::cout << "Primary neutrino candidate found in event "<< event.id().event() << " with "<< fAggregateNumberOfPFParticles << " PFPs." << std::endl;
            std::cout << "Vertex: (" << vtxX << ", " << vtxY << ", " << vtxZ << ")" << std::endl;
            std::cout << "Direction: (" << summary[1] << ", "<< summary[2] << ", "<< summary[3] << ")" << std::endl;
            std::cout << "Energy: " << summary[0] << " GeV"<< std::endl;
            std::cout << "PdgCode: " << pfp->PdgCode() << std::endl;

        }
        fRecoTree->Fill();
        if (fRecoPassSelectionCriterion)
            fAggregatePassSelectionCriterion = 1;
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

    // ------------------------------------------------------------------------
    //  Trigger candidate count (non-ground-shake)
    // ------------------------------------------------------------------------
    auto triggerCandHandle =
        event.getValidHandle<std::vector<dunedaq::trgdataformats::TriggerCandidateData>>(
            "triggerrawdecoder:daq");
    for (auto const& tc : *triggerCandHandle) {
        if (tc.type != dunedaq::trgdataformats::TriggerCandidateData::Type::kADCSimpleWindow)
            ++fAggregateTriggerCandidateCount;
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

    // Total number of hits in event
    fAggregateTotalNumberOfHits = static_cast<int>(allHitPtrs.size());

    // ------------------------------------------------------------------------
    //  Fill aggregate tree
    // ------------------------------------------------------------------------
    fAggregateTree->Fill();
}

// ============================================================================
//  Module registration macro
// ============================================================================
DEFINE_ART_MODULE(NeutrinoAna::FindNeutrinos)
