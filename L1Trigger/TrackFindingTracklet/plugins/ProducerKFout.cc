#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/EDPutToken.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Handle.h"

#include "L1Trigger/TrackTrigger/interface/Setup.h"
#include "L1Trigger/TrackerTFP/interface/DataFormats.h"
#include "L1Trigger/TrackTrigger/interface/L1TrackQuality.h"

#include <string>
#include <numeric>

using namespace std;
using namespace edm;
using namespace trackerTFP;
using namespace tt;

namespace trklet {

  /*! \class  trklet::ProducerKFout
   *  \brief  Converts KF output into TFP output
   *  A bit accurate emulation of the track transformation, the 
   *  eta routing and splitting of the 96-bit track words into 64-bit 
   *  packets. Also run is a bit accurate emulation of the track quality
   *  BDT, whose output is also added to the track word.
   *  \author Christopher Brown
   *  \date   2021, Aug
   */
  class ProducerKFout : public stream::EDProducer<> {
  public:
    explicit ProducerKFout(const ParameterSet&);
    ~ProducerKFout() override {}

  private:
    void beginRun(const Run&, const EventSetup&) override;
    void produce(Event&, const EventSetup&) override;
    void endJob() {}

    // ED input token of kf stubs
    EDGetTokenT<StreamsStub> edGetTokenStubs_;
    // ED input token of kf tracks
    EDGetTokenT<StreamsTrack> edGetTokenTracks_;
    // ED input token of kf input to kf output TTTrack map
    EDGetTokenT<TTTrackRefMap> edGetTokenTTTrackRefMap_;
    // ED output token for accepted kfout tracks
    EDPutTokenT<StreamsTrack> edPutTokenAccepted_;
    // ED output token for truncated kfout tracks
    EDPutTokenT<StreamsTrack> edPutTokenLost_;
    // Setup token
    ESGetToken<Setup, SetupRcd> esGetTokenSetup_;
    // DataFormats token
    ESGetToken<DataFormats, DataFormatsRcd> esGetTokenDataFormats_;
    // configuration
    ParameterSet iConfig_;
    // helper class to store configurations
    const Setup* setup_;
    // helper class to extract structured data from tt::Frames
    const DataFormats* dataFormats_;
    // Bins for dPhi/dZ use to create weight LUT
    vector<double> dPhiBins_;
    vector<double> dZBins_;

    std::unique_ptr<L1TrackQuality> trackQualityModel_;
    vector<int> tqBins_;
    double tqTanlScale_;
    double tqZ0Scale_;
    double ap_fixed_rescale = 32.0;

    // For convenience and keeping readable code, accessed many times
    int numWorkers_;

    int partialTrackWordBits_;

    // Helper function to convert floating chi2 to chi2 bin
    /*template <typename T>
    unsigned int digitise(const T& bins, double value, double factor) {
      unsigned int bin = 0;
      for (unsigned int i = 0; i < bins.size() - 1; i++) {
        if (value * factor > bins[i] && value * factor <= bins[i + 1])
          bin = i;
      }
      return bin;
    }*/
  };

  ProducerKFout::ProducerKFout(const ParameterSet& iConfig) : iConfig_(iConfig) {
    const string& labelKF = iConfig.getParameter<string>("LabelKF");
    const string& labelAS = iConfig.getParameter<string>("LabelAS");
    const string& branchStubs = iConfig.getParameter<string>("BranchAcceptedStubs");
    const string& branchTracks = iConfig.getParameter<string>("BranchAcceptedTracks");
    const string& branchLost = iConfig.getParameter<string>("BranchLostTracks");
    // book in- and output ED products
    edGetTokenStubs_ = consumes<StreamsStub>(InputTag(labelKF, branchStubs));
    edGetTokenTracks_ = consumes<StreamsTrack>(InputTag(labelKF, branchTracks));
    edGetTokenTTTrackRefMap_ = consumes<TTTrackRefMap>(InputTag(labelAS, branchTracks));
    edPutTokenAccepted_ = produces<StreamsTrack>(branchTracks);
    edPutTokenLost_ = produces<StreamsTrack>(branchLost);
    // book ES products
    esGetTokenSetup_ = esConsumes<Setup, SetupRcd, Transition::BeginRun>();
    esGetTokenDataFormats_ = esConsumes<DataFormats, DataFormatsRcd, Transition::BeginRun>();
    // initial ES products
    setup_ = nullptr;
    dataFormats_ = nullptr;

    trackQualityModel_ = std::make_unique<L1TrackQuality>(iConfig.getParameter<edm::ParameterSet>("TrackQualityPSet"));
    edm::ParameterSet trackQualityPSset = iConfig.getParameter<edm::ParameterSet>("TrackQualityPSet");
    tqBins_ = trackQualityPSset.getParameter<vector<int>>("tqemu_bins");
    tqTanlScale_ = trackQualityPSset.getParameter<double>("tqemu_TanlScale");
    tqZ0Scale_ = trackQualityPSset.getParameter<double>("tqemu_Z0Scale");
  }

  void ProducerKFout::beginRun(const Run& iRun, const EventSetup& iSetup) {
    // helper class to store configurations
    setup_ = &iSetup.getData(esGetTokenSetup_);
    // helper class to extract structured data from tt::Frames
    dataFormats_ = &iSetup.getData(esGetTokenDataFormats_);

    // Calculate 1/dz**2 and 1/dphi**2 bins for v0 and v1 weightings
    for (int i = 0; i < pow(2, dataFormats_->width(Variable::dPhi, Process::ht)) / pow(2, setup_->weightBinFraction());
         i++)
      dPhiBins_.push_back(
          pow(dataFormats_->base(Variable::dPhi, Process::ht) * (i + 1) * pow(2, setup_->weightBinFraction()), -2));

    for (int i = 0; i < pow(2, dataFormats_->width(Variable::dZ, Process::ht)) / pow(2, setup_->weightBinFraction());
         i++)
      dZBins_.push_back(
          pow(dataFormats_->base(Variable::dZ, Process::ht) * (i + 1) * pow(2, setup_->weightBinFraction()), -2));

    partialTrackWordBits_ = TTBV::S_ / 2;
    numWorkers_ = setup_->kfNumWorker();
  }

  // Helper function to convert floating chi2 to chi2 bin
  /*template <typename T>
  int ProducerKFout::digitise(const vector<T> Bins, T Value, T factor) {
    for (int i = 0; i < (int)Bins.size(); i++) {
      if (Value * factor > Bins[i] && Value * factor <= Bins[i + 1]) {
        return i;
      }
    }
    for (int i = 0; i < pow(2, dataFormats_->width(Variable::dZ, Process::kfin)) / pow(2, setup_->weightBinFraction());
         i++) {
      temp_dz =
          pow(dataFormats_->base(Variable::dZ, Process::kfin) * (i + 1) * pow(2, setup_->weightBinFraction()), -2);
      temp_dz = temp_dz * setup_->dzTruncation();
      temp_dz = std::ceil(temp_dz);
      dZBins_.push_back(temp_dz / setup_->dzTruncation());
    }
    numWorkers_ = setup_->kfNumWorker();
    partialTrackWordBits_ = TTBV::S_ / 2;
  }*/

  void ProducerKFout::produce(Event& iEvent, const EventSetup& iSetup) { /*
    // empty KFout product
    StreamsTrack accepted(setup_->numRegions() * setup_->tfpNumChannel());
    StreamsTrack lost(setup_->numRegions() * setup_->tfpNumChannel());
    // read in KF Product and produce KFout product
    if (setup_->configurationSupported()) {
      Handle<StreamsStub> handleStubs;
      iEvent.getByToken<StreamsStub>(edGetTokenStubs_, handleStubs);
      const StreamsStub& streamsStubs = *handleStubs.product();
      Handle<StreamsTrack> handleTracks;
      iEvent.getByToken<StreamsTrack>(edGetTokenTracks_, handleTracks);
      const StreamsTrack& streamsTracks = *handleTracks.product();
      Handle<TTTrackRefMap> handleTTTrackRefMap;
      iEvent.getByToken<TTTrackRefMap>(edGetTokenTTTrackRefMap_, handleTTTrackRefMap);
      const TTTrackRefMap& ttTrackRefMap = *handleTTTrackRefMap.product();
      // 18 Output Links (First Vector) each has a vector of tracks per event (second vector) each track is 3 32 bit TTBV partial tracks
      vector<vector<TTBV>> sortedPartialTracks(setup_->numRegions() * setup_->tfpNumChannel(), vector<TTBV>(0));

      TrackKFOutSAPtrCollectionss inTrackStreams;
      TrackKFOutSAPtrCollectionss outTrackStreams;

      // Setup empty collections for input tracks to be routed
      for (int iRegion = 0; iRegion < setup_->numRegions(); iRegion++) {
        TrackKFOutSAPtrCollections temp_collection;
        for (int iLink = 0; iLink < setup_->tfpNumChannel(); iLink++) {
          TrackKFOutSAPtrCollection temp;
          for (int iTrack = 0; iTrack < setup_->numFramesIO(); iTrack++)
            temp.emplace_back(std::make_shared<TrackKFOut>());
          temp_collection.push_back(temp);
        }
        outTrackStreams.push_back(temp_collection);
      }

      // Setup empty collections for output tracks from routing
      for (int iRegion = 0; iRegion < setup_->numRegions(); iRegion++) {
        TrackKFOutSAPtrCollections temp_collection;
        for (int iLink = 0; iLink < numWorkers_; iLink++) {
          TrackKFOutSAPtrCollection temp;
          for (int iTrack = 0; iTrack < setup_->numFramesIO(); iTrack++)
            temp.emplace_back(std::make_shared<TrackKFOut>());
          temp_collection.push_back(temp);
        }
        inTrackStreams.push_back(temp_collection);
      }

      StreamsTrack outputStreamsTracks(setup_->numRegions() * setup_->tfpNumChannel());

      // Setup containers for track quality
      float tempTQMVA = 0.0;
      // Due to ap_fixed implementation in CMSSW this 10,5 must be specified at compile time, TODO make this a changeable parameter
      std::vector<ap_fixed<10, 5>> trackQuality_inputs = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

      for (int iLink = 0; iLink < (int)streamsTracks.size(); iLink++) {
        for (int iTrack = 0; iTrack < (int)streamsTracks[iLink].size(); iTrack++) {
          const auto& track = streamsTracks[iLink].at(iTrack);
          TrackKF inTrack(track, dataFormats_);

          double temp_z0 = inTrack.zT() - ((inTrack.cot() * setup_->chosenRofZ()));

          // Correction to Phi calcuation depending if +ve/-ve phi sector
          const double baseSectorCorr = 0;//InTrack.sectorPhi() ? -setup_->baseSector() : setup_->baseSector();

          double temp_phi0 = InTrack.phiT() - ((InTrack.inv2R()) * setup_->chosenRofPhi()) + baseSectorCorr;

          double temp_tanL = inTrack.cotGlobal();

          TTBV hitPattern(0, setup_->numLayers());

          double tempchi2rphi = 0;
          double tempchi2rz = 0;

          int temp_nstub = 0;
          int temp_ninterior = 0;
          bool counter = false;

          for (int iStub = 0; iStub < setup_->numLayers() - 1; iStub++) {
            const auto& stub = streamsStubs[setup_->numLayers() * iLink + iStub].at(iTrack);
            StubKF inStub(stub, dataFormats_, iStub);

            if (!stub.first.isNonnull()) {
              if (counter)
                temp_ninterior += 1;
              continue;
            }

            counter = true;

            double tempv0 = dPhiBins_[(InStub.dPhi() / (dataFormats_->base(Variable::dPhi, Process::ht) *
                                                        pow(2, setup_->weightBinFraction())))];
            double tempv1 = dZBins_[(
                InStub.dZ() / (dataFormats_->base(Variable::dZ, Process::ht) * pow(2, setup_->weightBinFraction())))];

            double tempRphi = phiSquared * tempv0;
            double tempRz = zSquared * tempv1;

            tempchi2rphi += tempRphi;
            tempchi2rz += tempRz;
          }  // Iterate over track stubs

          // Create bit vectors for eacch output, including digitisation of chi2
          // TODO implement extraMVA, bendChi2, d0
          TTBV trackValid(1, TTTrack_TrackWord::TrackBitWidths::kValidSize, false);
          TTBV extraMVA(0, TTTrack_TrackWord::TrackBitWidths::kMVAOtherSize, false);
          TTBV bendChi2(0, TTTrack_TrackWord::TrackBitWidths::kBendChi2Size, false);
          TTBV chi2rphi(digitise(TTTrack_TrackWord::chi2RPhiBins, tempchi2rphi, (double)setup_->kfoutchi2rphiConv()),
                        TTTrack_TrackWord::TrackBitWidths::kChi2RPhiSize,
                        false);
          TTBV chi2rz(digitise(TTTrack_TrackWord::chi2RZBins, tempchi2rz, (double)setup_->kfoutchi2rzConv()),
                      TTTrack_TrackWord::TrackBitWidths::kChi2RZSize,
                      false);
          TTBV d0(0, TTTrack_TrackWord::TrackBitWidths::kD0Size, false);
          TTBV z0(
              temp_z0, dataFormats_->base(Variable::zT, Process::kf), TTTrack_TrackWord::TrackBitWidths::kZ0Size, true);
          TTBV tanL(temp_tanL,
                    dataFormats_->base(Variable::cot, Process::kf),
                    TTTrack_TrackWord::TrackBitWidths::kTanlSize,
                    true);
          TTBV phi0(temp_phi0,
                    dataFormats_->base(Variable::phiT, Process::kf),
                    TTTrack_TrackWord::TrackBitWidths::kPhiSize,
                    true);
          TTBV invR(-inTrack.inv2R(),
                    dataFormats_->base(Variable::inv2R, Process::kf),
                    TTTrack_TrackWord::TrackBitWidths::kRinvSize + 1,
                    true);
          invR.resize(TTTrack_TrackWord::TrackBitWidths::kRinvSize);

          int sortKey = 0;//(InTrack.sectorEta() < (int)(setup_->numSectorsEta() / 2)) ? 0 : 1;
          // Set correct bit to valid for track valid
          TrackKFOut temp_track(partialTrack1.set((partialTrackWordBits_ - 1)),
                                partialTrack2,
                                partialTrack3,
                                sortKey,
                                track,
                                iTrack,
                                iLink,
                                true);

          inTrackStreams[iLink / setup_->kfNumWorker()][iLink % setup_->kfNumWorker()][iTrack] =
              (std::make_shared<TrackKFOut>(temp_track));
        }  // Iterate over Tracks
      }    // Iterate over Links
           // Route Tracks in eta based on their sort key
      for (int iRegion = 0; iRegion < setup_->numRegions(); iRegion++) {
        int buffered_tracks[] = {0, 0};
        for (int iTrack = 0;
             iTrack < setup_->numFramesIO() * ((double)TTBV::S_ / TTTrack_TrackWord::TrackBitWidths::kTrackWordSize);
             iTrack++) {
          for (int iWorker = 0; iWorker < setup_->kfNumWorker(); iWorker++) {
            for (int iLink = 0; iLink < setup_->tfpNumChannel(); iLink++) {
              if ((inTrackStreams[iRegion][iWorker][iTrack]->sortKey() == iLink) &&
                  (inTrackStreams[iRegion][iWorker][iTrack]->dataValid() == true)) {
                outTrackStreams[iRegion][iLink][buffered_tracks[iLink]] = inTrackStreams[iRegion][iWorker][iTrack];
                buffered_tracks[iLink] = buffered_tracks[iLink] + 1;
              }
            }
          }
        }
      }

      // Pack output of router onto each link, with correct partial tracks in correct places
      for (int iRegion = 0; iRegion < setup_->numRegions(); iRegion++) {
        for (int iLink = 0; iLink < setup_->tfpNumChannel(); iLink++) {
          for (int iTrack = 0; iTrack < (int)outTrackStreams[iRegion][iLink].size(); iTrack++) {
            sortedPartialTracks[2 * iRegion + iLink].push_back(
                outTrackStreams[iRegion][iLink][iTrack]->PartialTrack1());
            sortedPartialTracks[2 * iRegion + iLink].push_back(
                outTrackStreams[iRegion][iLink][iTrack]->PartialTrack2());
            sortedPartialTracks[2 * iRegion + iLink].push_back(
                outTrackStreams[iRegion][iLink][iTrack]->PartialTrack3());
            outputStreamsTracks[2 * iRegion + iLink].emplace_back(outTrackStreams[iRegion][iLink][iTrack]->track());
          }
        }
      }
      // Fill products and match up tracks
      // store products
      const TTBV nullBitTrack(0, partialTrackWordBits_, false);
      for (int iLink = 0; iLink < (int)outputStreamsTracks.size(); iLink++) {
        // Iterate through partial tracks
        int numLinkTracks = (int)outputStreamsTracks[iLink].size();
        if (numLinkTracks == 0)
          continue;  // Don't fill links if no tracks
        if ((numLinkTracks % 2 != 0)) {
          sortedPartialTracks[iLink].push_back(nullBitTrack);  //Pad out final set of bits
          outputStreamsTracks[iLink].emplace_back(
              outputStreamsTracks[iLink][numLinkTracks++]);  //Pad out with final repeated track
        }                                                    //If there is an odd number of tracks
        for (int iTrack = 0; iTrack < (int)(sortedPartialTracks[iLink].size()); iTrack++) {
          if (iTrack % 2 != 1)  // Write to links every other partial track, 3 partial tracks per full TTTrack
            continue;
          TTTrackRef trackRef;
          for (auto& it : ttTrackRefMap) {  //Iterate through ttTrackRefMap to find TTTrackRef Key by a TTTrack Value
            if (it.second == outputStreamsTracks[iLink][(int)(iTrack - 1) / 3].first)
              trackRef = it.first;
          }
          if ((int)iTrack / 3 <=
              setup_->numFramesIO() * ((double)TTBV::S_ / TTTrack_TrackWord::TrackBitWidths::kTrackWordSize))
            accepted[iLink].emplace_back(
                std::make_pair(trackRef,
                               (sortedPartialTracks[iLink][iTrack - 1].slice(partialTrackWordBits_) +
                                sortedPartialTracks[iLink][iTrack].slice(partialTrackWordBits_))
                                   .bs()));
          else
            lost[iLink].emplace_back(
                std::make_pair(trackRef,
                               (sortedPartialTracks[iLink][iTrack - 1].slice(partialTrackWordBits_) +
                                sortedPartialTracks[iLink][iTrack].slice(partialTrackWordBits_))
                                   .bs()));
        }  //Iterate through sorted partial tracks
      }    // Iterate through links
    }      // Config Supported
    // store products
    iEvent.emplace(edPutTokenAccepted_, move(accepted));
    iEvent.emplace(edPutTokenLost_, move(lost));*/
  }
}  // namespace trklet

DEFINE_FWK_MODULE(trklet::ProducerKFout);
