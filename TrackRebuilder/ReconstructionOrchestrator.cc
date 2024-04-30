#ifndef RECONSTRUCTION_ORCHESTRATOR
#define RECONSTRUCTION_ORCHESTRATOR 1

//#include "../CCKaonProducer_module.h"
#include "ParticleDirectionFinder.h"
#include "TrackHitCollector.h"
#include "TrackRebuilder.h"
#include "AngularDistributionDrawer.h"


namespace kaon_reconstruction {

  class ReconstructionOrchestrator {
  public:
    ReconstructionOrchestrator() : directionFinder(), hitCollector() {}

    typedef std::vector<art::Ptr<recob::SpacePoint>> SPList;
    typedef std::vector<art::Ptr<recob::Hit>> HitList;

    void runReconstruction(SPList& sp_list, const std::map<art::Ptr<recob::SpacePoint>, art::Ptr<recob::Hit>>& spacepointToHitMap, const std::map<art::Ptr<recob::Hit>, art::Ptr<recob::SpacePoint>>& hitToSpacePointMap, const art::Ptr<recob::Track> k_track, const HitList& hits_from_track);

    const std::vector<HitList>& getHitLists();
    const std::vector<recob::Track>& getRebuildTrackList();


  private:
    ParticleDirectionFinder directionFinder;
    TrackHitCollector hitCollector;
    TrackRebuilder trackRebuilder;
    std::vector<recob::Track> rebuildTrackList;
    std::vector<HitList> trackHitLists;

    // Other members like unavailableHitList might be defined here if they are shared across methods
  };
  
  const std::vector<recob::Track>& ReconstructionOrchestrator::getRebuildTrackList()  { return rebuildTrackList; }

  const std::vector<ReconstructionOrchestrator::HitList>& ReconstructionOrchestrator::getHitLists() { return trackHitLists; }
  
  //void ReconstructionOrchestrator::runrecobnstruction(const SPList& sp_list, const recob::Track& k_track) {
  void ReconstructionOrchestrator::runReconstruction(SPList& sp_list, const std::map<art::Ptr<recob::SpacePoint>, art::Ptr<recob::Hit>>& spacepointToHitMap, const std::map<art::Ptr<recob::Hit>, art::Ptr<recob::SpacePoint>>& hitToSpacePointMap,  const art::Ptr<recob::Track> k_track, const HitList& hits_from_track) {
    
    // Container for peak direction vectors calculated by ParticleDirectionFinder
    std::vector<TVector3> peakDirectionVector;
    
    HitList unavailableHitList = hits_from_track;
    
    // Run the ParticleDirectionFinder
    auto status = directionFinder.Run(sp_list, k_track, unavailableHitList, peakDirectionVector);

    if (status != pandora::STATUS_CODE_SUCCESS) {
      std::cout << "ParticleDirectionFinder FAILED" << std::endl;
      return;
    }
    
    // Assuming you've obtained k_end and a peak direction from directionFinder
    // that you want to use in hitCollector. 
    if (peakDirectionVector.empty()) return;
    
    //loop over number of peaks
    for(const TVector3& peakDirection : peakDirectionVector){
      
      HitList trackHitList;
      
      // Run the TrackHitCollector

      auto collectorStatus = hitCollector.Run(directionFinder.get_k_end(), directionFinder.get_sp_list_roi(), peakDirection, unavailableHitList, trackHitList, spacepointToHitMap, hitToSpacePointMap);
      if (collectorStatus != pandora::STATUS_CODE_SUCCESS) {
	std::cout << "TrackHitCollector FAILED" << std::endl;
	continue;
      }
      
      //make run function with statuscode return
      //think about how to retrieve reco::track object define getrebuildtracj in this orchestrator?
      trackRebuilder.Run(trackHitList, *k_track, hitToSpacePointMap);
      rebuildTrackList.push_back( trackRebuilder.get_rebuild_reco_track() );
      //rebuildTrackList.push_back( trackRebuilder.track_rebuild(trackHitList, k_trac );
      trackHitLists.push_back( trackHitList );
    }
    
  }
  

} // namespace kaon_recontruction
#endif
