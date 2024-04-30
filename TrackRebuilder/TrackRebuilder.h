#ifndef TRACK_REBUILDER
#define TRACK_REBUILDER 1

#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "art/Persistency/Common/PtrMaker.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArObjects/LArPfoObjects.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "larpandora/LArPandoraInterface/LArPandoraOutput.h" 

#include <iostream>


using namespace pandora;

namespace kaon_reconstruction 
{
  class TrackRebuilder
  {
  public:

    TrackRebuilder();

    pandora::StatusCode Run(TrackHitCollector::HitList& track_hit_list, const recob::Track& primary_track, const std::map<art::Ptr<recob::Hit>, art::Ptr<recob::SpacePoint>>& hitToSpacePointMap);
    recob::Track get_rebuild_reco_track();

  private:

    recob::Track track_rebuild(TrackHitCollector::HitList& track_hit_list, const recob::Track& primary_track, recob::Track& rebuild_reco_track, const std::map<art::Ptr<recob::Hit>, art::Ptr<recob::SpacePoint>>&hitToSpacePointMap) const;
    recob::Track build_track(int track_id, lar_content::LArTrackStateVector& track_state_vector) const;
    recob::Track rebuild_reco_track;
    int m_rebuild_track_counter;
  };


}//namespace kaon_reconstruction

#endif
