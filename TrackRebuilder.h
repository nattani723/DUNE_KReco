/*
 *
 * Header file for track rebuilder
 * Mostly following LArPandoraTrackCreation 
 * under larpandora/larpandora/LArPandoraEventBuilding/LArPandoraTrackCreation_module.cc
 *
 *
 */

#ifndef TRACK_REBULIDER
#define TRACK_REBUILDER 1

#include "art/Framework/Core/EDProducer.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "art/Persistency/Common/PtrMaker.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArObjects/LArPfoObjects.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "larpandora/LArPandoraInterface/LArPandoraOutput.h" 

namespace kaon_reconstruction
{

  class TrackRebuilder
  {

  public:
    TrackRebuilder();

  private:

    recob::Track track_rebuild(HitList& track_hit_list, const recob::Track& track) const;

    recob::Track build_track(int track_id, lar_content::LArTrackStateVector& track_state_vector);

    const int m_rebuild_track_id = 1000;
    const unsigned int m_slidingFitHalfWindow = 20;

    float sliding_fit_pitch = TrackUtilities::get_wire_pitch();
  }

}// namespace kaon_reconstruction
