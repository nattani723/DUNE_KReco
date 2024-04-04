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

    /**
     * @brief run LArPandoraTrackCreation 
     *
     * @param track_hit_list: list of hits collected by TrackHitCollector
     * @param primary_track: in this case we assume K+ track
     * 
     */
    recob::Track track_rebuild(HitList& track_hit_list, const recob::Track& primary_track) const;

    /**
     * @brief build reco::Track object
     *
     * @param track_id: the ID of tracks in the event
     * @param track_state_vector
     * 
     */
    recob::Track build_track(int track_id, lar_content::LArTrackStateVector& track_state_vector);


    int m_rebuild_track_counter = 1000;
    const unsigned int m_sliding_fit_half_window = 20;


  }

}// namespace kaon_reconstruction
