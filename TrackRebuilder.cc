#include "TrackHitCollector.h"
#include "TrackRebuilder.h"

using namespace pandora;

namespace kaon_reconstruction 
{


  recob::Track TrackRebuilder::track_rebuild(HitList& track_hit_list, const recob::Track& rebuilt_track) const{


    const int m_rebuild_track_id = 1000;
    const unsigned int m_slidingFitHalfWindow = 20;

    float sliding_fit_pitch = TrackUtilities::get_wire_pitch();


    for (const auto& hit : hit_list) {
      auto sp = fHitsToSpacePoints_old.at(hit); // Assuming fHitsToSpacePoints_old maps hit to space points
      const auto& xyz = sp->XYZ(); // Assuming XYZ() returns a reference or a copy efficiently
      pandora::CartesianVector pandora_hit_position(xyz[0], xyz[1], xyz[2]);
      pandora_hit_positions.emplace_back(std::move(pandora_hit_position)); // Use std::move for efficiency
    }



  }

}//namespace kaon_reconstruction
