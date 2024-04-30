#include "TrackHitCollector.h"
#include "TrackRebuilder.h"

using namespace pandora;

namespace kaon_reconstruction 
{
  TrackRebuilder::TrackRebuilder() : 
    m_rebuild_track_counter(1000)
  {
  }

  pandora::StatusCode TrackRebuilder::Run(TrackHitCollector::HitList& track_hit_list, const recob::Track& primary_track, const std::map<art::Ptr<recob::Hit>, art::Ptr<recob::SpacePoint>>& hitToSpacePointMap)
  {
    this->track_rebuild(track_hit_list, primary_track, rebuild_reco_track, hitToSpacePointMap);
    return STATUS_CODE_SUCCESS;
  }

  recob::Track TrackRebuilder::get_rebuild_reco_track() {
    return rebuild_reco_track;
  }
  
  recob::Track TrackRebuilder::track_rebuild(TrackHitCollector::HitList& track_hit_list, const recob::Track& primary_track, recob::Track& rebuild_reco_track, const std::map<art::Ptr<recob::Hit>, art::Ptr<recob::SpacePoint>>& hitToSpacePointMap) const{
    
    
    const float wire_pitch_w = TrackUtilities::get_wire_pitch();

    pandora::CartesianPointVector pandora_hit_position_vec;


    for (const auto& hit : track_hit_list) {
      auto sp = hitToSpacePointMap.at(hit); // Assuming fHitsToSpacePoints_old maps hit to space points
      const auto& xyz = sp->XYZ();
      pandora::CartesianVector pandora_hit_position(xyz[0], xyz[1], xyz[2]);
      pandora_hit_position_vec.emplace_back(std::move(pandora_hit_position));
    }

    std::unique_ptr<std::vector<recob::Track>> output_tracks(new std::vector<recob::Track>);
    const pandora::CartesianVector vertex_position(primary_track.End().X(), primary_track.End().Y(), primary_track.End().Z());

    // these are dummy values to get GetSlidingFitTrajectory instance
    lar_content::LArTrackStateVector track_state_vector;
    pandora::IntVector index_vector;
    const unsigned int sliding_fit_half_window = TrackUtilities::m_sliding_fit_half_window;
 
    lar_content::LArPfoHelper::GetSlidingFitTrajectory(pandora_hit_position_vec,
						       vertex_position,
						       sliding_fit_half_window,
						       wire_pitch_w,
						       track_state_vector,
						       &index_vector);

    // the first parameter gives the id of tracks in the event
    // to avoid overriding exsisting tracks, m_rebuild_track_counter is set to 1000
    int track_id = m_rebuild_track_counter;
    output_tracks->emplace_back(this->build_track( track_id++, track_state_vector));

    // always expect to have single track from single track_hit_list
    rebuild_reco_track = output_tracks->at(0);
    std::cout << "rebuild_reco_track.Length(): " << rebuild_reco_track.Length() << std::endl;
    return rebuild_reco_track;

  }

  recob::Track TrackRebuilder::build_track(int track_id, lar_content::LArTrackStateVector& track_state_vector) const
  {
    if (track_state_vector.empty()) std::cout << "BuildTrack - No input trajectory points provided" << std::endl;
    
    recob::tracking::Positions_t xyz;
    recob::tracking::Momenta_t pxpypz;
    recob::TrackTrajectory::Flags_t flags;
    
    for (const lar_content::LArTrackState& trackState : track_state_vector) {
      
      xyz.emplace_back(recob::tracking::Point_t(trackState.GetPosition().GetX(),
						trackState.GetPosition().GetY(),
						trackState.GetPosition().GetZ()));
      pxpypz.emplace_back(recob::tracking::Vector_t(trackState.GetDirection().GetX(),
						    trackState.GetDirection().GetY(),
						    trackState.GetDirection().GetZ()));
      
      // Set flag NoPoint if point has bogus coordinates, otherwise use clean flag set
      if (std::fabs(trackState.GetPosition().GetX() - util::kBogusF) <
	  std::numeric_limits<float>::epsilon() &&
	  std::fabs(trackState.GetPosition().GetY() - util::kBogusF) <
	  std::numeric_limits<float>::epsilon() &&
	  std::fabs(trackState.GetPosition().GetZ() - util::kBogusF) <
	  std::numeric_limits<float>::epsilon()) {
	flags.emplace_back(recob::TrajectoryPointFlags(recob::TrajectoryPointFlags::InvalidHitIndex,
						       recob::TrajectoryPointFlagTraits::NoPoint));
      }
      else {
	flags.emplace_back(recob::TrajectoryPointFlags());
      }
    }
    
    // note from gc: eventually we should produce a TrackTrajectory, not a Track with empty covariance matrix and bogus chi2, etc.
    return recob::Track(
			recob::TrackTrajectory(std::move(xyz), std::move(pxpypz), std::move(flags), false),
			util::kBogusI,
			util::kBogusF,
			util::kBogusI,
			recob::tracking::SMatrixSym55(),
			recob::tracking::SMatrixSym55(),
			track_id);
    
  }
  
}//namespace kaon_reconstruction
