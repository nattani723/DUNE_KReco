#include "TrackHitCollector.h"

using namespace pandora;

namespace kaon_reconstruction
{
  
  TrackHitCollector::TrackHitCollector() :
      
    //m_peak_searh_region(15.),
    m_theta_bin_size(0.06),
    m_phi_bin_size(0.06),
    m_smoothing_window(1),
    m_peak_search_window(1),
    m_peak_open_angle(TMath::Pi() * 1/4),
    m_min_peak_height(0.4),
    m_max_num_peak(3)
  {
  }

  //------------------------------------------------------------------------------------------------------------------------------------------ 

  const double TrackHitCollocter::get_wire_pitch()
  {

    art::ServiceHandle<geo::Geometry> theGeometry;
    const unsigned int nWirePlanes(theGeometry->MaxPlanes());

    if (nWirePlanes > 3)
      throw cet::exception("LArPandoraTrackCreation") << " LArPandoraTrackCreation::produce --- More than three wire planes present ";

    if ((0 == theGeometry->Ncryostats()) || (0 == theGeometry->NTPC(geo::CryostatID{0})))
      throw cet::exception("LArPandoraTrackCreation") << " LArPandoraTrackCreation::produce --- unable to access first tpc in first cryostat ";

    std::unordered_set<geo::_plane_proj> planeSet;
    for (unsigned int iPlane = 0; iPlane < nWirePlanes; ++iPlane)
      (void) planeSet.insert(theGeometry->TPC().Plane(iPlane).View());

    if ((nWirePlanes != planeSet.size()) || !planeSet.count(geo::kU) || !planeSet.count(geo::kV) || (planeSet.count(geo::kW) && planeSet.count(geo::kY)))
      throw cet::exception("LArPandoraTrackCreation") << " LArPandoraTrackCreation::produce --- expect to find u and v views; if there is one further view, it must be w or y ";

    const bool useYPlane((nWirePlanes > 2) && planeSet.count(geo::kY));

    const float wirePitchU(theGeometry->WirePitch(geo::kU));
    const float wirePitchV(theGeometry->WirePitch(geo::kV));
    const float wirePitchW((nWirePlanes < 3) ? 0.5f * (wirePitchU + wirePitchV) : (useYPlane) ? theGeometry->WirePitch(geo::kY) : theGeometry->WirePitch(geo::kW));

    const float sliding_fit_pitch = wirePitchW;

  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void TrackHitCollector::find_track_hits(const SPList& sp_list, HitList& unavailable_hit_list, HitList& track_hit_list, const TVector3& k_end, TVector peak_direction) const
  {
  }

 
  //------------------------------------------------------------------------------------------------------------------------------------------    

  bool TrackHitCollector::collect_subsection_hits(const lar_content::ThreeDSlidingFitResult& extrapolated_fit. const TVector3& extrapolated_start_position, const TVector3& extrapolated_end_position, const TVector3& extrapolated_direction, const bool is_end_downstream, const SPList& sp_list, TVector3& running_fit_position_vector, pandora::CartesianPointVector& pandora_running_fit_position_vector, HitList& unavailable_hit_list, HitList& track_hit_list) const;


  //------------------------------------------------------------------------------------------------------------------------------------------

  bool TrackHitCollector::is_close_to_line(const TVector3& hit_position, const TVector3& line_start, const TVector3& line_direction, const double& distance_to_line) const;

  //------------------------------------------------------------------------------------------------------------------------------------------

  void TrackHitCollector::collect_connected_hits(HitList& collected_hit_list, const TVector3& extrapolated_start_position, const TVector3& extrapolated_direction, TVector3& running_fit_position_vector, pandora::CartesianPointVector& pandora_running_fit_position_vector, HitList& track_hit_list) const;

  //------------------------------------------------------------------------------------------------------------------------------------------

  double TrackHitCollector::get_closest_distance(const TVector3& hit_position, const vector<TVector3>& test_positions) const;

} // namespace kaon_reconstruction
