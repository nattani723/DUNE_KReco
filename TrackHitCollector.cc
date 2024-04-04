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

  void TrackHitCollector::find_track_hits(const SPList& sp_list, HitList& unavailable_hit_list, HitList& track_hit_list, const TVector3& k_end, TVector3& peak_direction) const
  {

    if (sp_list.empty()) return;

    // Use initial direction to find seed hits for a starting fit
    double highest_l = 0.;
    vector<TVector3> running_fit_position_vec;  
    pandora::CartesianPointVector pandora_running_fit_position_vec;


    // Iterate through each space point in the list
    for (const auto& sp : sp_list) {

      const TVector3 hit_position = sp->XYZ(); 
      pandora::CartesianVector pandora_hit_position(hit_position.X(), hit_position.Y(), hit_position.Z());
        
      // Calculate distance and alignment vectors relative to the endpoint and peak direction
      const TVector3 distance_vector = hit_position - k_end;
      const double l = peak_direction.Dot(distance_vector); // Projection length along peak direction
      const double t = peak_direction.Cross(distance_vector).Mag(); // Perpendicular distance to peak direction

      // Filter hits based on proximity to the peak direction and availability
      if ((l < m_growing_fit_initial_length) && (l > 0.) && (t < m_initial_fit_distance_to_line)) { 

	if(l>highest_l)
	  highest_l = l;

	auto corresponding_hit = fSpacePointsToHits.at(sp); // Retrieve corresponding hit (modify based on actual mapping)
	if (std::find(unavailable_hit_list.begin(), unavailable_hit_list.end(), corresponding_hit) == unavailable_hit_list.end())
	  track_hit_list.push_back(corresponding_hit); // Add hit to track hit list if NOT unavailable

	running_fit_position_vec.push_back(hit_position);
	pandora_running_fit_position_vec.push_back(pandora_hit_position);      

      }
    }

    // Require significant number of initial hits
    if (running_fit_position_vec.size() < min_initial_hits_found){
      // cout << "Requiring significant number of initial hits" << endl;
      shower_spine_hit_list.clear();
      return;
    }



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
