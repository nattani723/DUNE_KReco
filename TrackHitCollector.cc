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

    return sliding_fit_pitch;
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
      const TVector3 distance = hit_position - k_end;
      const double l = peak_direction.Dot(distance); // Projection length along peak direction
      const double t = peak_direction.Cross(distance).Mag(); // Perpendicular distance to peak direction

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


    // Perform a running fit to collect a pathway of hits
    unsigned int count = 0;
    bool hits_collected = true;
    bool is_end_downstream = false;
    if(initial_direction.Z() > 0.) is_end_downstream = true;

    TVector3 extrapolated_direction = peak_direction;
    TVector3 extrapolated_start_position = k_end;
    TVector3 extrapolated_end_position = extrapolated_start_position + (extrapolated_direction * highest_l); 
    pandora::CartesianVector pandora_extrapolated_end_position(extrapolated_end_position.X(), extrapolated_end_position.Y(), extrapolated_end_position.Z());

    const float sliding_fit_pitch = get_wire_pitch();

    // do we need to sort pandora_running_fit_position by vertex position?

    while (hits_collected){

      ++count;
      const int excess_hits_in_fit = pandora_running_fit_position_vec.size() - m_max_fitting_hits;

      // Remove furthest away hits
      if(excess_hits_in_fit>0){

	// sort running_fit_position_vec and pandora_fit_position_vec
	// is there a better way to do this?

	lar_content::LArConnectionPathwayHelper::SortByDistanceToPoint sort_end(pandora_extrapolated_end_position); 
	std::sort(pandora_running_fit_position_vec.begin(), pandora_running_fit_position_vec.end(), sort_end);

	for (int i = 0; i < excess_hits_in_fit; ++i){
	  pandora_running_fit_position_vec.erase(std::prev(pandora_running_fit_position_vec.end())); 
	}

	running_fit_position_vec.clear();
	TVector3 hit_position_tmp;
	for(auto pandora_running_fit_position : pandora_running_fit_position_vec){
	  hit_position_tmp.SetXYZ(pandora_running_fit_position.GetX(), pandora_running_fit_position.GetY(), pandora_running_fit_position.GetZ());
	  running_fit_position_vec.push_back(hit_position_tmp);
	}
	
      }

      const lar_content::ThreeDSlidingFitResult extrapolated_fit(&pandora_running_fit_position_vec, m_macro_sliding_fit_window, sliding_fit_pitch); 

      // apply nominal fit
      this->update_extrapolation(count, extrapolated_fit, extrapolated_start_position, extrapolated_end_position, extrapolated_direction, is_end_downstream);
      
      hits_collected = this->collect_subsection_hits(extrapolated_fit, extrapolated_start_position, extrapolated_end_position, extrapolated_direction, is_end_downstream, sp_list, running_fit_position_vec, pandora_running_fit_position_vec, unavailable_hit_list, track_hit_list, m_distance_to_line, m_hit_connection_distance);


      // If no hits found, fit by all collected hits
      if (!hits_collected){
	
	// fill pandora_running_fit_position_trackall_vector by storingg all positions of track hits collected so far
	pandora::CartesianPointVector pandora_running_fit_position_trackall_vec;

	for(size_t i_h=0; i_h<track_hit_list.size(); i_h++) { 
	  const TVector3 hit_position = fHitsToSpacePoints.at(track_hit_list[i_h])->XYZ(); 
	  const pandora::CartesianVector pandora_hit_position(hit_position.X(), hit_position.Y(), hit_position.Z());
	  pandora_running_fit_position_trackall_vector.push_back(pandora_hit_position);
	}

	const lar_content::ThreeDSlidingFitResult extrapolated_fit_trackall(&pandora_running_fit_position_trackall_vec, m_trackall_sliding_fit_window, sliding_fit_pitch);
	this->update_extrapolation(count, extrapolated_fit_trackall, extrapolated_start_position, extrapolated_end_position, extrapolated_direction, is_end_downstream, m_distance_to_line, m_hit_connection_distance);

      }


      // If no hits found, as a final effort, reduce the sliding fit window
      if (!hits_collected){
	
	const lar_content::ThreeDSlidingFitResult extrapolated_fit_micro(&pandora_running_fit_position_track_vec, m_high_resolution_sliding_fit_window, sliding_fit_pitch);
	this->update_extrapolation(count, extrapolated_fit_micro, extrapolated_start_position, extrapolated_end_position, extrapolated_direction, is_end_downstream, m_distance_to_line, m_hit_connection_distance);

      }

      // demand track has a significant number of collected hits
      if (track_hit_list.size() < m_hit_threshold_for_track) {
	track_hit_list.clear();
	return;
      }
    }


  }

 
  //------------------------------------------------------------------------------------------------------------------------------------------    

  //void TrackHitCollector::update_extrapolation(int count, const lar_content::ThreeDSlidingFitResult& extrapolated_fit, const TVector3& extrapolated_start_position, const TVector3& extrapolated_end_position, const TVector3& extrapolated_direction, const bool is_end_downstream, const SPList& sp_list, TVector3& running_fit_position_vector, pandora::CartesianPointVector& pandora_running_fit_position_vector, HitList& unavailable_hit_list, HitList& track_hit_list) const
void TrackHitCollector::update_extrapolation(int count, const lar_content::ThreeDSlidingFitResult& extrapolated_fit, const TVector3& extrapolated_start_position, const TVector3& extrapolated_end_position, const TVector3& extrapolated_direction, const bool is_end_downstream) const
  {

    /*
    extrapolated_start_position = count == 1 ? extrapolated_end_position
      : is_end_downstream ? extrapolated_fit.GetGlobalMaxLayerPosition()()
      : extrapolated_fit.GetGlobalMinLayerPosition()();
    
    extrapolated_direction = is_end_downstream ? extrapolated_fit.GetGlobalMaxLayerDirection()()
      : extrapolated_fit.GetGlobalMinLayerPosition() * (-1.f);
    */

    if (count == 1) {
      extrapolated_start_position = extrapolated_end_position;
    } else if (is_end_downstream) {
      extrapolated_start_position = extrapolated_fit.get_global_max_layer_position();
    } else {
      extrapolated_start_position = extrapolated_fit.get_global_min_layer_position();
    }

    if (is_end_downstream) {
      extrapolated_direction = extrapolated_fit.get_global_max_layer_direction();
    } else {
      extrapolated_direction = extrapolated_fit.get_global_min_layer_direction() * (-1.f);
    }
    
    extrapolated_end_position = extrapolated_start_position + (extrapolated_direction * m_growing_fit_segment_length);

    hits_collected = this->collect_subsection_hits(extrapolated_fit, extrapolated_start_position, extrapolated_end_position, extrapolated_direction, is_end_downstream, sp_list, running_fit_position_vector, pandora_running_fit_position_vector, unavailable_hit_list, track_hit_list);
    
  }



  //------------------------------------------------------------------------------------------------------------------------------------------    

bool TrackHitCollector::collect_subsection_hits(const lar_content::ThreeDSlidingFitResult& extrapolated_fit. const TVector3& extrapolated_start_position, const TVector3& extrapolated_end_position, const TVector3& extrapolated_direction, const bool is_end_downstream, const SPList& sp_list, TVector3& running_fit_position_vector, pandora::CartesianPointVector& pandora_running_fit_position_vector, HitList& unavailable_hit_list, HitList& track_hit_list, float& distance_to_line, float& hit_connection_distance)) const;


  //------------------------------------------------------------------------------------------------------------------------------------------

  bool TrackHitCollector::is_close_to_line(const TVector3& hit_position, const TVector3& line_start, const TVector3& line_direction, const double& distance_to_line) const;

  //------------------------------------------------------------------------------------------------------------------------------------------

  void TrackHitCollector::collect_connected_hits(HitList& collected_hit_list, const TVector3& extrapolated_start_position, const TVector3& extrapolated_direction, TVector3& running_fit_position_vector, pandora::CartesianPointVector& pandora_running_fit_position_vector, HitList& track_hit_list) const;

  //------------------------------------------------------------------------------------------------------------------------------------------

  double TrackHitCollector::get_closest_distance(const TVector3& hit_position, const vector<TVector3>& test_positions) const;

} // namespace kaon_reconstruction
