#include "TrackHitCollector.h"
//#include "../CCKaonProducer_module.h" 
#include "/exp/uboone/app/users/taniuchi/51_pandora/srcs/larpandoracontent/larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"


using namespace pandora;


namespace kaon_reconstruction
{

  float TrackUtilities::get_wire_pitch() 
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
 


  TrackHitCollector::TrackHitCollector() :
      
    m_hit_threshold_for_track(50),
    m_growing_fit_initial_length(10.f),
    m_initial_fit_distance_to_line(1.f),
    m_min_initial_hits_found(7),
    m_max_fitting_hits(15),
    m_local_sliding_fit_window(10),
    m_growing_fit_segment_length(5.f),
    m_distance_to_line(0.75),
    m_hit_connection_distance(0.75),
    m_trackall_sliding_fit_window(10),
    m_trackall_distance_to_line(0.75),
    m_trackall_hit_connection_distance(0.75),
    m_high_resolution_sliding_fit_window(5),
    m_high_resolution_distance_to_line(0.75),
    m_high_resolution_hit_connection_distance(0.75)
  {
  }

  pandora::StatusCode TrackHitCollector::Run(const TVector3& k_end, SPList& sp_list,
					     const TVector3& peak_direction, HitList& unavailable_hit_list, HitList& track_hit_list,
					     const std::map<art::Ptr<recob::SpacePoint>, art::Ptr<recob::Hit>>& spacepointToHitMap,
					     const std::map<art::Ptr<recob::Hit>, art::Ptr<recob::SpacePoint>>& hitToSpacePointMap)
  {

    this->find_track_hits(sp_list, unavailable_hit_list, track_hit_list, k_end, peak_direction, spacepointToHitMap, hitToSpacePointMap);

    // Demand that spine is significant, be lenient here as some have small stubs and a gap
    if (track_hit_list.size() < m_hit_threshold_for_track)
        return STATUS_CODE_NOT_FOUND;

    return STATUS_CODE_SUCCESS;

}


  void TrackHitCollector::find_track_hits(SPList& sp_list, HitList& unavailable_hit_list, HitList& track_hit_list, const TVector3& k_end, const TVector3& peak_direction, const std::map<art::Ptr<recob::SpacePoint>, art::Ptr<recob::Hit>>& spacepointToHitMap, const std::map<art::Ptr<recob::Hit>, art::Ptr<recob::SpacePoint>>& hitToSpacePointMap) const
  {

    if (sp_list.empty()) return;

    // Use initial direction to find seed hits for a starting fit
    double highest_l = 0.;
    std::vector<TVector3> running_fit_position_vec;  
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

	auto corresponding_hit = spacepointToHitMap.at(sp); // Retrieve corresponding hit (modify based on actual mapping)

	//skip hits in primary track
	if (std::find(unavailable_hit_list.begin(), unavailable_hit_list.end(), corresponding_hit) != unavailable_hit_list.end()) continue;
	  
	track_hit_list.push_back(corresponding_hit);	
	running_fit_position_vec.push_back(hit_position);
	pandora_running_fit_position_vec.push_back(pandora_hit_position);      

      }
    }

    // Require significant number of initial hits
    if (running_fit_position_vec.size() < m_min_initial_hits_found){
      std::cout << "Requiring significant number of initial hits" << std::endl;
      track_hit_list.clear();
      return;
    }

    // Perform a running fit to collect a pathway of hits
    unsigned int count = 0;
    bool hits_collected = true;
    bool is_end_downstream = false;
    //if(peak_direction.Z() > 0.) is_end_downstream = true;

    TVector3 extrapolated_direction = peak_direction;
    TVector3 extrapolated_start_position = k_end;
    TVector3 extrapolated_end_position = extrapolated_start_position + (extrapolated_direction * highest_l); 

    const float sliding_fit_pitch = TrackUtilities::get_wire_pitch();
    // do we need to sort pandora_running_fit_position by vertex position?
    while (hits_collected){

      ++count;

      const int excess_hits_in_fit = pandora_running_fit_position_vec.size() - m_max_fitting_hits;
      //std::cout << "excess_hits_in_fit: " << excess_hits_in_fit << std::endl;

      // sort running_fit_position_vec and pandora_fit_position_vec
      pandora::CartesianVector pandora_extrapolated_end_position(extrapolated_end_position.X(), extrapolated_end_position.Y(), extrapolated_end_position.Z());
      lar_content::LArConnectionPathwayHelper::SortByDistanceToPoint sort_end(pandora_extrapolated_end_position); 
      std::sort(pandora_running_fit_position_vec.begin(), pandora_running_fit_position_vec.end(), sort_end);

      // Remove furthest away hits
      if(excess_hits_in_fit>0){

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

      //test
      /*
	for(auto pandora_running_fit_position : pandora_running_fit_position_vec){
	  std::cout << "pandora_running_fit_position: " << pandora_running_fit_position.GetX() << " " << pandora_running_fit_position.GetY() << " " << pandora_running_fit_position.GetZ() << std::endl;
	}
      */


      // apply nominal fit
      const lar_content::ThreeDSlidingFitResult extrapolated_fit(&pandora_running_fit_position_vec, m_local_sliding_fit_window, sliding_fit_pitch); 
      is_end_downstream = this->determine_is_end_down_stream(peak_direction, extrapolated_fit);
      this->update_extrapolation(count, extrapolated_fit, extrapolated_start_position, extrapolated_end_position, extrapolated_direction, is_end_downstream);

      //std::cout << "count: " << count << ", peak_direction: " << extrapolated_direction.X() << " " << extrapolated_direction.Y() << " " << extrapolated_direction.Z() << std::endl;
      
      hits_collected = this->collect_subsection_hits(extrapolated_fit, extrapolated_start_position, extrapolated_end_position, extrapolated_direction, is_end_downstream, sp_list, running_fit_position_vec, pandora_running_fit_position_vec, unavailable_hit_list, track_hit_list, m_distance_to_line, m_hit_connection_distance, spacepointToHitMap, hitToSpacePointMap);


      // If no hits found, fit by all collected hits
      /*
      if (!hits_collected){

	// fill pandora_running_fit_position_trackall_vector by storingg all positions of track hits collected so far
	pandora::CartesianPointVector pandora_running_fit_position_trackall_vec;

	for(size_t i_h=0; i_h<track_hit_list.size(); i_h++) { 
	  const TVector3 hit_position = hitToSpacePointMap.at(track_hit_list[i_h])->XYZ(); 
	  const pandora::CartesianVector pandora_hit_position(hit_position.X(), hit_position.Y(), hit_position.Z());
	  pandora_running_fit_position_trackall_vec.push_back(pandora_hit_position);
	}

	const lar_content::ThreeDSlidingFitResult extrapolated_fit_trackall(&pandora_running_fit_position_trackall_vec, m_trackall_sliding_fit_window, sliding_fit_pitch);
	is_end_downstream = this->determine_is_end_down_stream(peak_direction, extrapolated_fit_trackall);
	this->update_extrapolation(count, extrapolated_fit_trackall, extrapolated_start_position, extrapolated_end_position, extrapolated_direction, is_end_downstream);
	
	hits_collected = this->collect_subsection_hits(extrapolated_fit_trackall, extrapolated_start_position, extrapolated_end_position, extrapolated_direction, is_end_downstream, sp_list, running_fit_position_vec, pandora_running_fit_position_vec, unavailable_hit_list, track_hit_list, m_trackall_distance_to_line, m_trackall_hit_connection_distance, spacepointToHitMap, hitToSpacePointMap);

      }
      */


      // If no hits found, as a final effort, reduce the sliding fit window
      if (!hits_collected){
	const lar_content::ThreeDSlidingFitResult extrapolated_fit_micro(&pandora_running_fit_position_vec, m_high_resolution_sliding_fit_window, sliding_fit_pitch);
	is_end_downstream = this->determine_is_end_down_stream(peak_direction, extrapolated_fit_micro);
	this->update_extrapolation(count, extrapolated_fit_micro, extrapolated_start_position, extrapolated_end_position, extrapolated_direction, is_end_downstream);

	hits_collected = this->collect_subsection_hits(extrapolated_fit_micro, extrapolated_start_position, extrapolated_end_position, extrapolated_direction, is_end_downstream, sp_list, running_fit_position_vec, pandora_running_fit_position_vec, unavailable_hit_list, track_hit_list, m_high_resolution_distance_to_line, m_high_resolution_hit_connection_distance, spacepointToHitMap, hitToSpacePointMap);

      }
	    
    }

  }

 
  //------------------------------------------------------------------------------------------------------------------------------------------    

  bool TrackHitCollector::determine_is_end_down_stream(const TVector3 peak_direction, const lar_content::ThreeDSlidingFitResult& extrapolated_fit) const
  {

    pandora::CartesianVector minLayerPosition = extrapolated_fit.GetGlobalMinLayerPosition();
    pandora::CartesianVector maxLayerPosition = extrapolated_fit.GetGlobalMaxLayerPosition();

    pandora::CartesianVector pandora_vector_from_min_to_max = maxLayerPosition - minLayerPosition;
    TVector3 vector_from_min_to_max(pandora_vector_from_min_to_max.GetX(), pandora_vector_from_min_to_max.GetY(), pandora_vector_from_min_to_max.GetZ());

    return (peak_direction.Dot(vector_from_min_to_max) > 0);

  }

  //------------------------------------------------------------------------------------------------------------------------------------------    

void TrackHitCollector::update_extrapolation(int count, const lar_content::ThreeDSlidingFitResult& extrapolated_fit, TVector3& extrapolated_start_position, TVector3& extrapolated_end_position, TVector3& extrapolated_direction, const bool is_end_downstream) const
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
      extrapolated_start_position.SetXYZ(extrapolated_fit.GetGlobalMaxLayerPosition().GetX(), extrapolated_fit.GetGlobalMaxLayerPosition().GetY(), extrapolated_fit.GetGlobalMaxLayerPosition().GetZ());
    } else {
      extrapolated_start_position.SetXYZ(extrapolated_fit.GetGlobalMinLayerPosition().GetX(), extrapolated_fit.GetGlobalMinLayerPosition().GetY(), extrapolated_fit.GetGlobalMinLayerPosition().GetZ());
    }

    if (is_end_downstream) {
      extrapolated_direction.SetXYZ(extrapolated_fit.GetGlobalMaxLayerDirection().GetX(), extrapolated_fit.GetGlobalMaxLayerDirection().GetY(), extrapolated_fit.GetGlobalMaxLayerDirection().GetZ());
    } else {
      extrapolated_direction.SetXYZ(extrapolated_fit.GetGlobalMinLayerDirection().GetX()*(-1.), extrapolated_fit.GetGlobalMinLayerDirection().GetY()*(-1.), extrapolated_fit.GetGlobalMinLayerDirection().GetZ()*(-1.));
    }
    
    extrapolated_end_position = extrapolated_start_position + (extrapolated_direction * m_growing_fit_segment_length);

    // hits_collected = this->collect_subsection_hits(extrapolated_fit, extrapolated_start_position, extrapolated_end_position, extrapolated_direction, is_end_downstream, sp_list, running_fit_position_vector, pandora_running_fit_position_vector, unavailable_hit_list, track_hit_list);
    
  }



  //------------------------------------------------------------------------------------------------------------------------------------------    

  bool TrackHitCollector::collect_subsection_hits(const lar_content::ThreeDSlidingFitResult& extrapolated_fit, const TVector3& extrapolated_start_position, const TVector3& extrapolated_end_position, const TVector3& extrapolated_direction, const bool is_end_downstream, SPList& sp_list, std::vector<TVector3>& running_fit_position_vec, pandora::CartesianPointVector& pandora_running_fit_position_vec, HitList& unavailable_hit_list, HitList& track_hit_list, const float& distance_to_line, const float& hit_connection_distance, const std::map<art::Ptr<recob::SpacePoint>, art::Ptr<recob::Hit>>& spacepointToHitMap, const std::map<art::Ptr<recob::Hit>, art::Ptr<recob::SpacePoint>>& hitToSpacePointMap) const
  {

    float extrapolated_start_l(0.f), extrapolated_start_t1(0.f), extrapolated_start_t2(0.f);
    float extrapolated_end_l(0.f), extrapolated_end_t1(0.f), extrapolated_end_t2(0.f);

    const pandora::CartesianVector pandora_extrapolated_start_position(extrapolated_start_position.X(), extrapolated_start_position.Y(), extrapolated_start_position.Z());
    const pandora::CartesianVector pandora_extrapolated_end_position(extrapolated_end_position.X(), extrapolated_end_position.Y(), extrapolated_end_position.Z());

    extrapolated_fit.GetLocalPosition(pandora_extrapolated_start_position, extrapolated_start_l, extrapolated_start_t1, extrapolated_start_t2);
    extrapolated_fit.GetLocalPosition(pandora_extrapolated_end_position, extrapolated_end_l, extrapolated_end_t1, extrapolated_end_t2); 


    HitList collected_hit_list;
    this->sort_sp_by_distance(sp_list, extrapolated_start_position);

    for(auto it_sp = sp_list.begin(); it_sp != sp_list.end(); ++it_sp) {

      auto hit = spacepointToHitMap.at(*it_sp);
      if(std::find(track_hit_list.begin(), track_hit_list.end(), hit) != track_hit_list.end()) continue;
      if(std::find(unavailable_hit_list.begin(), unavailable_hit_list.end(), hit) != unavailable_hit_list.end()) continue;


      float hit_l(0.f), hit_t1(0.f), hit_t2(0.f);
      const TVector3 hit_position = (*it_sp)->XYZ();
      pandora::CartesianVector pandora_hit_position(hit_position.X(), hit_position.Y(), hit_position.Z());
      extrapolated_fit.GetLocalPosition(pandora_hit_position, hit_l, hit_t1, hit_t2);
      
      // Assess whether hit is within section boundaries
      if(is_end_downstream && ((hit_l < extrapolated_start_l) || (hit_l > extrapolated_end_l))) continue;
      if(!is_end_downstream && ((hit_l > extrapolated_start_l) || (hit_l < extrapolated_end_l))) continue;

      //if(this->is_close_to_line(hit_position, extrapolated_start_position, extrapolated_direction, m_distance_to_line)) std::cout << "Assess whether hit is close to connecting line" << is_close_to_line(hit_position, extrapolated_start_position, extrapolated_direction, m_distance_to_line) << std::endl;

      if (this->is_close_to_line(hit_position, extrapolated_start_position, extrapolated_direction, m_distance_to_line))
	collected_hit_list.push_back(spacepointToHitMap.at(*it_sp));

    }

    // Now find a continuous path of collected hits
    const int n_initial_hits(track_hit_list.size()); 
    this->collect_connected_hits(collected_hit_list, extrapolated_start_position, extrapolated_direction, running_fit_position_vec, pandora_running_fit_position_vec, track_hit_list, hit_connection_distance, hitToSpacePointMap);
    const int n_final_hits(track_hit_list.size());

    //std::cout << "n_initial_hits: " << n_initial_hits << ", n_final_hits: " << n_final_hits << std::endl;

    return (n_final_hits != n_initial_hits); 

  }


  //------------------------------------------------------------------------------------------------------------------------------------------

  void TrackHitCollector::sort_sp_by_distance(SPList& sp_list, const TVector3& sort_position) const
  {

    // SP and hit positions are confusing, better rename to make it consistent 

    std::sort(sp_list.begin(), sp_list.end(),
              [&sort_position](const art::Ptr<recob::SpacePoint>& sp1, const art::Ptr<recob::SpacePoint>& sp2) {
		const TVector3 position1 = sp1->XYZ();
		const TVector3 position2 = sp2->XYZ();

		double distance1 = (position1 - sort_position).Mag();
		double distance2 = (position2 - sort_position).Mag();

		return distance1 < distance2;
              });
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  bool TrackHitCollector::is_close_to_line(const TVector3& hit_position, const TVector3& line_start, const TVector3& line_direction, const double& distance_to_line) const
  {
    const double transverse_distance_from_line = line_direction.Cross(hit_position - line_start).Mag();
    //std::cout << "transverse_distance_from_line is " << transverse_distance_from_line << ", and distance_to_line is " << distance_to_line << std::endl;
    if(transverse_distance_from_line <= distance_to_line) return true;
    return false;
  }
  
  //------------------------------------------------------------------------------------------------------------------------------------------
  
  void TrackHitCollector::collect_connected_hits(HitList& collected_hit_list, const TVector3& extrapolated_start_position, const TVector3& extrapolated_direction, std::vector<TVector3>& running_fit_position_vec, pandora::CartesianPointVector& pandora_running_fit_position_vec, HitList& track_hit_list, const float& hit_connection_distance, const std::map<art::Ptr<recob::Hit>, art::Ptr<recob::SpacePoint>>& hitToSpacePointMap) const
  {
    
    // Now add connected hits  

    bool found = true;
    
    while(found) {

      found = false;

      //for(size_t i_h=0; i_h<track_hit_list.size(); i_h++){ 
      for(size_t i_h=0; i_h<collected_hit_list.size(); i_h++){ 

	if(std::find(track_hit_list.begin(), track_hit_list.end(), collected_hit_list[i_h]) != track_hit_list.end()) continue;

	TVector3 hit_position;// = hitToSpacePointMap.at(collected_hit_list[i_h])->XYZ();

	if (hitToSpacePointMap.find(collected_hit_list[i_h]) != hitToSpacePointMap.end()) {
	  hit_position = hitToSpacePointMap.at(collected_hit_list[i_h])->XYZ();
	} else {
	  std::cerr << "Warning: Hit not found in map." << std::endl;
	  continue;  // Skip to the next iteration if appropriate
	}
	pandora::CartesianVector pandora_hit_position(hit_position.X(), hit_position.Y(), hit_position.Z());
	
	if( this->get_closest_distance(hit_position, running_fit_position_vec) > hit_connection_distance )
	  continue;
	
	found = true;
	running_fit_position_vec.push_back(hit_position);
	pandora_running_fit_position_vec.push_back(pandora_hit_position); 
	track_hit_list.push_back(collected_hit_list[i_h]);
	//std::cout << "collected hit: " << hit_position.X() << " "  << hit_position.Y() << " " << hit_position.Z() << std::endl;

      }      
      
    }

  }
  
  //------------------------------------------------------------------------------------------------------------------------------------------
  
  double TrackHitCollector::get_closest_distance(const TVector3& hit_position, const std::vector<TVector3>& test_positions) const
  {
    
    double closest_distance_squared(std::numeric_limits<double>::max());
    
    for(TVector3 test_position : test_positions){
      
      const double separation_squared = (test_position - hit_position).Mag2();
      if (separation_squared < closest_distance_squared) closest_distance_squared = separation_squared; 

    }
    return std::sqrt(closest_distance_squared);
  }
  
} // namespace kaon_reconstruction
