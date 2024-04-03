#include "ParticleDirectionFinder.h"

using namespace pandora;

namespace kaon_reconstruction
{
  
  ParticleDirectionFinder::ParticleDirectionFinder() :
      
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


  void ParticleDirectionFinder::collect_sp_in_roi(const SPList& sp_list, const TVector3& k_end, double& region_of_interest, SPList& sp_list_roi) const
  {

    for(auto it_sp = sp_list.begin(); it_sp != sp_list.end(); ++it_sp){

      const TVector3 hit_position = (*it_sp)->XYZ();
      const TVector3 displacement_vector = hit_position - k_end;

      if( displacement_vector.Mag() > region_of_interest) continue;

      sp_list_roi.push_back(*it_sp);

    }
    
  }

  //-----------------------------------------------------------------------------

  void ParticleDirectionFinder::fill_angular_distribution_map(const std::vector<art::Ptr<recob::SpacePoint>>& sp_list_roi, const TVector3& k_end, AngularDistributionMap3D& angular_distribution_map) const
  {

    const TVector3 x_axis(1.,0.,0.);

    for(auto it_sp = sp_list_roi.begin(); it_sp != sp_list_roi.end(); ++it_sp){

      const TVector3 hit_position = (*it_sp)->XYZ();
      const TVector3 displacement_vector = hit_position - k_end;

      double theta = distance_vector.Theta();
      double phi = distance_vector.Phi();
      int theta_factor = (int)(std::floor(theta / m_theta_nin_size));
      int phi_factor = (int)(std::floor(phi / m_phi_bin_size));

      if( angular_distribution_map.find(theta_factor) == angular_distribution_map.end() &&
	  angular_distribution_map[theta_factor].find(phi_factor) == angular_distribution_map[theta_factor].end() )
	angular_distribution_map[theta_factor][phi_factor] = TMath::Sin(theta); // weight by sintheta as this is 3d angular distribution
      else angular_distribution_map[theta_factor][phi_factor] += TMath::Sin(theta);

    }

  }

  //-----------------------------------------------------------------------------

  void ParticleDirectionFinder::smooth_angular_distribution_map(AngularDistributionMap3D& angular_distribution_map) const
  {

    //const int loop_min = (-1) * (m_smoothing_window - 1) / 2;
    //const int loop_max = (m_smoothing_window - 1) / 2;
    const int loop_min = (-1) * m_smoothing_window;
    const int loop_max = m_smoothing_window;

    AngularDistributionMap3D angular_distribution_map_temp( angular_distribution_map );
    angular_distribution_map_temp.clear();


    for(auto const& theta_entry : angular_distribution_map_temp) {

      const int current_bin_theta = theta_entry.first;

      for (const auto& phi_entry : theta_entry.second) {

	const int current_bin_phi = phi_entry.first;
	double total = 0;

	// Loop over the neighborhood of the current bin, defined by loop_min and loop_max
	for (int bin_offset_theta = loop_min; bin_offset_theta <= loop_max; ++bin_offset_theta) {
	  const int contributing_bin_theta = current_bin_theta + bin_offset_theta;

	  for (int bin_offset_phi = loop_min; bin_offset_phi <= loop_max; ++bin_offset_phi) {
	    const int contributing_bin_phi = current_bin_phi + bin_offset_phi;

	    // Check if the contributing bin exists in the temp map
	    // If the contributing bin exists, add its value to the total
	    // If not, add 0 to make entry for map

            if(angular_distribution_map_temp[ contributing_bin_theta ].find( contributing_bin_phi ) == angular_distribution_map_temp[ contributing_bin_theta ].end()) total += 0;
            else total += angular_distribution_map_temp[ contributing_bin_theta ].at( contributing_bin_phi );
	  }
	}

	int num_bins = (2 * m_smoothing_window + 1) * (2 * m_smoothing_window + 1);
	//int num_bins = smoothing_window * smoothing_window;
	angular_distribution_map[current_bin_theta][current_bin_phi] = total / static_cast<double>(num_bins);

      }

    }  

  }

  //-----------------------------------------------------------------------------

  void ParticleDirectionFinder::retrieve_peak_directions(const angular_distribution_map_3d& angular_distribution_map, std::map<double, TVector3, std::greater<>>& sort_peak_direction_map) const

  {

    //sort peaks by their height (high -> low)
    //std::map<double, TVector3, std::greater<>> sort_peak_direction_map;

    for(auto const& theta_entry : angular_distribution_map){

      const int bin_theta = theta_entry.first;

      for(auto const& phi_entry : theta_entry.second){

	const int bin_phi = phi_entry.first;
	const double bin_weight = phi_entry.second;

	bool is_peak = true; // Assume the current bin is a peak until proven otherwise

	// Check the neighbouring bins to see if any have a greater weight than the current bin
	for (int dTheta = -m_peak_search_window; is_peak && dTheta <= m_peak_search_window; ++dTheta){{
	  for (int dPhi = -m_peak_search_window; dPhi <= m_peak_search_window; ++dPhi) {
	    // Skip the current bin itself
	    if (dTheta == 0 && dPhi == 0) continue;

	    int neighbor_theta = bin_theta + dTheta;
	    int neighbor_phi = bin_phi + dPhi;

	    // Retrieve the weight of the neighbouring bin, if it exists
	    auto it_neighbor_theta = angular_distribution_map.find(neighbor_theta);
	    if (it_neighbor_theta != angular_distribution_map.end()) {
	      auto it_neighbor_phi = it_neighbor_theta->second.find(neighbor_phi);
	      if (it_neighbor_phi != it_neighbor_theta->second.end() && it_neighbor_phi->second > bin_weight) {
		// A neighboring bin has a greater weight, so the current bin cannot be a peak
		is_peak = false;
		break; // No need to check further neighbors
	      }
	    }
	  }
	}

	// If the current bin is determined to be a peak, record its weight and position
	  if (is_peak) {
	    //TVector2 peak_position(bin_theta, bin_phi);
	    TVector3 peak_direction.SetMagTheta(1, bin_theta*m_theta_bin_size, bin_phi*m_phi_bin_size);
	    sort_peak_direction_map[bin_weight] = peak_direction; // Peaks are sorted by weight
	  } 
	}
      }   
    }
  }

  //-----------------------------------------------------------------------------

  void ParticleDirectionFinder::refine_peak_directions(const std::map<double, TVector3, std::greater<>>& sort_peak_direction_map, vector<TVector3> &peak_direction_vector) const
  {

    if(sort_peak_direction_map.empty()) return;

    // Start with the highest-weighted peak
    const TVector3& highest_peak_direction = sort_peak_direction_map.begin()->second;
    peak_direction_vector.push_back(highest_peak_direction);

    // Iterate through remaining peaks to find those significantly different in direction

    for(const auto& entry : sort_peak_direction_map) {

      TVector3 current_peak_direction = entry.second;
      if(current_peak_direction == highest_peak_direction) continue; // Skip the highest peak itself

      double open_angle = current_peak_direction.Angle(highest_peak_direction);
      
      // Consider a peak significant if it's sufficiently angularly separated or too low
      if(open_angle > m_peak_open_angle && entry.first > m_min_peak_height)
	peak_direction_vector.push_back(current_peak_direction);

      if (best_peak_bins.size() >= m_max_num_peak) break;
	
    }

  }

} // namespace kaon_reconstruction


