#include "ParticleDirectionFinder.h"

using namespace pandora;

namespace kaon_reconstruction
{

  ParticleDirectionFinder::ParticleDirectionFinder() :

    m_peak_searh_region(15.),
    m_theta_bin_size(0.06),
    m_phi_bin_size(0.06),
    m_smoothing_window(1),
    m_peak_search_window(1)
  {
  }

  //-----------------------------------------------------------------------------

  void ParticleDirectionFinder::collect_sp_in_roi(const SPList& sp_list, const TVector3 k_end, SPList& sp_list_roi) const
  {

    for(auto it_sp = sp_list.begin(); it_sp != sp_list.end(); ++it_sp){

      const TVector3 hit_position = (*it_sp)->XYZ();
      const TVector3 displacement_vector = hit_position - k_end;

      if( displacement_vector.Mag() > m_peak_searh_region) continue;

      sp_list_roi.push_back(*it_sp);

    }
    
  }

  //-----------------------------------------------------------------------------

  void ParticleDirectionFinder::fill_angular_distribution_map(const std::vector<art::Ptr<recob::SpacePoint>>& sp_list_roi, const TVector3 k_end, AngularDistributionMap3D& angular_distribution_map) const
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

  void retrieve_peak_directions(const angular_distribution_map_3d& angular_distribution_map, std::vector<TVector2>& peak_direction_vectors) const

  {


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
	  TVector2 peak_position(bin_theta, bin_phi);
	  view_peak_map[bin_weight] = peak_position; // Peaks are sorted by weight
	} 

      }
    }
    
  }


  //-----------------------------------------------------------------------------

} // namespace kaon_reconstruction



  //------------------------------------------------------------------------------------------------------------------------------------------ 
 
  bool HitSplitAlg::isInsideROI( art::Ptr<recob::SpacePoint> &sp,
				TVector3 Kend_candidate)
  {
    const TVector3 hit_position = sp->XYZ();
    const TVector3 distance_vector = hit_position - Kend_candidate;
    if(distance_vector.Mag() < region_of_interest_hitcandidate) return true;
    else return false;
  }


  //------------------------------------------------------------------------------------------------------------------------------------------ 

  void HitSplitAlg::fillHistAngularDistributionMap3DCheat( std::vector<art::Ptr<recob::SpacePoint>>& sp_from_recoobj,
							   TVector3 Kend_candidate,
							   std::map<int, std::map<int, std::map<int, double>>> &angular_distribution_map_3D_cheat,
							   std::map<int, TH2D*> &h_angular_distribution_pfparticle_cheat_3D){
    
    const TVector3 Xaxis(1.,0.,0.);

    const simb::MCParticle *particletmp;
    if(!sp_from_recoobj.size()) return;

    for (auto spIter = sp_from_recoobj.begin(); spIter != sp_from_recoobj.end(); ++spIter) {

      const TVector3 hit_position = (*spIter)->XYZ();
      const TVector3 distance_vector = hit_position - Kend_candidate;
      //art::Ptr<recob::Hit> hit = fSpacePointsToHits.at(*spIter);
      art::Ptr<recob::Hit> hit = fSpacePointsToHits_old.at(*spIter);
      truthHitMatcher(hit, particletmp);

      if(distance_vector.Mag() > region_of_interest) continue;
      if(!particletmp) continue;

      double theta = distance_vector.Theta();
      double sintheta = TMath::Sin(theta);
      double phi = distance_vector.Phi();

      int theta_factor = (int)(std::floor(theta / thetaBinSize));
      int phi_factor = (int)(std::floor(phi / phiBinSize));
      int pdg = particletmp->PdgCode();

      if( (angular_distribution_map_3D_cheat.find(pdg) == angular_distribution_map_3D_cheat.end()) &&
	  (angular_distribution_map_3D_cheat[pdg].find(theta_factor) == angular_distribution_map_3D_cheat[pdg].end()) &&
	  (angular_distribution_map_3D_cheat[pdg][theta_factor].find(phi_factor) == angular_distribution_map_3D_cheat[pdg][theta_factor].end()) ){
	angular_distribution_map_3D_cheat[pdg][theta_factor][phi_factor] = sintheta;
      }
      else angular_distribution_map_3D_cheat[pdg][theta_factor][phi_factor] += sintheta;

    }
    
    //cout << "angular_distribution_map_3D_cheat.size(): " << angular_distribution_map_3D_cheat.size() << endl;
    if(angular_distribution_map_3D_cheat.size()==0) return;

    for (auto const& angular_distribution_map_3D : angular_distribution_map_3D_cheat){

      TH2D * h = new TH2D("", "", num_bin_theta, 0, 3.14, num_bin_phi, -3.14, 3.14);
      for (auto const& x : angular_distribution_map_3D.second){
	for (auto const& y : x.second)
	  h->Fill(thetaBinSize * x.first, phiBinSize * y.first, y.second);
      }
      //cout << angular_distribution_map_3D.first << " " << h->GetEntries() << endl;
      h_angular_distribution_pfparticle_cheat_3D[angular_distribution_map_3D.first] = h;
      //cout<< h_angular_distribution_pfparticle_cheat_3D.size() << endl;
      //h_angular_distribution_pfparticle_cheat_3D.push_back(h);
    }

  }


  //-----------------------------------------------------------------------------------------------------------------------------------------------------  



  //------------------------------------------------------------------------------------------------------------------------------------------  

  void HitSplitAlg::accumulateAngularDistributionMap3D(std::map<int, std::map<int, double>> &angular_distribution_map_3D_1,
                                                          std::map<int, std::map<int, double>> &angular_distribution_map_3D_2,
                                                          std::map<int, std::map<int, double>> &angular_distribution_map_3D){
 

    angular_distribution_map_3D = angular_distribution_map_3D_1;

    for( auto const& angular_distribution_map_3D_theta_2 : angular_distribution_map_3D_2){
      for( auto const& angular_distribution_map_3D_phi_2 : angular_distribution_map_3D_theta_2.second){
        if( angular_distribution_map_3D[angular_distribution_map_3D_theta_2.first].find( angular_distribution_map_3D_phi_2.first ) == angular_distribution_map_3D[angular_distribution_map_3D_theta_2.first].end() ) angular_distribution_map_3D[angular_distribution_map_3D_theta_2.first][angular_distribution_map_3D_phi_2.first] = angular_distribution_map_3D_phi_2.second;
        else{
          angular_distribution_map_3D[angular_distribution_map_3D_theta_2.first][angular_distribution_map_3D_phi_2.first] += angular_distribution_map_3D_phi_2.second;
          //cout << "ACCUMULATE "  << angular_distribution_map_3D_1[angular_distribution_map_3D_theta_2.first][angular_distribution_map_3D_phi_2.first] << " " << angular_distribution_map_3D_2[angular_distribution_map_3D_theta_2.first][angular_distribution_map_3D_phi_2.first] << " " << angular_distribution_map_3D[angular_distribution_map_3D_theta_2.first][angular_distribution_map_3D_phi_2.first] << endl;
        }
      }
    }
  }



  //------------------------------------------------------------------------------------------------------------------------------------------                                                                                                                                                                                

  void HitSplitAlg::obtainPeakVector3D(std::map<int, std::map<int, double>> &angular_distribution_map_3D,
                                          vector<bool>& v_trk_flg_peak,
                                          std::map<double, TVector2, std::greater<>>& view_peak_map,
                                          bool trk_flg){
    //vector<TVector2>& view_peak_vector,                                                                                                                                                                                                                                                                                     


  }


  //------------------------------------------------------------------------------------------------------------------------------------------                                                                                                                                                                                

  void HitSplitAlg::findBestAngularPeak3D(std::map<int, std::map<int, double>> &angular_distribution_map_3D,
                                             std::map<double, TVector2, std::greater<>>& view_peak_map,
                                             vector<TVector2> &best_peak_bins){ 

    if(view_peak_map.size()==0) return;
    const TVector2 best_peak_bin = view_peak_map.begin()->second;
    //if(view_peak_map.begin()->first < 1) return;
    best_peak_bins.push_back(best_peak_bin);
    //cout << "view_peak_map[0] is " << view_peak_map.begin()->first << " " << view_peak_map.begin()->second.X() << endl;
    //cout << "best_peak_bin is " << best_peak_bin.X()*thetaBinSize << " " << best_peak_bin.Y()*thetaBinSize << endl;

    TVector3 best_peak_dir;
    Double_t theta_best = best_peak_bin.X()*thetaBinSize;
    Double_t phi_best = best_peak_bin.Y()*phiBinSize;

    best_peak_dir.SetMagThetaPhi(1., theta_best, phi_best);
 
    cout << "view_peak_map.size(): " << view_peak_map.size() << endl;

    for (auto const& entry : view_peak_map){

      TVector2 theta_phi_bin = entry.second;
      TVector3 view_peak_dir;
      Double_t theta_view = theta_phi_bin.X()*thetaBinSize;
      Double_t phi_view = theta_phi_bin.Y()*phiBinSize;

      view_peak_dir.SetMagThetaPhi(1., theta_view, phi_view);
 
      if(!(theta_phi_bin.X() == best_peak_bin.X() && theta_phi_bin.Y() == best_peak_bin.Y())){

        double open_angle = view_peak_dir.Angle(best_peak_dir);

	//cout << "open_angle: " << open_angle << ", entry.first: " << entry.first << endl;
 
        //if(open_angle > TMath::Pi()*2/3 && entry.first>0.4) best_peak_bins.push_back(theta_phi_bin);
        if(open_angle > TMath::Pi()*1/4 && entry.first>0.4) best_peak_bins.push_back(theta_phi_bin);
 
      }
      if(best_peak_bins.size() >= 3) break; //just look second highest peak with large opening angle
    }

    /*
    for(auto const& entry : best_peak_bins){
      cout << "best peaks are: " << entry.X()*thetaBinSize << " " << entry.Y()*phiBinSize << endl;
    }
    */

  }
