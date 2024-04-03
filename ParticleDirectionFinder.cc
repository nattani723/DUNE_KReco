#include "ParticleDirectionFinder.h"

using namespace pandora;

namespace kaon_reconstruction
{

  ParticleDirectionFinder::ParticleDirectionFinder() :

    m_ThetaBinSize(0.06),
    m_PhiBinSize(0.06),
    m_SmoothingWindow(3)
  {
  }

  //-----------------------------------------------------------------------------

  void CollectSPsInROI(const std::vector<art::Ptr<recob::SpacePoint>>& SPList, const TVector3 KEnd, std::vector<art::Ptr<recob::SpacePoint>>& SPListROI) const
  {
  }

  //-----------------------------------------------------------------------------

  void FillAngularDistributionMap(const std::vector<art::Ptr<recob::SpacePoint>>& SPListROI, const TVector3 KEnd, AngularDistributionMap3D& AngularDistributionMap) const
  {
  }

  //-----------------------------------------------------------------------------

  void SmoothAngularDistributionMap(AngularDistributionMap3D& AngularDistributionMap) const
  {
  }

  //-----------------------------------------------------------------------------

  void RetrievePeakDirections(const AngularDistributionMap3D& AngularDistributionMap, std::vector<TVector2>& PeakDirectionVectors) const
  {
  }

  //-----------------------------------------------------------------------------

} // namespace kaon_reconstruction



  //------------------------------------------------------------------------------------------------------------------------------------------ 

  void HitSplitAlg::fillAngularDistributionMap3D( std::vector<art::Ptr<recob::SpacePoint>>& sp_from_recoobj,
                                                     TVector3 Kend_candidate,
                                                     std::map<int, std::map<int, double>> &angular_distribution_map_3D){

    const TVector3 Xaxis(1.,0.,0.);

    for (auto spIter = sp_from_recoobj.begin(); spIter != sp_from_recoobj.end(); ++spIter) {

      const TVector3 hit_position = (*spIter)->XYZ();
      //cout << hit_position.X() << endl;
      const TVector3 distance_vector = hit_position - Kend_candidate;

      if(distance_vector.Mag() > region_of_interest) continue;

      double theta = distance_vector.Theta();
      double sintheta = TMath::Sin(theta);
      double phi = distance_vector.Phi();

      int theta_factor = (int)(std::floor(theta / thetaBinSize));
      int phi_factor = (int)(std::floor(phi / phiBinSize));

      if( angular_distribution_map_3D.find(theta_factor) == angular_distribution_map_3D.end() &&
	  angular_distribution_map_3D[theta_factor].find(phi_factor) == angular_distribution_map_3D[theta_factor].end() )
	angular_distribution_map_3D[theta_factor][phi_factor] = sintheta; 
      else angular_distribution_map_3D[theta_factor][phi_factor] += sintheta;

      /*
     if(angular_distribution_map_3D.find(theta_factor) == angular_distribution_map_3D.end()){
        if(angular_distribution_map_3D[theta_factor].find(phi_factor) == angular_distribution_map_3D[theta_factor].end()){
          angular_distribution_map_3D[theta_factor][phi_factor] = sintheta;
        }
      }
      else angular_distribution_map_3D[theta_factor][phi_factor] += sintheta;
      */
    }

  }



  //------------------------------------------------------------------------------------------------------------------------------------------                                                                                                                                                                                

  void HitSplitAlg::fillAngularDistributionMap3D( std::vector<art::Ptr<recob::Hit>>& hits_from_pfparticle,
                                                     TVector3 Kend_candidate,
                                                     art::FindManyP<recob::SpacePoint>& spacepoint_per_hit,
                                                     std::map<int, std::map<int, double>> &angular_distribution_map_3D){

    //std::vector<art::Ptr<recob::SpacePoint>> spacepoint_vec;
    const TVector3 Xaxis(1.,0.,0.);

    
    //cout << "hits_from_pfparticle.size(): " << hits_from_pfparticle.size() << endl;
    for(size_t i_h=0; i_h<hits_from_pfparticle.size(); i_h++){
      //spacepoint_vec.clear();
      //cout << "hits_from_pfparticle[i_h].key(): " << hits_from_pfparticle[i_h].key() << endl;
      //cout << spacepoint_per_hit.at(hits_from_pfparticle[i_h].key())[0]->XYZ().X() << endl;
      auto spacepoint_vec = spacepoint_per_hit.at(hits_from_pfparticle[i_h].key());
      //cout << spacepoint_vec.size() << endl;
      if(!spacepoint_vec.size()) continue;
      const TVector3 hit_position = spacepoint_vec.at(0)->XYZ();
      //cout << hit_position.X() << endl;
      const TVector3 distance_vector = hit_position - Kend_candidate;

      if(distance_vector.Mag() > region_of_interest) continue;

      double theta = distance_vector.Theta();
      double sintheta = TMath::Sin(theta);
      double phi = distance_vector.Phi();

      int theta_factor = (int)(std::floor(theta / thetaBinSize));
      int phi_factor = (int)(std::floor(phi / phiBinSize));
      //cout << "theta: " <<  distance_vector.Theta() << ", factor: " << theta_factor << endl;                                                                                                                                                                                                                                
      //cout << "phi: " <<  distance_vector.Phi() << ", factor: " << phi_factor << endl;                                                                                                                                                                                                                         
     if(angular_distribution_map_3D.find(theta_factor) == angular_distribution_map_3D.end()){
        if(angular_distribution_map_3D[theta_factor].find(phi_factor) == angular_distribution_map_3D[theta_factor].end()){
          angular_distribution_map_3D[theta_factor][phi_factor] = sintheta;
        }
      }
      else angular_distribution_map_3D[theta_factor][phi_factor] += sintheta;
    }
    
  }


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


  void HitSplitAlg::smoothAngularDistributionMap3D(std::map<int, std::map<int, double>> &angular_distribution_map_3D){


    const int loop_min = (-1) * (smoothing_window - 1) / 2;
    const int loop_max = (smoothing_window - 1) / 2;
    std::map<int, std::map<int, double>> angular_distribution_map_3D_temp(angular_distribution_map_3D);
    angular_distribution_map_3D.clear();

    for (auto const& angular_distribution_theta : angular_distribution_map_3D_temp){
      const int current_bin_theta(angular_distribution_theta.first);

      for(auto const& angular_distribution_phi : angular_distribution_theta.second){
        const int current_bin_phi(angular_distribution_phi.first);
        int bin_count = 0;
        double total = 0;

        for (int bin_offset_theta = loop_min; bin_offset_theta <= loop_max; ++bin_offset_theta){
          const int contributing_bin_theta = current_bin_theta + bin_offset_theta;

          for (int bin_offset_phi = loop_min; bin_offset_phi <= loop_max; ++bin_offset_phi){
            ++bin_count;
            const int contributing_bin_phi = current_bin_phi + bin_offset_phi;

            if(angular_distribution_map_3D_temp[ contributing_bin_theta ].find( contributing_bin_phi ) == angular_distribution_map_3D_temp[ contributing_bin_theta ].end()) total += 0;
            else total += angular_distribution_map_3D_temp[ contributing_bin_theta ].at( contributing_bin_phi );

          }
        }
        angular_distribution_map_3D[ current_bin_theta ][ current_bin_phi ] = total/(double)bin_count;
      }
    }

  }


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


    for (auto const& angular_distribution_theta : angular_distribution_map_3D){

      for(auto const& angular_distribution_phi : angular_distribution_theta.second){

        const int bin_theta = angular_distribution_theta.first;
        const int bin_phi = angular_distribution_phi.first;
        const double bin_weight = angular_distribution_phi.second ;

        int preceding_bin_theta = bin_theta - 1;
        int preceding_bin_phi = bin_phi - 1;
        int following_bin_theta = bin_theta + 1;
        int following_bin_phi = bin_phi + 1;
        bool FoundBin = false;
        bool MaxWeight = true;

        for (int scan_bin_theta = preceding_bin_theta; scan_bin_theta <= following_bin_theta; ++scan_bin_theta){
          for (int scan_bin_phi = preceding_bin_phi; scan_bin_phi <= following_bin_phi; ++scan_bin_phi){

            if(scan_bin_theta == bin_theta && scan_bin_phi == bin_phi) continue;

            if(angular_distribution_map_3D[ scan_bin_theta ].find( scan_bin_phi ) == angular_distribution_map_3D[ scan_bin_theta ].end() ||
               (std::fabs(angular_distribution_map_3D[ scan_bin_theta ].at( scan_bin_phi ) - angular_distribution_map_3D[ bin_theta ].at( bin_phi )) > std::numeric_limits<double>::epsilon())){
              FoundBin = true;
            }else{
              FoundBin = false;
              break;
            }
          }
          if(FoundBin == false) break;
        }


        for (int scan_bin_theta = preceding_bin_theta; scan_bin_theta <= following_bin_theta; ++scan_bin_theta){
          for (int scan_bin_phi = preceding_bin_phi; scan_bin_phi <= following_bin_phi; ++scan_bin_phi){
            double scan_bin_weight;
            if( angular_distribution_map_3D[ scan_bin_theta ].find( scan_bin_phi ) == angular_distribution_map_3D[ scan_bin_theta ].end() ){
              scan_bin_weight = 0;
            }else{
              scan_bin_weight = angular_distribution_map_3D[scan_bin_theta].at(scan_bin_phi);
            }
            if(bin_weight < scan_bin_weight) MaxWeight = false;
          }
        }

        if(MaxWeight == false) continue;
	TVector2 peak_bins;
        peak_bins.Set((double)bin_theta, (double)bin_phi);

        view_peak_map[bin_weight] = peak_bins;
        v_trk_flg_peak.push_back(trk_flg);

      }
    }
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
