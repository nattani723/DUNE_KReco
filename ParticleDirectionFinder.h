/*
 * Header file for particle direction finder tool class.
 * Plots angular distribution of hits insider the region of interest,
 * and find peaks as candidates for daughter particle direction.
 */

#ifndef PARTICLE_DIRECTION_FINDER
#define PARTICLE_DIRECTION_FINDER 1

//#include "HitSplitAlg_module.h"

#include "TMath.h"
#include <algorithm>
#include <functional>
#include <TVector3.h>

#include "art/Persistency/Common/PtrMaker.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/PFParticle.h"

#include "Api/PandoraApi.h"
#include "Pandora/AlgorithmHeaders.h"
#include "larpandoracontent/LArContent.h"

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"


using namespace pandora;
//using namespace std;

namespace kaon_reconstruction
{

  class ParticleDirectionFinder
  {

  public:

    ParticleDirectionFinder();
    typedef std::vector<art::Ptr<recob::Hit>> HitList;
    typedef std::vector<art::Ptr<recob::SpacePoint>> SPList;
    typedef std::map<int, std::map<int, double>> AngularDistribution3DMap;

    pandora::StatusCode Run(const SPList& sp_list, const art::Ptr<recob::Track> primary_track, const  HitList& unavailable_hit_list, const std::map<art::Ptr<recob::SpacePoint>, art::Ptr<recob::Hit>>& spacepointToHitMap, std::vector<TVector3> &peak_direction_vector);
    float get_theta_bin_size() const;
    float get_phi_bin_size() const;

    SPList& get_sp_list_roi();
    const TVector3& get_k_end() const;
    const TVector3& get_k_vtx() const;

  private:

    /*
     * @brief  Collect spacepoints inside the region of interest (ROI)
     *
     * @param sp_list: spacepoint list of K+ daughter candidate hits
     * @param k_end: vector of K+ track candidate's end
     * @param sp_list_roi: spacepoint list inside the ROI
     *
     */

    //void collect_sp_in_roi(const SPList& sp_list, const TVector3& k_end, float& region_of_interest, SPList& sp_list_roi) const;

    void collect_sp_in_roi(const SPList& sp_list, const TVector3& k_end, float& region_of_interest, SPList& sp_list_roi, const HitList& unavailable_hit_list, const std::map<art::Ptr<recob::SpacePoint>, art::Ptr<recob::Hit>>& spacepointToHitMap) const;

    /*
     * @brief  Fill map of angular distribution for spacepoints inside the region of interest (ROI)
     *
     * @param sp_list_roi: spacepoint list inside the ROI
     * @param k_end: vector of K+ track candidate's end
     * @param angular_distribution_map: output map 
     *
     */

    //void FillAngularDistributionMap(const std::vector<art::Ptr<recob::SpacePoint>>& SPListROI, const TVector3 KEnd, AngularDistribution3DMap& AngularDistributionMap) const;
    void fill_angular_distribution_map(const SPList& sp_list_roi, const TVector3 k_end, AngularDistribution3DMap& angular_distribution_map) const;

    /*
     * @brief  Smooth out the angular distribution map
     *
     * @param angular_distribution_map: output map
   *
   */

    void smooth_angular_distribution_map(AngularDistribution3DMap& angular_distribution_map) const;

    /*
     * @brief  Obtain the directions from angular distribution peaks
     * 
     * @param angular_distribution_map: output map 
     * @param sort_peak_direction_map: vector of obtained peak directions sorted by peak height from high to low
     * 
     */

    //void retrieve_peak_directions(const angular_distribution_map_3d& angular_distribution_map, std::vector<TVector2>& peak_direction_vectors) const;
    void retrieve_peak_directions(const AngularDistribution3DMap& angular_distribution_map, std::map<double, TVector3, std::greater<>>& sort_peak_direction_map) const;

    /*
     * @brief  Remove some peak candidates depending on their height or opnening angles
     * 
     * @param sort_peak_direction_map: vector of obtained peak directions sorted by peak height from high to low
     * @param peak_direction_vector: vector of peak direction vector
     * 
     */

    void refine_peak_directions(const std::map<double, TVector3, std::greater<>>& sort_peak_direction_map, std::vector<TVector3> &peak_direction_vector) const;

    float m_region_of_interest;
    float m_peak_search_region;
    float m_theta_bin_size;
    float m_phi_bin_size;
    int m_smoothing_window;
    int m_peak_search_window;
    double m_peak_open_angle;
    double m_min_peak_height;
    size_t m_max_num_peak;

    SPList sp_list_roi;
    TVector3 k_end;
    TVector3 k_vtx;

  };


  //-------------------------------------


  ParticleDirectionFinder::ParticleDirectionFinder() :

    m_region_of_interest(100.),
    m_peak_search_region(20.),//30
    m_theta_bin_size(0.06),
    m_phi_bin_size(0.06),
    m_smoothing_window(1),
    m_peak_search_window(1),
    m_peak_open_angle(TMath::Pi() * 1/4),
    m_min_peak_height(0.4),
    m_max_num_peak(4)
  {
  }


  //------------------------------------------------------------------------------------------------------------------------------------------

  float ParticleDirectionFinder::get_theta_bin_size() const { 
    return m_theta_bin_size;
  }

  float ParticleDirectionFinder::get_phi_bin_size() const { 
    return m_phi_bin_size;
  }

  ParticleDirectionFinder::SPList& ParticleDirectionFinder::get_sp_list_roi() { 
    return sp_list_roi; 
  }

  const TVector3& ParticleDirectionFinder::get_k_end() const { 
    return k_end; 
  }

  const TVector3& ParticleDirectionFinder::get_k_vtx() const { 
    return k_vtx; 
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  pandora::StatusCode ParticleDirectionFinder::Run(const SPList& sp_list, const art::Ptr<recob::Track> primary_track,  const HitList& unavailable_hit_list, const std::map<art::Ptr<recob::SpacePoint>, art::Ptr<recob::Hit>>& spacepointToHitMap, std::vector<TVector3> &peak_direction_vector)
  {

    //get coordinates of k track end and store unavailable_hit_list
    k_end.SetXYZ(primary_track->End().x(), primary_track->End().y(), primary_track->End().z());
    k_vtx.SetXYZ(primary_track->Vertex().x(), primary_track->Vertex().y(), primary_track->Vertex().z());

    int nSP_roi_end=0;
    int nSP_roi_vtx=0;

    // get sp list inside region of interest
    this->collect_sp_in_roi(sp_list, k_end, m_region_of_interest, sp_list_roi, unavailable_hit_list, spacepointToHitMap);

    if (sp_list_roi.empty())
      return STATUS_CODE_NOT_FOUND;

    // get sp list for peak finder
    SPList sp_list_peak_search;
    this->collect_sp_in_roi(sp_list_roi, k_end, m_peak_search_region, sp_list_peak_search, unavailable_hit_list, spacepointToHitMap);
    nSP_roi_end = sp_list_peak_search.size();

    sp_list_peak_search.clear();
    this->collect_sp_in_roi(sp_list_roi, k_vtx, m_peak_search_region, sp_list_peak_search, unavailable_hit_list, spacepointToHitMap);
    nSP_roi_vtx = sp_list_peak_search.size();

    if(nSP_roi_vtx > nSP_roi_end) k_end = k_vtx;

    if (sp_list_peak_search.empty())
      return STATUS_CODE_NOT_FOUND;

    // Fill angular distribution map
    AngularDistribution3DMap angular_distribution_map;
    this->fill_angular_distribution_map(sp_list_roi, k_end, angular_distribution_map);

    if (angular_distribution_map.empty())
      return STATUS_CODE_NOT_FOUND;

    // Smooth angular decomposition map
    this->smooth_angular_distribution_map(angular_distribution_map);

    // Store peaks into map from highest to lowest
    std::map<double, TVector3, std::greater<>> sort_peak_direction_map;
    this->retrieve_peak_directions(angular_distribution_map, sort_peak_direction_map);

    if(sort_peak_direction_map.empty())
      return STATUS_CODE_NOT_FOUND;

    // Get vector of peak directions
    this->refine_peak_directions(sort_peak_direction_map, peak_direction_vector);

    if(peak_direction_vector.empty())
      return STATUS_CODE_NOT_FOUND;

    return STATUS_CODE_SUCCESS;
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void ParticleDirectionFinder::collect_sp_in_roi(const SPList& sp_list, const TVector3& k_end, float& region_of_interest, SPList& sp_list_roi, const HitList& unavailable_hit_list, const std::map<art::Ptr<recob::SpacePoint>, art::Ptr<recob::Hit>>& spacepointToHitMap) const
  {

    for(auto it_sp = sp_list.begin(); it_sp != sp_list.end(); ++it_sp){

      const TVector3 hit_position = (*it_sp)->XYZ();
      const TVector3 displacement_vector = hit_position - k_end;

      auto corresponding_hit = spacepointToHitMap.at(*it_sp);
      //skip hits in primary track 
      if (std::find(unavailable_hit_list.begin(), unavailable_hit_list.end(), corresponding_hit) != unavailable_hit_list.end()){
	//cout << "THIS IS A HIT IN PRIMARY TRACK: " << hit_position.X() << " " << hit_position.Y() << " " << hit_position.Z() << endl;
	continue;
      }

      if( displacement_vector.Mag() > region_of_interest) continue;
      sp_list_roi.push_back(*it_sp);

    }

  }

  //-----------------------------------------------------------------------------

  void ParticleDirectionFinder::fill_angular_distribution_map(const std::vector<art::Ptr<recob::SpacePoint>>& sp_list_roi, const TVector3 k_end, AngularDistribution3DMap& angular_distribution_map) const
  {

    //const TVector3 x_axis(1.,0.,0.);

    for (const auto& sp : sp_list_roi) {

      const TVector3 hit_position = sp->XYZ();
      const TVector3 displacement_vector = hit_position - k_end;

      double theta = displacement_vector.Theta();
      double phi = displacement_vector.Phi();
      int theta_factor = static_cast<int>(std::floor(theta / m_theta_bin_size));
      int phi_factor = static_cast<int>(std::floor(phi / m_phi_bin_size));

      // Using double brackets safely by checking or initializing correctly
      if (angular_distribution_map[theta_factor].find(phi_factor) == angular_distribution_map[theta_factor].end()) {
	// If phi_factor is not found under the current theta_factor, initialize it
	angular_distribution_map[theta_factor][phi_factor] = TMath::Sin(theta);
      } else {
	// If found, just add to the existing value
	angular_distribution_map[theta_factor][phi_factor] += TMath::Sin(theta);
      }

    }

    /*
    for(auto it_sp = sp_list_roi.begin(); it_sp != sp_list_roi.end(); ++it_sp){
      const TVector3 hit_position = (*it_sp)->XYZ();
      const TVector3 displacement_vector = hit_position - k_end;
      double theta = distance_vector.Theta();
      double phi = distance_vector.Phi();
      int theta_factor = (int)(std::floor(theta / m_theta_bin_size));
      int phi_factor = (int)(std::floor(phi / m_phi_bin_size));
      if( angular_distribution_map.find(theta_factor) == angular_distribution_map.end() &&
        angular_distribution_map[theta_factor].find(phi_factor) == angular_distribution_map[theta_factor].end() )
	angular_distribution_map[theta_factor][phi_factor] = TMath::Sin(theta); // weight by sintheta as this is 3d angular distribution
      else angular_distribution_map[theta_factor][phi_factor] += TMath::Sin(theta);
    }
    */

  }

  //-----------------------------------------------------------------------------

  void ParticleDirectionFinder::smooth_angular_distribution_map(AngularDistribution3DMap& angular_distribution_map) const
  {

    //const int loop_min = (-1) * (m_smoothing_window - 1) / 2;
    //const int loop_max = (m_smoothing_window - 1) / 2;
    const int loop_min = (-1) * m_smoothing_window;
    const int loop_max = m_smoothing_window;

    AngularDistribution3DMap angular_distribution_map_temp( angular_distribution_map );
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

  void ParticleDirectionFinder::retrieve_peak_directions(const AngularDistribution3DMap& angular_distribution_map, std::map<double, TVector3, std::greater<>>& sort_peak_direction_map) const

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
	    TVector3 peak_direction;
	    peak_direction.SetMagThetaPhi(1, bin_theta*m_theta_bin_size, bin_phi*m_phi_bin_size);
	    sort_peak_direction_map[bin_weight] = peak_direction; // Peaks are sorted by weight
	  } 
	}
      }   
    }
  }

  //-----------------------------------------------------------------------------

  void ParticleDirectionFinder::refine_peak_directions(const std::map<double, TVector3, std::greater<>>& sort_peak_direction_map, std::vector<TVector3> &peak_direction_vector) const
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

      if (peak_direction_vector.size() >= m_max_num_peak) break;

    }

  }



} // end of namespace kaon_reconstruction

#endif // #ifndef PARTICLE_DIRECTION_FINDER
