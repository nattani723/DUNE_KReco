/*

 * Header file for particle direction finder tool class.
 * Plots angular distribution of hits insider the region of interest,
 * and find peaks as candidates for daughter particle direction.

 */

#ifndef PARTICLE_DIRECTION_FINDER
#define PARTICLE_DIRECTION_FINDER 1

#include <bits/c++config.h>
#undef _GLIBCXX14_CONSTEXPR
#define _GLIBCXX14_CONSTEXPR

//#include "HitSplitAlg_module.h"
#include "TMath.h"
#include <algorithm>
#include <functional>

#include "Api/PandoraApi.h"
#include "Pandora/AlgorithmHeaders.h"
#include "larpandoracontent/LArContent.h"


//using namespace pandora;
//using namespace std;

namespace kaon_reconstruction
{

  class ParticleDirectionFinder
  {
    
  public:

    ParticleDirectionFinder();
    
    float get_theta_bin_size() const;
    
    float get_phi_bin_size() const;    

    
  private:
    
    typdef std::vector<art::Ptr<recob::SpacePoint>> SPList;
    typdef std::map<int, std::map<int, double>> AngularDistributionMap3D;

  /*
   * @brief  Collect spacepoints inside the region of interest (ROI)
   *
   * @param sp_list: spacepoint list of K+ daughter candidate hits
   * @param k_end: vector of K+ track candidate's end
   * @param sp_list_roi: spacepoint list inside the ROI
   *
   */

  //void collect_sp_in_roi(const SPList& sp_list, const TVector3 k_end, SPList& sp_list_roi) const;
  void collect_sp_in_roi(const SPList& sp_list, const TVector3& k_end, double& region_of_interest, SPList& sp_list_roi) const;

  /*
   * @brief  Fill map of angular distribution for spacepoints inside the region of interest (ROI)
   *
   * @param sp_list_roi: spacepoint list inside the ROI
   * @param k_end: vector of K+ track candidate's end
   * @param angular_distribution_map: output map 
   *
   */

  //void FillAngularDistributionMap(const std::vector<art::Ptr<recob::SpacePoint>>& SPListROI, const TVector3 KEnd, AngularDistributionMap3D& AngularDistributionMap) const;
  void fill_angular_distribution_map(const std::vector<art::Ptr<recob::SpacePoint>>& sp_list_roi, const TVector3 k_end, AngularDistributionMap3D& angular_distribution_map) const;

  /*
   * @brief  Smooth out the angular distribution map
   *
   * @param angular_distribution_map: output map
   *
   */

  void smooth_angular_distribution_map(AngularDistributionMap3D& angular_distribution_map) const;

  /*
   * @brief  Obtain the directions from angular distribution peaks
   * 
   * @param angular_distribution_map: output map 
   * @param sort_peak_direction_map: vector of obtained peak directions sorted by peak height from high to low
   * 
   */

  //void retrieve_peak_directions(const angular_distribution_map_3d& angular_distribution_map, std::vector<TVector2>& peak_direction_vectors) const;
  void retrieve_peak_directions(const angular_distribution_map_3d& angular_distribution_map, std::map<double, TVector3, std::greater<>>& sort_peak_direction_map) const;

  /*
   * @brief  Remove some peak candidates depending on their height or opnening angles
   * 
   * @param sort_peak_direction_map: vector of obtained peak directions sorted by peak height from high to low
   * @param peak_direction_vector: vector of peak direction vector
   * 
   */

  void refine_peak_directions(const std::map<double, TVector3, std::greater<>>& sort_peak_direction_map, vector<TVector3> &peak_direction_vector) const;


  //int m_peak_search_region;
  float m_theta_bin_size;
  float m_phi_bin_size;
  int m_smoothing_window;
  int m_peak_search_window;
  double m_peak_open_angle;
  double m_min_peak_height;
  int m_max_num_peak;

};


} // end of namespace kaon_reconstruction

#endif // #ifndef PARTICLE_DIRECTION_FINDER

