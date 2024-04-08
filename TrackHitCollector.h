/*
 *
 * Header file for track hit collector tool class.
 * Apply linear fit to candidate hits by using particle direction obtained by ParticleDirectionFinder
 *
 */

#ifndef TRACK_HIT_COLLECTOR
#define TRACK_HIT_COLLECTOR 1

#include "larpandoracontent/LArHelpers/LArHitWidthHelper.h"
#include "larpandoracontent/LArHelpers/LArConnectionPathwayHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"

namespace kaon_reconstruction 
{

  class TrackUtilities
  {

  public:
    /**
     * @brief Obtain wire pitch used for sliding linear fit
     *
     */
    static const float get_wire_pitch()
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

  }


  class TrackHitCollector
  {
    
  public:
    TrackHitCollector();

    pandora::StatusCode Run(const TVector3& k_end, const SPList& sp_list,
        const TVector3& peak_direction, HitList& unavailable_hit_list, HitList& track_hit_list);

  private:

    typdef std::vector<art::Ptr<recob::SpacePoint>> SPList;
    typdef std::vector<art::Ptr<recob::Hit>> HitList;


    //const double get_wire_pitch();

    /**
     * @brief Perform a running fit to collect the hits of the daughter track
     *
     * @param sp_list: spacepoint list of K+ daughter candidate hits
     * @param unavailable_hit_list: protected hits that cannot be collected
     * @param track_hit_list: the output list of daughter track hits
     * @param k_end: vector of K+ track candidate's end
     * @param peak_direction: vector of peak direction
     *
     */
    
    void find_track_hits(const SPList& sp_list, HitList& unavailable_hit_list, HitList& track_hit_list, const TVector3& k_end, TVector3& peak_direction) const;

    /**
     *  @brief  Reassign extrapolated direction and start/end positions
     *
     *  @param  count: n-th step path (how many times it has been looping)
     *  @param  extrapolated_fit: the fit to the collected hits on each step path
     *  @param  extrapolated_start_position: the track step path projection start position
     *  @param  extrapolated_end_position the track step path projection end position
     *  @param  extrapolated_direction: the track step path projection direction
     *  @param  is_end_downstream whether the shower direction is downstream (in Z) of the kaon track end
     *  @param  sp_list: spacepoint list of K+ daughter candidate hits 
     *  @param  running_fit_position_vector: the vector of the collected hit positions
     *  @param  pandora_running_fit_position_vector: CartesianPointVector of the collected hit positions
     *  @param  unavailable_hit_list: protected hits that cannot be collected
     *  @param  track_hit_list: the output list of daughter track hits 
     *
     */

    void update_extrapolation(int count, const lar_content::ThreeDSlidingFitResult& extrapolated_fit, const TVector3& extrapolated_start_position, const TVector3& extrapolated_end_position, const TVector3& extrapolated_direction, const bool is_end_downstream, const SPList& sp_list, TVector3& running_fit_position_vector, pandora::CartesianPointVector& pandora_running_fit_position_vector, HitList& unavailable_hit_list, HitList& track_hit_list) const;

    /**
     *  @brief  Reassign extrapolated direction and start/end positions
     *
     *  @param  count: n-th step path (how many times it has been looping)
     *  @param  extrapolated_fit: the fit to the collected hits on each step path
     *  @param  extrapolated_start_position: the track step path projection start position
     *  @param  extrapolated_end_position the track step path projection end position
     *  @param  extrapolated_direction: the track step path projection direction
     *  @param  is_end_downstream whether the shower direction is downstream (in Z) of the kaon track end
     *
     */

    void update_extrapolation(int count, const lar_content::ThreeDSlidingFitResult& extrapolated_fit, const TVector3& extrapolated_start_position, const TVector3& extrapolated_end_position, const TVector3& extrapolated_direction, const bool is_end_downstream) const;


    /**
     *  @brief  Perform a running fit step: collect hits which lie close to the track step path projection
     *
     *  @param  extrapolated_fit: the fit to the collected hits on each step path
     *  @param  extrapolated_start_position: the track step path projection start position
     *  @param  extrapolated_end_position the track step path projection end position
     *  @param  extrapolated_direction: the track step path projection direction
     *  @param  is_end_downstream whether the shower direction is downstream (in Z) of the kaon track end
     *  @param  sp_list: spacepoint list of K+ daughter candidate hits 
     *  @param  running_fit_position_vector: the vector of the collected hit positions
     *  @param  pandora_running_fit_position_vector: CartesianPointVector of the collected hit positions
     *  @param  unavailable_hit_list: protected hits that cannot be collected
     *  @param  track_hit_list: the output list of daughter track hits 
     *  @param  distance_to_line: the comparison distance for 'is close'  
     *  @param  hit_connection_distanc: separation between connected hits
     *
     *  @return whether any hits were collected in the running fit step
     */

    bool collect_subsection_hits(const lar_content::ThreeDSlidingFitResult& extrapolated_fit, const TVector3& extrapolated_start_position, const TVector3& extrapolated_end_position, const TVector3& extrapolated_direction, const bool is_end_downstream, const SPList& sp_list, TVector3& running_fit_position_vector, pandora::CartesianPointVector& pandora_running_fit_position_vector, HitList& unavailable_hit_list, HitList& track_hit_list, float& distance_to_line, float& hit_connection_distance) const;


    /**
     *  @brief  Sort the list of sp by distance from 2nd parameter
     *
     *  @param  sp_list: spacepoint list
     *  @param  sort_position: take the distance between sp in list and this position
     *
     */
    void TrackHitCollector::sort_sp_by_distance(SPList& sp_list, const TVector3& sort_position) cons;t

    /**
     *  @brief  Determine whether a hit lies close to the track projection
     *
     *  @param  hit_position: the hit position
     *  @param  line_start: the track projection start position
     *  @param  line_direction: the track projection direction
     *  @param  distance_to_line: the comparison distance for 'is close'
     *
     *  @return whether the hit is close to the track projection
     */

    bool is_close_to_line(const TVector3& hit_position, const TVector3& line_start, const TVector3& line_direction, const double& distance_to_line) const;

    /**
     *  @brief  Add the connecting hits to track
     *
     *  @param  collected_hits: the input list of close hits
     *  @param  extrapolated_start_position: the track projection start position
     *  @param  extrapolated_direction: the track projection direction
     *  @param  running_fit_position_vector: the vector of the collected hit positions
     *  @param  pandora_running_fit_position_vector: CartesianPointVector of the collected hit positions
     *  @param  track_hit_list: the list of collected track hits
     *  @param  distance_to_line: the comparison distance for 'is close'
     *
     */

    void collect_connected_hits(HitList& collected_hit_list, const TVector3& extrapolated_start_position, const TVector3& extrapolated_direction, TVector3& running_fit_position_vector, pandora::CartesianPointVector& pandora_running_fit_position_vector, HitList& track_hit_list, float& hit_connection_distance) const;

    /**
     *  @brief  Find the smallest distance between a position and a list of other positions
     *
     *  @param  position: the input position
     *  @param  test_positions the list of other positions
     *
     *  @return the closest distance
     */

    double get_closest_distance(const TVector3& hit_position, const vector<TVector3>& test_positions) const;


    unsigned int m_hit_threshold_for_track;           ///< The hit threshold for a significant track
    float m_growing_fit_initial_length;               ///< The first step distance
    float m_initial_fit_distance_to_line;              ///< The max. proximity to the track projection for collection in the first step
    unsigned int m_min_initial_hits_found;            ///< The min. number of hits collected in the first step for continuation
    unsigned int m_max_fitting_hits;                 ///< The number of hits to consider in the running fit
    unsigned int m_local_sliding_fit_window;          ///< The standard sliding fit window for track fits
    float m_growing_fit_segment_length;               ///< The standard step distance
    unsigned int m_high_resolution_sliding_fit_window; ///< The high resolution sliding fit window for track fits
    float m_distance_to_line;                        ///< The max. proximity to the spine projection for collection
    float m_hit_connection_distance;                 ///< The max. separation between connected hits
    float m_trackall_sliding_fit_window;

  } // end of class
  
} // namespace kaon_reconstruction 

#endif // #ifndef TRACK_HIT_COLLECTOR
