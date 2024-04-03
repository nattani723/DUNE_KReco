/*

 * Header file for particle direction finder tool class.
 * Plots angular distribution of hits insider the region of interest,
 * and find peaks as candidates for daughter particle direction.

 */

#ifndef PARTICLE_DIRECTION_FINDER
#define PARTICLE_DIRECTION_FINDER

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

#include "larpandoracontent/LArHelpers/LArHitWidthHelper.h"
#include "larpandoracontent/LArHelpers/LArConnectionPathwayHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"

//using namespace pandora;
//using namespace std;

namespace kaon_reconstruction
{

class ParticleDirectionFinder
{
 
public:
  /*
   * @brief Default constructor
   */
  ParticleDirectionFinder();

  //Main functionality

 private:

  typdef std::map<int, std::map<int, double>> AngularDistributionMap3D;

  /*
   * @brief  Collect spacepoints inside the region of interest (ROI)
   *
   * @param SPList: spacepoint list of K+ daughter candidate hits
   * @param KEnd: vector of K+ track candidate's end
   * @param SPListROI: spacepoint list inside the ROI
   *
   */

  void CollectSPsInROI(const std::vector<art::Ptr<recob::SpacePoint>>& SPList, const TVector3 KEnd, std::vector<art::Ptr<recob::SpacePoint>>& SPListROI) const;

  /*
   * @brief  Fill map of angular distribution for spacepoints inside the region of interest (ROI)
   *
   * @param SPListROI: spacepoint list inside the ROI
   * @param KEnd: vector of K+ track candidate's end
   * @param AngularDistributionMap: output map 
   *
   */

  void FillAngularDistributionMap(const std::vector<art::Ptr<recob::SpacePoint>>& SPListROI, const TVector3 KEnd, AngularDistributionMap3D& AngularDistributionMap) const;

  /*
   * @brief  Smooth out the angular distribution map
   *
   * @param AngularDistributionMap: output map
   *
   */

  void SmoothAngularDistributionMap(AngularDistributionMap3D& AngularDistributionMap) const;

  /*
   * @brief  Obtain the directions from angular distribution peaks
   * 
   * @param AngularDistributionMap: output map 
   * @param PeakDirectionVectors: vector of obtained peak directions
   * 
   */

  void RetrievePeakDirections(const AngularDistributionMap3D& AngularDistributionMap, std::vector<TVector2>& PeakDirectionVectors) const;

  float m_ThetaBinSize;
  float m_PhiBinSize;
  int m_SmoothingWindow;

};

} // end of namespace kaon_reconstruction

#endif // #ifndef PARTICLE_DIRECTION_FINDER





  //------------------------------------------------------------------------------------------------------------------------------------------                                                                                                                                                                                

  void HitSplitAlg::findShowerSpine3D(std::vector<art::Ptr<recob::SpacePoint>>& sp_from_pfparticle,
                                         std::vector<art::Ptr<recob::Hit>>& unavailable_hit_list,
                                         std::vector<std::vector<art::Ptr<recob::Hit>>>& shower_spine_hit_list_vector,
                                         TVector3 Kend_candidate,
                                         std::map<double, TVector2, std::greater<>>& view_peak_map,
                                         vector<TVector2> &best_peak_bins){

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

    if(best_peak_bins.size()==0 || sp_from_pfparticle.size()==0) return;
 
   for(auto& entry : best_peak_bins){

      double highest_l = 0.;
      vector<TVector3> running_fit_position_vec;
      pandora::CartesianPointVector pandora_running_fit_position_vec;
      TVector3 initial_direction;
      std::vector<art::Ptr<recob::Hit>> shower_spine_hit_list;

      double theta = entry.X()*thetaBinSize;
      double phi = entry.Y()*phiBinSize;
   
      initial_direction.SetMagThetaPhi(1, theta, phi);
      cout << "theta: " << theta << ", phi: " << phi << endl;
      cout << "initial_direction: " << initial_direction.X() << " " << initial_direction.Y() << " " << initial_direction.Z() << endl;


      cout << "sp_from_pfparticle.size(): " << sp_from_pfparticle.size() << endl;

      sort(sp_from_pfparticle.begin(), sp_from_pfparticle.end());
      sp_from_pfparticle.erase( unique( sp_from_pfparticle.begin(), sp_from_pfparticle.end() ), sp_from_pfparticle.end() );

      for (auto spIter = sp_from_pfparticle.begin(); spIter != sp_from_pfparticle.end(); ++spIter) {

	const TVector3 hit_position = (*spIter)->XYZ();
	const pandora::CartesianVector pandora_hit_position(hit_position.X(), hit_position.Y(), hit_position.Z());
	const TVector3 distance_vector = hit_position - Kend_candidate;
	const double l = initial_direction.Dot(distance_vector);
	const double t = initial_direction.Cross(distance_vector).Mag();

	//if(spIter != sp_from_pfparticle.begin() )
	//if( (*spIter)->XYZ() == (*prev(spIter))->XYZ() ) continue;

	//cout << "l: " << l << ", t: " << t << ", x: "  << hit_position.X() << ", y: " << hit_position.Y() << ", z: " << hit_position.Z() << endl;

	//cout << "growing_fit_initial_length" << growing_fit_initial_length << endl;

	if( (l<growing_fit_initial_length) && (l>0.) && (t<initial_fit_distance_to_line) ){  

	  if(l>highest_l) highest_l = l;
	  //if(std::find(unavailable_hit_list.begin(), unavailable_hit_list.end(), fSpacePointsToHits.at(*spIter)) == unavailable_hit_list.end()){
	  if(std::find(unavailable_hit_list.begin(), unavailable_hit_list.end(), fSpacePointsToHits_old.at(*spIter)) == unavailable_hit_list.end()){
	    //if(std::find(unavailable_hit_list.begin(), unavailable_hit_list.end(), fSpacePointsToHits.at(*spIter)) == unavailable_hit_list.end()){
 	    //shower_spine_hit_list.push_back(fSpacePointsToHits.at(*spIter));
 	    shower_spine_hit_list.push_back(fSpacePointsToHits_old.at(*spIter));
	  }

	  running_fit_position_vec.push_back(hit_position);
	  pandora_running_fit_position_vec.push_back(pandora_hit_position);
	  //cout << "HITS passed the first selection" << endl;
	  //cout << "l: " << l << ", t: " << t << endl;
	}
      }

      //cout << "running_fit_position_vec.size(): " << running_fit_position_vec.size() << endl;
      //cout << "pandora_running_fit_position_vec.size(): " << pandora_running_fit_position_vec.size() << endl;

	// Require significant number of initial hits                                                                                                        
                     
	if (running_fit_position_vec.size() < min_initial_hits_found){
	  cout << "Requiring significant number of initial hits" << endl;
	  shower_spine_hit_list.clear();
	  shower_spine_hit_list_vector.push_back(shower_spine_hit_list);
	  continue;
	}

   
	// Perform a running fit to follow pathway
	
	bool is_end_downstream = false;
	if(initial_direction.Z() > 0.) is_end_downstream = true;
	TVector3 extrapolated_direction;
     
	extrapolated_direction.SetMagThetaPhi(1, theta, phi);
	TVector3 extrapolated_start_position = Kend_candidate;
	TVector3 extrapolated_end_position = extrapolated_start_position + (extrapolated_direction * highest_l);

	if( find(pandora_running_fit_position_vec.begin(),  pandora_running_fit_position_vec.end(), pandora_running_fit_position_vec[0] ) != pandora_running_fit_position_vec.end()){
	  //cout << "this hit is suspecious" << endl;
	  //continue;
	}

	const pandora::CartesianVector pandora_vertex_position(extrapolated_start_position.X(), extrapolated_start_position.Y(), extrapolated_start_position.Z());
	lar_content::LArConnectionPathwayHelper::SortByDistanceToPoint vtx_distance(pandora_vertex_position);
	std::sort(pandora_running_fit_position_vec.begin(), pandora_running_fit_position_vec.end(), vtx_distance);
	TVector3 closest_hit;
	closest_hit.SetXYZ(pandora_running_fit_position_vec[0].GetX(), pandora_running_fit_position_vec[0].GetY(), pandora_running_fit_position_vec[0].GetZ());
	//cout << "closest distance from Kend_candidate is " << (Kend_candidate - closest_hit).Mag() << endl;

	/*
	if( (Kend_candidate - closest_hit).Mag()>10 ){//temporary
	  cout << "called (Kend_candidate - closest_hit).Mag()>10" << endl;
	  shower_spine_hit_list.clear();
	  shower_spine_hit_list_vector.push_back(shower_spine_hit_list);
	  continue;
	}
	*/
   
	unsigned int count = 0;
	bool hits_collected = true;

	while (hits_collected){

	  ++count;

	  const int excess_hits_in_fit = pandora_running_fit_position_vec.size() - max_fitting_hits;
	  cout << "excess_hits_in_fit: " << excess_hits_in_fit << endl;

          const pandora::CartesianVector pandora_extrapolated_end_position(extrapolated_end_position.X(), extrapolated_end_position.Y(), extrapolated_end_position.Z());
	  lar_content::LArConnectionPathwayHelper::SortByDistanceToPoint smpl(pandora_extrapolated_end_position);
	  std::sort(pandora_running_fit_position_vec.begin(), pandora_running_fit_position_vec.end(), smpl);

	  if(excess_hits_in_fit>0){

	    // Remove furthest away hits
	    //cout << "before loop pandora_running_fit_position_vec.size(): " << pandora_running_fit_position_vec.size() << endl;
	    for (int i = 0; i < excess_hits_in_fit; ++i){
	      //running_fit_position_vec.erase(std::prev(running_fit_position_vec.end()));
	      pandora_running_fit_position_vec.erase(std::prev(pandora_running_fit_position_vec.end()));
	    }

	    running_fit_position_vec.clear();
	    TVector3 hit_position_tmp;
	    for(auto pandora_running_fit_position : pandora_running_fit_position_vec){
	      hit_position_tmp.SetXYZ(pandora_running_fit_position.GetX(), pandora_running_fit_position.GetY(), pandora_running_fit_position.GetZ());
	      running_fit_position_vec.push_back(hit_position_tmp);
	    }
	  }
	  //cout << "after loop pandora_running_fit_position_vec.size(): " << pandora_running_fit_position_vec.size() << endl;                                                                                                                                           
	  distance_to_line = distance_to_line_default;
	  hit_connection_distance = hit_connection_distance_default;
	  
	  
	  for(auto pandora_running_fit_position : pandora_running_fit_position_vec){
	    cout << pandora_running_fit_position.GetX() << " " << pandora_running_fit_position.GetY() << " " << pandora_running_fit_position.GetZ() << endl;
	  }
	  
	  
	  
	  const lar_content::ThreeDSlidingFitResult extrapolated_fit(&pandora_running_fit_position_vec, macro_sliding_fit_window, sliding_fit_pitch);

	  //cout << "GetGlobalMaxLayerPosition(): " << extrapolated_fit.GetGlobalMaxLayerPosition().GetX() << " " << extrapolated_fit.GetGlobalMaxLayerPosition().GetY() << " " << extrapolated_fit.GetGlobalMaxLayerPosition().GetZ() << endl;
	  //cout << "GetGlobalMinLayerPosition(): " << extrapolated_fit.GetGlobalMinLayerPosition().GetX() << " " << extrapolated_fit.GetGlobalMinLayerPosition().GetY() << " " << extrapolated_fit.GetGlobalMinLayerPosition().GetZ() << endl;

	  cout << "is_end_downstream: " << is_end_downstream << endl;
	  if(count == 1){
	    extrapolated_start_position = extrapolated_end_position;
	  }else{
	    if(is_end_downstream){
	      extrapolated_start_position.SetXYZ(extrapolated_fit.GetGlobalMaxLayerPosition().GetX(), extrapolated_fit.GetGlobalMaxLayerPosition().GetY(), extrapolated_fit.GetGlobalMaxLayerPosition().GetZ());
	      //cout << "GetGlobalFitPosition(): " << extrapolated_fit.GetGlobalFitPosition().GetX() << " " << endl;
	    }else{
	      extrapolated_start_position.SetXYZ(extrapolated_fit.GetGlobalMinLayerPosition().GetX(), extrapolated_fit.GetGlobalMinLayerPosition().GetY(), extrapolated_fit.GetGlobalMinLayerPosition().GetZ());
	    }
	  }

	  if(is_end_downstream){
	    extrapolated_direction.SetXYZ(extrapolated_fit.GetGlobalMaxLayerDirection().GetX(), extrapolated_fit.GetGlobalMaxLayerDirection().GetY(), extrapolated_fit.GetGlobalMaxLayerDirection().GetZ());
	    if(initial_direction.Angle(extrapolated_direction) > TMath::Pi()*2/3){
	      cout << "THIS IS CALLED" << endl;
	      extrapolated_start_position.SetXYZ(extrapolated_fit.GetGlobalMinLayerPosition().GetX(), extrapolated_fit.GetGlobalMinLayerPosition().GetY(), extrapolated_fit.GetGlobalMinLayerPosition().GetZ());
	      extrapolated_direction *= -1.;
	      //is_end_downstream = false;
	    }

	  }else{
	    extrapolated_direction.SetXYZ(extrapolated_fit.GetGlobalMinLayerDirection().GetX()*(-1.), extrapolated_fit.GetGlobalMinLayerDirection().GetY()*(-1.), extrapolated_fit.GetGlobalMinLayerDirection().GetZ()*(-1.));
	    if(initial_direction.Angle(extrapolated_direction) > TMath::Pi()*2/3){
	      cout << "THIS IS CALLED" << endl;
	      extrapolated_start_position.SetXYZ(extrapolated_fit.GetGlobalMaxLayerPosition().GetX(), extrapolated_fit.GetGlobalMaxLayerPosition().GetY(), extrapolated_fit.GetGlobalMaxLayerPosition().GetZ());
	      extrapolated_direction *= -1.;
	      //is_end_downstream = true;
	    }
	  }

	  extrapolated_end_position = extrapolated_start_position + (extrapolated_direction * growing_fit_segment_length);

	  //cout << "initial_direction.Angle(extrapolated_direction): " << initial_direction.Angle(extrapolated_direction) << endl;
	  cout << "extrapolated_direction for count " << count << " is " << extrapolated_direction.X() << " " << extrapolated_direction.Y() << " "  << extrapolated_direction.Z() << endl;

	  hits_collected = this->collectSubsectionHits(extrapolated_fit, extrapolated_start_position, extrapolated_end_position, extrapolated_direction, is_end_downstream, sp_from_pfparticle, running_fit_position_vec, pandora_running_fit_position_vec, unavailable_hit_list, shower_spine_hit_list);

	  //cout << "are hits collected? " << hits_collected << endl;

	  ///// as a final effort, fit with all spine hit list                                                                                                
   
	  if(!hits_collected){

	    distance_to_line = distance_to_line_spineall;
	    hit_connection_distance = hit_connection_distance_spineall;

	    pandora::CartesianPointVector pandora_running_fit_position_spineall_vec;

	    for(size_t i_h=0; i_h<shower_spine_hit_list.size(); i_h++){
	      /*
	      spacepoint_vec.clear();
	      spacepoint_vec = spacepoint_per_hit.at(shower_spine_hit_list[i_h].key());
	      if(spacepoint_vec.size()!=1) continue;
	      */
	      //const TVector3 hit_position = fHitsToSpacePoints.at(shower_spine_hit_list[i_h])->XYZ();
	      const TVector3 hit_position = fHitsToSpacePoints_old.at(shower_spine_hit_list[i_h])->XYZ();

	      //const TVector3 hit_position = spacepoint_vec[0]->XYZ();
	      const pandora::CartesianVector pandora_hit_position(hit_position.X(), hit_position.Y(), hit_position.Z());
	      pandora_running_fit_position_spineall_vec.push_back(pandora_hit_position);
	    }

	    const lar_content::ThreeDSlidingFitResult spineall_extrapolated_fit(&pandora_running_fit_position_spineall_vec, 20, sliding_fit_pitch);

	    if(count == 1)
	      {
		extrapolated_start_position = extrapolated_start_position;
	      }else{
	      if(is_end_downstream){
		extrapolated_start_position.SetXYZ(spineall_extrapolated_fit.GetGlobalMaxLayerPosition().GetX(), spineall_extrapolated_fit.GetGlobalMaxLayerPosition().GetY(), spineall_extrapolated_fit.GetGlobalMaxLayerPosition().GetZ());
	      }else{
		extrapolated_start_position.SetXYZ(spineall_extrapolated_fit.GetGlobalMinLayerPosition().GetX(), spineall_extrapolated_fit.GetGlobalMinLayerPosition().GetY(), spineall_extrapolated_fit.GetGlobalMinLayerPosition().GetZ());
	      }
	    }

            if(is_end_downstream){
              extrapolated_direction.SetXYZ(spineall_extrapolated_fit.GetGlobalMaxLayerDirection().GetX(), spineall_extrapolated_fit.GetGlobalMaxLayerDirection().GetY(), spineall_extrapolated_fit.GetGlobalMaxLayerDirection().GetZ());
              if(initial_direction.Angle(extrapolated_direction) > TMath::Pi()*2/3){
                extrapolated_direction *= -1.;
                is_end_downstream = false;
              }
            }else{
              extrapolated_direction.SetXYZ(spineall_extrapolated_fit.GetGlobalMinLayerDirection().GetX()*(-1.), spineall_extrapolated_fit.GetGlobalMinLayerDirection().GetY()*(-1.), spineall_extrapolated_fit.GetGlobalMinLayerDirection().GetZ()*(-1.));
              if(initial_direction.Angle(extrapolated_direction) > TMath::Pi()*2/3){
                extrapolated_direction *= -1.;
                is_end_downstream = true;
              }
            }

            extrapolated_end_position = extrapolated_start_position + (extrapolated_direction * growing_fit_segment_length);
	    //cout << "extrapolated_direction for count " << count << " is " << extrapolated_direction.X() << " " << extrapolated_direction.Y() << " "  << extrapolated_direction.Z() << endl;

            hits_collected = this->collectSubsectionHits(spineall_extrapolated_fit, extrapolated_start_position, extrapolated_end_position, extrapolated_direction, is_end_downstream, sp_from_pfparticle, running_fit_position_vec, pandora_running_fit_position_vec, unavailable_hit_list, shower_spine_hit_list);

	  }


	  // As a final effort, reduce the sliding fit window                                                                                                                                                                                                                                                        

	  if (!hits_collected)
	    {

	      const int excess_hits_in_fit_micro = pandora_running_fit_position_vec.size() - max_fitting_hits_micro;
	      //cout << "excess_hits_in_fit_micro: " << excess_hits_in_fit_micro << endl;

	      pandora::CartesianPointVector pandora_running_fit_position_micro_vec;
	      pandora_running_fit_position_micro_vec = pandora_running_fit_position_vec;

	      if(excess_hits_in_fit_micro>0){
		// Remove furthest away hits                                                                                                                                                                                                                                                                                
		for (int i = 0; i < excess_hits_in_fit_micro; ++i){
		  pandora_running_fit_position_micro_vec.erase(std::prev(pandora_running_fit_position_micro_vec.end()));
		}
	      }

	      distance_to_line = distance_to_line_micro;
	      hit_connection_distance = hit_connection_distance_micro;

	      const lar_content::ThreeDSlidingFitResult micro_extrapolated_fit(&pandora_running_fit_position_micro_vec, micro_sliding_fit_window, sliding_fit_pitch);
       

	      if(count == 1)
		{
		  extrapolated_start_position = extrapolated_start_position;
		}else{
		if(is_end_downstream){
		  extrapolated_start_position.SetXYZ(micro_extrapolated_fit.GetGlobalMaxLayerPosition().GetX(), micro_extrapolated_fit.GetGlobalMaxLayerPosition().GetY(), micro_extrapolated_fit.GetGlobalMaxLayerPosition().GetZ());
		}else{
		  extrapolated_start_position.SetXYZ(micro_extrapolated_fit.GetGlobalMinLayerPosition().GetX(), micro_extrapolated_fit.GetGlobalMinLayerPosition().GetY(), micro_extrapolated_fit.GetGlobalMinLayerPosition().GetZ());
		}
	      }

	      if(is_end_downstream){
		extrapolated_direction.SetXYZ(micro_extrapolated_fit.GetGlobalMaxLayerDirection().GetX(), micro_extrapolated_fit.GetGlobalMaxLayerDirection().GetY(), micro_extrapolated_fit.GetGlobalMaxLayerDirection().GetZ());
		if(initial_direction.Angle(extrapolated_direction) > TMath::Pi()*2/3){
		  extrapolated_direction *= -1.;
		  is_end_downstream = false;
		}

	      }else{
		extrapolated_direction.SetXYZ(micro_extrapolated_fit.GetGlobalMinLayerDirection().GetX()*(-1.), micro_extrapolated_fit.GetGlobalMinLayerDirection().GetY()*(-1.), micro_extrapolated_fit.GetGlobalMinLayerDirection().GetZ()*(-1.));
		if(initial_direction.Angle(extrapolated_direction) > TMath::Pi()*2/3){
		  extrapolated_direction *= -1.;
		  is_end_downstream = true;
		}

	      }

	      extrapolated_end_position = extrapolated_start_position + (extrapolated_direction * growing_fit_segment_length);
	      //cout << "extrapolated_direction for count " << count << " is " << extrapolated_direction.X() << " " << extrapolated_direction.Y() << " "  << extrapolated_direction.Z() << endl;

	      hits_collected = this->collectSubsectionHits(micro_extrapolated_fit, extrapolated_start_position, extrapolated_end_position, extrapolated_direction, is_end_downstream, sp_from_pfparticle, running_fit_position_vec, pandora_running_fit_position_vec, unavailable_hit_list, shower_spine_hit_list);
	    }

	}

	//if(shower_spine_hit_list.size()==0) cout << "shower_spine_hit_list.size() is 0" << endl;
	//cout << "the number of hits in shower spine is " << shower_spine_hit_list.size() << endl;
	if(shower_spine_hit_list.size()<30){
	  cout << "shower_spine_hit_list.size() is " << shower_spine_hit_list.size() << ", Requiring significant number of final hits" << endl;
	  shower_spine_hit_list.clear();
	  shower_spine_hit_list_vector.push_back(shower_spine_hit_list);
	  continue;
	}
	//cout << "shower_spine_hit_list.size() before pushing back to the vector is " << shower_spine_hit_list.size() << endl;
	shower_spine_hit_list_vector.push_back(shower_spine_hit_list);

      }
    }


    //----------------------------------------------------------------------------------------------------------------                                                                                                                                                                                                          


    bool HitSplitAlg::collectSubsectionHits(const lar_content::ThreeDSlidingFitResult &extrapolated_fit,
					       const TVector3 &extrapolated_start_position,
					       const TVector3 &extrapolated_end_position,
					       const TVector3 &extrapolated_direction,
					       const bool is_end_downstream,
					       std::vector<art::Ptr<recob::SpacePoint>>& sp_from_pfparticle,
					       vector<TVector3> &running_fit_position_vec,
					       pandora::CartesianPointVector &pandora_running_fit_position_vec,
					       std::vector<art::Ptr<recob::Hit>>& unavailable_hit_list,
					       std::vector<art::Ptr<recob::Hit>>& shower_spine_hit_list)
    {

      float extrapolated_start_l = 0.;
      float extrapolated_start_t = 0.;
      float extrapolated_start_t2 = 0.;

      pandora::CartesianVector pandora_extrapolated_start_position(extrapolated_start_position.X(), extrapolated_start_position.Y(), extrapolated_start_position.Z());
      pandora::CartesianVector pandora_extrapolated_end_position(extrapolated_end_position.X(), extrapolated_end_position.Y(), extrapolated_end_position.Z());

      extrapolated_fit.GetLocalPosition(pandora_extrapolated_start_position, extrapolated_start_l, extrapolated_start_t, extrapolated_start_t2);

      float extrapolated_end_l = 0.;
      float extrapolated_end_t = 0.;
      float extrapolated_end_t2 = 0.;
      extrapolated_fit.GetLocalPosition(pandora_extrapolated_end_position, extrapolated_end_l, extrapolated_end_t, extrapolated_end_t2);

      //cout << "extrapolated_start_position: " << extrapolated_start_position.X() << " " << extrapolated_start_position.Y() << " " << extrapolated_start_position.Z() << endl;
      //cout << "extrapolated_end_position: " << extrapolated_end_position.X() << " " << extrapolated_end_position.Y() << " " << extrapolated_end_position.Z() << endl;
      //cout << "extrapolated_start_l: " << extrapolated_start_l << endl;;
      //cout << "extrapolated_end_l: " << extrapolated_end_l << endl;;

      std::vector<art::Ptr<recob::Hit>> collected_hits;   
      //std::vector<art::Ptr<recob::SpacePoint>> spacepoint_vec;

      pandora::CartesianPointVector pandora_hits_from_pfparticle_position_vec;

      cout << "sp_from_pfparticle.size(): " << sp_from_pfparticle.size() << endl;
      for (auto spIter = sp_from_pfparticle.begin(); spIter != sp_from_pfparticle.end(); ++spIter) {
      /*
      for(size_t i_h=0; i_h<hits_from_pfparticle.size(); i_h++){
	spacepoint_vec.clear();
	spacepoint_vec = spacepoint_per_hit.at(hits_from_pfparticle[i_h].key());
	if(spacepoint_vec.size()==0) continue;
      */
	const TVector3 hit_position = (*spIter)->XYZ();
	//const TVector3 hit_position = spacepoint_vec[0]->XYZ();
	pandora::CartesianVector pandora_hit_from_pfparticle_position(hit_position.X(), hit_position.Y(), hit_position.Z());
	pandora_hits_from_pfparticle_position_vec.push_back(pandora_hit_from_pfparticle_position);
      }
      lar_content::LArConnectionPathwayHelper::SortByDistanceToPoint smpl(pandora_extrapolated_start_position);
      std::sort(pandora_hits_from_pfparticle_position_vec.begin(), pandora_hits_from_pfparticle_position_vec.end(), smpl);
   
      //std::vector<art::Ptr<recob::Hit>> hits_from_pfparticle_sort;
      std::vector<art::Ptr<recob::SpacePoint>> sp_from_pfparticle_sort;
      for(size_t i_p=0; i_p<pandora_hits_from_pfparticle_position_vec.size(); i_p++){
	TVector3 sort_position;
	sort_position.SetXYZ(pandora_hits_from_pfparticle_position_vec[i_p].GetX(), pandora_hits_from_pfparticle_position_vec[i_p].GetY(), pandora_hits_from_pfparticle_position_vec[i_p].GetZ());

	for (auto spIter = sp_from_pfparticle.begin(); spIter != sp_from_pfparticle.end(); ++spIter) {
	  const TVector3 hit_position = (*spIter)->XYZ();
	/*
	for(size_t i_h=0; i_h<hits_from_pfparticle.size(); i_h++){
	  spacepoint_vec.clear();
	  spacepoint_vec = spacepoint_per_hit.at(hits_from_pfparticle[i_h].key());
	  const TVector3 hit_position = spacepoint_vec[0]->XYZ();
	*/
	  if(sort_position == hit_position){
	    sp_from_pfparticle_sort.push_back(*spIter);
	    //hits_from_pfparticle_sort.push_back(hits_from_pfparticle[i_h]);
	    break;
	  }
	}
      }

      //for(size_t i_h=0; i_h<hits_from_pfparticle_sort.size(); i_h++){
      cout << "sp_from_pfparticle_sort.size(): " << sp_from_pfparticle_sort.size() << endl;;
      for (auto spIter = sp_from_pfparticle_sort.begin(); spIter != sp_from_pfparticle_sort.end(); ++spIter) {
	//spacepoint_vec.clear();

	//if(std::find(shower_spine_hit_list.begin(), shower_spine_hit_list.end(), hits_from_pfparticle_sort[i_h]) != shower_spine_hit_list.end()) continue;
	//if(std::find(unavailable_hit_list.begin(), unavailable_hit_list.end(), hits_from_pfparticle_sort[i_h]) != unavailable_hit_list.end()) continue;
   
	/*
	spacepoint_vec = spacepoint_per_hit.at(hits_from_pfparticle_sort[i_h].key());
	if(spacepoint_vec.size()!=1) continue;
	const TVector3 hit_position = spacepoint_vec[0]->XYZ();
	*/
	const TVector3 hit_position = (*spIter)->XYZ();
	pandora::CartesianVector pandora_hit_position(hit_position.X(), hit_position.Y(), hit_position.Z());

	float hit_l = 0.;
	float hit_t = 0.;
	float hit_t2 = 0.;
	extrapolated_fit.GetLocalPosition(pandora_hit_position, hit_l, hit_t, hit_t2);

	// Assess whether hit is within section boundaries                                                                                                                                                                                                                                        
	//cout << "hit_l: " << hit_l << ", extrapolated_start_l: " << extrapolated_start_l << ", extrapolated_end_l: " << extrapolated_end_l << endl;
	if(is_end_downstream && ((hit_l < extrapolated_start_l) || (hit_l > extrapolated_end_l))) continue;
	if(!is_end_downstream && ((hit_l > extrapolated_start_l) || (hit_l < extrapolated_end_l))) continue;
	//cout << "hit_l: " << hit_l << endl;

	// Assess whether hit is close to connecting line - taking account hit width if necessary                                                                                                                                                                                                                         
	//cout << hit_position.X() << " " << hit_position.Y() << " " << hit_position.Z() << endl;
	//cout << "Assess whether hit is close to connecting line" << this->isCloseToLine(hit_position, extrapolated_start_position, extrapolated_direction, distance_to_line)  << endl;

	//this->isCloseToLine(hit_position, extrapolated_start_position, extrapolated_direction, distance_to_line);

	if (this->isCloseToLine(hit_position, extrapolated_start_position, extrapolated_direction, distance_to_line)){
	  //cout << "called collected_hits.push_back" << endl;
	  //collected_hits.push_back(fSpacePointsToHits.at(*spIter));
	  collected_hits.push_back(fSpacePointsToHits_old.at(*spIter));
	}
      }

      const int n_initial_hits(shower_spine_hit_list.size());
      this->collectConnectedHits(collected_hits, extrapolated_start_position, extrapolated_direction, running_fit_position_vec, pandora_running_fit_position_vec, shower_spine_hit_list);
      const int n_final_hits(shower_spine_hit_list.size());

      //cout << "n_initial_hits: " << n_initial_hits << ", n_final_hits: " << n_final_hits << endl;
      //cout << "n_final_hits != n_initial_hits " << std::boolalpha(n_final_hits != n_initial_hits) << endl;
      return (n_final_hits != n_initial_hits);
    }


    //------------------------------------------------------------------------------------------------------------------------------------------  


    void HitSplitAlg::collectConnectedHits(std::vector<art::Ptr<recob::Hit>>& collected_hits,
					      const TVector3 &extrapolated_start_position,
					      const TVector3 &extrapolated_direction,
					      vector<TVector3> &running_fit_position_vec,
					      pandora::CartesianPointVector &pandora_running_fit_position_vec,
					      std::vector<art::Ptr<recob::Hit>>& shower_spine_hit_list)
    {
      // Now add connected hits

      std::vector<art::Ptr<recob::SpacePoint>> spacepoint_vec;
      bool found = true;

      while(found)
	{
	  found = false;

	  vector<TVector3> shower_spine_hit_list_position;
	  //std::vector<art::Ptr<recob::SpacePoint>> spacepoint_vec_2;

	  for(size_t i_h=0; i_h<shower_spine_hit_list.size(); i_h++){
	    //const TVector3 hit2_position = fHitsToSpacePoints.at(shower_spine_hit_list[i_h])->XYZ();
	    const TVector3 hit2_position = fHitsToSpacePoints_old.at(shower_spine_hit_list[i_h])->XYZ();

	    //spacepoint_vec_2 = spacepoint_per_hit.at(shower_spine_hit_list[i_h].key());
	    //const TVector3 hit2_position = spacepoint_vec_2[0]->XYZ();
	    shower_spine_hit_list_position.push_back(hit2_position);
	  }

	  cout << "collected_hits.size(): " << collected_hits.size() << endl;
	  for(size_t i_h=0; i_h<collected_hits.size(); i_h++){


	    //if(std::find(shower_spine_hit_list.begin(), shower_spine_hit_list.end(), collected_hits[i_h]) != shower_spine_hit_list.end()) cout << "THIS HIT IS ALREADY INSIDE SHOWER SPINE\n";
	    if(std::find(shower_spine_hit_list.begin(), shower_spine_hit_list.end(), collected_hits[i_h]) != shower_spine_hit_list.end()) continue;


	    /*
	    spacepoint_vec.clear();
	    spacepoint_vec = spacepoint_per_hit.at(collected_hits[i_h].key());
	    TVector3 hit_position = spacepoint_vec[0]->XYZ();
	    */
	    //TVector3 hit_position = fHitsToSpacePoints.at(collected_hits[i_h])->XYZ();
	    TVector3 hit_position = fHitsToSpacePoints_old.at(collected_hits[i_h])->XYZ();
	    pandora::CartesianVector pandora_hit_position(hit_position.X(), hit_position.Y(), hit_position.Z());

	    cout << "closest distance: " << this->getClosestDistance(hit_position, running_fit_position_vec) << " " << hit_connection_distance << endl;

	    if(this->getClosestDistance(hit_position, running_fit_position_vec) > hit_connection_distance){
	      continue;
	    }

	    found = true;
	    running_fit_position_vec.push_back(hit_position);
	    pandora_running_fit_position_vec.push_back(pandora_hit_position);
	    //cout << "calling shower_spine_hit_list.push_back" << endl;
	    shower_spine_hit_list.push_back(collected_hits[i_h]);
	  }
	}

    }


    //------------------------------------------------------------------------------------------------------------------------------------------              


    //spacepoints

  //------------------------------------------------------------------------------------------------------------------------------------------                                                                                                                                                                                

    //------------------------------------------------------------------------------------------------------------------------------------------                                                                                                                                                                                


    bool HitSplitAlg::isCloseToLine(const TVector3 &hit_position,
				       const TVector3 &line_start,
				       const TVector3 &line_direction,
				       const double distance_to_line)
    {
      const double transverse_distance_from_line = line_direction.Cross(hit_position - line_start).Mag();
      //cout << "transverse_distance_from_line is " << transverse_distance_from_line << ", and distance_to_line is " << distance_to_line << endl;
      if(transverse_distance_from_line > distance_to_line) return false;
      return true;
    }


    //------------------------------------------------------------------------------------------------------------------------------------------                                                                                                                                                                                

    double HitSplitAlg::getClosestDistance(const TVector3 &position,
					      vector<TVector3> &test_positions)
    {
      double closest_distance_sqaured(std::numeric_limits<double>::max());

      for(TVector3 test_position : test_positions){
	const double separation_squared = (test_position - position).Mag2();
	if (separation_squared < closest_distance_sqaured) closest_distance_sqaured = separation_squared;
      }

      return std::sqrt(closest_distance_sqaured);
    }


    //------------------------------------------------------------------------------------------------------------------------------------------                                                                                                                                                                   
    void HitSplitAlg::drawHistAngularDistributionMap3DCheat( std::map<int, TH2D*> &h_angular_distribution_pfparticle_cheat_3D,
						TString outfile_name,
						TCanvas* &c){
      
      TFile outfile(outfile_name, "update");
      auto hs = new THStack("hs", "");
      TLegend * leg = new TLegend(0.7, 0.7, 0.9, 0.9, "");
      c->SetFillStyle(1001);
      
      if(h_angular_distribution_pfparticle_cheat_3D.size()==0) return;
      
      std::map<int, TH2D*>::iterator it;
      //cout << "h_angular_distribution_pfparticle_cheat_3D.size(): " << h_angular_distribution_pfparticle_cheat_3D.size() << endl;

      for(it=h_angular_distribution_pfparticle_cheat_3D.begin(); it!=h_angular_distribution_pfparticle_cheat_3D.end(); it++){

	it->second->SetFillStyle(1001);
	//cout << "it->first: " << it->first << ", it->second->GetEntries(): " << it->second->GetEntries() << endl;

	if(it->first == 321){
	  it->second->SetFillColor(kBlue);
	  //leg->AddEntry(it->second, "K+");
	}
	else if(it->first == -13){
	  it->second->SetFillColor(kCyan);
	  //leg->AddEntry(it->second, "mu+");
	}
	else if(it->first == 211){
	  it->second->SetFillColor(kMagenta);
	  //leg->AddEntry(it->second, "pi+");
	}
	else if(it->first == -11 || it->first == 11){
	  it->second->SetFillColor(kGreen+2);
	  //leg->AddEntry(it->second, "shower");
	}
	else if(it->first == 2212 || it->first == 2112){
	  it->second->SetFillColor(kRed);
	  //leg->AddEntry(it->second, "nucleon");
	}
	else{
	  it->second->SetFillColor(kBlack);
	}
	hs->Add(it->second);

      }

      /*
      for(auto const& h_org: h_angular_distribution_pfparticle_cheat_3D){

	TH2D h = h_org;	
	h.second.SetFillStyle(1001);
	
	if(h.first == 321) h.second.SetFillColor(kBlue);
	else if(h.first == -13) h.second.SetFillColor(kCyan);
	else if(h.first == 211) h.second.SetFillColor(kMagenta);
	else if(h.first == -11 || h.first == 11) h.second.SetFillColor(kGreen+2);
	else if(h.first == 2212 || 2112) h.second.SetFillColor(kRed);
	else h.second.SetFillColor(kBlack);
	hs->Add(h.second);
	
      }
      */
      hs->Draw("hist""lego3 0");
      //hs->Draw("hist nostack");                                                                                                                                                                                                                                                                                      
      leg->Draw();
      c->Write();
      
      delete hs;
      delete leg;
      
      return;
      
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------   
   
#endif
