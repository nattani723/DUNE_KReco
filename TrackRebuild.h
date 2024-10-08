#ifndef TRACK_REBUILDER
#define TRACK_REBUILDER

//#include "HitSplitAlg_module.h"


#include "art/Framework/Core/EDProducer.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "art/Persistency/Common/PtrMaker.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArObjects/LArPfoObjects.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "larpandora/LArPandoraInterface/LArPandoraOutput.h" 

#include <iostream>

using namespace pandora;
using namespace std;
namespace HitSplitAlg_module
{
  

  recob::Track HitSplitAlg::buildTrack(int id,
				       lar_content::LArTrackStateVector& trackStateVector)
    {

      if (trackStateVector.empty()) std::cout << "BuildTrack - No input trajectory points provided" << std::endl;

      recob::tracking::Positions_t xyz;
      recob::tracking::Momenta_t pxpypz;
      recob::TrackTrajectory::Flags_t flags;

      for (const lar_content::LArTrackState& trackState : trackStateVector) {

        xyz.emplace_back(recob::tracking::Point_t(trackState.GetPosition().GetX(),
                                                  trackState.GetPosition().GetY(),
                                                  trackState.GetPosition().GetZ()));
        pxpypz.emplace_back(recob::tracking::Vector_t(trackState.GetDirection().GetX(),
                                                      trackState.GetDirection().GetY(),
                                                      trackState.GetDirection().GetZ()));

        // Set flag NoPoint if point has bogus coordinates, otherwise use clean flag set                                   
  
        if (std::fabs(trackState.GetPosition().GetX() - util::kBogusF) <
            std::numeric_limits<float>::epsilon() &&
            std::fabs(trackState.GetPosition().GetY() - util::kBogusF) <
            std::numeric_limits<float>::epsilon() &&
            std::fabs(trackState.GetPosition().GetZ() - util::kBogusF) <
            std::numeric_limits<float>::epsilon()) {
          flags.emplace_back(recob::TrajectoryPointFlags(recob::TrajectoryPointFlags::InvalidHitIndex,
                                                         recob::TrajectoryPointFlagTraits::NoPoint));
        }
        else {
          flags.emplace_back(recob::TrajectoryPointFlags());
        }
      }

      // note from gc: eventually we should produce a TrackTrajectory, not a Track with empty covariance matrix and bogus chi2, etc.                    
  
      return recob::Track(
                          recob::TrackTrajectory(std::move(xyz), std::move(pxpypz), std::move(flags), false),
                          util::kBogusI,
                          util::kBogusF,
                          util::kBogusI,
                          recob::tracking::SMatrixSym55(),
                          recob::tracking::SMatrixSym55(),
                          id);

    }
  



  recob::Track HitSplitAlg::trackRebuid(std::vector<art::Ptr<recob::Hit>>& hit_list,
					const recob::Track& track)
    {
      
      //cout << "aaa" << endl;

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


      int trackCounter=1000;
      unsigned  int m_slidingFitHalfWindow = 20;

      pandora::CartesianPointVector pandora_hit_positions;

      cout << "hit_list.size(): " << hit_list.size() << endl;

      for(size_t i_h=0; i_h<hit_list.size(); i_h++){

	//art::Ptr<recob::SpacePoint> sp = fHitsToSpacePoints.at(hit_list[i_h]);
	art::Ptr<recob::SpacePoint> sp = fHitsToSpacePoints_old.at(hit_list[i_h]);
	pandora::CartesianVector pandora_hit_position(sp->XYZ()[0], sp->XYZ()[1], sp->XYZ()[2]);
	//cout << "sp->XYZ()[0]: " << sp->XYZ()[0] << endl;
	pandora_hit_positions.emplace_back(pandora_hit_position);
	//}
      }

      std::unique_ptr<std::vector<recob::Track>> outputTracks(new std::vector<recob::Track>);

      const pandora::CartesianVector vertexPosition(track.End().X(), track.End().Y(),track.End().Z());

      lar_content::LArTrackStateVector trackStateVector;
      pandora::IntVector indexVector;

      lar_content::LArPfoHelper::GetSlidingFitTrajectory(pandora_hit_positions,
							 vertexPosition,
							 m_slidingFitHalfWindow,
							 wirePitchW,
							 trackStateVector,
							 &indexVector);

      cout << "trackStateVector.size(): " << trackStateVector.size() << endl;
      cout << "hit_list.size(): " << hit_list.size() << endl;

      outputTracks->emplace_back(buildTrack(trackCounter++, trackStateVector));
      cout << "Rebuilt track length: " <<  outputTracks->at(0).Length() << endl;

      recob::Track reco_track = outputTracks->at(0);
      return reco_track;

    }
  //}



  recob::Track HitSplitAlg::trackRebuid(std::vector<art::Ptr<recob::Hit>>& hit_list,
					art::FindManyP<recob::SpacePoint>& spacepoint_per_hit,
					const recob::Track& track)
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

      int trackCounter=1000;
      unsigned  int m_slidingFitHalfWindow = 20;

      std::vector<art::Ptr<recob::SpacePoint>> spacepoint_vec;
      pandora::CartesianPointVector pandora_hit_positions;

      for(size_t i_h=0; i_h<hit_list.size(); i_h++){

	spacepoint_vec.clear();
	spacepoint_vec = spacepoint_per_hit.at(hit_list[i_h].key());
	if(spacepoint_vec.size()!=1) continue;

	for(auto spacepoint : spacepoint_vec){
	  pandora::CartesianVector pandora_hit_position(spacepoint->XYZ()[0], spacepoint->XYZ()[1], spacepoint->XYZ()[2]);
	  pandora_hit_positions.emplace_back(pandora_hit_position);
	}

      }

      std::unique_ptr<std::vector<recob::Track>> outputTracks(new std::vector<recob::Track>);

      const pandora::CartesianVector vertexPosition(track.End().X(), track.End().Y(),track.End().Z());

      lar_content::LArTrackStateVector trackStateVector;
      pandora::IntVector indexVector;

      lar_content::LArPfoHelper::GetSlidingFitTrajectory(pandora_hit_positions,
							 vertexPosition,
							 m_slidingFitHalfWindow,
							 wirePitchW,
							 trackStateVector,
							 &indexVector);

      cout << "trackStateVector.size(): " << trackStateVector.size() << endl;
      cout << "hit_list.size(): " << hit_list.size() << endl;

      outputTracks->emplace_back(buildTrack(trackCounter++, trackStateVector));
      cout << "Rebuilt track length: " <<  outputTracks->at(0).Length() << endl;

      recob::Track reco_track = outputTracks->at(0);
      return reco_track;

    }
}

#endif
