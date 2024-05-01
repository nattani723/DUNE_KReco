#ifndef HitSplitAlg_Module
#define HitSplitAlg_Module


#include "art/Framework/Core/EDProducer.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "art/Persistency/Common/PtrMaker.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArObjects/LArPfoObjects.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "larpandora/LArPandoraInterface/LArPandoraOutput.h" 

// LArSoft includes
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardata/ArtDataHelper/MVAReader.h"
// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"

// ROOT includes 
#include "TDirectory.h" 

#include "art_root_io/TFileService.h"
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/View.h"
#include "art/Utilities/make_tool.h"

#include "canvas/Utilities/InputTag.h"
#include "canvas/Utilities/ensurePointer.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Persistency/Common/TriggerResults.h"
#include "canvas/Utilities/InputTag.h"

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/ParameterSetRegistry.h"

#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"

#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "larreco/RecoAlg/TrajectoryMCSFitter.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"

#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h"  
#include "larcoreobj/SummaryData/POTSummary.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

#include "larcorealg/Geometry/GeometryCore.h" 

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/FlashMatch.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/MCBase/MCStep.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/Simulation/GeneratedParticleInfo.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/BeamInfo.h"
#include "lardataobj/RawData/TriggerData.h"
#include "larcoreobj/SummaryData/POTSummary.h"



#include "larevt/SpaceCharge/SpaceCharge.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

#include "larsim/EventWeight/Base/MCEventWeight.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/simb.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/GTruth.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "usimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include "Pandora/PdgTable.h"

#include "TCanvas.h"
#include "TTree.h"
#include "TMath.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphDelaunay.h"
#include "TRandom3.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TH1.h"
#include "TFile.h"
#include "TTimeStamp.h"

#include <iostream>
#include <TH2D.h>
#include <TH3D.h>
#include <TH1D.h>
#include <TFile.h>
#include <TROOT.h>
#include <TChain.h>
#include <TTree.h>
#include <TVectorT.h>
#include <THStack.h>
#include <TPDF.h>
#include <TLegend.h>
#include <vector>
#include <iomanip>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <string.h>
#include <array>
#include <map>
#include <cstddef>
#include <cstring>
#include <iterator>
#include <string>
#include <numeric>
#include <algorithm>
#include <functional>
#include <typeinfo>
#include <memory>
#include <chrono>
#include <sys/stat.h>


#define MAX_TRACKS 40
#define MAX_GEN_PARTS 100
#define MAX_MC_PARTS 2500
#define MAX_MC_TRAJ_PTS 100
#define MAX_FLASHES 200
#define MAX_SHOWERS 100
#define MVA_LENGTH 4
#define MAX_HITS 1000
#define MAX_CALO_PTS 500

using namespace pandora;
using namespace std;

namespace HitSplitAlg_module{
class HitSplitAlg : public art::EDAnalyzer {
public:
    
  explicit HitSplitAlg(fhicl::ParameterSet const& pset);
  virtual ~HitSplitAlg();

  void beginJob();
  void endJob();
  void beginRun(const art::Run& run);
  void analyze(const art::Event& evt);

  void reconfigure(fhicl::ParameterSet const& pset);
  void PIDAcal( std::vector<const anab::Calorimetry*> cal, std::vector<double> &PIDA );
  void Process(const art::Event& evt, bool &isFiducial);
  void truthHitMatcher( art::Ptr<recob::Hit> hit, const simb::MCParticle *&MCparticle);
  void truthMatcher( std::vector<art::Ptr<recob::Hit>> all_hits, std::vector<art::Ptr<recob::Hit>> track_hits, const simb::MCParticle *&MCparticle, double &Efrac, double &Ecomplet);
  double truthLength( const simb::MCParticle *MCparticle );
  bool insideFV(double vertex[4]);

    
  void fillAngularDistributionMapCheat(std::vector<art::Ptr<recob::Hit>>& hits_from_pfparticle,
				       TVector3 Kend_candidate,
				       art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData>& particles_per_hit,
				       art::FindManyP<recob::SpacePoint>& spacepoint_per_hit,
				       std::map<int, std::map<int, double>>& angular_distribution_mrgid_map,
				       std::map<int,int> &mrgidpdg,
				       int& currentMergedId);

  void fillAngularDistributionMapCheat3D(std::vector<art::Ptr<recob::Hit>>& hits_from_pfparticle,
					 std::vector<art::Ptr<recob::Hit>>& true_pi0_hit_list,
					 TVector3 Kend_candidate,
					 art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData>& particles_per_hit,
					 art::FindManyP<recob::SpacePoint>& spacepoint_per_hit,
					 std::map<int, std::map<int, std::map<int, double>>>& angular_distribution_mrgid_map_3D,
					 std::map<int,int> &mrgidpdg,
					 std::map<int,TVector3> &mrgidmom,
					 int& currentMergedId);

  void fillAngularDistributionMap(std::vector<art::Ptr<recob::Hit>>& hits_from_pfparticle,
				  TVector3 Kend_candidate,
				  art::FindManyP<recob::SpacePoint>& spacepoint_per_hit,
				  std::map<int, double>& angular_distribution_map);


  void fillAngularDistributionMap3D(std::vector<art::Ptr<recob::Hit>>& hits_from_pfparticle,
				    TVector3 Kend_candidate,
				    art::FindManyP<recob::SpacePoint>& spacepoint_per_hit,
				    std::map<int, std::map<int, double>> &angular_distribution_map_3D);

  void fillAngularDistributionMap3D(std::vector<art::Ptr<recob::SpacePoint>>& sp_from_recoobj,
				    TVector3 Kend_candidate,
				    std::map<int, std::map<int, double>> &angular_distribution_map_3D);

  int fillHistAngularDistributionMapCheat( std::map<int, std::map<int, double>>& angular_distribution_mrgid_map,
					   std::map<int,int>& mrgidpdg,
					   std::vector<TH1D*>& h_angular_distribution_pfparticle_cheat,
					   vector<int>& v_pdg);


  int fillHistAngularDistributionMapCheat3D( std::map<int, std::map<int, std::map<int, double>>>& angular_distribution_mrgid_map_3D,
					     std::map<int,int>& mrgidpdg,
					     std::vector<TH2D*>& h_angular_distribution_pfparticle_cheat_3D,
					     vector<int>& v_pdg);

  void fillHistAngularDistributionMap3DCheat( std::vector<art::Ptr<recob::SpacePoint>>& sp_from_recoobj,
					      TVector3 Kend_candidate,
					      std::map<int, std::map<int, std::map<int, double>>> &angular_distribution_map_3D_cheat,
					      std::map<int, TH2D*> &h_angular_distribution_pfparticle_cheat_3D);

  void drawHistAngularDistributionMap3DCheat( std::map<int, TH2D*> &h_angular_distribution_pfparticle_cheat_3D,
					      TString outfile_name,
					      TCanvas* &c);

  int fillHistAngularDistributionMap( std::map<int, double>& angular_distribution_map,
				      vector<TH1D*>& h_angular_distribution_pfparticle,
				      vector<bool>& v_trk_flg,
				      bool trk_flg);

  int fillHistAngularDistributionMap3D( std::map<int, std::map<int, double>>& angular_distribution_map_3D,
					vector<TH2D*>& h_angular_distribution_pfparticle_3D,
					vector<bool>& v_trk_flg,
					bool trk_flg);

  int fillHistAngularDistributionMapSurface( std::map<int, std::map<int, double>>& angular_distribution_map_3D,
					     vector<TH2D*>& h_angular_distribution_pfparticle_surface,
					     vector<bool>& v_trk_flg,
					     bool trk_flg);

  int fillHistAngularDistributionMapSphere( std::map<int, std::map<int, double>>& angular_distribution_map_3D,
					    vector<TH3D*>& h_angular_distribution_pfparticle_sphere,
					    vector<bool>& v_trk_flg,
					    bool trk_flg);

  void smoothAngularDistributionMap(std::map<int, double> &angular_distribution_map);

  void smoothAngularDistributionMapCheat3D(std::map<int, std::map<int, std::map<int, double>>> &angular_distribution_mrgid_map_3D);

  void smoothAngularDistributionMap3D(std::map<int, std::map<int, double>> &angular_distribution_map_3D);

  void accumulateAngularDistributionMap3D(std::map<int, std::map<int, double>> &angular_distribution_map_3D_1,
					  std::map<int, std::map<int, double>> &angular_distribution_map_3D_2,
					  std::map<int, std::map<int, double>> &angular_distribution_map_3D);

  void obtainPeakVectorCheat3D(std::map<int, std::map<int, std::map<int, double>>> &angular_distribution_mrgid_map_3D,
			       std::map<int,int>& mrgidpdg,
			       std::map<int, std::map<double, TVector2, std::greater<>>>& view_peak_map_cheat,
			       vector<int>& v_pdg_peak);

  void obtainPeakVector3D(std::map<int, std::map<int, double>> &angular_distribution_map_3D,
			  vector<bool>& v_trk_flg,
			  std::map<double, TVector2, std::greater<>>& view_peak_map,
			  bool trk_flg);

  void findBestAngularPeakCheat3D( std::map<int, std::map<double, TVector2, std::greater<>>>& view_peak_map_cheat,
				   std::map<int, vector<TVector2>> &best_peak_bins_cheat);


  void findBestAngularPeak3D(std::map<int, std::map<int, double>> &angular_distribution_map_3D,
			     std::map<double, TVector2, std::greater<>>& view_peak_vector,
			     vector<TVector2> &best_peak_bins);

  void findShowerSpine3D(std::vector<art::Ptr<recob::SpacePoint>>& sp_from_pfparticle,
			 std::vector<art::Ptr<recob::Hit>>& unavailable_hit_list,
			 std::vector<std::vector<art::Ptr<recob::Hit>>>& shower_spine_hit_list_vector,
			 TVector3 Kend_candidate,
			 std::map<double, TVector2, std::greater<>>& view_peak_map,
			 vector<TVector2> &best_peak_bins);


  void findShowerSpine3D(std::vector<art::Ptr<recob::Hit>>& hits_from_pfparticle,
			 art::FindManyP<recob::SpacePoint>& spacepoint_per_hit,
			 std::vector<art::Ptr<recob::Hit>>& unavailable_hit_list,
			 std::vector<std::vector<art::Ptr<recob::Hit>>>& shower_spine_hit_list_vector,
			 TVector3 Kend_candidate,
			 std::map<double, TVector2, std::greater<>>& view_peak_map,
			 vector<TVector2> &best_peak_bins);


  void findShowerSpine3D(std::vector<art::Ptr<recob::Hit>>& hits_from_pfparticle,
			 art::FindManyP<recob::SpacePoint>& spacepoint_per_hit,
			 std::vector<art::Ptr<recob::Hit>>& unavailable_hit_list,
			 std::vector<art::Ptr<recob::Hit>>& shower_spine_hit_list,
			 TVector3 Kend_candidate,
			 std::map<double, TVector2, std::greater<>>& view_peak_map,
			 vector<TVector2> &best_peak_bins);

  void findShowerSpine3D(std::vector<art::Ptr<recob::Hit>>& hits_from_pfparticle,
			 art::FindManyP<recob::SpacePoint>& spacepoint_per_hit,
			 std::vector<art::Ptr<recob::Hit>>& unavailable_hit_list,
			 std::vector<art::Ptr<recob::Hit>>& shower_spine_hit_list,
			 TVector3 Kend_candidate,
			 std::map<double, TVector2, std::greater<>>& view_peak_map,
			 vector<TVector2> &best_peak_bins,
			 double true_length);
 
  bool collectSubsectionHits(const lar_content::ThreeDSlidingFitResult &extrapolated_fit,
			     const TVector3 &extrapolated_start_position,
			     const TVector3 &extrapolated_end_position,
			     const TVector3 &extrapolated_direction,
			     const bool is_end_downstream,
			     std::vector<art::Ptr<recob::SpacePoint>>& sp_from_pfparticle,
			     vector<TVector3> &running_fit_position_vec,
			     pandora::CartesianPointVector &pandora_running_fit_position_vec,
			     std::vector<art::Ptr<recob::Hit>>& unavailable_hit_list,
			     std::vector<art::Ptr<recob::Hit>>& shower_spine_hit_list);


  bool collectSubsectionHits(const lar_content::ThreeDSlidingFitResult &extrapolated_fit,
			     const TVector3 &extrapolated_start_position,
			     const TVector3 &extrapolated_end_position,
			     const TVector3 &extrapolated_direction,
			     const bool is_end_downstream,
			     std::vector<art::Ptr<recob::Hit>>& hits_from_pfparticle,
			     art::FindManyP<recob::SpacePoint>& spacepoint_per_hit,
			     vector<TVector3> &running_fit_position_vec,
			     pandora::CartesianPointVector &pandora_running_fit_position_vec,
			     std::vector<art::Ptr<recob::Hit>>& unavailable_hit_list,
			     std::vector<art::Ptr<recob::Hit>>& shower_spine_hit_list);


  bool collectSubsectionHits(const lar_content::ThreeDSlidingFitResult &extrapolated_fit,
			     const TVector3 &extrapolated_start_position,
			     const TVector3 &extrapolated_end_position,
			     const TVector3 &extrapolated_direction,
			     const bool is_end_downstream,
			     std::vector<art::Ptr<recob::Hit>>& hits_from_pfparticle,
			     art::FindManyP<recob::SpacePoint>& spacepoint_per_hit,
			     vector<TVector3> &running_fit_position_vec,
			     pandora::CartesianPointVector &pandora_running_fit_position_vec,
			     std::vector<art::Ptr<recob::Hit>>& unavailable_hit_list,
			     std::vector<art::Ptr<recob::Hit>>& shower_spine_hit_list,
			     TVector3 Kend_candidate,
			     double true_length);

  void collectConnectedHits(std::vector<art::Ptr<recob::Hit>>& collected_hits,
			    const TVector3 &extrapolated_start_position,
			    const TVector3 &extrapolated_direction,
			    vector<TVector3> &running_fit_position_vec,
			    pandora::CartesianPointVector &pandora_running_fit_position_vec,
			    std::vector<art::Ptr<recob::Hit>>& shower_spine_hit_list);

  void collectConnectedHits(std::vector<art::Ptr<recob::Hit>>& collected_hits,
			    art::FindManyP<recob::SpacePoint>& spacepoint_per_hit,
			    const TVector3 &extrapolated_start_position,
			    const TVector3 &extrapolated_direction,
			    vector<TVector3> &running_fit_position_vec,
			    pandora::CartesianPointVector &pandora_running_fit_position_vec,
			    std::vector<art::Ptr<recob::Hit>>& shower_spine_hit_list);

  TVector3 getClosestPointToLine3D(const TVector3 &extrapolated_start_position,
				   const TVector3 &extrapolated_direction,
				   art::Ptr<recob::Hit>& collected_hit,
				   const TVector3 &hit_position);

  bool isCloseToLine(const TVector3 &hit_position,
		     const TVector3 &line_start,
		     const TVector3 &line_direction,
		     const double distance_to_line);

  bool isInsideROI( art::Ptr<recob::SpacePoint> &sp,
		    TVector3 Kend_candidate);

  double getClosestDistance(art::Ptr<recob::Hit>& collected_hit,
			    art::Ptr<recob::Hit>& shower_spine_hit,
			    art::FindManyP<recob::SpacePoint>& spacepoint_per_hit);

  double getClosestDistance(art::Ptr<recob::Hit>& collected_hit,
			    std::vector<art::Ptr<recob::Hit>>& shower_spine_hit_list,
			    art::FindManyP<recob::SpacePoint>& spacepoint_per_hit);

  double getClosestDistance(const TVector3 &position,
			    vector<TVector3> &test_positions);

  void obtainLength(std::vector<art::Ptr<recob::Hit>>& shower_spine_hit_list,
		    art::FindManyP<recob::SpacePoint>& spacepoint_per_hit,
		    TVector3 Kend_candidate);

  void obtainLongitudinalDecomposition(std::vector<art::Ptr<recob::Hit>>& shower_spine_hit_list,
				       art::FindManyP<recob::SpacePoint>& spacepoint_per_hit);

  int drawHistAngularDistributionMapCheat( std::vector<TH1D*> &h_angular_distribution_pfparticle_cheat,
					   std::vector<int> &v_pdg,
					   TString outfile_name,
					   TCanvas*& c);

  int drawHistAngularDistributionMapCheat3D( std::vector<TH2D*> &h_angular_distribution_pfparticle_cheat_3D,
					     std::vector<int> &v_pdg,
					     TString outfile_name,
					     TCanvas*& c);

  int drawHistAngularDistributionMap( std::vector<TH1D*>& h_angular_distribution_pfparticle,
				      std::vector<bool>& v_trk_flg,
				      TString outfile_name,
				      TCanvas*& c);

  int drawHistAngularDistributionMap3D( std::vector<TH2D*>& h_angular_distribution_pfparticle_3D,
					std::vector<bool>& v_trk_flg,
					TString outfile_name,
					TCanvas*& c);

  int drawHistAngularDistributionMapSurface( std::vector<TH2D*> &h_angular_distribution_pfparticle_surface,
					     std::vector<bool> &v_trk_flg,
					     TString outfile_name,
					     TCanvas* &c);

  int drawHistAngularDistributionMapSphere( std::vector<TH3D*> &h_angular_distribution_pfparticle_sphere,
					    std::vector<bool> &v_trk_flg,
					    TString outfile_name,
					    TCanvas* &c);

  void fillTrueMatching(std::vector<art::Ptr<recob::Hit>>& hits_from_track,
			art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData>& particles_per_hit,
			int track_i=-1,
			int daughter_i=-1);
  void fillTrueMatching_sh(std::vector<art::Ptr<recob::Hit>>& hits_from_shower,
			   art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData>& particles_per_hit,
			   int track_i=-1,
			   int daughter_i=-1);
  void fillTrackMatching(std::vector<art::Ptr<recob::Hit>>& hits_from_track,
			 art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData>& particles_per_hit,
			 int track_i=-1,
			 int daughter_i=-1);

  void fillShowerMatching(std::vector<art::Ptr<recob::Hit>>& hits_from_shower,
			  art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData>& particles_per_hit,
			  int track_i=-1,
			  int daughter_i=-1);

  recob::Track trackRebuid(std::vector<art::Ptr<recob::Hit>>& hit_list,
			   art::FindManyP<recob::SpacePoint>& spacepoint_per_hit,
			   const recob::Track& track);

  recob::Track trackRebuid(std::vector<art::Ptr<recob::Hit>>& hit_list,
			   const recob::Track& track);

  recob::Track buildTrack(int id,
                          lar_content::LArTrackStateVector& trackStateVector);
  



 private:
  void reset();

 private:

  // the parameters we'll read from the .fcl
  std::string fHitModuleLabel; //uB: fHitsModuleLabel
  std::string fHitModuleLabelOLD;
  std::string fHitSPAssns;
  std::string fTrackModuleLabel;
  std::string fPidModuleLabel; //uB: fPIDLabel
  std::string fCaloModuleLabel; //uB: fCalorimetryModuleLabel
  std::string fOpFlashModuleLabel;
  std::string fShowerModuleLabel;
  std::string fPointIDModuleLabel;
  std::string fNNetModuleLabel;
  std::string fPFParticleModuleLabel; //uB: fPFParticleLabel
  std::string fSpacePointproducer;
  std::string fMCgenieLabel;
  std::string fHitTrackAssns;
  //std::string fTracksToSpacePoints;
  double      fExponentConstant;
  double fPIDA_endPoint;
  doublefMaxPIDAValue;
  doublefMinPIDAValue;
  boolfSaveMCTree;
  TCanvas * c = new TCanvas("c", "c", 800, 600);

  std::vector<art::Ptr<recob::Hit>> fHits;
  std::vector<art::Ptr<recob::Track>> fTracks;
  std::vector<art::Ptr<recob::Shower>> fShowers;
  std::vector<art::Ptr<recob::SpacePoint>> fSpacePoints;
  std::map<art::Ptr<recob::Track>, std::vector<art::Ptr<recob::Hit>>> fTracksToHits; 
  std::map<art::Ptr<recob::Track>, std::vector<art::Ptr<recob::SpacePoint>>> fTracksToSpacePoints;
  std::map<art::Ptr<recob::Shower>, std::vector<art::Ptr<recob::Hit>>> fShowersToHits;
  std::map<art::Ptr<recob::Shower>, std::vector<art::Ptr<recob::SpacePoint>>> fShowersToSpacePoints;
  std::map<art::Ptr<recob::Hit>, art::Ptr<recob::SpacePoint>> fHitsToSpacePoints;
  std::map<art::Ptr<recob::Hit>, art::Ptr<recob::SpacePoint>> fHitsToSpacePoints_old;
  std::map<art::Ptr<recob::SpacePoint>, art::Ptr<recob::Hit>> fSpacePointsToHits;
  std::map<art::Ptr<recob::SpacePoint>, art::Ptr<recob::Hit>> fSpacePointsToHits_old;



  TTree *fEventTree;

  // Event
  int Event;
  int Run;
  int SubRun;

  //MC truth
  double MC_Ev;
  int    MC_nuPDG;
  double MC_Q2;
  double MC_hit_nucleon;
  int    MC_cc;
  int    MCgenie_npart;
  int    MCgenie_id[MAX_GEN_PARTS];
  int    MCgenie_pdg[MAX_GEN_PARTS];
  int    MCgenie_mother[MAX_GEN_PARTS];
  int    MCgenie_statusCode[MAX_GEN_PARTS];
  int    MCgenie_fate[MAX_GEN_PARTS];
  double MCgenie_startMomentum[MAX_GEN_PARTS][4];
  double MCgenie_endMomentum[MAX_GEN_PARTS][4];

  double MC_vertex[4];
  int    MC_npart;
  int    MC_id[MAX_MC_PARTS];
  int    MC_pdg[MAX_MC_PARTS];
  int    MC_mother[MAX_MC_PARTS];
  double MC_startXYZT[MAX_MC_PARTS][4];
  double MC_endXYZT[MAX_MC_PARTS][4];
  int    MC_npoints[MAX_MC_PARTS];
  double MC_xyz[MAX_MC_PARTS][MAX_MC_TRAJ_PTS][3];
  double MC_startMomentum[MAX_MC_PARTS][4];
  double MC_endMomentum[MAX_MC_PARTS][4];
  double MC_truthlength[MAX_MC_PARTS];
  double MC_Prange[MAX_MC_PARTS];
  int    MC_statusCode[MAX_MC_PARTS];
  double MC_kaontruthlength;
  double MC_dautruthlength;
  int Kaon_BR;
  std::vector<std::string> MC_process;
  std::vector<std::string> MC_Endprocess;

  int    n_vertices;
  double vertex[MAX_TRACKS][4];
  int    vtx_ID[MAX_TRACKS];
  int    n_decayVtx;
  double decayVtx[MAX_TRACKS][3];
  int    n_recoTracks;
  int    n_recoDauTracks[MAX_TRACKS];
  int    n_recoRebDauTracks[MAX_TRACKS];
  int    track_isContained[MAX_TRACKS];
  int    track_ID[MAX_TRACKS];
  int    vtxID_trk[MAX_TRACKS][10];
  double track_vtxDir[MAX_TRACKS][3];
  double track_vtx[MAX_TRACKS][4];
  double track_end[MAX_TRACKS][4];
  double track_length[MAX_TRACKS];
  double dautrack_length[MAX_TRACKS][10];
  double rebdautracktrue_length[MAX_TRACKS];
  double rebdautracktruedir_length[MAX_TRACKS];
  double rebdautrack_length[MAX_TRACKS][10];
  double dautrack_pdg[MAX_TRACKS][10];
  double rebdautrack_pdg[MAX_TRACKS][10];
  double best_peak_theta[MAX_TRACKS][10];
  double best_peak_phi[MAX_TRACKS][10];
  double best_peak_theta_true[MAX_TRACKS];
  double best_peak_phi_true[MAX_TRACKS];

  double track_dir_vtx[MAX_TRACKS][4];
  double track_PIDA[MAX_TRACKS][3];
  int   track_PID_pdg[MAX_TRACKS][3];
  double track_KE[MAX_TRACKS][3];
  int    track_bestplane[MAX_TRACKS];
  double track_Prange[MAX_TRACKS];
  double track_Efrac[MAX_TRACKS];
  double track_complet[MAX_TRACKS];
  int    track_mcID[MAX_TRACKS];
  int    track_mcPDG[MAX_TRACKS];
  int    dautrack_mcPDG[MAX_TRACKS][10];
  int    n_track_points[MAX_TRACKS];
  double track_point_xyz[MAX_TRACKS][MAX_CALO_PTS][3];
  int    n_cal_points[MAX_TRACKS];
  double track_dQ_dx[MAX_TRACKS][MAX_CALO_PTS];
  double track_dE_dx[MAX_TRACKS][MAX_CALO_PTS];
  double track_range[MAX_TRACKS][MAX_CALO_PTS];
  double track_pitch[MAX_TRACKS][MAX_CALO_PTS];
  int    n_cal_points_byplane[MAX_TRACKS][3];
  double track_dQ_dx_byplane[MAX_TRACKS][3][MAX_CALO_PTS];
  double track_dE_dx_byplane[MAX_TRACKS][3][MAX_CALO_PTS];
  double track_range_byplane[MAX_TRACKS][3][MAX_CALO_PTS];
  double track_pitch_byplane[MAX_TRACKS][3][MAX_CALO_PTS];
  double track_calo_xyz_byplane[MAX_TRACKS][3][MAX_CALO_PTS][3]; // storing xyz coordinates of all calo points

  int     n_recoHits;
  int     hit_channel[MAX_HITS];
  int     hit_tpc[MAX_HITS];
  int     hit_plane[MAX_HITS];
  int     hit_wire[MAX_HITS];
  double  hit_peakT[MAX_HITS];
  double  hit_charge[MAX_HITS];
  double  hit_ph[MAX_HITS];
  double  hit_startT[MAX_HITS];
  double  hit_endT[MAX_HITS];
  double  hit_rms[MAX_HITS];
  double  hit_electrons[MAX_HITS];


  double Em_ch;
  double Em_e;
  double trk_e;
  double Emichel_e;

  int    n_recoShowers;
  double sh_direction_X[MAX_SHOWERS];
  double sh_direction_Y[MAX_SHOWERS];
  double sh_direction_Z[MAX_SHOWERS];
  double sh_start_X[MAX_SHOWERS];
  double sh_start_Y[MAX_SHOWERS];
  double sh_start_Z[MAX_SHOWERS];
  double sh_energy[MAX_SHOWERS][3];
  double sh_MIPenergy[MAX_SHOWERS][3];
  double sh_dEdx[MAX_SHOWERS][3];
  int    sh_bestplane[MAX_SHOWERS];
  double sh_length[MAX_SHOWERS];


  int    n_flashes;
  double flash_time[MAX_FLASHES];
  double flash_pe[MAX_FLASHES];
  double flash_ycenter[MAX_FLASHES];
  double flash_zcenter[MAX_FLASHES];
  double flash_ywidth[MAX_FLASHES];
  double flash_zwidth[MAX_FLASHES];
  double flash_timewidth[MAX_FLASHES];
  double flash_abstime[MAX_FLASHES];
  int    flash_frame[MAX_FLASHES];
  double flash_PE_ndk[MAX_FLASHES];
  double flash_PE_Ar39[MAX_FLASHES];

  double fFidVolCutX;
  double fFidVolCutY;
  double fFidVolCutZ;

  double fFidVolXmin;
  double fFidVolXmax;
  double fFidVolYmin;
  double fFidVolYmax;
  double fFidVolZmin;
  double fFidVolZmax;

  double fPidValue;
  unsigned int    fView;

  bool fSaveHits;

  detinfo::DetectorClocksData const clockdata = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
  detinfo::DetectorPropertiesData const detprop = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob(clockdata);
  //auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);

  double XDriftVelocity      = detprop.DriftVelocity()*1e-3; //cm/ns
  double WindowSize          = detprop.NumberTimeSamples() * clockdata.TPCClock().TickPeriod() * 1e3;

  /*
    detinfo::DetectorProperties const *detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    detinfo::DetectorClocks const *ts = lar::providerFrom<detinfo::DetectorClocksService>();
    double XDriftVelocity = detprop->DriftVelocity()*1e-3; //cm/ns
    double WindowSize     = detprop->NumberTimeSamples() * ts->TPCClock().TickPeriod() * 1e3;
  */

  art::ServiceHandle<geo::Geometry> geom;
  calo::CalorimetryAlg fCalorimetryAlg;

  /*
  const double region_of_interest = 20;
  const double theta0XZBinSize = 0.005;
  const int smoothing_window = 3;

  const double thetaBinSize = 0.06;
  const double costhetaBinSize = 0.06;
  const double phiBinSize = 0.06;

  const int num_bin = (int)(6.28/theta0XZBinSize);
  const int num_bin_theta = (int)(3.14/thetaBinSize);
  const int num_bin_phi = (int)(6.28/phiBinSize);

  const double growing_fit_initial_length = 10.;
  const double initial_fit_distance_to_line = 1;

  long unsigned int min_initial_hits_found = 5;
  //const int min_initial_hits_found = 5;
  const int min_initial_hits = 15;
  const int max_fitting_hits_micro = 15;

  const int macro_sliding_fit_window =10;
  const int micro_sliding_fit_window = 5;
  const double growing_fit_segment_length = 5.;

  //const double distance_to_line;
  //const double hit_connection_distance;

  const double distance_to_line_default = 0.75;
  const double hit_connection_distance_default = 0.75;

  const double distance_to_line_spineall = 3;                                                             
                                      
  const double hit_connection_distance_spineall = 5;

  const double distance_to_line_micro = 3;
  const double hit_connection_distance_micro = 5;

  double distance_to_line = 1;
  double hit_connection_distance = 1;
  double max_fitting_hits = 15;
  */
  
};// class HitSplitAlg
 DEFINE_ART_MODULE(HitSplitAlg);
}
#endif
