#include "ReconstructionOrchestrator.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/Geometry.h"
#include "lardata/ArtDataHelper/HitCreator.h"


#include "art/Persistency/Common/PtrMaker.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArObjects/LArPfoObjects.h"
#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "larpandora/LArPandoraInterface/LArPandoraOutput.h" 

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/ArtDataHelper/MVAReader.h"

#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "larsim/EventWeight/Base/MCEventWeight.h"

#include "lardataobj/RecoBase/TrackHitMeta.h"
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

#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"
#include "larcoreobj/SummaryData/POTSummary.h"

#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h"  

#include "larcorealg/Geometry/GeometryCore.h" 

#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "larreco/RecoAlg/TrajectoryMCSFitter.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"




// Framework includes

#include "art_root_io/TFileService.h"
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/View.h"
#include "art/Utilities/make_tool.h"

#include "canvas/Utilities/InputTag.h"
#include "canvas/Utilities/ensurePointer.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
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

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/simb.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/GTruth.h"

#include "Pandora/PdgTable.h"

// ROOT includes 
#include "TDirectory.h" 

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
#include <TVector3.h>
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

#include <memory>

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

namespace kaon_reconstruction{
  
  class RunningProducerwKAlg : public art::EDProducer {
    
  public:
    explicit RunningProducerwKAlg(fhicl::ParameterSet const & p);
    virtual ~RunningProducerwKAlg();
    RunningProducerwKAlg(RunningProducerwKAlg const &) = delete;
    RunningProducerwKAlg(RunningProducerwKAlg &&) = delete;
    RunningProducerwKAlg & operator = (RunningProducerwKAlg const &) = delete;
    RunningProducerwKAlg & operator = (RunningProducerwKAlg &&) = delete;
    
    // Required functions.
    void produce(art::Event & event) override;
    void reconfigure(fhicl::ParameterSet const& pset);
    double CalculateTrackDauDistance(art::Ptr<recob::Track>& track, art::Ptr<recob::Track>& dau_track);
    void CalculateTrackDauDistance(art::Ptr<recob::Track>& track, art::Ptr<recob::Track>& dau_track, double& track_dau_distance, bool& is_vtx, bool& is_end);
    
    simb::MCParticle const* truthMatchTrack(std::vector<art::Ptr<recob::Hit>>& hit_list, art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData>& particles_per_hit);
    simb::MCParticle const* truthMatchHit(art::Ptr<recob::Hit>& hit, art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData>& particles_per_hit);
    
    
  private:
    
    // The module name of the raw digits we're reading in
    
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

    double fFidVolCutX;
    double fFidVolCutY;
    double fFidVolCutZ;

    double fFidVolXmin;
    double fFidVolXmax;
    double fFidVolYmin;
    double fFidVolYmax;
    double fFidVolZmin;
    double fFidVolZmax;

    double      fExponentConstant;
    double 	fPIDA_endPoint;
    double	fMaxPIDAValue;
    double	fMinPIDAValue;
    bool	fSaveMCTree;
    
    double fPidValue;
    unsigned int    fView;
    int    n_recoTracks;
    int    n_recoShowers;
    
    //detinfo::DetectorProperties const * theDetector ;
    //detinfo::DetectorClocks    const *  detClocks   ;

    std::vector<art::Ptr<recob::Hit>> fHits;
    std::vector<art::Ptr<recob::Track>> fTracks;
    std::vector<art::Ptr<recob::Shower>> fShowers;
    std::vector<art::Ptr<recob::SpacePoint>> fSpacePoints;
    std::map<art::Ptr<recob::Hit>, art::Ptr<recob::SpacePoint>> fHitsToSpacePoints;
    std::map<art::Ptr<recob::SpacePoint>, art::Ptr<recob::Hit>> fSpacePointsToHits;

  };
  

  RunningProducerwKAlg::RunningProducerwKAlg(fhicl::ParameterSet const& parameterSet)
    : EDProducer(parameterSet)
  {
    reconfigure(parameterSet);

    produces< std::vector<recob::Track> >(); 
    produces< art::Assns<recob::Track, recob::Hit> >();

    //theDetector = lar::providerFrom<detinfo::DetectorPropertiesService>();
    //detClocks   = lar::providerFrom<detinfo::DetectorClocksService>();
  }

  RunningProducerwKAlg::~RunningProducerwKAlg(){
    //destructor 
  }

  void RunningProducerwKAlg::reconfigure(fhicl::ParameterSet const& p)
  {

    fHitModuleLabel      = p.get<std::string>("HitModuleLabel");
    fHitModuleLabelOLD      = p.get<std::string>("HitModuleLabelOLD");
    fHitSPAssns          = p.get<std::string>("HitSPAssns");
    fTrackModuleLabel    = p.get<std::string>("TrackModuleLabel");
    fPidModuleLabel    = p.get<std::string>("PidModuleLabel");
    fCaloModuleLabel     = p.get<std::string>("CaloModuleLabel");
    fOpFlashModuleLabel  = p.get<std::string>("OpFlashModuleLabel");
    fShowerModuleLabel = p.get<std::string>("ShowerModuleLabel");
    fPointIDModuleLabel  = p.get<std::string>("PointIDModuleLabel");
    fNNetModuleLabel     = p.get<std::string>("NNetModuleLabel");
    fMCgenieLabel = p.get<std::string>("MCgenieLabel");
    fSpacePointproducer  = p.get<std::string>("SpacePointproducer");
    fHitTrackAssns = p.get<std::string>("HitTrackAssns");
    fMaxPIDAValue  = p.get<double>("MaxPIDAValue");
    fMinPIDAValue  = p.get<double>("MinPIDAValue");
    fPIDA_endPoint = p.get<double>("PIDA_endPoint");
    fView                = p.get<double>("View");
    fPidValue = p.get<double>("PidValue");
    fFidVolCutX          = p.get<double>("FidVolCutX");
    fFidVolCutY          = p.get<double>("FidVolCutY");
    fFidVolCutZ          = p.get<double>("FidVolCutZ");
  }

  

void RunningProducerwKAlg::produce(art::Event & event)
{

  std::vector< art::Ptr<simb::MCParticle> > ptList;

  std::unique_ptr< std::vector<recob::Track> > anaTrackCollection(new std::vector<recob::Track>);
  std::unique_ptr<art::Assns<recob::Track, recob::Hit>>  anaTrackHitAssociations(new art::Assns<recob::Track, recob::Hit>);


  //save GENIE stuff for atmos

  art::Handle<std::vector<simb::MCTruth>> MCtruthHandle;
  std::vector<art::Ptr<simb::MCTruth>> MCtruthlist;
  if(event.getByLabel(fMCgenieLabel, MCtruthHandle))
    art::fill_ptr_vector(MCtruthlist, MCtruthHandle);
  
  //---- Reco stuff ----//

  // get tracks
  art::Handle< std::vector<recob::Track> > trackListHandle;
  std::vector<art::Ptr<recob::Track>> tracklist;
  if( event.getByLabel(fTrackModuleLabel, trackListHandle))
    art::fill_ptr_vector(tracklist, trackListHandle);
  n_recoTracks = tracklist.size();
  if( n_recoTracks > MAX_TRACKS || n_recoTracks == 0) {
    n_recoTracks = 0; // make sure we don't save any fake data
    event.put(std::move(anaTrackCollection));
    event.put(std::move(anaTrackHitAssociations));
    return;
  }

  // get hits
  art::Handle<std::vector<recob::Hit>> HitHandle;
  std::vector<art::Ptr<recob::Hit>> all_hits;
  if(event.getByLabel(fHitModuleLabel,HitHandle))
    art::fill_ptr_vector(all_hits, HitHandle);
  
  // get spacepoints
  art::Handle<std::vector<recob::SpacePoint>> SpacePointHandle;
  std::vector<art::Ptr<recob::SpacePoint>> SpacePointlist;
  if(event.getByLabel(fSpacePointproducer,SpacePointHandle))
    art::fill_ptr_vector(SpacePointlist, SpacePointHandle);
  
  // get showers
  art::Handle<std::vector<recob::Shower>> showerHandle;
  std::vector<art::Ptr<recob::Shower>> showerlist;
  if(event.getByLabel(fShowerModuleLabel,showerHandle))
    art::fill_ptr_vector(showerlist, showerHandle);
  n_recoShowers= showerlist.size();
  
  
    // get vtx related to track
  art::FindManyP<recob::Vertex> trk_from_vtx(trackListHandle,event,"pmtrack");
  
  art::Handle< std::vector<recob::Vertex> > vtxListHandle;
  std::vector<art::Ptr<recob::Vertex>> vtxlist;
  if(event.getByLabel("pmtrack", vtxListHandle))
    art::fill_ptr_vector(vtxlist, vtxListHandle);
  
  // get hits and calo of tracks
  art::FindManyP<recob::Hit> track_hits(trackListHandle, event, fHitTrackAssns);
  art::Handle<std::vector<recob::Hit>> hitsHandle;
  event.getByLabel(fHitModuleLabel, hitsHandle);

  fHits.clear();
  fSpacePoints.clear();
  fTracks.clear();
  fShowers.clear();
  fSpacePointsToHits.clear();
  fHitsToSpacePoints.clear();

  for (unsigned int iHit = 0; iHit < hitsHandle->size(); ++iHit) {
    const art::Ptr<recob::Hit> hit(hitsHandle, iHit);
    fHits.push_back(hit);
  }
  cout << "fHits.size(): " << fHits.size() << endl;

  art::Handle<std::vector<recob::Track>> tracksHandle;  
  event.getByLabel(fTrackModuleLabel, tracksHandle);
  for (unsigned int iTrack = 0; iTrack < tracksHandle->size(); ++iTrack) {
    const art::Ptr<recob::Track> track(tracksHandle, iTrack);
    fTracks.push_back(track); 
  }
  cout << "fTracks.size(): " << fTracks.size() << endl;
  
  art::Handle<std::vector<recob::Shower>> showersHandle; 
  event.getByLabel(fShowerModuleLabel, showersHandle); 
  for (unsigned int iShower = 0; iShower < showersHandle->size(); ++iShower) {
    const art::Ptr<recob::Shower> shower(showersHandle, iShower);
    fShowers.push_back(shower);
  }
  cout << "fShowers.size(): " << fShowers.size() << endl;
  
  art::Handle<std::vector<recob::SpacePoint>> spHandle;
  event.getByLabel(fSpacePointproducer, spHandle);
  for (unsigned int iSP = 0; iSP < spHandle->size(); ++iSP) {
    const art::Ptr<recob::SpacePoint> spacePoint(spHandle, iSP);
    fSpacePoints.push_back(spacePoint);
  }
  cout << "fSpacePoints.size(): " << fSpacePoints.size() << endl;
  
  art::FindOneP<recob::Hit> findSPToHits(fSpacePoints, event, fSpacePointproducer); 
  art::FindManyP<recob::Hit> findTracksToHits(fTracks, event, fTrackModuleLabel);
  art::FindManyP<recob::Hit> findShowersToHits(fShowers, event, fShowerModuleLabel);

  auto hit_handle = event.getValidHandle<std::vector<recob::Hit>>("hitfd");
  std::vector<art::Ptr<recob::Hit>> hit_list;
  art::fill_ptr_vector(hit_list, hit_handle);
  auto findHitToSPs = art::FindManyP<recob::SpacePoint>(hit_handle, event, fHitSPAssns); 
  

  for (unsigned int iSP = 0; iSP < fSpacePoints.size(); ++iSP) { 

    const art::Ptr<recob::SpacePoint> spacePoint = fSpacePoints.at(iSP);
    const art::Ptr<recob::Hit> hit = findSPToHits.at(iSP);
    fSpacePointsToHits[spacePoint] = hit;
    fHitsToSpacePoints[hit] = spacePoint;
    
  }


  std::vector<art::Ptr<recob::Hit>> unused_hits;	
  std::vector<int> unused_tracks;

  // loop reco tracks
  for(int i=0; i<n_recoTracks; ++i) {
    
    art::Ptr<recob::Track> track = tracklist[i];
    
    std::vector<art::Ptr<recob::Hit>> hits_from_track = track_hits.at(i);

    bool is_vtx = false;
    bool is_end = false;

    for(int j=0; j<n_recoTracks; ++j) {

      art::Ptr<recob::Track> dau_track = tracklist[j];

      std::vector<art::Ptr<recob::Hit>> hits_from_dau_track = track_hits.at(j);

      if(dau_track->ID() == track->ID()) continue;

      //double track_dau_distance = CalculateTrackDauDistance(track, dau_track);
      double track_dau_distance = 999;
      bool is_vtx_tmp = false;
      bool is_end_tmp = false;

      CalculateTrackDauDistance(track, dau_track, track_dau_distance, is_vtx_tmp, is_end_tmp);

      if(track_dau_distance<5){

	if(is_vtx == false && is_vtx_tmp == true) is_vtx = is_vtx_tmp; //has a common vertex at vtx of tracklist[i]
	if(is_end == false && is_end_tmp == true) is_end = is_end_tmp; //has a common vertex at end of tracklist[i]
	/*
	if(std::find(unused_tracks.begin(), unused_tracks.end(), i) == unused_tracks.end()){
	  unused_tracks.push_back(i);
	  unused_hits.insert(unused_hits.end(), hits_from_track.begin(), hits_from_track.end());
	}
	if(std::find(unused_tracks.begin(), unused_tracks.end(), j) == unused_tracks.end()){
	  unused_tracks.push_back(j);
	  unused_hits.insert(unused_hits.end(), hits_from_dau_track.begin(), hits_from_dau_track.end());
	}
	*/

      }
    }

    std::vector<recob::Track> rebuildTrackList;
    std::vector<std::vector<art::Ptr<recob::Hit>>> trackHitLists;

    if(is_vtx&&is_end) continue;
    if(!is_vtx){

      ReconstructionOrchestrator orchestrator_vtx;
      orchestrator_vtx.runReconstruction(SpacePointlist, fSpacePointsToHits, fHitsToSpacePoints, track, hits_from_track, "vtx");
      std::vector<recob::Track> rebuildTrackList_vtx = orchestrator_vtx.getRebuildTrackList();
      std::vector<std::vector<art::Ptr<recob::Hit>>> trackHitLists_vtx = orchestrator_vtx.getHitLists();

      rebuildTrackList.insert(rebuildTrackList.end(), rebuildTrackList_vtx.begin(), rebuildTrackList_vtx.end());
      trackHitLists.insert(trackHitLists.end(), trackHitLists_vtx.begin(), trackHitLists_vtx.end());

    }
    if(!is_end){

      ReconstructionOrchestrator orchestrator_end;
      orchestrator_end.runReconstruction(SpacePointlist, fSpacePointsToHits, fHitsToSpacePoints, track, hits_from_track, "end");
      std::vector<recob::Track> rebuildTrackList_end = orchestrator_end.getRebuildTrackList(); 
      std::vector<std::vector<art::Ptr<recob::Hit>>> trackHitLists_end = orchestrator_end.getHitLists();

      rebuildTrackList.insert(rebuildTrackList.end(), rebuildTrackList_end.begin(), rebuildTrackList_end.end());
      trackHitLists.insert(trackHitLists.end(), trackHitLists_end.begin(), trackHitLists_end.end());

    }

    for(unsigned int j=0; j < rebuildTrackList.size(); j++) {
      
      if(trackHitLists[j].empty()) continue;
      
      anaTrackCollection->push_back(rebuildTrackList[j]);
      std::vector<art::Ptr<recob::Hit>> hits_from_track_rebuild = trackHitLists[j];
      
      lar_pandora::HitVector anaHitCollection_rebuild_tmp;
      for(auto hitptr : hits_from_track_rebuild){
	anaHitCollection_rebuild_tmp.push_back(hitptr);
      }
      
      util::CreateAssn(*this, event, *(anaTrackCollection.get()), anaHitCollection_rebuild_tmp, *(anaTrackHitAssociations.get()));
      
    }

  }

  /*
  for(int i=0; i<n_recoTracks; ++i) { 

    if(std::find(unused_tracks.begin(), unused_tracks.end(), i) != unused_tracks.end()) continue;

    art::Ptr<recob::Track> track = tracklist[i];
    std::vector<art::Ptr<recob::Hit>> hits_from_track = track_hits.at(i);

    std::vector<art::Ptr<recob::Hit>> unused_hits_tmp = unused_hits;
    unused_hits_tmp.insert(unused_hits_tmp.end(), hits_from_track.begin(), hits_from_track.end());

    ReconstructionOrchestrator orchestrator;
    //orchestrator.runReconstruction(SpacePointlist, fSpacePointsToHits, fHitsToSpacePoints, track, hits_from_track);
    orchestrator.runReconstruction(SpacePointlist, fSpacePointsToHits, fHitsToSpacePoints, track, unused_hits_tmp);
    std::cout << "escaped reconstruction" << std::endl;

    std::vector<recob::Track> rebuildTrackList = orchestrator.getRebuildTrackList();
    std::vector<std::vector<art::Ptr<recob::Hit>>> trackHitLists = orchestrator.getHitLists();

      for(unsigned int j=0; j < rebuildTrackList.size(); j++) {
	
	if(trackHitLists[j].empty()) continue;

	anaTrackCollection->push_back(rebuildTrackList[j]);
	std::vector<art::Ptr<recob::Hit>> hits_from_track_rebuild = trackHitLists[j];

	lar_pandora::HitVector anaHitCollection_rebuild_tmp;
	for(auto hitptr : hits_from_track_rebuild){
	  anaHitCollection_rebuild_tmp.push_back(hitptr);
	}

	util::CreateAssn(*this, event, *(anaTrackCollection.get()), anaHitCollection_rebuild_tmp, *(anaTrackHitAssociations.get()));
	
      }
  }
*/

  event.put(std::move(anaTrackCollection));
  event.put(std::move(anaTrackHitAssociations));
}
  
  double RunningProducerwKAlg::CalculateTrackDauDistance(art::Ptr<recob::Track>& track, art::Ptr<recob::Track>& dau_track) {
    double track_dau_distance_ev = TMath::Sqrt((track->End().x()-dau_track->Vertex().x())*(track->End().x()-dau_track->Vertex().x()) +
					       (track->End().y()-dau_track->Vertex().y())*(track->End().y()-dau_track->Vertex().y()) +
					       (track->End().z()-dau_track->Vertex().z())*(track->End().z()-dau_track->Vertex().z()));
    
    double track_dau_distance_vv = TMath::Sqrt((track->Vertex().x()-dau_track->Vertex().x())*(track->Vertex().x()-dau_track->Vertex().x()) +
					       (track->Vertex().y()-dau_track->Vertex().y())*(track->Vertex().y()-dau_track->Vertex().y()) +
					       (track->Vertex().z()-dau_track->Vertex().z())*(track->Vertex().z()-dau_track->Vertex().z()));
    
    double track_dau_distance_ve = TMath::Sqrt((track->Vertex().x()-dau_track->End().x())*(track->Vertex().x()-dau_track->End().x()) +
					       (track->Vertex().y()-dau_track->End().y())*(track->Vertex().y()-dau_track->End().y()) +
					       (track->Vertex().z()-dau_track->End().z())*(track->Vertex().z()-dau_track->End().z()));
    
    double track_dau_distance_ee = TMath::Sqrt((track->End().x()-dau_track->End().x())*(track->End().x()-dau_track->End().x()) +
					       (track->End().y()-dau_track->End().y())*(track->End().y()-dau_track->End().y()) +
					       (track->End().z()-dau_track->End().z())*(track->End().z()-dau_track->End().z()));
   
    return std::min({track_dau_distance_ev, track_dau_distance_vv, track_dau_distance_ve, track_dau_distance_ee}); 
    //double track_dau_distance = std::min({track_dau_distance_ev, track_dau_distance_vv, track_dau_distance_ve, track_dau_distance_ee});
  }

  void RunningProducerwKAlg::CalculateTrackDauDistance(art::Ptr<recob::Track>& track, art::Ptr<recob::Track>& dau_track, double& track_dau_distance, bool& is_vtx, bool& is_end) {

    double track_dau_distance_ev = TMath::Sqrt((track->End().x()-dau_track->Vertex().x())*(track->End().x()-dau_track->Vertex().x()) +
					       (track->End().y()-dau_track->Vertex().y())*(track->End().y()-dau_track->Vertex().y()) +
					       (track->End().z()-dau_track->Vertex().z())*(track->End().z()-dau_track->Vertex().z()));
    
    double track_dau_distance_vv = TMath::Sqrt((track->Vertex().x()-dau_track->Vertex().x())*(track->Vertex().x()-dau_track->Vertex().x()) +
					       (track->Vertex().y()-dau_track->Vertex().y())*(track->Vertex().y()-dau_track->Vertex().y()) +
					       (track->Vertex().z()-dau_track->Vertex().z())*(track->Vertex().z()-dau_track->Vertex().z()));
    
    double track_dau_distance_ve = TMath::Sqrt((track->Vertex().x()-dau_track->End().x())*(track->Vertex().x()-dau_track->End().x()) +
					       (track->Vertex().y()-dau_track->End().y())*(track->Vertex().y()-dau_track->End().y()) +
					       (track->Vertex().z()-dau_track->End().z())*(track->Vertex().z()-dau_track->End().z()));
    
    double track_dau_distance_ee = TMath::Sqrt((track->End().x()-dau_track->End().x())*(track->End().x()-dau_track->End().x()) +
					       (track->End().y()-dau_track->End().y())*(track->End().y()-dau_track->End().y()) +
					       (track->End().z()-dau_track->End().z())*(track->End().z()-dau_track->End().z()));
   
    track_dau_distance = std::min({track_dau_distance_ev, track_dau_distance_vv, track_dau_distance_ve, track_dau_distance_ee}); 
    
    if(track_dau_distance == track_dau_distance_vv || track_dau_distance == track_dau_distance_ve)
      is_vtx = true;

    if(track_dau_distance == track_dau_distance_ee || track_dau_distance == track_dau_distance_ev)
      is_end = true;

  }


  DEFINE_ART_MODULE(RunningProducerwKAlg)
}

