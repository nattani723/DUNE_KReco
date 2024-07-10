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
  
  class RunningKRcoAlg : public art::EDProducer {
    
  public:
    explicit RunningKRcoAlg(fhicl::ParameterSet const & p);
    virtual ~RunningKRcoAlg();
    RunningKRcoAlg(RunningKRcoAlg const &) = delete;
    RunningKRcoAlg(RunningKRcoAlg &&) = delete;
    RunningKRcoAlg & operator = (RunningKRcoAlg const &) = delete;
    RunningKRcoAlg & operator = (RunningKRcoAlg &&) = delete;
    
    // Required functions.
    void produce(art::Event & event) override;
    void reconfigure(fhicl::ParameterSet const& pset);
    
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
  

  RunningKRcoAlg::RunningKRcoAlg(fhicl::ParameterSet const& parameterSet)
    : EDProducer(parameterSet)
  {
    reconfigure(parameterSet);

    produces< std::vector<recob::Track> >(); 
    produces< art::Assns<recob::Track, recob::Hit> >();

    //theDetector = lar::providerFrom<detinfo::DetectorPropertiesService>();
    //detClocks   = lar::providerFrom<detinfo::DetectorClocksService>();
  }

  RunningKRcoAlg::~RunningKRcoAlg(){
    //destructor 
  }

  void RunningKRcoAlg::reconfigure(fhicl::ParameterSet const& p)
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

  

void RunningKRcoAlg::produce(art::Event & event)
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


  // loop reco tracks
  for(int i=0; i<n_recoTracks; ++i) {
    
    art::Ptr<recob::Track> track = tracklist[i];
    
    TVector3 end(track->End().x(),track->End().y(),track->End().z());
    
    for(int j=0; j<n_recoTracks; ++j) {

      std::vector<art::Ptr<recob::Hit>> hits_from_track = track_hits.at(i);

      ReconstructionOrchestrator orchestrator;
      orchestrator.runReconstruction(SpacePointlist, fSpacePointsToHits, fHitsToSpacePoints, track, hits_from_track);
      std::cout << "escaped reconstruction" << std::endl;

      std::vector<recob::Track> rebuildTrackList = orchestrator.getRebuildTrackList();
      std::vector<std::vector<art::Ptr<recob::Hit>>> trackHitLists = orchestrator.getHitLists();

      for(unsigned int j=0; j < rebuildTrackList.size(); j++) {
	
	const std::vector<TVector3> peakDirectionVector = orchestrator.getPeakDirectionList();

	/*
	if(trackHitLists[j].empty()){
	  event.put(std::move(anaTrackCollection));
	  event.put(std::move(anaTrackHitAssociations));
	  continue;
	}
	*/

	TVector3 pos(rebuildTrackList[j].Vertex().X(),rebuildTrackList[j].Vertex().Y(),rebuildTrackList[j].Vertex().Z());

	double track_dau_distance=TMath::Sqrt((end.X()-pos.X())*(end.X()-pos.X()) +
					      (end.Y()-pos.Y())*(end.Y()-pos.Y()) +
					      (end.Z()-pos.Z())*(end.Z()-pos.Z()));

	/*
	if (track_dau_distance>10){
	  event.put(std::move(anaTrackCollection));
	  event.put(std::move(anaTrackHitAssociations));	
	  continue;
	}
	*/
	std::cout << track_dau_distance << std::endl;

	anaTrackCollection->push_back(rebuildTrackList[j]);
	
	std::vector<art::Ptr<recob::Hit>> hits_from_track_rebuild = trackHitLists[j];
	
	lar_pandora::HitVector anaHitCollection_rebuild_tmp;
	for(auto hitptr : hits_from_track_rebuild){
	  anaHitCollection_rebuild_tmp.push_back(hitptr);
	}
	
	util::CreateAssn(*this, event, *(anaTrackCollection.get()), anaHitCollection_rebuild_tmp, *(anaTrackHitAssociations.get()));
	
      }
    }
    std::cout << i << std::endl;
  }
    event.put(std::move(anaTrackCollection));
    event.put(std::move(anaTrackHitAssociations));

}

DEFINE_ART_MODULE(RunningKRcoAlg)
}
