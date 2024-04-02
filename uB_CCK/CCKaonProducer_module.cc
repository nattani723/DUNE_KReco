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

 
#ifdef __CLING__
#pragma link C++ class std::vector < std::vector<Float_t> >+;
#pragma link C++ class std::vector < std::vector< std::vector<Float_t> > >+;
#endif


using namespace std;
namespace Kaon_Analyzer{

  CCKaonProducer::CCKaonProducer(fhicl::ParameterSet const& pset) :
    EDProducer(pset),
    fHitsModuleLabel          (pset.get< std::string >("HitsModuleLabel","gaushit")),
    fLArG4ModuleLabel         (pset.get< std::string >("LArG4ModuleLabel","largeant")),
    fGenieGenModuleLabel      (pset.get< std::string >("GenieGenModuleLabel","generator")),
    fTrackModuleLabel         (pset.get< std::string >("TrackModuleLabel","pandora")),
    fShowerModuleLabel        (pset.get< std::string >("ShowerModuleLabel","pandora")),
    fCalorimetryModuleLabel   (pset.get< std::string >("CalorimetryModuleLabel","pandoracaliSCE")),
    fPIDLabel                 (pset.get< std::string >("PIDLabel","pandoracalipidSCE")),
    fHitTruthAssns            (pset.get< std::string >("HitTruthAssn","gaushitTruthMatch")),
    fHitTrackAssns            (pset.get< std::string >("HitTrackAssn","pandora")),
    fHitShowerAssns            (pset.get< std::string >("HitShowerAssn","pandora")),
    m_pfp_producer            (pset.get< std::string >("pfp_producer","pandora")),
    fPFParticleLabel          (pset.get< std::string >("PFParticleLabel", "pandora")),
    fSpacePointproducer       (pset.get< std::string >("SpacePointproducer", "pandora")),
    fPandoraLabel             (pset.get< std::string >("PandoraLabel", "pandora")),
    isMC                      (pset.get< bool >("IsMC",true))

  {

    theDetector = lar::providerFrom<detinfo::DetectorPropertiesService>();
    detClocks   = lar::providerFrom<detinfo::DetectorClocksService>();
    SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
    geom = lar::providerFrom<geo::Geometry>();

    produces< std::vector<recob::Track> >();
    produces< art::Assns<recob::Track, recob::Hit> >();
 
  }

  CCKaonProducer::~CCKaonProducer()
  {
    //destructor                                                                                                                                                                                                                                                                                                             
  }


  void CCKaonProducer::beginJob()
  {
    //art::ServiceHandle<art::TFileService> tfs;                                                                                                                                                                                                                                                                             
    //fEventTree->Branch();                                                                                                                                                                                                                                                                                                  
  }


  void CCKaonProducer::produce(art::Event& evt)
  {
    art::Handle< std::vector<simb::MCTruth> > mctruths;
    std::vector< art::Ptr<simb::MCParticle> > ptList;

    std::unique_ptr< std::vector<recob::Track> > anaTrackCollection(new std::vector<recob::Track> );
    std::unique_ptr<art::Assns<recob::Track, recob::Hit>>  anaTrackHitAssociations(new art::Assns<recob::Track, recob::Hit>);

    // get MCTruth                                                                                                                                                                                                                                                                                                        
    evt.getByLabel("generator", mctruths);
    if (mctruths->size()!=1) {
      //return;                                                                                                                                                                                                                                                                                                              
      //continue;                                                                                                                                                                                                                                                                                                            
    }
    simb::MCTruth mctruth = mctruths->at(0);


    // get MCParticles                                                                                                                                                                                                                                                                                                 
    art::Handle< std::vector<simb::MCParticle> > mcParticleHandle;
    if (evt.getByLabel(fLArG4ModuleLabel, mcParticleHandle)){
      art::fill_ptr_vector(ptList, mcParticleHandle);
    }


    // Collect all recontructed particles                                                                                                                                                                                                                                                                                    
    art::Handle<std::vector<recob::PFParticle>> pfparticles;
    std::vector< art::Ptr<recob::PFParticle> > pfpVector;
    if(evt.getByLabel(m_pfp_producer, pfparticles)){
      art::fill_ptr_vector(pfpVector, pfparticles);
    }

    if (pfparticles->size()==0) {
      std::cout << "No PFParticles found" << std::endl;
 
      evt.put(std::move(anaTrackCollection));
      evt.put(std::move(anaTrackHitAssociations));
      return;
    }


    // Get PFParticle associations                                                                                                                                                                                                                                                              
    art::FindManyP<anab::T0> pfp_muon_assn(pfparticles, evt, "NuCCproducer");
    if(!pfp_muon_assn.isValid()){
      evt.put(std::move(anaTrackCollection));
      evt.put(std::move(anaTrackHitAssociations));
      return;
    }

    art::FindManyP<recob::Track> pfparticleTrackAssn(pfparticles, evt, "pandora");
    if(!pfparticleTrackAssn.isValid()){
      evt.put(std::move(anaTrackCollection));
      evt.put(std::move(anaTrackHitAssociations));
      return;
    }

    art::FindManyP<recob::Vertex> pfparticleVertexAssn(pfparticles, evt, "pandora");
    if(!pfparticleVertexAssn.isValid()){
      evt.put(std::move(anaTrackCollection));
      evt.put(std::move(anaTrackHitAssociations));
      return;
    }


    // Find recontructed neutrino (there should be one)                                                                                                                                                                                                                                                                      
    lar_pandora::PFParticleVector pfneutrinos(0);
    for (unsigned int i=0; i<pfparticles->size(); ++i) {

      art::Ptr<recob::PFParticle> pfparticle(pfparticles,i);
      if (pfparticle->IsPrimary() && pfparticle->PdgCode()==14) {
        pfneutrinos.push_back(pfparticle);
      }

    }

    if (pfneutrinos.size() != 1) {
      evt.put(std::move(anaTrackCollection));
      evt.put(std::move(anaTrackHitAssociations));
      return;
    }

    art::Ptr<recob::PFParticle> pfnu = pfneutrinos.front();


    // Find CC muon and daughters                                                                                                                                                                                                                                                                                            
    lar_pandora::PFParticleVector pfmuons(0);
    vector<int> reco_nu_daughters_id(0);

    for (unsigned int i=0; i<pfparticles->size(); ++i) {

      art::Ptr<recob::PFParticle> pfparticle(pfparticles,i);

      // look at particles with neutrino parent and one associated track                                                                                                                                                                                                                                                     
      if (pfparticle->Parent()==pfnu->Self() && pfparticleTrackAssn.at(i).size()==1) {

	art::Ptr<recob::Track> track = pfparticleTrackAssn.at(i).front();
        reco_nu_daughters_id.push_back(track.key());

        // CC muon has a T0 associated                                                                                                                                                                                                                                                                                       
        if (pfp_muon_assn.at(i).size()==1) {
          pfmuons.push_back(pfparticle);
        }
      }
    }

    reco_nu_ndaughters = reco_nu_daughters_id.size();

    if (pfmuons.size()!=1) {
      evt.put(std::move(anaTrackCollection));
      evt.put(std::move(anaTrackHitAssociations));
      return;
    }

    art::Ptr<recob::PFParticle> pfmuon = pfmuons.front();
    art::Ptr<recob::Track> trkmuon = pfparticleTrackAssn.at(pfmuon.key()).front();


    // get hits                                                                                                                                                                                                                                                                                                              
    art::Handle< std::vector<recob::Hit> > hitListHandle;
    std::vector<art::Ptr<recob::Hit> > hitlist;
    if(evt.getByLabel(fHitsModuleLabel,hitListHandle)){
      art::fill_ptr_vector(hitlist, hitListHandle);
    }

    // get tracks                                                                                                                                                                                                                                                                                                            
    art::Handle< std::vector<recob::Track> > trackListHandle;
    std::vector<art::Ptr<recob::Track> > tracklist;
    if(evt.getByLabel(fTrackModuleLabel,trackListHandle)) {
      art::fill_ptr_vector(tracklist, trackListHandle);
    }

    // get showers                                                                                                                                                                                                                                                                                                           
    art::Handle< std::vector<recob::Shower> > showerListHandle;
    std::vector<art::Ptr<recob::Shower> > showerlist;
    if(evt.getByLabel(fShowerModuleLabel,showerListHandle)) {
      art::fill_ptr_vector(showerlist, showerListHandle);
    }

    // get clusters
    art::Handle< std::vector<recob::Cluster> > clusterHandle;
    std::vector< art::Ptr<recob::Cluster> > clusterVector;
    if(evt.getByLabel(fPandoraLabel,clusterHandle)){
      art::fill_ptr_vector(clusterVector,clusterHandle);
    }

    // get spacepoints
    art::Handle< std::vector<recob::SpacePoint> > spacepointHandle;
    std::vector< art::Ptr<recob::SpacePoint> > spacepointVector;
    if(evt.getByLabel(fPandoraLabel,spacepointHandle)){
      art::fill_ptr_vector(spacepointVector,spacepointHandle);
    }


    art::FindManyP<recob::Hit> hits_from_tracks(trackListHandle, evt, fHitTrackAssns);
    art::FindManyP<recob::Hit> hits_from_showers(showerListHandle, evt, fHitShowerAssns);
    if(!hits_from_tracks.isValid()){
      evt.put(std::move(anaTrackCollection));
      evt.put(std::move(anaTrackHitAssociations));
      return;
    }

    art::FindManyP<recob::Cluster> clusters_per_pfparticle(pfparticles, evt, fPandoraLabel);
    art::FindManyP<recob::Hit> hits_per_cluster(clusterHandle, evt, fPandoraLabel);
    std::map<art::Ptr<recob::PFParticle>,  std::vector<art::Ptr<recob::Cluster>> > pfParticleToClustersMap;
    std::map<art::Ptr<recob::Cluster>,  std::vector<art::Ptr<recob::Hit>> > clusterToHitsMap;

    art::FindManyP<recob::SpacePoint> spacepoints_per_pfparticle(pfparticles, evt, fSpacePointproducer);
    art::FindManyP<recob::Hit> hits_per_spacepoint(spacepointHandle, evt, fSpacePointproducer);
    std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::SpacePoint>> > pfParticleToSpacePointsMap;
    std::map<art::Ptr<recob::SpacePoint>,  std::vector<art::Ptr<recob::Hit>> > spacepointToHitsMap;

    std::map<art::Ptr<recob::PFParticle>,  std::vector<art::Ptr<recob::Hit>> > pfParticleToHitsMap;
    std::map<art::Ptr<recob::PFParticle>,  std::vector<art::Ptr<recob::Hit>> > pfParticleToHitsFromSpacePointsMap;
    std::map<art::Ptr<recob::PFParticle>, int> pfParticleToNHitsFromSpacePoints;

    if( !clusters_per_pfparticle.isValid() || !hits_per_cluster.isValid() || !spacepoints_per_pfparticle.isValid() || !hits_per_spacepoint.isValid()  ){
      evt.put(std::move(anaTrackCollection));
      evt.put(std::move(anaTrackHitAssociations));
      return;
    }


    const art::PtrMaker<recob::Track> makeTrackPtr(evt);

    for (unsigned int i=0; i<pfpVector.size(); ++i) {
      auto pfp = pfpVector[i];
      pfParticleToClustersMap[pfp] = clusters_per_pfparticle.at(pfp.key());
      pfParticleToSpacePointsMap[pfp] = spacepoints_per_pfparticle.at(pfp.key());
    }

    for(size_t i=0; i< clusterVector.size(); ++i){
      auto cluster = clusterVector[i];
      clusterToHitsMap[cluster] = hits_per_cluster.at(cluster.key());
    }

    for(size_t i=0; i< spacepointVector.size(); ++i){
      auto spacepoint = spacepointVector[i];
      spacepointToHitsMap[spacepoint] = hits_per_spacepoint.at(spacepoint.key());
    }


    for (unsigned int i=0; i<pfpVector.size(); ++i) {
      auto pfp = pfpVector[i];
      std::vector<art::Ptr<recob::Cluster>> clusters_vec  = pfParticleToClustersMap[pfp];
      std::vector<art::Ptr<recob::SpacePoint>> spacepoints_vec  = pfParticleToSpacePointsMap[pfp];
      std::vector<art::Ptr<recob::Hit>> hits_for_pfp = {};

      lar_pandora::HitSet hitsInParticleSet;

      for (art::Ptr<recob::SpacePoint> spacepoint: spacepoints_vec){
	std::vector<art::Ptr<recob::Hit>> hits_vec =  spacepointToHitsMap[spacepoint];
        hits_for_pfp.insert( hits_for_pfp.end(), hits_vec.begin(), hits_vec.end() );
 
        for(auto hit : hits_vec){
          (void) hitsInParticleSet.insert(hit);
        }

      }
      pfParticleToNHitsFromSpacePoints[pfp] = hits_for_pfp.size();
      pfParticleToHitsFromSpacePointsMap[pfp] = hits_for_pfp;
 
      cout << "hits_for_pfp size before cluster " << hits_for_pfp.size() << endl;

      for (art::Ptr<recob::Cluster> cluster: clusters_vec){
	std::vector<art::Ptr<recob::Hit>> hits_vec =  clusterToHitsMap[cluster];

        for(auto hit : hits_vec){
          if(hitsInParticleSet.count(hit) == 0) hits_for_pfp.push_back(hit);
        }
      }
      cout << "hits_for_pfp size after cluster " << hits_for_pfp.size() <<endl; 
      pfParticleToHitsMap[pfp] = hits_for_pfp;
    }



    // find track multiplicity around vertex                                                                                                                                                                                                                                                                                 
    int NTracks=tracklist.size();
    int NShowers=showerlist.size();

    for (int i=0; i<reco_nu_ndaughters; i++) {

      vector<bool> v_trk_flg_peak;
      std::map<double, TVector2, std::greater<>> view_peak_map;
      vector<TVector2> best_peak_bins;

      std::map<int, std::map<int, double>> angular_distribution_map_3D_track;
      std::map<int, std::map<int, double>> angular_distribution_map_3D_shower;
      std::map<int, std::map<int, double>> angular_distribution_map_3D_pfparticle;

      std::vector<art::Ptr<recob::Hit>> unavailable_hit_list;
      std::vector<art::Ptr<recob::Hit>> shower_spine_hit_list;
      std::vector<std::vector<art::Ptr<recob::Hit>>> shower_spine_hit_list_vector;
      std::vector<art::Ptr<recob::Hit>> hits_from_reco;


      art::Ptr<recob::Track> ptrack(trackListHandle,reco_nu_daughters_id[i]);
      const recob::Track& track = *ptrack;
      lar_pandora::HitVector hitsInParticle;

      unsigned int hitsFromSpacePointsSize = 0;

      for (unsigned int i=0; i<pfparticles->size(); ++i) {
   
	art::Ptr<recob::PFParticle> pfparticle(pfparticles,i);
	if (pfparticle->Parent()==pfnu->Self() && pfparticleTrackAssn.at(i).size()==1){
	  art::Ptr<recob::Track> ppfptrack = pfparticleTrackAssn.at(i).front();

	  if(ppfptrack.key() == ptrack.key()){
	    hitsInParticle = pfParticleToHitsMap[pfparticle];
	    hitsFromSpacePointsSize = pfParticleToNHitsFromSpacePoints[pfparticle];
   
	    cout << "hitsInParticle.size()" << hitsInParticle.size() << endl;
	    break;
	  }
	}
      }

      std::vector<art::Ptr<recob::Hit>> hits_from_track = hits_from_tracks.at(ptrack.key());
      art::FindManyP<recob::SpacePoint> spacepoint_per_hit_test(hitListHandle, evt, fSpacePointproducer);

      std::vector<art::Ptr<recob::SpacePoint>> spacepoint_vec_test;

      for(size_t i_h=0; i_h<hits_from_track.size(); i_h++){
	spacepoint_vec_test.clear();
	spacepoint_vec_test = spacepoint_per_hit_test.at(hits_from_track[i_h].key());

      }

      cout << "hits_from_track.size(): " << hits_from_track.size() << endl;
      cout << "hitsFromSpacePointsSize: " << hitsFromSpacePointsSize << endl;

      lar_pandora::HitVector anaHitCollection_tmp;
      anaHitCollection_tmp.clear();

      for(auto hitptr : hits_from_track){   
	anaHitCollection_tmp.push_back(hitptr);
      }


      art::Ptr<recob::Track> pTrack(makeTrackPtr(anaTrackCollection->size() - 1));

      // ATTN: metadata added with index from space points if available, null for others                                                                                                                                                                                                                                   
      for (unsigned int hitIndex = 0; hitIndex < hitsInParticle.size(); hitIndex++){
        const art::Ptr<recob::Hit> phit(hitsInParticle.at(hitIndex));
        const int index((hitIndex < hitsFromSpacePointsSize) ? hitIndex : std::numeric_limits<int>::max());
	recob::TrackHitMeta metadata(index, -std::numeric_limits<double>::max());
      }


      // skip cc muon track                                                                                                                                                                                                                                                                         
      if (ptrack.key()==trkmuon.key()) continue;

      // check track start and end                                                                                                                                                                                                                                                                                         
      TVector3 pos(track.Vertex().X(),track.Vertex().Y(),track.Vertex().Z());
      TVector3 end(track.End().X(),track.End().Y(),track.End().Z());


      int ndaughters = 0;
      int ndaughters_sh = 0;

      cout << "before NTracks loop" << endl;
      cout << "Ntracks: " << NTracks << endl;

      for (int j=0; j<NTracks; j++) {

	art::Ptr<recob::Track> ptrack_dau(trackListHandle,j);
	const recob::Track& track_dau = *ptrack_dau;

	// skip all primary tracks                                                                                                                                                                                                                                                                        
	bool skip = false;
	for (int k=0; k<reco_nu_ndaughters; k++) {
	  if (int(ptrack_dau.key())==reco_nu_daughters_id[k]) {
	    skip=true;
	    break;
	  }
	}
	if (skip) continue;

	TVector3 pos2(track_dau.Vertex().X(),track_dau.Vertex().Y(),track_dau.Vertex().Z());

	double track_dau_distance=TMath::Sqrt((end.X()-pos2.X())*(end.X()-pos2.X()) +
					      (end.Y()-pos2.Y())*(end.Y()-pos2.Y()) +
					      (end.Z()-pos2.Z())*(end.Z()-pos2.Z()));

	cout << "distance to vertex: " << track_dau_distance << endl;

	if (track_dau_distance<40) { 

	  TVector3 end2(track_dau.End().X(),track_dau.End().Y(),track_dau.End().Z());

	  //if(true_kaon_end_process == 1 || true_kaon_end_process == 1){
	  //if(reco_track_true_pdg[ntracks] == 321){

	  std::vector<art::Ptr<recob::Hit>> hits_from_dau_track = hits_from_tracks.at(ptrack_dau.key());
	  art::FindManyP<recob::SpacePoint> spacepoint_per_hit(hitListHandle, evt, fSpacePointproducer);
	  hits_from_reco.insert(hits_from_reco.end(), hits_from_dau_track.begin(), hits_from_dau_track.end());

	  //fillAngularDistributionMap(hits_from_track, end, spacepoint_per_hit, angular_distribution_map_track);                                                                                                                                                        
	  fillAngularDistributionMap3D(hits_from_dau_track, end, spacepoint_per_hit, angular_distribution_map_3D_track);
	  ndaughters++;

	}//gap between K+ end and daughter                                                                                                                                                                                                   

      }//NTrack loop                            


      for (int j_s=0; j_s<NShowers; j_s++) {
	art::Ptr<recob::Shower> pshower_dau(showerListHandle,j_s);
	const recob::Shower& shower_dau = *pshower_dau;

	// skip all primary showers 
	bool skip2 = false;
	for (int k_sh=0; k_sh<reco_nu_ndaughters; k_sh++) {
	  if (int(pshower_dau.key())==reco_nu_daughters_id[k_sh]) {
	    skip2=true;
	    break;
	  }
	}
	if (skip2) continue;
	TVector3 pos2_sh(shower_dau.ShowerStart().X(),shower_dau.ShowerStart().Y(),shower_dau.ShowerStart().Z());

	double shower_dau_distance=TMath::Sqrt((end.X()-pos2_sh.X())*(end.X()-pos2_sh.X()) +
					       (end.Y()-pos2_sh.Y())*(end.Y()-pos2_sh.Y()) +
					       (end.Z()-pos2_sh.Z())*(end.Z()-pos2_sh.Z()));

	if (shower_dau_distance<40) {

	  if(ndaughters_sh>20) break;

	  std::vector<art::Ptr<recob::Hit>> hits_from_shower = hits_from_showers.at(pshower_dau.key());
	  art::FindManyP<recob::SpacePoint> spacepoint_per_hit(hitListHandle, evt, fSpacePointproducer);
	  hits_from_reco.insert(hits_from_reco.end(), hits_from_shower.begin(), hits_from_shower.end());


	  //fillTrueMatching_sh(hits_from_shower, particles_per_hit, ntracks, ndaughters_sh);
	  //fillShowerMatching(hits_from_shower, particles_per_hit, ntracks, ndaughters_sh);

	  //if(true_kaon_end_process == 1 || true_kaon_end_process == 1){
	  //if(reco_track_true_pdg[ntracks] == 321){

	  //fillAngularDistributionMap(hits_from_shower, end, spacepoint_per_hit, angular_distribution_map_shower);    
	  fillAngularDistributionMap3D(hits_from_shower, end, spacepoint_per_hit, angular_distribution_map_3D_shower);

	  //fillAngularDistributionMapCheat(hits_from_shower, end, particles_per_hit, spacepoint_per_hit, angular_distribution_mrgid_map, mrgidpdg, currentMergedId); 
	  //fillAngularDistributionMapCheat3D(hits_from_shower, true_pi0_hit_list, end, particles_per_hit, spacepoint_per_hit, angular_distribution_mrgid_map_3D, mrgidpdg_3D, mrgidmom, currentMergedId_3D);
	  //}
	  //}

	  ndaughters_sh++;
	}
      }


      smoothAngularDistributionMap3D(angular_distribution_map_3D_track);
      smoothAngularDistributionMap3D(angular_distribution_map_3D_shower);

      accumulateAngularDistributionMap3D(angular_distribution_map_3D_track, angular_distribution_map_3D_shower, angular_distribution_map_3D_pfparticle);

      //obtainPeakVectorCheat3D(angular_distribution_mrgid_map_3D, mrgidpdg_3D, view_peak_map_cheat, v_pdg_peak);                                                                                                                       
      obtainPeakVector3D(angular_distribution_map_3D_track, v_trk_flg_peak, view_peak_map, true);
      obtainPeakVector3D(angular_distribution_map_3D_shower, v_trk_flg_peak, view_peak_map, false);

      //findBestAngularPeakCheat3D(view_peak_map_cheat, best_peak_bins_cheat);                                                                                                                                                              
      findBestAngularPeak3D(angular_distribution_map_3D_pfparticle, view_peak_map, best_peak_bins);

      art::FindManyP<recob::SpacePoint> spacepoint_per_hit(hitListHandle, evt, fSpacePointproducer);
      findShowerSpine3D(hits_from_reco, spacepoint_per_hit, unavailable_hit_list, shower_spine_hit_list_vector, end, view_peak_map, best_peak_bins);



      cout << "best_peak_bins.size(): " << best_peak_bins.size() << endl;
      cout << "shower_spine_hit_list_vector.size(): " << shower_spine_hit_list_vector.size() << endl;

      if(best_peak_bins.size()!=0){
	int i=0;
	for(auto const& entry : best_peak_bins){

	  cout << "entry.X(): " << entry.X() << endl;

	  if(shower_spine_hit_list_vector.size()==0) continue;
	  cout << "shower_spine_hit_list_vector[i]: " << shower_spine_hit_list_vector[i].size() << endl;
	  if(shower_spine_hit_list_vector[i].size()<10 || shower_spine_hit_list_vector[i].size() > 10000) continue;
	  cout << "shower_spine_hit_list_vector[i]: " << shower_spine_hit_list_vector[i].size() << endl;
	  //trkln = trackRebuid(shower_spine_hit_list_vector[i], spacepoint_per_hit, track);

	  recob::Track reco_track = trackRebuid(shower_spine_hit_list_vector[i], spacepoint_per_hit, track);

	  cout << "adding reco_track" << endl;
	  anaTrackCollection->push_back(reco_track);

	  std::vector<art::Ptr<recob::Hit>> hits_from_track_rebuild = shower_spine_hit_list_vector[i];

	  art::FindManyP<recob::SpacePoint> spacepoint_per_hit_reco(hitListHandle, evt, fSpacePointproducer);
	  std::vector<art::Ptr<recob::SpacePoint>> spacepoint_vec_reco;
	  unsigned int hitsFromSpacePointsRecoSize = 0;
	  for(size_t i_h=0; i_h<hits_from_track_rebuild.size(); i_h++){
	    spacepoint_vec_reco.clear();
	    spacepoint_vec_reco = spacepoint_per_hit_reco.at(hits_from_track_rebuild[i_h].key());
	    hitsFromSpacePointsRecoSize += spacepoint_vec_reco.size();
	  }

	  art::Ptr<recob::Track> pTrackdau(makeTrackPtr(anaTrackCollection->size() - 1));

	  lar_pandora::HitVector anaHitCollection_rebuild_tmp;
	  for(auto hitptr : hits_from_track_rebuild){
	    anaHitCollection_rebuild_tmp.push_back(hitptr);
	  }
	  util::CreateAssn(*this, evt, *(anaTrackCollection.get()), anaHitCollection_rebuild_tmp, *(anaTrackHitAssociations.get()));

	  for (unsigned int hitIndex = 0; hitIndex < hits_from_track_rebuild.size(); hitIndex++){
	    const art::Ptr<recob::Hit> phit(hits_from_track_rebuild.at(hitIndex));
	    const int index((hitIndex < hitsFromSpacePointsRecoSize) ? hitIndex : std::numeric_limits<int>::max());
	    recob::TrackHitMeta metadata(index, -std::numeric_limits<double>::max());
	  }
	  i++;
	}
      }

    }
    evt.put(std::move(anaTrackCollection));
    evt.put(std::move(anaTrackHitAssociations));
  }
}
