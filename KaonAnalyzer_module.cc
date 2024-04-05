//Module analyzer
//Ana module for nucleon decay and atmospheric analysis
//Ana TTree contains MC truth and reconstruction info
//ahiguera@central.uh.edu


//header file
#include "HitSplitAlg_module.h"
#include "HitSplitAlg.h"
#include "TrackRebuild.h"


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
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"

#ifdef __CLING__
#pragma link C++ class std::vector < std::vector<Float_t> >+;
#pragma link C++ class std::vector < std::vector< std::vector<Float_t> > >+;
#endif

using namespace std;

//========================================================================

namespace HitSplitAlg_module{

HitSplitAlg::HitSplitAlg(fhicl::ParameterSet const& parameterSet)
    : EDAnalyzer(parameterSet), fCalorimetryAlg(parameterSet.get< fhicl::ParameterSet >("CalorimetryAlg"))

{
    reconfigure(parameterSet);
}
//========================================================================
HitSplitAlg::~HitSplitAlg(){
  //destructor
}
//========================================================================
void HitSplitAlg::reconfigure(fhicl::ParameterSet const& p)
{

    fHitModuleLabel      = p.get<std::string>("HitModuleLabel");
    fHitModuleLabelOLD      = p.get<std::string>("HitModuleLabelOLD");
    fHitSPAssns          = p.get<std::string>("HitSPAssns");
    fTrackModuleLabel    = p.get<std::string>("TrackModuleLabel");
    fPidModuleLabel    = p.get<std::string>("PidModuleLabel");
    fCaloModuleLabel     = p.get<std::string>("CaloModuleLabel");
    fOpFlashModuleLabel  = p.get<std::string>("OpFlashModuleLabel");
    fShowerModuleLabel	 = p.get<std::string>("ShowerModuleLabel");
    fPointIDModuleLabel  = p.get<std::string>("PointIDModuleLabel");
    fNNetModuleLabel     = p.get<std::string>("NNetModuleLabel");
    fMCgenieLabel	 = p.get<std::string>("MCgenieLabel");
    fSpacePointproducer  = p.get<std::string>("SpacePointproducer");
    fHitTrackAssns = p.get<std::string>("HitTrackAssns");
    fSaveHits		 = p.get<bool>("SaveHits");
    fExponentConstant 	 = p.get<double>("ExponentConstant");
    fMaxPIDAValue 	 = p.get<double>("MaxPIDAValue");
    fMinPIDAValue 	 = p.get<double>("MinPIDAValue");
    fPIDA_endPoint	 = p.get<double>("PIDA_endPoint");
    fView                = p.get<double>("View");
    fPidValue		 = p.get<double>("PidValue");
    fFidVolCutX          = p.get<double>("FidVolCutX");
    fFidVolCutY          = p.get<double>("FidVolCutY");
    fFidVolCutZ          = p.get<double>("FidVolCutZ");
}
//========================================================================
void HitSplitAlg::beginJob(){
  cout<<"job begin..."<<endl;
  // Get geometry.
  auto const* geo = lar::providerFrom<geo::Geometry>();
  //art::ServiceHandle<geo::Geometry> geo;
  // Define histogram boundaries (cm).
  // For now only draw cryostat=0.
  double minx = 1e9;
  double maxx = -1e9;
  double miny = 1e9;
  double maxy = -1e9;
  double minz = 1e9;
  double maxz = -1e9;
  //for (size_t i = 0; i<geo->NTPC(); ++i){
  for (geo::TPCID const& TPCID: geo->Iterate<geo::TPCID>()) { 
    //double local[3] = {0.,0.,0.};
    double world[3] = {0.,0.,0.};
    //const geo::TPCGeo &tpc = geo->TPC(TPCID);
    if (minx>world[0]-geo->DetHalfWidth(TPCID))
      minx = world[0]-geo->DetHalfWidth(TPCID);
    if (maxx<world[0]+geo->DetHalfWidth(TPCID))
      maxx = world[0]+geo->DetHalfWidth(TPCID);
    if (miny>world[1]-geo->DetHalfHeight(TPCID))
      miny = world[1]-geo->DetHalfHeight(TPCID);
    if (maxy<world[1]+geo->DetHalfHeight(TPCID))
      maxy = world[1]+geo->DetHalfHeight(TPCID);
    if (minz>world[2]-geo->DetLength(TPCID)/2.)
      minz = world[2]-geo->DetLength(TPCID)/2.;
    if (maxz<world[2]+geo->DetLength(TPCID)/2.)
      maxz = world[2]+geo->DetLength(TPCID)/2.;
 
    /*
    tpc.LocalToWorld(local,world);
    if (minx>world[0]-geo->DetHalfWidth(i))
      minx = world[0]-geo->DetHalfWidth(i);
    if (maxx<world[0]+geo->DetHalfWidth(i))
      maxx = world[0]+geo->DetHalfWidth(i);
    if (miny>world[1]-geo->DetHalfHeight(i))
      miny = world[1]-geo->DetHalfHeight(i);
    if (maxy<world[1]+geo->DetHalfHeight(i))
      maxy = world[1]+geo->DetHalfHeight(i);
    if (minz>world[2]-geo->DetLength(i)/2.)
      minz = world[2]-geo->DetLength(i)/2.;
    if (maxz<world[2]+geo->DetLength(i)/2.)
      maxz = world[2]+geo->DetLength(i)/2.;
    */
  }

  //if (std::abs(nuvtxx_truth)<360-50&&
  //  std::abs(nuvtxy_truth)<600-50&&
  //  nuvtxz_truth>50&&nuvtxz_truth<1394-150){
  minx = -360.;
  maxx = +360.;
  miny = -600.;
  maxy = +600.;
  maxz = 0.;
  maxz = 1394.;

  fFidVolXmin = minx + fFidVolCutX;
  fFidVolXmax = maxx - fFidVolCutX;
  fFidVolYmin = miny + fFidVolCutY;
  fFidVolYmax = maxy - fFidVolCutY;
  fFidVolZmin = minz + fFidVolCutZ;
  fFidVolZmax = maxz - fFidVolCutZ;

  /*
  cout<<"Fiducial volume:"<<"\n"
	   <<fFidVolXmin<<"\t< x <\t"<<fFidVolXmax<<"\n"
	   <<fFidVolYmin<<"\t< y <\t"<<fFidVolYmax<<"\n"
	   <<fFidVolZmin<<"\t< z <\t"<<fFidVolZmax<<"\n";
  cout << "THIS IS MODULE FOR HITSPLIT ALG" << endl;
  */

  art::ServiceHandle<art::TFileService> tfs;
  fEventTree = tfs->make<TTree>("Event", "Event Tree from Sim & Reco");

  fEventTree->Branch("eventNo", &Event);
  fEventTree->Branch("runNo", &Run);
  fEventTree->Branch("subRunNo", &SubRun);

  //GENIE info & list of particles
  fEventTree->Branch("MC_Ev", &MC_Ev);
  fEventTree->Branch("MC_cc", &MC_cc);
  fEventTree->Branch("MC_Q2", &MC_Q2);
  fEventTree->Branch("MC_nuPDG", &MC_nuPDG);
  fEventTree->Branch("MC_hit_nucleon", &MC_hit_nucleon);
  fEventTree->Branch("mcgenie_npart", &MCgenie_npart);  // number of particles
  fEventTree->Branch("mcgenie_id", MCgenie_id, "mcgenie_id[mcgenie_npart]/I");
  fEventTree->Branch("mcgenie_fate", MCgenie_fate, "mcgenie_fate[mcgenie_npart]/I");
  fEventTree->Branch("mcgenie_statusCode", MCgenie_statusCode, "mcgenie_statusCode[mcgenie_npart]/I");
  fEventTree->Branch("mcgenie_pdg", MCgenie_pdg, "mcgenie_pdg[mcgenie_npart]/I");
  fEventTree->Branch("mcgenie_mother", MCgenie_mother, "mcgenie_mother[mcgenie_npart]/I");
  fEventTree->Branch("mcgenie_startMomentum", MCgenie_startMomentum, "mcgenie_startMomentum[mcgenie_npart][4]/D");
  fEventTree->Branch("mcgenie_endMomentum", MCgenie_endMomentum, "mcgenie_endMomentum[mcgenie_npart][4]/D");

  //Geant4 list of particles
  fEventTree->Branch("mc_vertex", MC_vertex, "mc_vertex[4]/D");
  fEventTree->Branch("mc_npart", &MC_npart);  // number of particles
  fEventTree->Branch("mc_id", MC_id, "mc_id[mc_npart]/I");
  fEventTree->Branch("mc_pdg", MC_pdg, "mc_pdg[mc_npart]/I");
  fEventTree->Branch("mc_statusCode", MC_statusCode, "mc_statusCode[mc_npart]/I");
  fEventTree->Branch("mc_mother", MC_mother, "mc_mother[mc_npart]/I");
  fEventTree->Branch("mc_startXYZT", MC_startXYZT, "mc_startXYZT[mc_npart][4]/D");
  fEventTree->Branch("mc_endXYZT", MC_endXYZT, "mc_endXYZT[mc_npart][4]/D");
  fEventTree->Branch("mc_npoints", MC_npoints, "mc_npoints[mc_npart]/I");
  fEventTree->Branch("mc_xyz", MC_xyz, "mc_xyz[mc_npart][100][3]/D");
  fEventTree->Branch("mc_startMomentum", MC_startMomentum, "mc_startMomentum[mc_npart][4]/D");
  fEventTree->Branch("mc_endMomentum", MC_endMomentum, "mc_endMomentum[mc_npart][4]/D");
  fEventTree->Branch("mc_Prange", MC_Prange, "mc_Prange[mc_npart]/D");
  fEventTree->Branch("mc_truthlength", MC_truthlength, "mc_truthlength[mc_npart]/D");
  fEventTree->Branch("mc_kaontruthlength", &MC_kaontruthlength);
  fEventTree->Branch("mc_dautruthlength", &MC_dautruthlength);
  fEventTree->Branch("kaon_br", &Kaon_BR);
  fEventTree->Branch("mc_process", &MC_process);
  fEventTree->Branch("mc_Endprocess", &MC_Endprocess);

  fEventTree->Branch("n_vertices", &n_vertices);
  fEventTree->Branch("vertex", vertex,"vertex[n_vertices][4]/D");
  fEventTree->Branch("vtx_ID", vtx_ID,"vtx_ID[n_vertices]/I");
  fEventTree->Branch("n_reco_tracks", &n_recoTracks);
  fEventTree->Branch("n_reco_dautracks", n_recoDauTracks, "n_recoDauTracks[n_reco_tracks]/D");
  fEventTree->Branch("n_reco_rebdautracks", n_recoRebDauTracks, "n_recoRebDauTracks[n_reco_tracks]/D");
  fEventTree->Branch("n_decayVtx", &n_decayVtx);
  fEventTree->Branch("decayVtx", decayVtx,"decayVtx[n_decayVtx][3]/D");  //vertices found using decayID point Alg
  fEventTree->Branch("vtxID_trk", vtxID_trk,"vtxID_trk[n_reco_tracks][10]/I"); //track-vertex association
  fEventTree->Branch("track_vtx", track_vtx,"track_vtx[n_reco_tracks][4]/D");
  fEventTree->Branch("track_vtxDir", track_vtxDir,"track_vtxDir[n_reco_tracks][3]/D");
  fEventTree->Branch("track_end", track_end,"track_end[n_reco_tracks][4]/D");
  fEventTree->Branch("track_isContained", track_isContained,"track_isContained[n_reco_tracks]/I");
  fEventTree->Branch("track_ID", track_ID,"track_ID[n_reco_tracks]/I");
  fEventTree->Branch("track_length", track_length,"track_length[n_reco_tracks]/D");
  fEventTree->Branch("dautrack_length", dautrack_length,"dautrack_length[n_reco_tracks][10]/D");
  fEventTree->Branch("rebdautracktrue_length", rebdautracktrue_length,"rebdautracktrue_length[n_reco_tracks]/D");
  fEventTree->Branch("rebdautracktruedir_length", rebdautracktruedir_length,"rebdautracktruedir_length[n_reco_tracks]/D");
  fEventTree->Branch("rebdautrack_length", rebdautrack_length,"rebdautrack_length[n_reco_tracks][10]/D");
  fEventTree->Branch("dautrack_pdg", dautrack_pdg,"dautrack_pdg[n_reco_tracks][10]/D");
  fEventTree->Branch("rebdautrack_pdg", rebdautrack_pdg,"rebdautrack_pdg[n_reco_tracks][10]/D");
  fEventTree->Branch("best_peak_theta", best_peak_theta,"best_peak_theta[n_reco_tracks][10]/D");
  fEventTree->Branch("best_peak_phi", best_peak_phi,"best_peak_phi[n_reco_tracks][10]/D");
  fEventTree->Branch("best_peak_theta_true", best_peak_theta_true,"best_peak_theta_true[n_reco_tracks]/D");
  fEventTree->Branch("best_peak_phi_true", best_peak_phi_true,"best_peak_phi_true[n_reco_tracks]/D");
  fEventTree->Branch("track_PIDA", track_PIDA,"track_PIDA[n_reco_tracks][3]/D");
  fEventTree->Branch("track_PID_pdg", track_PID_pdg,"track_PID_pdg[n_reco_tracks][3]/I");
  fEventTree->Branch("track_KE", track_KE,"track_KE[n_reco_tracks][3]/D");
  fEventTree->Branch("track_Prange", track_Prange,"track_Prange[n_reco_tracks]/D");
  fEventTree->Branch("track_bestplane", track_bestplane,"track_bestplane[n_reco_tracks]/I");
  fEventTree->Branch("n_track_points", n_track_points,"n_track_points[n_reco_tracks]/I");
  fEventTree->Branch("track_point_xyz", track_point_xyz,Form("track_point_xyz[n_reco_tracks][%d][3]/D",MAX_CALO_PTS));

  fEventTree->Branch("n_cal_points", n_cal_points,"n_cal_points[n_reco_tracks]/I");
  fEventTree->Branch("track_dQ_dx", track_dQ_dx,"track_dQ_dx[n_reco_tracks][500]/D");
  fEventTree->Branch("track_dE_dx", track_dE_dx,"track_dE_dx[n_reco_tracks][500]/D");
  fEventTree->Branch("track_range", track_range,"track_range[n_reco_tracks][500]/D");
  fEventTree->Branch("track_pitch", track_pitch,"track_pitch[n_reco_tracks][500]/D");

  fEventTree->Branch("n_cal_points_byplane", n_cal_points_byplane,"n_cal_points_byplane[n_reco_tracks][3]/I");
  fEventTree->Branch("track_dQ_dx_byplane", track_dQ_dx_byplane,"track_dQ_dx_byplane[n_reco_tracks][500][3]/D");
  fEventTree->Branch("track_dE_dx_byplane", track_dE_dx_byplane,"track_dE_dx_byplane[n_reco_tracks][500][3]/D");
  fEventTree->Branch("track_range_byplane", track_range_byplane,"track_range_byplane[n_reco_tracks][500][3]/D");
  fEventTree->Branch("track_pitch_byplane", track_pitch_byplane,"track_pitch_byplane[n_reco_tracks][500][3]/D");
  fEventTree->Branch("track_calo_xyz_byplane", track_calo_xyz_byplane,"track_calo_xyz_byplane[n_reco_tracks][3][500][3]/D");

  fEventTree->Branch("track_complet", track_complet,"track_complet[n_reco_tracks]/D");  //track quality variable (completeness)
  fEventTree->Branch("track_Efrac", track_Efrac,"track_Efrac[n_reco_tracks]/D");        //track quality variable (purity)
  fEventTree->Branch("track_mcID", track_mcID,"track_mcID[n_reco_tracks]/I");           //true MC ID for a given track
  fEventTree->Branch("track_mcPDG", track_mcPDG,"track_mcPDG[n_reco_tracks]/I");        //true MC PDG for a given track
  fEventTree->Branch("dautrack_mcPDG", dautrack_mcPDG,"dautrack_mcPDG[n_reco_tracks][10]/I");        //true MC PDG for a given daughter track

  fEventTree->Branch("Em_ch", &Em_ch);
  fEventTree->Branch("Em_e", &Em_e);
  fEventTree->Branch("trk_e", &trk_e);
  fEventTree->Branch("Emichel_e", &Emichel_e);
  if (fSaveHits){
  fEventTree->Branch("n_recoHits", &n_recoHits);
  fEventTree->Branch("hit_channel", &hit_channel, "hit_channel[n_recoHits]/I");
  fEventTree->Branch("hit_tpc", &hit_tpc, "hit_tpc[n_recoHits]/I");
  fEventTree->Branch("hit_plane", &hit_plane, "hit_plane[n_recoHits]/I");
  fEventTree->Branch("hit_wire", &hit_wire, "hit_wire[n_recoHits]/I");
  fEventTree->Branch("hit_peakT", &hit_peakT, "hit_peakT[n_recoHits]/D");
  fEventTree->Branch("hit_charge", &hit_charge, "hit_charge[n_recoHits]/D");
  fEventTree->Branch("hit_ph", &hit_ph, "hit_ph[n_recoHits]/D");
  fEventTree->Branch("hit_endT", &hit_endT, "hit_endT[n_recoHits]/D");
  fEventTree->Branch("hit_rms", &hit_rms, "hit_rms[n_recoHits]/D");
  fEventTree->Branch("hit_electrons", &hit_electrons, "hit_electrons[n_recoHits]/D");
  }

  fEventTree->Branch("n_showers", &n_recoShowers);
  fEventTree->Branch("sh_direction_X", &sh_direction_X, "sh_direction_X[n_showers]/D");
  fEventTree->Branch("sh_direction_Y", &sh_direction_Y, "sh_direction_Y[n_showers]/D");
  fEventTree->Branch("sh_direction_Z", &sh_direction_Z, "sh_direction_Z[n_showers]/D");
  fEventTree->Branch("sh_start_X", &sh_start_X, "sh_start_X[n_showers]/D");
  fEventTree->Branch("sh_start_Y", &sh_start_Y, "sh_start_Y[n_showers]/D");
  fEventTree->Branch("sh_start_Z", &sh_start_Z, "sh_start_Z[n_showers]/D");
  fEventTree->Branch("sh_energy", &sh_energy, "sh_energy[n_showers][3]/D");
  fEventTree->Branch("sh_MIPenergy", &sh_MIPenergy, "sh_MIPenergy[n_showers][3]/D");
  fEventTree->Branch("sh_dEdx", &sh_dEdx, "sh_dEdx[n_showers][3]/D");
  fEventTree->Branch("sh_bestplane", &sh_bestplane, "sh_bestplane[n_showers]/I");
  fEventTree->Branch("sh_length", &sh_length, "sh_length[n_showers]/D");


  fEventTree->Branch("n_flashes", &n_flashes);
  fEventTree->Branch("flash_time", &flash_time,"flash_time[n_flashes]/D");
  fEventTree->Branch("flash_pe", &flash_pe,"flash_pe[n_flashes]/D");
  fEventTree->Branch("flash_ycenter", &flash_ycenter,"flash_ycenter[n_flashes]/D");
  fEventTree->Branch("flash_zcenter", &flash_zcenter,"flash_zcenter[n_flashes]/D");
  fEventTree->Branch("flash_ywidth", &flash_ywidth,"flash_ywidth[n_flashes]/D");
  fEventTree->Branch("flash_zwidth", &flash_zwidth,"flash_zwidth[n_flashes]/D");
  fEventTree->Branch("flash_timewidth", &flash_timewidth, "flash_timewidth[n_flashes]/D");
  fEventTree->Branch("flash_abstime", &flash_abstime, "flash_abstime[n_flashes]/D");
  fEventTree->Branch("flash_frame", &flash_frame, "flash_frame[n_flashes]/I");

  fEventTree->Branch("flash_PE_ndk", &flash_PE_ndk, "flash_PE_ndk[n_flashes]/D");
  fEventTree->Branch("flash_PE_Ar39", &flash_PE_Ar39, "flash_PE_Ar39[n_flashes]/D");



}
//========================================================================
void HitSplitAlg::endJob(){
}
//========================================================================
void HitSplitAlg::beginRun(const art::Run& /*run*/){
  mf::LogInfo("HitSplitAlg")<<"begin run..."<<endl;
}
//========================================================================
void HitSplitAlg::analyze( const art::Event& event ){
    if (event.isRealData()) return;
    reset(); // reset data to be stored in the tree

    Event  = event.id().event();
    Run    = event.run();
    SubRun = event.subRun();
    bool isFiducial = false;
    Process(event, isFiducial);
    //if(isFiducial) fEventTree->Fill();
    fEventTree->Fill();

}
//========================================================================
void HitSplitAlg::Process( const art::Event& event, bool &isFiducial){

    //save GENIE stuf for atmos
    art::Handle<std::vector<simb::MCTruth>> MCtruthHandle;
    std::vector<art::Ptr<simb::MCTruth>> MCtruthlist;
    if(event.getByLabel(fMCgenieLabel, MCtruthHandle))
      art::fill_ptr_vector(MCtruthlist, MCtruthHandle);

    //For now assume that there is only one neutrino interaction...
    art::Ptr<simb::MCTruth> MCtruth;

    if( MCtruthlist.size()>0 ){
      MCtruth = MCtruthlist[0];
      if( MCtruth->NeutrinoSet() ){
        simb::MCNeutrino nu = MCtruth->GetNeutrino();
        if( nu.CCNC() == 0 ) MC_cc = 1;
        else if ( nu.CCNC() == 1 ) MC_cc = 0;
        simb::MCParticle neutrino = nu.Nu();
        MC_nuPDG = nu.Nu().PdgCode();
        const TLorentzVector& nu_momentum = nu.Nu().Momentum(0);
        double MC_incoming_P[4];
        nu_momentum.GetXYZT(MC_incoming_P);
        MC_Ev = nu_momentum[3];
        MC_Q2 = nu.QSqr();
        MC_hit_nucleon = nu.HitNuc();
      }
      MCgenie_npart = MCtruth->NParticles();
      if (MCgenie_npart > MAX_GEN_PARTS) MCgenie_npart = MAX_GEN_PARTS;
      for( int i =0; i<MCgenie_npart; ++i ){
         simb::MCParticle particle = MCtruth->GetParticle(i);
         MCgenie_id[i] = particle.TrackId();
         MCgenie_pdg[i] = particle.PdgCode();
         MCgenie_mother[i] = particle.Mother();
         MCgenie_statusCode[i] =particle.StatusCode();
         MCgenie_fate[i] = particle.Rescatter();
         const TLorentzVector& momentumStart = particle.Momentum(0);
         const TLorentzVector& momentumEnd   = particle.EndMomentum();
         //!Save the true vertex as the vertex using primaries ...hmmm do you have another suggestion?
         momentumStart.GetXYZT(MCgenie_startMomentum[i]);
         momentumEnd.GetXYZT(MCgenie_endMomentum[i]);

      } 
    }

    //art::ServiceHandle<cheat::BackTrackerService> bt;
    art::ServiceHandle<cheat::ParticleInventoryService> part_inv;
    const sim::ParticleList& plist = part_inv->ParticleList();
    simb::MCParticle *particle=0;
    int i=0; // particle index
    MC_npart = plist.size();

    auto particleHandle = event.getValidHandle<vector<simb::MCParticle>> ("largeant");
    map< int, const simb::MCParticle* > particleMap;
    int fSimTrackID =0.; 
    for(auto const& particle : (*particleHandle))
      {
	fSimTrackID = particle.TrackId();
	particleMap[fSimTrackID] = &particle;
      }
    
    if( MC_npart > MAX_MC_PARTS )
	return;
    for( sim::ParticleList::const_iterator ipar = plist.begin(); ipar!=plist.end(); ++ipar){
       particle = ipar->second;
       MC_id[i] = particle->TrackId();
       MC_pdg[i] = particle->PdgCode();
       MC_mother[i] = particle->Mother();
       MC_process.push_back(particle->Process());
       MC_Endprocess.push_back(particle->EndProcess());
       MC_statusCode[i] = particle->StatusCode();
       const TLorentzVector& positionStart = particle->Position(0);
       const TLorentzVector& positionEnd   = particle->EndPosition();
       const TLorentzVector& momentumStart = particle->Momentum(0);
       const TLorentzVector& momentumEnd   = particle->EndMomentum();
       //!Save the true vertex as the vertex using primaries ...hmmm do you have another suggestion?
       if( particle->Mother() == 0 ) positionStart.GetXYZT(MC_vertex);
       positionStart.GetXYZT(MC_startXYZT[i]);
       positionEnd.GetXYZT(MC_endXYZT[i]);
       momentumStart.GetXYZT(MC_startMomentum[i]);
       momentumEnd.GetXYZT(MC_endMomentum[i]);
       MC_truthlength[i] = truthLength(particle);

       if (MC_pdg[i] == 321) {
	 MC_kaontruthlength = MC_truthlength[i];

	 for (int iDaughter=0; iDaughter < particle->NumberDaughters(); iDaughter++){
	   auto d_search = particleMap.find(particle->Daughter(iDaughter)); 
	   auto const & daughter = *((*d_search).second);

	   if(daughter.PdgCode()==-13){
	     Kaon_BR=1;
	     MC_dautruthlength = truthLength(&daughter);
	   }
	   if(daughter.PdgCode()==211){
	     Kaon_BR=2;
	     MC_dautruthlength = truthLength(&daughter);
	   }
	   else if(Kaon_BR<0) Kaon_BR=0;
	 }
       }

       // copy all trajectory point positions
       int tot_pts = 0;
       double* trj_xyz = (double*) MC_xyz[i];
       for (auto ipt = particle->Trajectory().begin(); ipt != particle->Trajectory().end(); ipt++) {
	 if (tot_pts >= MAX_MC_TRAJ_PTS)
	   break;
	 ipt->first.Vect().GetXYZ(trj_xyz+tot_pts*3);
	 tot_pts++;
       }
       MC_npoints[i] = tot_pts;

       //MC_Prange[i] = myPrange(MC_truthlength[i]);
       ++i; //paticle index

    }



    MC_process.resize(i);
    MC_Endprocess.resize(i);
    //cout << "insideFV: " << insideFV(MC_vertex) << endl;
    isFiducial = insideFV( MC_vertex );
    if( !isFiducial ) return;
  
    //========================================================================
    //========================================================================
    // Reco  stuff
    //========================================================================
    //========================================================================

    //this is a track base analysis so it must be at least one track
    // get tracks
    art::Handle< std::vector<recob::Track> > trackListHandle;
    std::vector<art::Ptr<recob::Track>> tracklist;
    if( event.getByLabel(fTrackModuleLabel, trackListHandle))
      art::fill_ptr_vector(tracklist, trackListHandle);
    n_recoTracks = tracklist.size();
    if( n_recoTracks > MAX_TRACKS || n_recoTracks == 0) {
	n_recoTracks = 0; // make sure we don't save any fake data
	return;
    }

    // get vtx related to track?
    art::FindManyP<recob::Vertex> trk_from_vtx(trackListHandle,event,"pmtrack");

    art::Handle< std::vector<recob::Vertex> > vtxListHandle;
    std::vector<art::Ptr<recob::Vertex>> vtxlist;
    if(event.getByLabel("pmtrack", vtxListHandle))
      art::fill_ptr_vector(vtxlist, vtxListHandle);

    n_vertices = vtxlist.size();
    if( n_vertices != 0 )
    for( int i =0; i<n_vertices; ++i){
       double tmp_vtx[3] ={-999.0,-999.0,-999.0};
       vtx_ID[i] = vtxlist[i]->ID();
       vtxlist[i]->XYZ(tmp_vtx);
       for( int j=0; j<3; ++j) vertex[i][j]=tmp_vtx[j];
    }

    // get hits and calo of tracks
    art::FindManyP<recob::Hit> track_hits(trackListHandle, event, fHitTrackAssns);

    art::FindMany<anab::Calorimetry>  reco_cal(trackListHandle, event, fCaloModuleLabel);
    trkf::TrackMomentumCalculator trackP;

    art::FindMany<anab::ParticleID> reco_PID(trackListHandle, event, fPidModuleLabel);

    fHits.clear();
    fSpacePoints.clear();
    fTracks.clear();
    fShowers.clear();
    fSpacePointsToHits.clear();
    fSpacePointsToHits_old.clear();
    fHitsToSpacePoints.clear();
    fHitsToSpacePoints_old.clear();
    fTracksToHits.clear();
    fTracksToSpacePoints.clear();
    fShowersToHits.clear(); 
    fShowersToSpacePoints.clear();

    art::Handle<std::vector<recob::Hit>> hitsHandle;
    event.getByLabel(fHitModuleLabel, hitsHandle);

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

    //art::FindOneP<recob::Hit> findSPToHits(fSpacePoints, event, fSpacePointproducer); 
    art::FindOneP<recob::Hit> findSPToHits(fSpacePoints, event, fSpacePointproducer); 
    art::FindManyP<recob::Hit> findTracksToHits(fTracks, event, fTrackModuleLabel);
    art::FindManyP<recob::Hit> findShowersToHits(fShowers, event, fShowerModuleLabel);
    //art::FindManyP<recob::SpacePoint> findHitToSPs_old(fHits, event, fHitModuleLabelOLD);
    //art::FindOneP<recob::SpacePoint> findHitToSPs(fHits, event, "pandora");

    //art::FindManyP<recob::SpacePoint> findHitToSPs(fHits, event, fHitSPAssns);


    auto hit_handle = event.getValidHandle<std::vector<recob::Hit>>("hitfd");
    std::vector<art::Ptr<recob::Hit>> hit_list;
    art::fill_ptr_vector(hit_list, hit_handle);
    // auto hit_to_sp = art::FindManyP<recob::SpacePoint>(hit_handle, event, fHitToSpacePointLabel);
    auto findHitToSPs = art::FindManyP<recob::SpacePoint>(hit_handle, event, fHitSPAssns); 

    for (unsigned int iSP = 0; iSP < fSpacePoints.size(); ++iSP) { 
      const art::Ptr<recob::SpacePoint> spacePoint = fSpacePoints.at(iSP);
      const art::Ptr<recob::Hit> hit = findSPToHits.at(iSP); 
      fSpacePointsToHits_old[spacePoint] = hit;
      fHitsToSpacePoints_old[hit] = spacePoint; 
    }
    
    /*
    for (unsigned int iHit = 0; iHit < fHits.size(); iHit++){
      const art::Ptr<recob::Hit> hit = fHits.at(iHit);
      const std::vector<art::Ptr<recob::SpacePoint> > spacePoints = findHitToSPs_old.at(iHit);

      if(!spacePoints.empty()){
	//for (unsigned int iSP = 0; iSP < spacePoints.size(); iSP++){ 
	  //const art::Ptr<recob::SpacePoint> spacePoint = spacePoints.at(iSP);
	  //fSpacePointsToHits[spacePoint] = hit;
	//}
	if(!fHitsToSpacePoints_old.count(hit))
	  fHitsToSpacePoints_old[hit] = spacePoints[0];
      }
    }
    */

    for (unsigned int iHit = 0; iHit < hit_list.size(); iHit++) {
      //const art::Ptr<recob::SpacePoint> spacePoints = find
      //const art::Ptr<recob::SpacePoint> spacePoints = findHitToSPs.at(iHit);
      
      const art::Ptr<recob::Hit> hit = hit_list.at(iHit);
      const std::vector<art::Ptr<recob::SpacePoint> > spacePoints = findHitToSPs.at(iHit);
      //auto const hit = hit_list.at(iHit);
      //auto spacePoints = findHitToSPs.at(iHit);

      //cout << spacePoints << endl;
      /*
      if(!fHitsToSpacePoints.count(hit)){
	  cout << "adding this hit" << endl;
	  cout << spacePoints->XYZ()[0] << endl;
	  fHitsToSpacePoints[hit] = spacePoints;
      }
      */
	
      if(!spacePoints.empty()){

	for (unsigned int iSP = 0; iSP < spacePoints.size(); iSP++){
	  const art::Ptr<recob::SpacePoint> spacePoint = spacePoints.at(iSP);
	  fSpacePointsToHits[spacePoint] = hit;
	}
	
	if(!fHitsToSpacePoints.count(hit)){
	  //cout << "adding this hit" << endl;
	  //cout << spacePoints[0]->XYZ()[0] << endl;
	  fHitsToSpacePoints[hit] = spacePoints[0];
	}
	//else if(fHitsToSpacePoints.count(hit)>0) cout << "this hit is already stored in the map" << endl;
      }
      
    }

    for (unsigned int iTrack = 0; iTrack < fTracks.size(); ++iTrack) {
      const art::Ptr<recob::Track> track = fTracks.at(iTrack); 
      const std::vector<art::Ptr<recob::Hit>> trackHits = findTracksToHits.at(iTrack);
      cout << "findTracksToHits.at(iTrack): " << findTracksToHits.at(iTrack).size() << endl;

      for (unsigned int iHit = 0; iHit < trackHits.size(); ++iHit) {
	const art::Ptr<recob::Hit> hit = trackHits.at(iHit);
	fTracksToHits[track].push_back(hit);
	////if (fHitsToSpacePoints.count(hit)) {
	//cout << "fHitsToSpacePoints.count(hit): " << fHitsToSpacePoints_old.count(hit) << endl;
	//fTracksToSpacePoints[track].push_back(fHitsToSpacePoints.at(hit));
	if (fHitsToSpacePoints_old.count(hit)) {
	  fTracksToSpacePoints[track].push_back(fHitsToSpacePoints_old.at(hit));
	  //cout << fHitsToSpacePoints_old.at(hit)->XYZ()[0] << endl;
	}
      }
      cout << "fTracksToSpacePoints[track].size(): " << fTracksToSpacePoints[track].size() << endl;
    }

    for (unsigned int iShower = 0; iShower < fShowers.size(); ++iShower) {
      const art::Ptr<recob::Shower> shower = fShowers.at(iShower);
      const std::vector<art::Ptr<recob::Hit>> showerHits = findShowersToHits.at(iShower);

      for (unsigned int iHit = 0; iHit < showerHits.size(); ++iHit) {
	const art::Ptr<recob::Hit> hit = showerHits.at(iHit);
	fShowersToHits[shower].push_back(hit);
	if (fHitsToSpacePoints_old.count(hit)) {
	  fShowersToSpacePoints[shower].push_back(fHitsToSpacePoints_old.at(hit));
	}
      }
      cout << "fShowersToSpacePoints[shower].size(): " << fShowersToSpacePoints[shower].size() << endl;
    }

    /*
    cout << "fTracksToSpacePoints.at(0).size(): " << fTracksToSpacePoints.at(fTracks.at(0)).size() << endl;
    const std::vector<art::Ptr<recob::SpacePoint>>& sp = fTracksToSpacePoints.at(fTracks.at(0));
    for (auto spIter = sp.begin(); spIter != sp.end(); ++spIter) {
      TVector3 point = (*spIter)->XYZ();
      cout << "point.X(): " << point.X() << endl;
    }
    */

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
    art::FindManyP<recob::Hit> shower_hits(showerHandle, event, fShowerModuleLabel);
    cout << "shower_hits " << shower_hits.size() << endl;

    cout << "n_recoShowers: " << n_recoShowers << ", n_recoTracks: " << n_recoTracks << endl;
    cout << "tracklist.size(): " << tracklist.size() << ", showerlist.size(): " << showerlist.size() << endl;
    for(int i=0; i<n_recoTracks; ++i){
      cout << "i: " << i << ", track_hits[i].size(): " << track_hits.at(i).size() << endl;
    }
    for(int i=0; i<n_recoShowers; ++i){
      cout << "i: " << i << ", shower_hits[i].size(): " << shower_hits.at(i).size() << endl;
    }

    // loop reco tracks
    for(int i=0; i<n_recoTracks; ++i) {
      //int nK_primary = 0;

      int ndau_tracks = 0;
      //int ndau_showers = 0;
      n_recoDauTracks[i] = 0;
      n_recoRebDauTracks[i] = 0;

      vector<bool> v_trk_flg_peak;
      vector<TVector2> best_peak_bins;
      vector<TVector2> best_peak_bins_truepi;
      vector<TVector2> best_peak_bins_truemu;
      std::map<double, TVector2, std::greater<>> view_peak_map;
      std::map<double, TVector2, std::greater<>> view_peak_map_truepi;
      std::map<double, TVector2, std::greater<>> view_peak_map_truemu;

      std::map<int, std::map<int, double>> angular_distribution_map_3D_track;
      std::map<int, std::map<int, double>> angular_distribution_map_3D_shower;
      std::map<int, std::map<int, double>> angular_distribution_map_3D_pfparticle;
      std::map<int, std::map<int, double>> angular_distribution_map_3D_truepi;
      std::map<int, std::map<int, double>> angular_distribution_map_3D_truemu;
      std::vector<art::Ptr<recob::Hit>> unavailable_hit_list;
      std::vector<art::Ptr<recob::Hit>> shower_spine_hit_list;
      std::vector<std::vector<art::Ptr<recob::Hit>>> shower_spine_hit_list_vector;
      std::vector<std::vector<art::Ptr<recob::Hit>>> shower_spine_hit_list_vector_truepi;
      std::vector<std::vector<art::Ptr<recob::Hit>>> shower_spine_hit_list_vector_truemu;
      std::vector<std::vector<art::Ptr<recob::Hit>>> shower_spine_hit_list_vector_truepidir;
      std::vector<std::vector<art::Ptr<recob::Hit>>> shower_spine_hit_list_vector_truemudir;
      std::vector<art::Ptr<recob::Hit>> hits_from_reco_obj;
      std::vector<art::Ptr<recob::SpacePoint>> sp_from_recoobj;
      std::vector<art::Ptr<recob::SpacePoint>> sp_from_truemu;
      std::vector<art::Ptr<recob::SpacePoint>> sp_from_truepi;
      std::vector<art::Ptr<recob::SpacePoint>> sp_from_allhits;
      std::vector<art::Ptr<recob::SpacePoint>> sp_from_allhits_old;

      std::map<int, std::map<int, std::map<int, double>>> angular_distribution_map_3D_cheat;
      std::map<int, TH2D*> h_angular_distribution_pfparticle_cheat_3D;
      //TCanvas * c = new TCanvas("c", "c", 800, 600);

       art::Ptr<recob::Track> track = tracklist[i];
       /*
       const std::vector<art::Ptr<recob::SpacePoint>>& sp = fTracksToSpacePoints.at(track);
       for (auto spIter = sp.begin(); spIter != sp.end(); ++spIter) {
	 TVector3 point = (*spIter)->XYZ();
	 cout << "X(): " << point.X() << endl;
       }
       */
       //vtx associations
       std::vector<art::Ptr<recob::Vertex>> vtxs = trk_from_vtx.at(i);
       for( size_t j=0; j<vtxs.size(); ++j){
          art::Ptr<recob::Vertex> vtx = vtxs[j];
          vtxID_trk[i][j] = vtx->ID();
       }
       track_length[i] = track->Length();
       const TVector3 tmp_track_vtx(track->Vertex().x(),
				    track->Vertex().y(),
				    track->Vertex().z());
       const TVector3 tmp_track_end(track->End().x(),
				    track->End().y(),
				    track->End().z());
       track_vtx[i][0] =tmp_track_vtx[0];
       track_vtx[i][1] =tmp_track_vtx[1];
       track_vtx[i][2] =tmp_track_vtx[2];
       track_vtx[i][3] = -999.0;

       track_end[i][0] =tmp_track_end[0];
       track_end[i][1] =tmp_track_end[1];
       track_end[i][2] =tmp_track_end[2];
       track_end[i][3] = -999.0;

       const TVector3 tmp_vtx_dir(track->VertexDirection().x(),
				  track->VertexDirection().y(),
				  track->VertexDirection().Z());
       track_vtxDir[i][0] = tmp_vtx_dir[0];
       track_vtxDir[i][1] = tmp_vtx_dir[1];
       track_vtxDir[i][2] = tmp_vtx_dir[2];

       /* import from uB code*/
       const recob::Track& atrack = *track;

       TVector3 end(track->End().x(),
		    track->End().y(),
		    track->End().z());
       
       
       // add some conditions to remove NON K+-like tracks?
       //make the list of candidate hits for daughter reco track 
       //for (unsigned int ihit = 0; ihit < fHits.size(); ++ihit){
       //if(!fHits.size()) continue;
    
       /*
       for (unsigned int iHit = 0; iHit < fHits.size(); ++iHit) {
	 const art::Ptr<recob::Hit> hit = fHits.at(iHit); 
	 //cout << "WireID().Wire: " << hit->WireID().Wire << endl;
	 if(!fHitsToSpacePoints.count(hit)) continue;
	 art::Ptr<recob::SpacePoint> sp = fHitsToSpacePoints.at(hit);
	 if(!isInsideROI(sp, end)) continue;
       }
*/
       
       cout << "fHitsToSpacePoints.size(): " << fHitsToSpacePoints.size() << endl;
       //for(auto const& hit : fHits){
       for(auto const& hit : hit_list){
	 //if(!fHitsToSpacePoints.count(hit)) continue;
	 if(fHitsToSpacePoints.find(hit) == fHitsToSpacePoints.end()) continue;
	 //cout << "passed !fHitsToSpacePoints.count(hit)" << endl;
	 art::Ptr<recob::SpacePoint> sp = fHitsToSpacePoints.at(hit);
	 //cout << "get sp" << endl;
	 //the hit is inside ROI
	 if(!isInsideROI(sp, end)) continue;
	 //cout << "passed ROI" << endl;
	 //remove the hits consisting the primary track

	 if( std::find(track_hits.at(i).begin(), track_hits.at(i).end(), hit) != track_hits.at(i).end()) continue;

	 //cout << sp->XYZ()[0] << ' '<< sp->XYZ()[1] << ' ' << sp->XYZ()[2] << endl; 
	 //cout << "passed primary track filter" << endl;
	 sp_from_allhits.push_back(sp);
       }
       cout << "sp_from_allhits.size(): " << sp_from_allhits.size() << endl;


       cout << '\n';

       for(auto & sp : fSpacePoints){ 
	 //check sp is stored inside SP-Hit map
	 if(fSpacePointsToHits_old.find(sp) == fSpacePointsToHits_old.end()) continue;

	 //the hit is inside ROI
	 if(!isInsideROI(sp, end)) continue;

	 //skip if sp is a part of primary track
	 std::vector<art::Ptr<recob::SpacePoint>> sp_from_primary = fTracksToSpacePoints[track];
	 if( std::find(sp_from_primary.begin(), sp_from_primary.end(), sp) != sp_from_primary.end()) continue;

	 //cout << sp->XYZ()[0] << ' '<< sp->XYZ()[1] << ' ' << sp->XYZ()[2] << endl; 

	 sp_from_allhits_old.push_back(sp);
       }
       cout << "sp_from_allhits_old.size(): " << sp_from_allhits_old.size() << endl;

       const simb::MCParticle *particletmp;
       //if(!sp_from_allhits.size()) continue;
       if(!sp_from_allhits_old.size()) continue;
       //for(auto const& sp : sp_from_allhits){
       for(auto const& sp : sp_from_allhits_old){
	 //art::Ptr<recob::Hit> hit = fSpacePointsToHits.at(sp);
	 art::Ptr<recob::Hit> hit = fSpacePointsToHits_old.at(sp);
	 truthHitMatcher(hit, particletmp);
	     if(!particletmp) continue;
	     //if(particle->PdgCode() == 321) ++nK_dautrk;
	     if(particletmp->PdgCode() == 211) sp_from_truepi.push_back(sp);
	     else if(particletmp->PdgCode() == -13) sp_from_truemu.push_back(sp);
       }       

       //fillAngularDistributionMap3D(sp_from_allhits, end, angular_distribution_map_3D_track); 
       fillAngularDistributionMap3D(sp_from_allhits_old, end, angular_distribution_map_3D_track); 


       // truth match       
       double tmpEfrac = 0;
       const simb::MCParticle *particle;
       double tmpComplet = 0;
       std::vector<art::Ptr<recob::Hit>> all_trackHits = track_hits.at(i);
       truthMatcher( all_hits,  all_trackHits, particle, tmpEfrac, tmpComplet );
       if(!particle) continue;
       track_mcID[i] = particle->TrackId();
       track_mcPDG[i] = particle->PdgCode();
       track_Efrac[i] = tmpEfrac;
       track_complet[i] = tmpComplet;
       if(track_mcPDG[i]!=321) continue;
       /*
       else{
	 //cout << "THIS IS K+ TRACK" << endl;

	 std::vector<art::Ptr<recob::SpacePoint>>& sp_from_track = fTracksToSpacePoints.at(fTracks.at(i));

	 const simb::MCParticle *particletmp;
	 if(!sp_from_track.size()) continue;
	 for(auto const& sp : sp_from_track){
	   art::Ptr<recob::Hit> hit = fSpacePointsToHits.at(sp);
	   truthHitMatcher(hit, particletmp);
	   if(!particletmp) continue;
	   //if(particletmp->PdgCode() == 321) ++nK_primary;
	   ++nK_primary;
	   //cout << "PDG of this hit consisting K+ track is " << particletmp->PdgCode() << endl;                                          
 	 }

       }
       */

       // loop over daughter track candidates
       //int nK_dautrk = 0;
	 cout << "Kaon_BR: " << Kaon_BR << ", MC_dautruthlength: " << MC_dautruthlength << endl;

       for(int j=0; j<n_recoTracks; ++j) {
	 //nK_dautrk = 0;

	 art::Ptr<recob::Track> dau_track = tracklist[j];
	 //const recob::Track& adau_track = *dau_track;

	 // exclude the track itself
	 //cout << "track->ID(): " << track->ID() << ", dau_track->ID(): " << dau_track->ID() << endl;
	 if(dau_track->ID() == track->ID()) continue;
	 
	 double track_dau_distance = TMath::Sqrt((track->End().x()-dau_track->Vertex().x())*(track->End().x()-dau_track->Vertex().x()) +
						 (track->End().y()-dau_track->Vertex().y())*(track->End().y()-dau_track->Vertex().y()) +
						 (track->End().z()-dau_track->Vertex().z())*(track->End().z()-dau_track->Vertex().z()));

	 //cout << "candidate of daughter track" << endl;

	 if(track_dau_distance<10){
	   dautrack_length[i][n_recoDauTracks[i]] = dau_track->Length();
	   n_recoDauTracks[i]++;
	   cout << "this event has daughter track" << endl;
	   cout << "i: " << i << ", j: " << j << ", dau_track->Length(): " << dau_track->Length() << endl;

	   //track_mcID[i] = particle->TrackId(); 
	   dautrack_mcPDG[i][n_recoDauTracks[i]] = track_mcPDG[j];
	   cout << "and its PDG is " << track_mcPDG[j] << endl;
	   //cout << "hits_from_dau_track: " << track_hits.at(j).size() << endl;
	   //cout << "sp_from_dau_track.size(): " << fTracksToSpacePoints.at(fTracks.at(j)).size() << endl;

	   const simb::MCParticle *mcparticle_dau;
	   std::map<int,int> hits_pdg_map_dau;
	   for(auto const& hit : track_hits.at(j)){
	     truthHitMatcher(hit, mcparticle_dau);
	     if(!mcparticle_dau) continue;
	     //cout << "mcparticle->PdgCode(): " << mcparticle->PdgCode() << endl;
	     if(hits_pdg_map_dau.find(mcparticle_dau->PdgCode()) == hits_pdg_map_dau.end()) hits_pdg_map_dau[mcparticle_dau->PdgCode()] = 1;
	     else hits_pdg_map_dau[mcparticle_dau->PdgCode()] += 1;
	   }

	   vector<pair<int, int>> v;
	   if(hits_pdg_map_dau.size()){
	     for (map<int, int>::iterator it = hits_pdg_map_dau.begin(); it != hits_pdg_map_dau.end(); it++)
	       v.push_back({ it->second, it->first });
	     
	     sort(v.rbegin(), v.rend());
	     dautrack_pdg[i][n_recoDauTracks[i]] = v[0].second;
	     cout << "dautrack_pdg[i][n_recoDauTracks[i]]: " << v[0].second << endl;
	   }

	 }

	 if(track_dau_distance<20){// consider as a daughter object
	   //get hits and spacepoints
	   std::vector<art::Ptr<recob::Hit>> hits_from_dau_track = track_hits.at(j);
	     //cout << "hits_from_dau_track: " << hits_from_dau_track.size() << endl;

	   std::vector<art::Ptr<recob::SpacePoint>>& sp_from_dau_track = fTracksToSpacePoints.at(fTracks.at(j));
	   cout << "sp_from_dau_track.size(): " << sp_from_dau_track.size() << endl;
	   //fillAngularDistributionMap3D(sp_from_dau_track, end, angular_distribution_map_3D_track);


	   /*
	   const simb::MCParticle *particle;
	   if(sp_from_dau_track.size()) cout << "take sps from candidate recoobj" << endl;
	   if(!sp_from_dau_track.size()) continue;
	   for(auto const& sp : sp_from_dau_track){
	     art::Ptr<recob::Hit> hit = fSpacePointsToHits.at(sp);
	     truthHitMatcher(hit, particle);
	     if(!particle) continue;
	     //if(particle->PdgCode() == 321) ++nK_dautrk;
	     //++nK_dautrk;
	     if(particle->PdgCode() == 211) sp_from_truepi.push_back(sp);
	     else if(particle->PdgCode() == -13) sp_from_truemu.push_back(sp);
	     //cout << "PDG of this daughter track hit is " << particle->PdgCode() << endl;
	   }
	   */

	   //art::FindManyP<recob::SpacePoint> spacepoint_per_hit(HitHandle, event, fSpacePointproducer);
	   //art::FindOneP<recob::SpacePoint> spacepoint_per_hit(HitHandle, event, fSpacePointproducer);
	   //cout << "spacepoint_per_hit.size() " << spacepoint_per_hit.size() << endl;
	   sp_from_recoobj.insert(sp_from_recoobj.end(), sp_from_dau_track.begin(), sp_from_dau_track.end());
	   //hits_from_reco_obj.insert(hits_from_reco_obj.end(), hits_from_dau_track.begin(), hits_from_dau_track.end());
	   //fillAngularDistributionMap3D(hits_from_dau_track, end, spacepoint_per_hit, angular_distribution_map_3D_track);
	   //fillAngularDistributionMap3D(hits_from_dau_track, end, spacepoint_per_track, angular_distribution_map_3D_track);
	   ndau_tracks++;
	 }

       }// end of track loop

       // loop over daughter shower candidates
       //int nK_daushw =0;

       /*
       for(int j=0; j<n_recoShowers; ++j) {
	 //nK_daushw = 0;

	 art::Ptr<recob::Shower> dau_shower = showerlist[j];
	 //const recob::Shower& adau_shower = *dau_shower;

	 double shower_dau_distance = TMath::Sqrt((track->End().x()-dau_shower->ShowerStart().X())*(track->End().x()-dau_shower->ShowerStart().X()) +
						  (track->End().y()-dau_shower->ShowerStart().Y())*(track->End().y()-dau_shower->ShowerStart().Y()) + 
						  (track->End().z()-dau_shower->ShowerStart().Z())*(track->End().z()-dau_shower->ShowerStart().Z()) );
	 if (shower_dau_distance<20){
	   if(ndau_showers>20) break;
	   //get hits and spacepoints
	   std::vector<art::Ptr<recob::Hit>> hits_from_dau_shower = shower_hits.at(j);

	   std::vector<art::Ptr<recob::SpacePoint>>& sp_from_dau_shower = fShowersToSpacePoints.at(fShowers.at(j));
	   //skip all showers
	   //fillAngularDistributionMap3D(sp_from_dau_shower, end, angular_distribution_map_3D_shower);

	   const simb::MCParticle *particle;
	   if(!sp_from_dau_shower.size()) continue;
	   for(auto const& sp : sp_from_dau_shower){
	     art::Ptr<recob::Hit> hit = fSpacePointsToHits.at(sp);
	     truthHitMatcher(hit, particle);
	     if(!particle) continue;
	     //if(particle->PdgCode() == 321) ++nK_daushw;
	     //++nK_daushw;
	     //skip
	     //if(particle->PdgCode() == 211) sp_from_truepi.push_back(sp);
	     //else if(particle->PdgCode() == -13) sp_from_truemu.push_back(sp);
	     //cout << "PDG of this daughter shower hit is " << particle->PdgCode() << endl;
	     //else cout << "PDG of this hit is " << particle->PdgCode() << endl;
	   }
	   //skip all showers
	   //sp_from_recoobj.insert(sp_from_recoobj.end(), sp_from_dau_shower.begin(), sp_from_dau_shower.end()); 

	   // art::FindManyP<recob::SpacePoint> spacepoint_per_hit(HitHandle, event, fSpacePointproducer);
	   //hits_from_reco_obj.insert(hits_from_reco_obj.end(), hits_from_dau_shower.begin(), hits_from_dau_shower.end());

	   //fillAngularDistributionMap3D(hits_from_dau_shower, end, spacepoint_per_hit, angular_distribution_map_3D_shower);
	   ndau_showers++;
	 }

       }// end of shower loop
       */
    
       /*
       cout << "track->ID(): " << track->ID() << endl;
       const simb::MCParticle *particletmp;
       for(auto const& sp : sp_from_recoobj){
	 art::Ptr<recob::Hit> hit = fSpacePointsToHits.at(sp);
	 truthHitMatcher(hit, particletmp);
	 if(!particletmp) continue;
	 cout << "PDG of this hit is " << particletmp->PdgCode() << endl;
       }
       */
       


       /*BEGIN for truth information*/

       cout << "sp_from_truepi.size(): " << sp_from_truepi.size() << endl;
       cout << "sp_from_truemu.size(): " << sp_from_truemu.size() << endl;

       fillAngularDistributionMap3D(sp_from_truepi, end, angular_distribution_map_3D_truepi);
       fillAngularDistributionMap3D(sp_from_truemu, end, angular_distribution_map_3D_truemu);

       //cout << "angular_distribution_map_3D_truepi.size(): " << angular_distribution_map_3D_truepi.size() << endl;
       //cout << "angular_distribution_map_3D_truemu.size(): " << angular_distribution_map_3D_truemu.size() << endl;

       smoothAngularDistributionMap3D(angular_distribution_map_3D_truepi);
       smoothAngularDistributionMap3D(angular_distribution_map_3D_truemu);

       //cout << "angular_distribution_map_3D_truepi.size(): " << angular_distribution_map_3D_truepi.size() << endl;
       //cout << "angular_distribution_map_3D_truemu.size(): " << angular_distribution_map_3D_truemu.size() << endl;

       obtainPeakVector3D(angular_distribution_map_3D_truepi, v_trk_flg_peak, view_peak_map_truepi, true);
       obtainPeakVector3D(angular_distribution_map_3D_truemu, v_trk_flg_peak, view_peak_map_truemu, false);

       findBestAngularPeak3D(angular_distribution_map_3D_truepi, view_peak_map_truepi, best_peak_bins_truepi);
       findBestAngularPeak3D(angular_distribution_map_3D_truemu, view_peak_map_truemu, best_peak_bins_truemu);

       cout << "calling findShowerSpine3D with view_peak_map_truepi" << endl;
       findShowerSpine3D(sp_from_truepi,  unavailable_hit_list, shower_spine_hit_list_vector_truepi, end, view_peak_map_truepi, best_peak_bins_truepi);
       findShowerSpine3D(sp_from_truemu,  unavailable_hit_list, shower_spine_hit_list_vector_truemu, end, view_peak_map_truemu, best_peak_bins_truemu);
       findShowerSpine3D(sp_from_recoobj,  unavailable_hit_list, shower_spine_hit_list_vector_truepidir, end, view_peak_map, best_peak_bins_truepi);
       findShowerSpine3D(sp_from_recoobj,  unavailable_hit_list, shower_spine_hit_list_vector_truemudir, end, view_peak_map, best_peak_bins_truemu);

       if(shower_spine_hit_list_vector_truepi.size()>0 && shower_spine_hit_list_vector_truepi[0].size()>0){
	 recob::Track reco_track_truepi = trackRebuid(shower_spine_hit_list_vector_truepi[0], atrack);
	 cout << "this event has true pi" << endl;
	 rebdautracktrue_length[i] = reco_track_truepi.Length();
       }

       if(shower_spine_hit_list_vector_truemu.size()>0 && shower_spine_hit_list_vector_truemu[0].size()>0){
	 recob::Track reco_track_truemu = trackRebuid(shower_spine_hit_list_vector_truemu[0], atrack);
	 cout << "this event has true mu" << endl; 
	 rebdautracktrue_length[i] = reco_track_truemu.Length();
       }

       if(shower_spine_hit_list_vector_truepidir.size()>0 && shower_spine_hit_list_vector_truepidir[0].size()>0){
	 recob::Track reco_track_truepidir = trackRebuid(shower_spine_hit_list_vector_truepidir[0], atrack);
	 rebdautracktruedir_length[i] = reco_track_truepidir.Length();
       }

       if(shower_spine_hit_list_vector_truemudir.size()>0 && shower_spine_hit_list_vector_truemudir[0].size()>0){
	 recob::Track reco_track_truemudir = trackRebuid(shower_spine_hit_list_vector_truemudir[0], atrack);
	 rebdautracktruedir_length[i] = reco_track_truemudir.Length();
       }

       /*
       if(shower_spine_hit_list_vector_truepi.size()>1 && shower_spine_hit_list_vector_truepi[1].size()>0){
	 recob::Track reco_track_truepi = trackRebuid(shower_spine_hit_list_vector_truepi[1], atrack);
	 cout << "this event has true pi 1" << endl;
       }

       if(shower_spine_hit_list_vector_truemu.size()>1 && shower_spine_hit_list_vector_truemu[1].size()>0){
	 recob::Track reco_track_truemu = trackRebuid(shower_spine_hit_list_vector_truemu[1], atrack);
	 cout << "this event has true mu 1" << endl; 
       }
       */


       /* all cheated */
       /*
       std::vector<art::Ptr<recob::Hit>> hit_from_truepi;
       std::vector<art::Ptr<recob::Hit>> hit_from_truemu;

       for(auto const &sp : sp_from_truepi){
	 if(fSpacePointsToHits.at(sp)) hit_from_truepi.push_back(fSpacePointsToHits.at(sp));
       }
       for(auto const &sp : sp_from_truemu){
	 if(fSpacePointsToHits.at(sp)) hit_from_truemu.push_back(fSpacePointsToHits.at(sp));
       }

       cout << "all cheated\n" << endl;
       if(hit_from_truepi.size()) recob::Track reco_track_cheatpi = trackRebuid(hit_from_truepi, atrack);
       if(hit_from_truemu.size()) recob::Track reco_track_cheatmu = trackRebuid(hit_from_truemu, atrack);
       */
       
       /*
       fillAngularDistributionMap3D(sp_from_truepi, end, angular_distribution_map_3D_truepi);
       fillAngularDistributionMap3D(sp_from_truemu, end, angular_distribution_map_3D_truemu);

       obtainPeakVector3D(angular_distribution_map_3D_truepi, v_trk_flg_peak, view_peak_map_truepi, true);
       obtainPeakVector3D(angular_distribution_map_3D_truemu, v_trk_flg_peak, view_peak_map_truemu, false);

       findBestAngularPeak3D(angular_distribution_map_3D_truepi, view_peak_map_truepi, best_peak_bins_truepi);
       findBestAngularPeak3D(angular_distribution_map_3D_truemu, view_peak_map_truemu, best_peak_bins_truemu);

       findShowerSpine3D(sp_from_truepi,  unavailable_hit_list, shower_spine_hit_list_vector_truepi, end, view_peak_map, best_peak_bins_truepi);
       findShowerSpine3D(sp_from_truemu,  unavailable_hit_list, shower_spine_hit_list_vector_truemu, end, view_peak_map, best_peak_bins_truemu);
       findShowerSpine3D(sp_from_recoobj,  unavailable_hit_list, shower_spine_hit_list_vector_truepidir, end, view_peak_map, best_peak_bins_truepi);
       findShowerSpine3D(sp_from_recoobj,  unavailable_hit_list, shower_spine_hit_list_vector_truemudir, end, view_peak_map, best_peak_bins_truemu);

       if(shower_spine_hit_list_vector_truepidir.size()>0 && shower_spine_hit_list_vector_truepidir[0].size()>0){
	 recob::Track reco_track_truepidir = trackRebuid(shower_spine_hit_list_vector_truepidir[0], atrack);
	 rebdautracktruedir_length[i] = reco_track_truepidir.Length();
       }
       if(shower_spine_hit_list_vector_truemudir.size()>0 && shower_spine_hit_list_vector_truemudir[0].size()>0){
	 recob::Track reco_track_truemudir = trackRebuid(shower_spine_hit_list_vector_truemudir[0], atrack);
	 rebdautracktruedir_length[i] = reco_track_truemudir.Length();
       }
       if(shower_spine_hit_list_vector_truepi.size()>0 && shower_spine_hit_list_vector_truepi[0].size()>0){
	 recob::Track reco_track_truepi = trackRebuid(shower_spine_hit_list_vector_truepi[0], atrack);
	 rebdautracktrue_length[i] = reco_track_truepi.Length();
       }
       if(shower_spine_hit_list_vector_truemu.size()>0 && shower_spine_hit_list_vector_truemu[0].size()>0){
	 recob::Track reco_track_truemu = trackRebuid(shower_spine_hit_list_vector_truemu[0], atrack);
	 rebdautracktrue_length[i] = reco_track_truemu.Length();
       }
       if(best_peak_bins_truepi.size()){
	 best_peak_theta_true[i] = best_peak_bins_truepi[0].X() * thetaBinSize; 
	 best_peak_phi_true[i] = best_peak_bins_truepi[0].Y() * thetaBinSize; 
       }
       if(best_peak_bins_truemu.size()){
	 best_peak_theta_true[i] = best_peak_bins_truemu[0].X() * thetaBinSize; 
	 best_peak_phi_true[i] = best_peak_bins_truemu[0].Y() * thetaBinSize; 
       }
       */
       /*END for truth information*/

       smoothAngularDistributionMap3D(angular_distribution_map_3D_track);
       //smoothAngularDistributionMap3D(angular_distribution_map_3D_shower);
       //accumulateAngularDistributionMap3D(angular_distribution_map_3D_track, angular_distribution_map_3D_shower, angular_distribution_map_3D_pfparticle);

       obtainPeakVector3D(angular_distribution_map_3D_track, v_trk_flg_peak, view_peak_map, true);
       //obtainPeakVector3D(angular_distribution_map_3D_shower, v_trk_flg_peak, view_peak_map, false);
       ////obtainPeakVector3D(angular_distribution_map_3D_pfparticle, v_trk_flg_peak, view_peak_map, true);

       findBestAngularPeak3D(angular_distribution_map_3D_track, view_peak_map, best_peak_bins);
       cout << "best_peak_bins.size(): " << best_peak_bins.size() << endl;
       //findBestAngularPeak3D(angular_distribution_map_3D_pfparticle, view_peak_map, best_peak_bins);
       //findShowerSpine3D(sp_from_recoobj,  unavailable_hit_list, shower_spine_hit_list_vector, end, view_peak_map, best_peak_bins);
       cout << "calling findShowerSpine3D with view_peak_map" << endl;  
       //findShowerSpine3D(sp_from_allhits,  unavailable_hit_list, shower_spine_hit_list_vector, end, view_peak_map, best_peak_bins);
       findShowerSpine3D(sp_from_allhits_old,  unavailable_hit_list, shower_spine_hit_list_vector, end, view_peak_map, best_peak_bins);
       cout << "n_recoDauTracks[i]: " << n_recoDauTracks[i] << endl;
       cout << "shower_spine_hit_list_vector.size(): " << shower_spine_hit_list_vector.size() << endl;


       //draw histos
       //cout << "call fillHistAngularDistributionMap3DCheat" << endl;
       fillHistAngularDistributionMap3DCheat(sp_from_allhits_old, end, angular_distribution_map_3D_cheat, h_angular_distribution_pfparticle_cheat_3D);
       //fillHistAngularDistributionMap3DCheat(sp_from_recoobj, end, angular_distribution_map_3D_cheat, h_angular_distribution_pfparticle_cheat_3D);
       //cout << "call drawHistAngularDistributionMap3DCheat" << endl;
       drawHistAngularDistributionMap3DCheat(h_angular_distribution_pfparticle_cheat_3D, "cheat_angle_distribution.root", c);


       /*
       for(auto const& shower_spine_hit_list : shower_spine_hit_list_vector){

	 const simb::MCParticle *mcparticle;
	 std::map<int,int> hits_pdg_map;
	 for(auto const& hit : shower_spine_hit_list){
	   truthHitMatcher(hit, mcparticle);
	   if(hits_pdg_map.find(mcparticle->PdgCode()) == hits_pdg_map.end()) hits_pdg_map[mcparticle->PdgCode()] = 1;
	   else hits_pdg_map[mcparticle->PdgCode()] += 1;
	 }
	 
	 vector<pair<int, int>> v;
	 for (map<int, int>::iterator it = hits_pdg_map.begin(); it != hits_pdg_map.end(); it++) {
	   v.push_back({ it->second, it->first });
	 }
	 sort(v.rbegin(), v.rend());
	 rebdautrack_pdg[i][n_recoRebDauTracks[i]] = v[0].second;

       }
       */
       /*
       cout << "run loop over hit PDG map" << endl;
       for(auto const& pair : v){
	 cout << pair.first << " " << pair.second << endl;
       }
       */

       //art::FindManyP<recob::SpacePoint> spacepoint_per_hit(HitHandle, event, fSpacePointproducer);
       //findShowerSpine3D(hits_from_reco_obj, spacepoint_per_hit, unavailable_hit_list, shower_spine_hit_list_vector, end, view_peak_map, best_peak_bins);

       //cout << "best_peak_bins.size(): " << best_peak_bins.size() << endl;
       //cout << "shower_spine_hit_list_vector.size(): " << shower_spine_hit_list_vector.size() << endl;
       // cout << "shower_spine_hit_list_vector[0].size(): " << shower_spine_hit_list_vector[0].size() << endl;
       //cout << "shower_spine_hit_list_vector[1].size(): " << shower_spine_hit_list_vector[1].size() << endl;

       if(best_peak_bins.size()){
	 //int k=0;
	 //for(auto const& entry : best_peak_bins){
	 for(long unsigned int k=0; k<best_peak_bins.size(); k++){

	   best_peak_theta[i][k] = best_peak_bins[k].X() * thetaBinSize;
	   best_peak_phi[i][k] = best_peak_bins[k].Y() * thetaBinSize;

	   // rebuild a track
	   if(shower_spine_hit_list_vector.size()==0) continue;
	   //cout << "Before running trackRebuid" << endl;
	   cout << "shower_spine_hit_list_vector[k].size(): " << shower_spine_hit_list_vector[k].size() << endl;
	   if(!shower_spine_hit_list_vector[k].size()) continue;
	     recob::Track reco_track = trackRebuid(shower_spine_hit_list_vector[k], atrack);
	     cout << "this event has rebuilt daughter track" << endl;
	     //cout << "reco_track.Length(): " << reco_track.Length() << endl;
	   rebdautrack_length[i][n_recoRebDauTracks[i]] = reco_track.Length();
	   //cout << i << " " << n_recoRebDauTracks[i] << " " << rebdautrack_length[i][n_recoRebDauTracks[i]] << endl;

	   const simb::MCParticle *mcparticle;
	   std::map<int,int> hits_pdg_map;
	   for(auto const& hit : shower_spine_hit_list_vector[k]){
	     truthHitMatcher(hit, mcparticle);
	     if(!mcparticle) continue;
	     //cout << "mcparticle->PdgCode(): " << mcparticle->PdgCode() << endl;
	     if(hits_pdg_map.find(mcparticle->PdgCode()) == hits_pdg_map.end()) hits_pdg_map[mcparticle->PdgCode()] = 1;
	     else hits_pdg_map[mcparticle->PdgCode()] += 1;
	   }

	   vector<pair<int, int>> v;
	   if(hits_pdg_map.size()){
	     for (map<int, int>::iterator it = hits_pdg_map.begin(); it != hits_pdg_map.end(); it++) {
	       v.push_back({ it->second, it->first });
	       cout << it->second << " " << it->first << endl;
	     }
	     sort(v.rbegin(), v.rend());
	     rebdautrack_pdg[i][n_recoRebDauTracks[i]] = v[0].second;
	     cout << "rebdautrack_pdg[i][n_recoRebDauTracks[i]]: " << v[0].second << endl;
	   }

	   n_recoRebDauTracks[i]++;
	   // rebuild associations
	   /*
	   cout << "adding reco_track" << endl;
	   anaTrackCollection->push_back(reco_track);

	   std::vector<art::Ptr<recob::Hit>> hits_from_track_rebuild = shower_spine_hit_list_vector[k];
	   std::vector<art::Ptr<recob::SpacePoint>> spacepoint_vec_reco;
	   unsigned int hitsFromSpacePointsRecoSize = 0;

	   for(size_t k_h=0; k_h<hits_from_track_rebuild.size(); k_h++){
	     spacepoint_vec_reco.clear();
	     spacepoint_vec_reco = spacepoint_per_hit_reco.at(hits_from_track_rebuild[k_h].key());
	     hitsFromSpacePointsRecoSize += spacepoint_vec_reco.size();
	   }
	   art::Ptr<recob::Track> pTrackdau(makeTrackPtr(anaTrackCollection->size() - 1));
	   lar_pandora::HitVector anaHitCollection_rebuild_tmp;
	   or(auto hitptr : hits_from_track_rebuild){
	     anaHitCollection_rebuild_tmp.push_back(hitptr);
	   }
	   util::CreateAssn(*this, evt, *(anaTrackCollection.get()), anaHitCollection_rebuild_tmp, *(anaTrackHitAssociations.get()));

	   for (unsigned int hitIndex = 0; hitIndex < hits_from_track_rebuild.size(); hitIndex++){
	     const art::Ptr<recob::Hit> phit(hits_from_track_rebuild.at(hitIndex));
	     const int index((hitIndex < hitsFromSpacePointsRecoSize) ? hitIndex : std::numeric_limits<int>::max());
	     recob::TrackHitMeta metadata(index, -std::numeric_limits<double>::max());
	   }
	   */
	   //k++;
	 }
       }
       /* end of imported bits*/

       // add track points
       n_track_points[i] = track->NPoints();
       for (int ipt=0; ipt < n_track_points[i]; ipt++) {
	   // FIXME: get the track point coordinates!
	   recob::tracking::Point_t pos_point = track->TrajectoryPoint(ipt).position;
	   TVector3 pos(pos_point.x(),pos_point.y(),pos_point.z());
	   pos.GetXYZ(track_point_xyz[i][ipt]);
       }

       track_ID[i] = track->ID();
       track_Prange[i] = trackP.GetTrackMomentum(track_length[i],13);
       double trk_end[4] ={tmp_track_end[0],tmp_track_end[1],tmp_track_end[2],-999};
       bool track_isInside = insideFV( trk_end );
       //check if the track ends within the FV
       if( track_isInside ) track_isContained[i] =1;
       else track_isContained[i] =0;
       //calculate PID
       std::vector<const anab::Calorimetry*> trk_cal = reco_cal.at(i);
       std::vector<const anab::ParticleID*> trk_pid = reco_PID.at(i);
       //std::vector AlgScoresVec = trk_pid.at(0)->ParticleIDAlgScores();

       int plane0=0; 
       int plane1=0;
       int plane2=0;

          for (size_t ipid = 0; ipid < trk_pid.size(); ++ipid){
            if (!trk_pid[ipid]->PlaneID().isValid) continue;
            int planenum = trk_pid[ipid]->PlaneID().Plane;
            if (planenum<0||planenum>2) continue;

            auto pidScore = trk_pid[ipid]->ParticleIDAlgScores();
            for(auto pScore: pidScore){
              //double chi2value = pScore.fValue;
	      
	      if(pScore.fAssumedPdg != 0){
		if(planenum==0) plane0 = pScore.fNdf;
		if(planenum==1) plane1 = pScore.fNdf;
		if(planenum==2) plane2 = pScore.fNdf;
	      }
      	      track_PID_pdg[i][ipid] = pScore.fAssumedPdg;


              // PIDA is always the last one and ndf there is -9999
	      /*
              if(pScore.fAssumedPdg != 0) TrackerData.trkpidndf[iTrk][planenum] = pScore.fNdf; // This value is the same for each particle type, but different in each plane
              switch(pScore.fAssumedPdg){
                case 2212:
                  TrackerData.trkpidchipr[iTrk][planenum] = chi2value;
                  break;
                case 321:
                  TrackerData.trkpidchika[iTrk][planenum] = chi2value;
                  break;
                case 211:
                  TrackerData.trkpidchipi[iTrk][planenum] = chi2value;
                  break;
                case 13:
                  TrackerData.trkpidchimu[iTrk][planenum] = chi2value;
                  break;
                case 0:
                  TrackerData.trkpidpida[iTrk][planenum] = chi2value;
                  break;	      
	    }
	      */
	    }
	  }


	  //int plane0 =   trk_pid[0]->Ndf();
       //int plane1 =   trk_pid[1]->Ndf();
       //int plane2 =   trk_pid[2]->Ndf();
       int best_plane =-1;
       int most_ndf = std::max({plane0, plane1, plane2});
       //for( size_t p =0; p<3; p++) if( most_ndf == trk_pid[p]->Ndf() ) best_plane = p;
       if( most_ndf == plane0 ) best_plane = 0; 
       if( most_ndf == plane1 ) best_plane = 1; 
       if( most_ndf == plane2 ) best_plane = 2; 
       //cout << "best_plane is " << best_plane << endl;

       //track_PID_pdg[i][0] = trk_pid[0]->Pdg();
       //track_PID_pdg[i][1] = trk_pid[1]->Pdg();
       //track_PID_pdg[i][2] = trk_pid[2]->Pdg();
       track_KE[i][0] = trk_cal[0]->KineticEnergy();
       track_KE[i][1] = trk_cal[1]->KineticEnergy();
       track_KE[i][2] = trk_cal[2]->KineticEnergy();

       std::vector<double> PIDAval;
       std::vector<double> chi2;
       PIDAcal( trk_cal, PIDAval);
       int idx =0;
       for(auto const& val : PIDAval){
	 track_PIDA[i][idx] = val;
	 idx ++;
       }
       track_bestplane[i] = best_plane;

       //save dE/dx & dQ/dx
       for (unsigned iplane = 0; iplane < 3; iplane++ ) {
	 int npts = trk_cal[iplane]->dEdx().size();
	 if (npts > MAX_CALO_PTS)
	   npts = MAX_CALO_PTS;
	 if (best_plane >= 0 && iplane == (unsigned)best_plane) {
	   n_cal_points[i] = npts;
	   // if (npts < 0)
	   //   n_cal_points[i] = 0;
	 }

	 n_cal_points_byplane[i][iplane] = npts;

	 const vector<float> &dqdx = trk_cal[iplane]->dQdx();
	 const vector<float> &dedx = trk_cal[iplane]->dEdx();
	 const vector<float> &resr = trk_cal[iplane]->ResidualRange();
	 const vector<float> &pitch = trk_cal[iplane]->TrkPitchVec();
	 auto &xyz = trk_cal[iplane]->XYZ();

	 if (best_plane >= 0 && iplane == (unsigned)best_plane) {
	   std::copy(dqdx.begin(), dqdx.end(), track_dQ_dx[i]);
	   std::copy(dedx.begin(), dedx.end(), track_dE_dx[i]);
	   std::copy(resr.begin(), resr.end(), track_range[i]);
	   std::copy(pitch.begin(), pitch.end(), track_pitch[i]);
	 }
	 std::copy(dqdx.begin(), dqdx.end(), track_dQ_dx_byplane[i][iplane]);
	 std::copy(dedx.begin(), dedx.end(), track_dE_dx_byplane[i][iplane]);
	 std::copy(resr.begin(), resr.end(), track_range_byplane[i][iplane]);
	 std::copy(pitch.begin(), pitch.end(), track_pitch_byplane[i][iplane]);

	 // save calo point's XYZ coords
	 double* coords = (double*)track_calo_xyz_byplane[i][iplane];
	 for(int j = 0; j < npts; j++) {
	   (coords+j*3)[0] = xyz[j].X();
	   (coords+j*3)[1] = xyz[j].Y();
	   (coords+j*3)[2] = xyz[j].Z();
	 }
       }


       //truth matcher
       /*
       double tmpEfrac = 0;
       const simb::MCParticle *particle;
       double tmpComplet = 0;
       std::vector<art::Ptr<recob::Hit>> all_trackHits = track_hits.at(i);
       truthMatcher( all_hits,  all_trackHits, particle, tmpEfrac, tmpComplet );
       if(!particle) continue;
       track_mcID[i] = particle->TrackId();
       track_mcPDG[i] = particle->PdgCode();
       track_Efrac[i] = tmpEfrac;
       track_complet[i] = tmpComplet;
       cout << "track_mcPDG[i]: " << track_mcPDG[i] << endl;
       */

    } // end of track loop


    //CNN dacayID vertex
    art::Handle<std::vector<recob::Vertex>> dcy_vtxHandle;
    std::vector<art::Ptr<recob::Vertex>> dcy_vtxlist;
    if(event.getByLabel(fPointIDModuleLabel,dcy_vtxHandle))
      art::fill_ptr_vector(dcy_vtxlist, dcy_vtxHandle);
    n_decayVtx= dcy_vtxlist.size();

    //art::FindManyP<recob::Track> decay_tracklist(dcy_vtxHandle, event, fPointIDModuleLabel);
    if( n_decayVtx !=0 )
    for( int i=0; i< n_decayVtx; ++i){
       double tmp_vtx[3] ={-999.0,-999.0,-999.0};
       dcy_vtxlist[i]->XYZ(tmp_vtx);
       for( int j=0; j<3; ++j) decayVtx[i][j]=tmp_vtx[j];
       //std::vector<art::Ptr<recob::Track>>  decay_track = decay_tracklist.at(i);
       //cout<<"how many tracks? "<<decay_track.size()<<endl;
    }


    //Showers... for background rejection
    /*
    art::Handle<std::vector<recob::Shower>> showerHandle;
    std::vector<art::Ptr<recob::Shower>> showerlist;
    if(event.getByLabel(fShowerModuleLabel,showerHandle))
      art::fill_ptr_vector(showerlist, showerHandle);
    n_recoShowers= showerlist.size();
    */
    if( n_recoShowers != 0 )
    for(int i=0; i<n_recoShowers && i< MAX_SHOWERS; ++i){
       art::Ptr<recob::Shower> shower = showerlist[i];
       sh_direction_X[i] = shower->Direction().X();
       sh_direction_Y[i] = shower->Direction().Y();
       sh_direction_Z[i] = shower->Direction().Z();
       sh_start_X[i] = shower->ShowerStart().X();
       sh_start_Y[i] = shower->ShowerStart().Y();
       sh_start_Z[i] = shower->ShowerStart().Z();
       sh_bestplane[i] = shower->best_plane();
       sh_length[i] = shower->Length();
       for( size_t j =0; j<shower->Energy().size(); j ++) sh_energy[i][j] = shower->Energy()[j];
       for( size_t j =0; j<shower->MIPEnergy().size(); j++) sh_MIPenergy[i][j] = shower->MIPEnergy()[j];
       for( size_t j =0; j<shower->dEdx().size(); j++) sh_dEdx[i][j] = shower->dEdx()[j];
    }

    if( fSaveHits ){
      //Hits
      n_recoHits= all_hits.size();
      if( n_recoHits != 0 )
      for(int i = 0; i < n_recoHits && i < MAX_HITS ; ++i){//loop over hits
         hit_channel[i] = all_hits[i]->Channel();
         hit_tpc[i]   = all_hits[i]->WireID().TPC;
         hit_plane[i]   = all_hits[i]->WireID().Plane;
         hit_wire[i]    = all_hits[i]->WireID().Wire;
         hit_peakT[i]   = all_hits[i]->PeakTime();
         hit_charge[i]  = all_hits[i]->Integral();
         hit_ph[i]  = all_hits[i]->PeakAmplitude();
         hit_startT[i] = all_hits[i]->PeakTimeMinusRMS();
         hit_endT[i] = all_hits[i]->PeakTimePlusRMS();
         hit_rms[i] = all_hits[i]->RMS();
         //hit_electrons[i]  = fCalorimetryAlg.ElectronsFromADCArea( all_hits[i]->Integral(), all_hits[i]->WireID().Plane) * fCalorimetryAlg.LifetimeCorrection( all_hits[i]->PeakTime() );
      }
    }

}
//========================================================================
void HitSplitAlg::PIDAcal( std::vector<const anab::Calorimetry*> cal, std::vector<double> &PIDA){
  std::vector<double> pida_vec_v0;
  std::vector<double> pida_vec_v1;
  std::vector<double> pida_vec_v2;

  int used_points[3] ={0,0,0};
  double tmp_pida[3];
  for( unsigned j =0; j<3; ++j){
     for( unsigned i =0; i<cal[j]->dEdx().size(); ++i ) { // loop through hits on each plane
        if( cal[j]->ResidualRange()[i] < fPIDA_endPoint ) { // Only want PIDA for last x cm
          tmp_pida[j] = cal[j]->dEdx()[i]* pow(cal[j]->ResidualRange()[i], fExponentConstant );
          if(fMinPIDAValue > tmp_pida[j] || tmp_pida[j] > fMaxPIDAValue) continue;
          if( j ==  0 )pida_vec_v0.push_back(tmp_pida[j]);
          if( j ==  1 )pida_vec_v1.push_back(tmp_pida[j]);
          if( j ==  2 )pida_vec_v2.push_back(tmp_pida[j]);
          used_points[j] ++;
        } // If ResRange < x cm
     }// Loop over hits on each plane
  }


  //for each view calculate PIDA median value
  std::sort(pida_vec_v0.begin(), pida_vec_v0.end() );
  int size_v0 = pida_vec_v0.size();
  double median_v0 = -999;
  if( size_v0 > 0 ) median_v0 = size_v0 % 2 ? pida_vec_v0[size_v0 / 2] : (pida_vec_v0[size_v0 / 2 - 1] + pida_vec_v0[size_v0 / 2]) / 2;

  std::sort(pida_vec_v1.begin(), pida_vec_v1.end() );
  int size_v1 = pida_vec_v1.size();
  double median_v1 = -999;
  if( size_v1 > 0 ) median_v1 = size_v1 % 2 ? pida_vec_v1[size_v1 / 2] : (pida_vec_v1[size_v1 / 2 - 1] + pida_vec_v1[size_v1 / 2]) / 2;

  std::sort(pida_vec_v2.begin(), pida_vec_v2.end() );
  int size_v2 = pida_vec_v2.size();
  double median_v2 = -999;
  if( size_v2 > 0 ) median_v2 = size_v2 % 2 ? pida_vec_v2[size_v2 / 2] : (pida_vec_v2[size_v2 / 2 - 1] + pida_vec_v2[size_v2 / 2]) / 2;

  PIDA.push_back(median_v0);
  PIDA.push_back(median_v1);
  PIDA.push_back(median_v2);

}

//========================================================================
void HitSplitAlg::truthHitMatcher( art::Ptr<recob::Hit> hit, const simb::MCParticle *&MCparticle){

    art::ServiceHandle<cheat::BackTrackerService> bt;
    art::ServiceHandle<cheat::ParticleInventoryService> part_inv;
    //detinfo::DetectorClocksData const& clockData;

    std::map<int,double> trkID_E;
    std::vector<sim::TrackIDE> TrackIDs = bt->HitToTrackIDEs(clockdata, hit);

    if(!TrackIDs.size()){
      MCparticle = 0;
      return;
    }
    for(size_t k = 0; k < TrackIDs.size(); k++){
      //cout << TrackIDs[k].trackID << " " << TrackIDs[k].energy << endl;
      trkID_E[TrackIDs[k].trackID] += TrackIDs[k].energy;
    }

    double E_em =0.0;
    double max_E = -999.0;
    double total_E = 0.0;
    int TrackID = -999;
    // amount of energy deposited by the particle that deposited more energy... tomato potato... blabla
    //!if the collection of hits have more than one particle associate save the particle w/ the highest energy deposition
    //!since we are looking for muons/pions/protons this should be enough
    if( !trkID_E.size() ) {
      MCparticle = 0;
      return; //Ghost track???
    }
    for(std::map<int,double>::iterator ii = trkID_E.begin(); ii!=trkID_E.end(); ++ii){
       total_E += ii->second;
       if((ii->second)>max_E){
         max_E = ii->second;
         TrackID = ii->first;
         if( TrackID < 0 ) E_em += ii->second;
       }
    }
    // consider the most energetic hit to estimate PDG etc
    MCparticle = part_inv->TrackIdToParticle_P(TrackID);
}


//========================================================================
void HitSplitAlg::truthMatcher( std::vector<art::Ptr<recob::Hit>> all_hits, std::vector<art::Ptr<recob::Hit>> track_hits, const simb::MCParticle *&MCparticle, double &Efrac, double &Ecomplet){

    art::ServiceHandle<cheat::BackTrackerService> bt;
    art::ServiceHandle<cheat::ParticleInventoryService> part_inv;
    std::map<int,double> trkID_E;
    for(size_t j = 0; j < track_hits.size(); ++j){
      art::Ptr<recob::Hit> hit = track_hits[j];
      std::vector<sim::TrackIDE> TrackIDs = bt->HitToTrackIDEs(clockdata, hit);
      for(size_t k = 0; k < TrackIDs.size(); k++){
	trkID_E[TrackIDs[k].trackID] += TrackIDs[k].energy;
      }      
    }
    double E_em =0.0;
    double max_E = -999.0;
    double total_E = 0.0;
    int TrackID = -999;
    double partial_E =0.0; // amount of energy deposited by the particle that deposited more energy... tomato potato... blabla
    //!if the collection of hits have more than one particle associate save the particle w/ the highest energy deposition
    //!since we are looking for muons/pions/protons this should be enough
    if( !trkID_E.size() ) {
      MCparticle = 0;
      return; //Ghost track???
    }
    for(std::map<int,double>::iterator ii = trkID_E.begin(); ii!=trkID_E.end(); ++ii){
       total_E += ii->second;
       if((ii->second)>max_E){
         partial_E = ii->second;
         max_E = ii->second;
         TrackID = ii->first;
         if( TrackID < 0 ) E_em += ii->second;
       }
    }
    // consider the most energetic hit to estimate PDG etc
    MCparticle = part_inv->TrackIdToParticle_P(TrackID);

    //In the current simulation, we do not save EM Shower daughters in GEANT. But we do save the energy deposition in TrackIDEs. If the energy deposition is from a particle that is the daughter of
    //an EM particle, the negative of the parent track ID is saved in TrackIDE for the daughter particle
    //we don't want to track gammas or any other EM activity
    if( TrackID < 0 ) return;

    //Efrac = (partial_E+E_em)/total_E;
    Efrac = (partial_E)/total_E;

    //completeness
    double totenergy =0;
    for(size_t k = 0; k < all_hits.size(); ++k){
      //art::Ptr<recob::Hit> hit = all_hits[k];
       /*
       std::vector<sim::TrackIDE> TrackIDs = bt->HitToTrackIDEs(hit);
       for(size_t l = 0; l < TrackIDs.size(); ++l){
          if(TrackIDs[l].trackID==TrackID) totenergy += TrackIDs[l].energy;
       }
       */
    }
    Ecomplet = partial_E/totenergy;

}
//========================================================================

double HitSplitAlg::truthLength( const simb::MCParticle *MCparticle ){
   //calculate the truth length considering only the part that is inside the TPC
   //Base on a peace of code from dune/TrackingAna/TrackingEfficiency_module.cc

   if( !MCparticle ) return -999.0;
   double TPCLength = 0.0;
   int numberTrajectoryPoints = MCparticle->NumberTrajectoryPoints();
   std::vector<double> TPCLengthHits(numberTrajectoryPoints, 0);
   int FirstHit=0, LastHit=0;
   bool BeenInVolume = false;

   for(int MCHit=0; MCHit < numberTrajectoryPoints; ++MCHit) {
     //const TLorentzVector& tmpPosition= MCparticle->Position(MCHit);
     auto const tmpPosition = geo::vect::toPoint(MCparticle->Position(MCHit).Vect());
      //double const tmpPosArray[]={tmpPosition[0],tmpPosition[1],tmpPosition[2]};
      if (MCHit!=0) TPCLengthHits[MCHit] = sqrt( pow( (MCparticle->Vx(MCHit-1)-MCparticle->Vx(MCHit)),2)+ pow( (MCparticle->Vy(MCHit-1)-MCparticle->Vy(MCHit)),2)+ pow( (MCparticle->Vz(MCHit-1)-MCparticle->Vz(MCHit)),2));

      geo::TPCID tpcid = geom->FindTPCAtPosition(tmpPosition);
      //geo::TPCID tpcid = geom->FindTPCAtPosition(tmpPosArray);
      if(tpcid.isValid) {
        // -- Check if hit is within drift window...
        //geo::CryostatGeo const& cryo = geom->Cryostat(geo::CryostatID{0});
	//geo::CryostatGeo const& cryo = geom->Cryostat(tpcid.Cryostat); 
        //geo::TPCGeo      const& tpc  = cryo.TPC(tpcid.TPC);
        //double XPlanePosition      = tpc.fPlaneLocation(0)[0];
	geo::TPCGeo      const& tpc  = geom->TPC(tpcid);
	double XPlanePosition      = tpc.Plane(0).GetCenter().X();
        double DriftTimeCorrection = fabs( tmpPosition.X() - XPlanePosition ) / XDriftVelocity;
        double TimeAtPlane         = MCparticle->T() + DriftTimeCorrection;
        if( TimeAtPlane < detprop.TimeOffsetZ() || TimeAtPlane > detprop.TimeOffsetZ() + WindowSize ) continue;
        LastHit = MCHit;
        if( !BeenInVolume ) {
	  BeenInVolume = true;
          FirstHit = MCHit;
	}
      }
   }
   for (int Hit = FirstHit+1; Hit <= LastHit; ++Hit ) TPCLength += TPCLengthHits[Hit];

   return TPCLength;
}

//========================================================================
//========================================================================
bool HitSplitAlg::insideFV( double vertex[4]){

     double x = vertex[0];
     double y = vertex[1];
     double z = vertex[2];

     cout << "FidVolXmin: " << fFidVolXmin << ", fFidVolXmax: " << fFidVolXmax << ", X: " << x << endl;
     cout << "FidVolYmin: " << fFidVolYmin << ", fFidVolYmax: " << fFidVolYmax << ", Y: " << y << endl;
     cout << "FidVolZmin: " << fFidVolZmin << ", fFidVolZmax: " << fFidVolZmax << ", Z: " << z << endl;

     if (x>fFidVolXmin && x<fFidVolXmax&&
	 y>fFidVolYmin && y<fFidVolYmax&&
	 z>fFidVolZmin && z<fFidVolZmax)
       return true;
     else
       return false;
}

void HitSplitAlg::reset() {
    Event = 0;
    Run = 0;
    SubRun = 0;

    //MC truth
    MC_Ev = 0;
    MC_nuPDG = 0;
    MC_Q2 = 0;
    MC_hit_nucleon = 0;
    MC_cc = 0;
    MCgenie_npart = 0;
    std::memset( MCgenie_id, 0, sizeof(MCgenie_id) );
    std::memset( MCgenie_pdg, 0, sizeof(MCgenie_pdg) );
    std::memset( MCgenie_mother, 0, sizeof(MCgenie_mother) );
    std::memset( MCgenie_statusCode, 0, sizeof(MCgenie_statusCode) );
    std::memset( MCgenie_fate, 0, sizeof(MCgenie_fate) );
    std::memset( MCgenie_startMomentum, 0, sizeof(MCgenie_startMomentum) );
    std::memset( MCgenie_endMomentum, 0, sizeof(MCgenie_endMomentum) );

    std::memset( MC_vertex, 0, sizeof(MC_vertex) );
    MC_npart = 0;
    std::memset( MC_id, 0, sizeof(MC_id) );
    std::memset( MC_pdg, 0, sizeof(MC_pdg) );
    std::memset( MC_mother, 0, sizeof(MC_mother) );
    std::memset( MC_startXYZT, 0, sizeof(MC_startXYZT) );
    std::memset( MC_endXYZT, 0, sizeof(MC_endXYZT) );
    std::memset( MC_startMomentum, 0, sizeof(MC_startMomentum) );
    std::memset( MC_endMomentum, 0, sizeof(MC_endMomentum) );
    std::memset( MC_truthlength, 0, sizeof(MC_truthlength) );
    std::memset( MC_Prange, 0, sizeof(MC_Prange) );
    std::memset( MC_statusCode, 0, sizeof(MC_statusCode) );
    MC_kaontruthlength = 0;
    MC_dautruthlength = 0;
    Kaon_BR = -1;
    MC_process.clear();
    MC_Endprocess.clear();

    n_vertices = 0;
    std::memset( vertex, 0, sizeof(vertex) );
    std::memset( vtx_ID, 0, sizeof(vtx_ID) );
    n_decayVtx = 0;
    std::memset( decayVtx, 0, sizeof(decayVtx) );
    n_recoTracks = 0;
    std::memset( n_recoDauTracks, 0, sizeof(n_recoDauTracks) );
    std::memset( n_recoRebDauTracks, 0, sizeof(n_recoRebDauTracks) );
    std::memset( track_isContained, 0, sizeof(track_isContained) );
    std::memset( track_ID, 0, sizeof(track_ID) );
    std::memset( vtxID_trk, 0, sizeof(vtxID_trk) );
    std::memset( track_vtxDir, 0, sizeof(track_vtxDir) );
    std::memset( track_vtx, 0, sizeof(track_vtx) );
    std::memset( track_end, 0, sizeof(track_end) );
    std::memset( track_length, 0, sizeof(track_length) );
    std::memset( dautrack_length, 0, sizeof(dautrack_length) );
    std::memset( rebdautracktrue_length, 0, sizeof(rebdautracktrue_length) );
    std::memset( rebdautracktruedir_length, 0, sizeof(rebdautracktruedir_length) );
    std::memset( rebdautrack_length, 0, sizeof(rebdautrack_length) );
    std::memset( dautrack_pdg, 0, sizeof(dautrack_pdg) );
    std::memset( rebdautrack_pdg, 0, sizeof(rebdautrack_pdg) );
    std::memset( best_peak_theta, 0, sizeof(best_peak_theta) );
    std::memset( best_peak_phi, 0, sizeof(best_peak_phi) );
    std::memset( best_peak_theta_true, 0, sizeof(best_peak_theta_true) );
    std::memset( best_peak_phi_true, 0, sizeof(best_peak_phi_true) );
    std::memset( track_dir_vtx, 0, sizeof(track_dir_vtx) );
    std::memset( track_PIDA, 0, sizeof(track_PIDA) );
    std::memset( track_PID_pdg, 0, sizeof(track_PID_pdg) );
    std::memset( track_KE, 0, sizeof(track_KE) );
    std::memset( track_bestplane, 0, sizeof(track_bestplane) );
    std::memset( track_Prange, 0, sizeof(track_Prange) );
    std::memset( track_Efrac, 0, sizeof(track_Efrac) );
    std::memset( track_complet, 0, sizeof(track_complet) );
    std::memset( track_mcID, 0, sizeof(track_mcID) );
    std::memset( track_mcPDG, 0, sizeof(track_mcPDG) );
    std::memset( dautrack_mcPDG, 0, sizeof(dautrack_mcPDG) );
    std::memset( n_track_points, 0, sizeof(n_track_points) );
    std::memset( track_point_xyz, 0, sizeof(track_point_xyz) );
    std::memset( n_cal_points, 0, sizeof(n_cal_points) );
    std::memset( track_dQ_dx, 0, sizeof(track_dQ_dx) );
    std::memset( track_dE_dx, 0, sizeof(track_dE_dx) );
    std::memset( track_range, 0, sizeof(track_range) );
    std::memset( track_pitch, 0, sizeof(track_pitch) );
    std::memset( n_cal_points_byplane, 0, sizeof(n_cal_points_byplane) );
    std::memset( track_dQ_dx_byplane, 0, sizeof(track_dQ_dx_byplane) );
    std::memset( track_dE_dx_byplane, 0, sizeof(track_dE_dx_byplane) );
    std::memset( track_range_byplane, 0, sizeof(track_range_byplane) );
    std::memset( track_pitch_byplane, 0, sizeof(track_pitch_byplane) );

    n_recoHits = 0;
    std::memset( hit_channel, 0, sizeof(hit_channel) );
    std::memset( hit_tpc, 0, sizeof(hit_tpc) );
    std::memset( hit_plane, 0, sizeof(hit_plane) );
    std::memset( hit_wire, 0, sizeof(hit_wire) );
    std::memset( hit_peakT, 0, sizeof(hit_peakT) );
    std::memset( hit_charge, 0, sizeof(hit_charge) );
    std::memset( hit_ph, 0, sizeof(hit_ph) );
    std::memset( hit_startT, 0, sizeof(hit_startT) );
    std::memset( hit_endT, 0, sizeof(hit_endT) );
    std::memset( hit_rms, 0, sizeof(hit_rms) );
    std::memset( hit_electrons, 0, sizeof(hit_electrons) );


    Em_ch = 0;
    Em_e = 0;
    trk_e = 0;
    Emichel_e = 0;

    n_recoShowers = 0;
    std::memset( sh_direction_X, 0, sizeof(sh_direction_X) );
    std::memset( sh_direction_Y, 0, sizeof(sh_direction_Y) );
    std::memset( sh_direction_Z, 0, sizeof(sh_direction_Z) );
    std::memset( sh_start_X, 0, sizeof(sh_start_X) );
    std::memset( sh_start_Y, 0, sizeof(sh_start_Y) );
    std::memset( sh_start_Z, 0, sizeof(sh_start_Z) );
    std::memset( sh_energy, 0, sizeof(sh_energy) );
    std::memset( sh_MIPenergy, 0, sizeof(sh_MIPenergy) );
    std::memset( sh_dEdx, 0, sizeof(sh_dEdx) );
    std::memset( sh_bestplane, 0, sizeof(sh_bestplane) );
    std::memset( sh_length, 0, sizeof(sh_length) );


    n_flashes = 0;
    std::memset( flash_time, 0, sizeof(flash_time) );
    std::memset( flash_pe, 0, sizeof(flash_pe) );
    std::memset( flash_ycenter, 0, sizeof(flash_ycenter) );
    std::memset( flash_zcenter, 0, sizeof(flash_zcenter) );
    std::memset( flash_ywidth, 0, sizeof(flash_ywidth) );
    std::memset( flash_zwidth, 0, sizeof(flash_zwidth) );
    std::memset( flash_timewidth, 0, sizeof(flash_timewidth) );
    std::memset( flash_abstime, 0, sizeof(flash_abstime) );
    std::memset( flash_frame, 0, sizeof(flash_frame) );
    std::memset( flash_PE_ndk, 0, sizeof(flash_PE_ndk) );
    std::memset( flash_PE_Ar39, 0, sizeof(flash_PE_Ar39) );
}



}
