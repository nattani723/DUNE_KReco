//Module analyzer
//Ana module for nucleon decay and atmospheric analysis
//Ana TTree contains MC truth and reconstruction info
//ahiguera@central.uh.edu


//header file
#include "RunningAnalyzer_module.h"
#include "ReconstructionOrchestrator.h"


#ifdef __CLING__
#pragma link C++ class std::vector < std::vector<Float_t> >+;
#pragma link C++ class std::vector < std::vector< std::vector<Float_t> > >+;
#endif

using namespace std;

//========================================================================

namespace kaon_reconstruction{

  RunningAnalyzer::RunningAnalyzer(fhicl::ParameterSet const& parameterSet)
    : EDAnalyzer(parameterSet), fCalorimetryAlg(parameterSet.get< fhicl::ParameterSet >("CalorimetryAlg"))

  {
    reconfigure(parameterSet);
  }
  //========================================================================
  RunningAnalyzer::~RunningAnalyzer(){
    //destructor
  }
  //========================================================================
  void RunningAnalyzer::reconfigure(fhicl::ParameterSet const& p)
  {

    fHitModuleLabel      = p.get<std::string>("HitModuleLabel");
    fHitModuleLabelOLD      = p.get<std::string>("HitModuleLabelOLD");
    fHitSPAssns          = p.get<std::string>("HitSPAssns");
    fTrackModuleLabel    = p.get<std::string>("TrackModuleLabel");
    fTrackRecoAlgModuleLabel    = p.get<std::string>("TrackRecoAlgModuleLabel");
    fPidModuleLabel    = p.get<std::string>("PidModuleLabel");
    fPidRecoAlgModuleLabel    = p.get<std::string>("PidRecoAlgModuleLabel");
    fCaloModuleLabel     = p.get<std::string>("CaloModuleLabel");
    fCaloRecoAlgModuleLabel     = p.get<std::string>("CaloRecoAlgModuleLabel");
    fOpFlashModuleLabel  = p.get<std::string>("OpFlashModuleLabel");
    fShowerModuleLabel = p.get<std::string>("ShowerModuleLabel");
    fPointIDModuleLabel  = p.get<std::string>("PointIDModuleLabel");
    fNNetModuleLabel     = p.get<std::string>("NNetModuleLabel");
    fMCgenieLabel = p.get<std::string>("MCgenieLabel");
    fSpacePointproducer  = p.get<std::string>("SpacePointproducer");
    fHitTrackAssns = p.get<std::string>("HitTrackAssns");
    fHitTrackRecoAlgAssns = p.get<std::string>("HitTrackRecoAlgAssns");
    fSaveHits = p.get<bool>("SaveHits");
    fExponentConstant  = p.get<double>("ExponentConstant");
    fMaxPIDAValue  = p.get<double>("MaxPIDAValue");
    fMinPIDAValue  = p.get<double>("MinPIDAValue");
    fPIDA_endPoint = p.get<double>("PIDA_endPoint");
    fView                = p.get<double>("View");
    fPidValue = p.get<double>("PidValue");
    fFidVolCutX          = p.get<double>("FidVolCutX");
    fFidVolCutY          = p.get<double>("FidVolCutY");
    fFidVolCutZ          = p.get<double>("FidVolCutZ");
  }
  //========================================================================
  void RunningAnalyzer::beginJob(){
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
 
    }

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
    fEventTree->Branch("n_reco_dautracks", n_recoDauTracks, "n_recoDauTracks[n_reco_tracks]/I");
    fEventTree->Branch("n_decayVtx", &n_decayVtx);
    fEventTree->Branch("decayVtx", decayVtx,"decayVtx[n_decayVtx][3]/D");  //vertices found using decayID point Alg
    fEventTree->Branch("vtxID_trk", vtxID_trk,"vtxID_trk[n_reco_tracks][10]/I"); //track-vertex association
    fEventTree->Branch("track_vtxDir", track_vtxDir,"track_vtxDir[n_reco_tracks][3]/D");
    fEventTree->Branch("track_vtx", track_vtx,"track_vtx[n_reco_tracks][4]/D");
    fEventTree->Branch("track_end", track_end,"track_end[n_reco_tracks][4]/D");
    fEventTree->Branch("track_isContained", track_isContained,"track_isContained[n_reco_tracks]/I");
    fEventTree->Branch("track_ID", track_ID,"track_ID[n_reco_tracks]/I");
    fEventTree->Branch("track_length", track_length,"track_length[n_reco_tracks]/D");
    fEventTree->Branch("dau_track_length", dau_track_length,"dau_track_length[n_reco_tracks][10]/D");
    fEventTree->Branch("dau_track_distance", dau_track_distance,"dau_track_distance[n_reco_tracks][10]/D");
    fEventTree->Branch("dau_track_pdg", dau_track_pdg,"dau_track_pdg[n_reco_tracks][10]/D");
    fEventTree->Branch("dau_track_mcPDG", dau_track_mcPDG,"dau_track_mcPDG[n_reco_tracks][10]/I");        //true MC PDG for a given daughter track
    fEventTree->Branch("dau_track_vtx", dau_track_vtx,"dau_track_vtx[n_reco_tracks][10][4]/D");
    fEventTree->Branch("dau_track_end", dau_track_end,"dau_track_end[n_reco_tracks][10][4]/D");


    fEventTree->Branch("n_reco_tracks_RecoAlg", &n_recoTracks_RecoAlg);
    fEventTree->Branch("n_reco_dautracks_RecoAlg", n_recoDauTracks_RecoAlg, "n_recoDauTracks_RecoAlg[n_reco_tracks]/I");
    fEventTree->Branch("dau_track_length_RecoAlg", dau_track_length_RecoAlg,"dau_track_length_RecoAlg[n_reco_tracks][10]/D");
    fEventTree->Branch("dau_track_distance_RecoAlg", dau_track_distance_RecoAlg,"dau_track_distance_RecoAlg[n_reco_tracks][10]/D");
    fEventTree->Branch("dau_track_pdg_RecoAlg", dau_track_pdg_RecoAlg,"dau_track_pdg_RecoAlg[n_reco_tracks][10]/D");
    fEventTree->Branch("dau_track_mcPDG_RecoAlg", dau_track_mcPDG_RecoAlg,"dau_track_mcPDG_RecoAlg[n_reco_tracks][10]/I");
    fEventTree->Branch("dau_track_vtx_RecoAlg", dau_track_vtx_RecoAlg,"dau_track_vtx_RecoAlg[n_reco_tracks][10][4]/D");
    fEventTree->Branch("dau_track_end_RecoAlg", dau_track_end_RecoAlg,"dau_track_end_RecoAlg[n_reco_tracks][10][4]/D");

    fEventTree->Branch("n_reco_rebdautracks", n_recoRebDauTracks, "n_recoRebDauTracks[n_reco_tracks]/I");
    fEventTree->Branch("rebdautrack_distance", rebdautrack_distance,"rebdautrack_distance[n_reco_tracks][10]/D");
    fEventTree->Branch("rebdautracktrue_length", rebdautracktrue_length,"rebdautracktrue_length[n_reco_tracks]/D");
    fEventTree->Branch("rebdautracktruedir_length", rebdautracktruedir_length,"rebdautracktruedir_length[n_reco_tracks]/D");
    fEventTree->Branch("rebdautrack_length", rebdautrack_length,"rebdautrack_length[n_reco_tracks][10]/D");
    fEventTree->Branch("rebdautrack_pdg", rebdautrack_pdg,"rebdautrack_pdg[n_reco_tracks][10]/D");

    fEventTree->Branch("best_peak_x", best_peak_x,"best_peak_x[n_reco_tracks][10]/D");
    fEventTree->Branch("best_peak_y", best_peak_y,"best_peak_y[n_reco_tracks][10]/D");
    fEventTree->Branch("best_peak_z", best_peak_z,"best_peak_z[n_reco_tracks][10]/D");
    fEventTree->Branch("best_peak_x_true", best_peak_x_true,"best_peak_x_true[n_reco_tracks]/D");  
    fEventTree->Branch("best_peak_y_true", best_peak_y_true,"best_peak_y_true[n_reco_tracks]/D");  
    fEventTree->Branch("best_peak_z_true", best_peak_z_true,"best_peak_z_true[n_reco_tracks]/D");  

    fEventTree->Branch("track_PIDA", track_PIDA,"track_PIDA[n_reco_tracks][3]/D");//
    fEventTree->Branch("track_PID_pdg", track_PID_pdg,"track_PID_pdg[n_reco_tracks][3]/I");
    fEventTree->Branch("track_KE", track_KE,"track_KE[n_reco_tracks][3]/D");
    fEventTree->Branch("track_Prange", track_Prange,"track_Prange[n_reco_tracks]/D");
    fEventTree->Branch("track_bestplane", track_bestplane,"track_bestplane[n_reco_tracks]/I");
    fEventTree->Branch("n_track_points", n_track_points,"n_track_points[n_reco_tracks]/I");
    fEventTree->Branch("track_point_xyz", track_point_xyz,"track_point_xyz[n_reco_tracks][500][3]/D");

    fEventTree->Branch("n_cal_points", n_cal_points,"n_cal_points[n_reco_tracks]/I");
    fEventTree->Branch("track_dQ_dx", track_dQ_dx,"track_dQ_dx[n_reco_tracks][500]/D");
    fEventTree->Branch("track_dE_dx", track_dE_dx,"track_dE_dx[n_reco_tracks][500]/D");
    fEventTree->Branch("track_range", track_range,"track_range[n_reco_tracks][500]/D");
    fEventTree->Branch("track_pitch", track_pitch,"track_pitch[n_reco_tracks][500]/D");

    fEventTree->Branch("n_cal_points_byplane", n_cal_points_byplane,"n_cal_points_byplane[n_reco_tracks][3]/I");
    fEventTree->Branch("track_dQ_dx_byplane", track_dQ_dx_byplane,"track_dQ_dx_byplane[n_reco_tracks][3][500]/D");
    fEventTree->Branch("track_dE_dx_byplane", track_dE_dx_byplane,"track_dE_dx_byplane[n_reco_tracks][3][500]/D");
    fEventTree->Branch("track_range_byplane", track_range_byplane,"track_range_byplane[n_reco_tracks][3][500]/D");
    fEventTree->Branch("track_pitch_byplane", track_pitch_byplane,"track_pitch_byplane[n_reco_tracks][3][500]/D");
    fEventTree->Branch("track_calo_xyz_byplane", track_calo_xyz_byplane,"track_calo_xyz_byplane[n_reco_tracks][3][500][3]/D");

    fEventTree->Branch("track_chi_pr", track_chi_pr, "track_chi_pr[n_reco_tracks][3]/D");
    fEventTree->Branch("track_chi_ka", track_chi_ka, "track_chi_ka[n_reco_tracks][3]/D");
    fEventTree->Branch("track_chi_pi", track_chi_pi, "track_chi_pi[n_reco_tracks][3]/D");
    fEventTree->Branch("track_chi_mu", track_chi_mu, "track_chi_mu[n_reco_tracks][3]/D");
    fEventTree->Branch("track_pida", track_pida, "track_pida[n_reco_tracks][3]/D");
    fEventTree->Branch("track_pidndf", track_pidndf, "track_pidndf[n_reco_tracks][3]/D");


    fEventTree->Branch("dau_track_complet", dau_track_complet,"dau_track_complet[n_reco_tracks][10]/D");  //track quality variable (completeness)
    fEventTree->Branch("dau_track_Efrac", dau_track_Efrac,"dau_track_Efrac[n_reco_tracks][10]/D");        //track quality variable (purity)
    fEventTree->Branch("dau_track_mcID", dau_track_mcID,"dau_track_mcID[n_reco_tracks][10]/I");           //true MC ID for a given track
    //fEventTree->Branch("dau_track_mcPDG", dau_track_mcPDG,"dau_track_mcPDG[n_reco_tracks][10]/I");        //true MC PDG for a given track

    fEventTree->Branch("dau_track_PIDA", dau_track_PIDA,"dau_track_PIDA[n_reco_tracks][10][3]/D");
    fEventTree->Branch("dau_track_PID_pdg", dau_track_PID_pdg,"dau_track_PID_pdg[n_reco_tracks][10][3]/I");
    fEventTree->Branch("dau_track_KE", dau_track_KE,"dau_track_KE[n_reco_tracks][10][3]/D");
    fEventTree->Branch("dau_track_bestplane", dau_track_bestplane,"dau_track_bestplane[n_reco_tracks][10]/I");
    fEventTree->Branch("n_dau_track_points", n_dau_track_points,"n_dau_track_points[n_reco_tracks][10]/I");
    fEventTree->Branch("n_dau_cal_points", n_dau_cal_points,"n_dau_cal_points[n_reco_tracks][10]/I");
    fEventTree->Branch("n_dau_cal_points_byplane", n_dau_cal_points_byplane,"n_dau_cal_points_byplane[n_reco_tracks][10][3]/I");
    fEventTree->Branch("dau_track_dQ_dx", dau_track_dQ_dx,"dau_track_dQ_dx[n_reco_tracks][10][500]/D");
    fEventTree->Branch("dau_track_dE_dx", dau_track_dE_dx,"dau_track_dE_dx[n_reco_tracks][10][500]/D");
    fEventTree->Branch("dau_track_range", dau_track_range,"dau_track_range[n_reco_tracks][10][500]/D");
    fEventTree->Branch("dau_track_pitch", dau_track_pitch,"dau_track_pitch[n_reco_tracks][10][500]/D");
    fEventTree->Branch("dau_track_dQ_dx_byplane", dau_track_dQ_dx_byplane,"dau_track_dQ_dx_byplane[n_reco_tracks][10][3][500]/D");
    fEventTree->Branch("dau_track_dE_dx_byplane", dau_track_dE_dx_byplane,"dau_track_dE_dx_byplane[n_reco_tracks][10][3][500]/D");
    fEventTree->Branch("dau_track_range_byplane", dau_track_range_byplane,"dau_track_range_byplane[n_reco_tracks][10][3][500]/D");
    fEventTree->Branch("dau_track_pitch_byplane", dau_track_pitch_byplane,"dau_track_pitch_byplane[n_reco_tracks][10][3][500]/D");
    fEventTree->Branch("dau_track_calo_xyz_byplane", dau_track_calo_xyz_byplane,"dau_track_calo_xyz_byplane[n_reco_tracks][10][3][500][3]/D");

    fEventTree->Branch("dau_track_chi_pr", dau_track_chi_pr, "dau_track_chi_pr[n_reco_tracks][10][3]/D");
    fEventTree->Branch("dau_track_chi_ka", dau_track_chi_ka, "dau_track_chi_ka[n_reco_tracks][10][3]/D");
    fEventTree->Branch("dau_track_chi_pi", dau_track_chi_pi, "dau_track_chi_pi[n_reco_tracks][10][3]/D");
    fEventTree->Branch("dau_track_chi_mu", dau_track_chi_mu, "dau_track_chi_mu[n_reco_tracks][10][3]/D");
    fEventTree->Branch("dau_track_pida", dau_track_pida, "dau_track_pida[n_reco_tracks][10][3]/D");
    fEventTree->Branch("dau_track_pidndf", dau_track_pidndf, "dau_track_pidndf[n_reco_tracks][10][3]/D");


    fEventTree->Branch("dau_track_complet_RecoAlg", dau_track_complet_RecoAlg,"dau_track_complet_RecoAlg[n_reco_tracks][10]/D");  //track quality variable (completeness)
    fEventTree->Branch("dau_track_Efrac_RecoAlg", dau_track_Efrac_RecoAlg,"dau_track_Efrac_RecoAlg[n_reco_tracks][10]/D");        //track quality variable (purity)
    fEventTree->Branch("dau_track_mcID_RecoAlg", dau_track_mcID_RecoAlg,"dau_track_mcID_RecoAlg[n_reco_tracks][10]/I");           //true MC ID for a given track
    //fEventTree->Branch("dau_track_mcPDG_RecoAlg", dau_track_mcPDG_RecoAlg,"dau_track_mcPDG_RecoAlg[n_reco_tracks][10]/I");        //true MC PDG for a given track

    fEventTree->Branch("dau_track_PIDA_RecoAlg", dau_track_PIDA_RecoAlg,"dau_track_PIDA_RecoAlg[n_reco_tracks][10][3]/D");
    fEventTree->Branch("dau_track_PID_pdg_RecoAlg", dau_track_PID_pdg_RecoAlg,"dau_track_PID_pdg_RecoAlg[n_reco_tracks][10][3]/I");
    fEventTree->Branch("dau_track_KE_RecoAlg", dau_track_KE_RecoAlg,"dau_track_KE_RecoAlg[n_reco_tracks][10][3]/D");
    fEventTree->Branch("dau_track_bestplane_RecoAlg", dau_track_bestplane_RecoAlg,"dau_track_bestplane_RecoAlg[n_reco_tracks][10]/I");
    fEventTree->Branch("n_dau_track_points_RecoAlg", n_dau_track_points_RecoAlg,"n_dau_track_points_RecoAlg[n_reco_tracks][10]/I");
    fEventTree->Branch("n_dau_cal_points_RecoAlg", n_dau_cal_points_RecoAlg,"n_dau_cal_points_RecoAlg[n_reco_tracks][10]/I");
    fEventTree->Branch("n_dau_cal_points_byplane_RecoAlg", n_dau_cal_points_byplane_RecoAlg,"n_dau_cal_points_byplane_RecoAlg[n_reco_tracks][10][3]/I");
    fEventTree->Branch("dau_track_dQ_dx_RecoAlg", dau_track_dQ_dx_RecoAlg,"dau_track_dQ_dx_RecoAlg[n_reco_tracks][10][500]/D");
    fEventTree->Branch("dau_track_dE_dx_RecoAlg", dau_track_dE_dx_RecoAlg,"dau_track_dE_dx_RecoAlg[n_reco_tracks][10][500]/D");
    fEventTree->Branch("dau_track_range_RecoAlg", dau_track_range_RecoAlg,"dau_track_range_RecoAlg[n_reco_tracks][10][500]/D");
    fEventTree->Branch("dau_track_pitch_RecoAlg", dau_track_pitch_RecoAlg,"dau_track_pitch_RecoAlg[n_reco_tracks][10][500]/D");
    fEventTree->Branch("dau_track_dQ_dx_byplane_RecoAlg", dau_track_dQ_dx_byplane_RecoAlg,"dau_track_dQ_dx_byplane_RecoAlg[n_reco_tracks][10][3][500]/D");
    fEventTree->Branch("dau_track_dE_dx_byplane_RecoAlg", dau_track_dE_dx_byplane_RecoAlg,"dau_track_dE_dx_byplane_RecoAlg[n_reco_tracks][10][3][500]/D");
    fEventTree->Branch("dau_track_range_byplane_RecoAlg", dau_track_range_byplane_RecoAlg,"dau_track_range_byplane_RecoAlg[n_reco_tracks][10][3][500]/D");
    fEventTree->Branch("dau_track_pitch_byplane_RecoAlg", dau_track_pitch_byplane_RecoAlg,"dau_track_pitch_byplane_RecoAlg[n_reco_tracks][10][3][500]/D");
    fEventTree->Branch("dau_track_calo_xyz_byplane_RecoAlg", dau_track_calo_xyz_byplane_RecoAlg,"dau_track_calo_xyz_byplane_RecoAlg[n_reco_tracks][10][3][500][3]/D");
    fEventTree->Branch("dau_track_chi_pr_RecoAlg", dau_track_chi_pr_RecoAlg, "dau_track_chi_pr_RecoAlg[n_reco_tracks][10][3]/D");
    fEventTree->Branch("dau_track_chi_ka_RecoAlg", dau_track_chi_ka_RecoAlg, "dau_track_chi_ka_RecoAlg[n_reco_tracks][10][3]/D");
    fEventTree->Branch("dau_track_chi_pi_RecoAlg", dau_track_chi_pi_RecoAlg, "dau_track_chi_pi_RecoAlg[n_reco_tracks][10][3]/D");
    fEventTree->Branch("dau_track_chi_mu_RecoAlg", dau_track_chi_mu_RecoAlg, "dau_track_chi_mu_RecoAlg[n_reco_tracks][10][3]/D");
    fEventTree->Branch("dau_track_pida_RecoAlg", dau_track_pida_RecoAlg, "dau_track_pida_RecoAlg[n_reco_tracks][10][3]/D");
    fEventTree->Branch("dau_track_pidndf_RecoAlg", dau_track_pidndf_RecoAlg, "dau_track_pidndf_RecoAlg[n_reco_tracks][10][3]/D");

    fEventTree->Branch("track_complet", track_complet,"track_complet[n_reco_tracks]/D");  //track quality variable (completeness)
    fEventTree->Branch("track_Efrac", track_Efrac,"track_Efrac[n_reco_tracks]/D");        //track quality variable (purity)
    fEventTree->Branch("track_mcID", track_mcID,"track_mcID[n_reco_tracks]/I");           //true MC ID for a given track
    fEventTree->Branch("track_mcPDG", track_mcPDG,"track_mcPDG[n_reco_tracks]/I");        //true MC PDG for a given track

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
  void RunningAnalyzer::endJob(){
  }
  //========================================================================
  void RunningAnalyzer::beginRun(const art::Run& /*run*/){
    mf::LogInfo("RunningAnalyzer")<<"begin run..."<<endl;
  }
  //========================================================================
  void RunningAnalyzer::analyze( const art::Event& event ){
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
  void RunningAnalyzer::Process( const art::Event& event, bool &isFiducial){

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

    // get tracks with recoalg
    art::Handle< std::vector<recob::Track> > trackRecoAlgListHandle;
    std::vector<art::Ptr<recob::Track>> trackRecoAlglist;
    if( event.getByLabel(fTrackRecoAlgModuleLabel, trackRecoAlgListHandle))
      art::fill_ptr_vector(trackRecoAlglist, trackRecoAlgListHandle);
    n_recoTracks_RecoAlg = trackRecoAlglist.size();
    if( n_recoTracks_RecoAlg > MAX_TRACKS || n_recoTracks_RecoAlg == 0) {
      n_recoTracks_RecoAlg = 0; // make sure we don't save any fake data
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

    // get hits and calo of tracks with recoclg
    art::FindManyP<recob::Hit> track_hits_RecoAlg(trackRecoAlgListHandle, event, fHitTrackRecoAlgAssns);
    art::FindMany<anab::Calorimetry>  reco_cal_RecoAlg(trackRecoAlgListHandle, event, fCaloRecoAlgModuleLabel);
    trkf::TrackMomentumCalculator trackPRecoAlg;
    art::FindMany<anab::ParticleID> reco_PID_RecoAlg(trackRecoAlgListHandle, event, fPidRecoAlgModuleLabel);

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


    // loop reco tracks
    for(int i=0; i<n_recoTracks; ++i) {

      n_recoDauTracks[i] = 0;
      n_recoDauTracks_RecoAlg[i] = 0;

      art::Ptr<recob::Track> track = tracklist[i];

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

      TVector3 end(track->End().x(),
		   track->End().y(),
		   track->End().z());

      //cout << "track_vtx: "  << track_vtx[i][0] << " " << track_vtx[i][1] << " " << track_vtx[i][2]  << endl;
      //cout << "track_end: " << track_end[i][0] << " " << track_end[i][1] << " " << track_end[i][2] << endl;

      //truth matcher
      
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
      //cout << "track_mcPDG[i]: " << track_mcPDG[i] << endl;
      //if(track_mcPDG[i]!=321) continue;

      for(int j=0; j<n_recoTracks; ++j) {
	
	art::Ptr<recob::Track> dau_track = tracklist[j];
	
	if(dau_track->ID() == track->ID()) continue;
	
	double track_dau_distance = TMath::Sqrt((track->End().x()-dau_track->Vertex().x())*(track->End().x()-dau_track->Vertex().x()) +
						(track->End().y()-dau_track->Vertex().y())*(track->End().y()-dau_track->Vertex().y()) +
						(track->End().z()-dau_track->Vertex().z())*(track->End().z()-dau_track->Vertex().z()));

	if(track_dau_distance<15){

	  dau_track_length[i][n_recoDauTracks[i]] = dau_track->Length();
	  //dau_track_mcPDG[i][n_recoDauTracks[i]] = track_mcPDG[n_recoDauTracks[i]];
	  dau_track_distance[i][n_recoDauTracks[i]] = track_dau_distance;

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
	    dau_track_pdg[i][n_recoDauTracks[i]] = v[0].second;
	  }

	  //calculate PID
	  std::vector<const anab::Calorimetry*> trk_cal = reco_cal.at(j);
	  std::vector<const anab::ParticleID*> trk_pid = reco_PID.at(j);
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
	      dau_track_PID_pdg[i][n_recoDauTracks[i]][ipid] = pScore.fAssumedPdg;

	      for(auto pScore: pidScore){	    
		
		double chi2value = pScore.fValue;
		
		// PIDA is always the last one and ndf there is -9999  
		if(pScore.fAssumedPdg != 0){
		  dau_track_pidndf[i][n_recoDauTracks[i]][planenum] = pScore.fNdf; // This value is the same for each particle type, but different in each plane
		}
		switch(pScore.fAssumedPdg){
		case 2212:
		  dau_track_chi_pr[i][n_recoDauTracks[i]][planenum] = chi2value;
		  break;
		case 321:
		  dau_track_chi_ka[i][n_recoDauTracks[i]][planenum] = chi2value;
		  break;
		case 211:
		  dau_track_chi_pi[i][n_recoDauTracks[i]][planenum] = chi2value;
		  break;
		case 13:
		  dau_track_chi_mu[i][n_recoDauTracks[i]][planenum] = chi2value;
		  break;
		case 0:
		  dau_track_pida[i][n_recoDauTracks[i]][planenum] = chi2value; 
		  break;      
		}
	      }	  

	    }
	  }

	  int best_plane =-1;
	  int most_ndf = std::max({plane0, plane1, plane2});
	  //for( size_t p =0; p<3; p++) if( most_ndf == trk_pid[p]->Ndf() ) best_plane = p;
	  if( most_ndf == plane0 ) best_plane = 0; 
	  if( most_ndf == plane1 ) best_plane = 1; 
	  if( most_ndf == plane2 ) best_plane = 2; 
	  //cout << "best_plane is " << best_plane << endl;
	  
	  dau_track_KE[i][n_recoDauTracks[i]][0] = trk_cal[0]->KineticEnergy();
	  dau_track_KE[i][n_recoDauTracks[i]][1] = trk_cal[1]->KineticEnergy();
	  dau_track_KE[i][n_recoDauTracks[i]][2] = trk_cal[2]->KineticEnergy();
	  
	  std::vector<double> PIDAval;
	  std::vector<double> chi2;
	  PIDAcal( trk_cal, PIDAval);
	  int idx=0;
	  for(auto const& val : PIDAval){
	    dau_track_PIDA[i][n_recoDauTracks[i]][idx] = val;
	    idx++;
	  }
	  dau_track_bestplane[i][n_recoDauTracks[i]] = best_plane;
	  
	  //save dE/dx & dQ/dx
	  for (unsigned iplane = 0; iplane < 3; iplane++ ) {
	    int npts = trk_cal[iplane]->dEdx().size();
	    if (npts > MAX_CALO_PTS)
	      npts = MAX_CALO_PTS;
	    if (best_plane >= 0 && iplane == (unsigned)best_plane) {
	      n_dau_cal_points[i][n_recoDauTracks[i]] = npts;
	    }

	    n_dau_cal_points_byplane[i][n_recoDauTracks[i]][iplane] = npts;
	    
	    const vector<float> &dqdx = trk_cal[iplane]->dQdx();
	    const vector<float> &dedx = trk_cal[iplane]->dEdx();
	    const vector<float> &resr = trk_cal[iplane]->ResidualRange();
	    const vector<float> &pitch = trk_cal[iplane]->TrkPitchVec();
	    //auto &xyz = trk_cal[iplane]->XYZ();
	    
	    size_t n_dqdx  = std::min(dqdx.size(),  static_cast<size_t>(MAX_CALO_PTS));
	    size_t n_dedx  = std::min(dedx.size(),  static_cast<size_t>(MAX_CALO_PTS));
	    size_t n_resr  = std::min(resr.size(),  static_cast<size_t>(MAX_CALO_PTS));
	    size_t n_pitch = std::min(pitch.size(), static_cast<size_t>(MAX_CALO_PTS));

	    if (best_plane >= 0 && iplane == (unsigned)best_plane) {
	      /*
	      std::copy(dqdx.begin(), dqdx.end(), dau_track_dQ_dx[i][n_recoDauTracks[i]]);
	      std::copy(dedx.begin(), dedx.end(), dau_track_dE_dx[i][n_recoDauTracks[i]]);
	      std::copy(resr.begin(), resr.end(), dau_track_range[i][n_recoDauTracks[i]]);
	      std::copy(pitch.begin(), pitch.end(), dau_track_pitch[i][n_recoDauTracks[i]]);
	      */

	      std::copy(dqdx.begin(),  dqdx.begin()  + n_dqdx,  dau_track_dQ_dx[i][n_recoDauTracks[i]]);
	      std::copy(dedx.begin(),  dedx.begin()  + n_dedx,  dau_track_dE_dx[i][n_recoDauTracks[i]]);
	      std::copy(resr.begin(),  resr.begin()  + n_resr,  dau_track_range[i][n_recoDauTracks[i]]);
	      std::copy(pitch.begin(), pitch.begin() + n_pitch, dau_track_pitch[i][n_recoDauTracks[i]]);
	    }
	    /*
	    std::copy(dqdx.begin(), dqdx.end(), dau_track_dQ_dx_byplane[i][n_recoDauTracks[i]][iplane]);
	    std::copy(dedx.begin(), dedx.end(), dau_track_dE_dx_byplane[i][n_recoDauTracks[i]][iplane]);
	    std::copy(resr.begin(), resr.end(), dau_track_range_byplane[i][n_recoDauTracks[i]][iplane]);
	    std::copy(pitch.begin(), pitch.end(), dau_track_pitch_byplane[i][n_recoDauTracks[i]][iplane]);
	    */
	    std::copy(dqdx.begin(),  dqdx.begin()  + n_dqdx,  dau_track_dQ_dx_byplane[i][n_recoDauTracks[i]][iplane]);
	    std::copy(dedx.begin(),  dedx.begin()  + n_dedx,  dau_track_dE_dx_byplane[i][n_recoDauTracks[i]][iplane]);
	    std::copy(resr.begin(),  resr.begin()  + n_resr,  dau_track_range_byplane[i][n_recoDauTracks[i]][iplane]);
	    std::copy(pitch.begin(), pitch.begin() + n_pitch, dau_track_pitch_byplane[i][n_recoDauTracks[i]][iplane]);
	    
	    // save calo point's XYZ coords
	    /*
	    double* coords = (double*)track_calo_xyz_byplane[n_recoDauTracks[i]][iplane];
	    for(int k = 0; k < npts; j++) {
	      (coords+k*3)[0] = xyz[k].X();
	      (coords+k*3)[1] = xyz[k].Y();
	      (coords+k*3)[2] = xyz[k].Z();
	    }
	    */
	  }

	  //truth matcher
      
	  double dau_tmpEfrac = 0;
	  const simb::MCParticle *dau_particle;
	  double dau_tmpComplet = 0;
	  std::vector<art::Ptr<recob::Hit>> dau_all_trackHits = track_hits.at(i);
	  truthMatcher( all_hits,  dau_all_trackHits, dau_particle, dau_tmpEfrac, dau_tmpComplet );
	  if(!dau_particle) continue;
	  dau_track_mcID[i][n_recoDauTracks[i]] = dau_particle->TrackId();
	  dau_track_mcPDG[i][n_recoDauTracks[i]] = dau_particle->PdgCode();
	  dau_track_Efrac[i][n_recoDauTracks[i]] = dau_tmpEfrac;
	  dau_track_complet[i][n_recoDauTracks[i]] = dau_tmpComplet;
	  
	  //dau_track_mcPDG[i][n_recoDauTracks[i]] = dau_track_mcPDG_[j];

	  n_recoDauTracks[i]++;
	}
      }     


      // loop for tracks reconstrcted by the algorithm
      std::cout << "n_recoTracks_RecoAlg: " << n_recoTracks_RecoAlg << std::endl;
      for(int j=0; j<n_recoTracks_RecoAlg; ++j) {
	
	art::Ptr<recob::Track> dau_track = trackRecoAlglist[j];
	
	if(dau_track->ID() == track->ID()) continue;
	
	double track_dau_distance = TMath::Sqrt((track->End().x()-dau_track->Vertex().x())*(track->End().x()-dau_track->Vertex().x()) +
						(track->End().y()-dau_track->Vertex().y())*(track->End().y()-dau_track->Vertex().y()) +
						(track->End().z()-dau_track->Vertex().z())*(track->End().z()-dau_track->Vertex().z()));
	
	if(track_dau_distance<15){
	  
	  std::cout << "dau_track_length_RecoAlg: " << dau_track->Length() << std::endl;
	  dau_track_length_RecoAlg[i][n_recoDauTracks_RecoAlg[i]] = dau_track->Length();
	  dau_track_distance_RecoAlg[i][n_recoDauTracks_RecoAlg[i]] = track_dau_distance;

	  const simb::MCParticle *mcparticle_dau;
	  std::map<int,int> hits_pdg_map_dau;
	  for(auto const& hit : track_hits_RecoAlg.at(j)){
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
	    dau_track_pdg_RecoAlg[i][n_recoDauTracks_RecoAlg[i]] = v[0].second;
	  }

	  //calculate PID
	  std::vector<const anab::Calorimetry*> trk_cal_RecoAlg = reco_cal_RecoAlg.at(j);
	  std::vector<const anab::ParticleID*> trk_pid_RecoAlg = reco_PID_RecoAlg.at(j);
	  //std::vector AlgScoresVec = trk_pid.at(0)->ParticleIDAlgScores();
	  
	  int plane0=0; 
	  int plane1=0;
	  int plane2=0;

	  for (size_t ipid = 0; ipid < trk_pid_RecoAlg.size(); ++ipid){
	    if (!trk_pid_RecoAlg[ipid]->PlaneID().isValid) continue;
	    int planenum = trk_pid_RecoAlg[ipid]->PlaneID().Plane;
	    if (planenum<0||planenum>2) continue;
	    
	    auto pidScore = trk_pid_RecoAlg[ipid]->ParticleIDAlgScores();
	    for(auto pScore: pidScore){
	      //double chi2value = pScore.fValue;
	      
	      if(pScore.fAssumedPdg != 0){
		if(planenum==0) plane0 = pScore.fNdf;
		if(planenum==1) plane1 = pScore.fNdf;
		if(planenum==2) plane2 = pScore.fNdf;
	      }
	      dau_track_PID_pdg_RecoAlg[i][n_recoDauTracks_RecoAlg[i]][ipid] = pScore.fAssumedPdg;

	      for(auto pScore: pidScore){	    
		
		double chi2value = pScore.fValue;
		
		// PIDA is always the last one and ndf there is -9999  
		if(pScore.fAssumedPdg != 0){
		  dau_track_pidndf_RecoAlg[i][n_recoDauTracks_RecoAlg[i]][planenum] = pScore.fNdf; // This value is the same for each particle type, but different in each plane
		}
		switch(pScore.fAssumedPdg){
		case 2212:
		  dau_track_chi_pr_RecoAlg[i][n_recoDauTracks_RecoAlg[i]][planenum] = chi2value;
		  break;
		case 321:
		  dau_track_chi_ka_RecoAlg[i][n_recoDauTracks_RecoAlg[i]][planenum] = chi2value;
		  break;
		case 211:
		  dau_track_chi_pi_RecoAlg[i][n_recoDauTracks_RecoAlg[i]][planenum] = chi2value;
		  break;
		case 13:
		  dau_track_chi_mu_RecoAlg[i][n_recoDauTracks_RecoAlg[i]][planenum] = chi2value;
		  break;
		case 0:
		  dau_track_pida_RecoAlg[i][n_recoDauTracks_RecoAlg[i]][planenum] = chi2value; 
		  break;      
		}
	      }	  

	    }
	  }

	  int best_plane =-1;
	  int most_ndf = std::max({plane0, plane1, plane2});
	  //for( size_t p =0; p<3; p++) if( most_ndf == trk_pid[p]->Ndf() ) best_plane = p;
	  if( most_ndf == plane0 ) best_plane = 0; 
	  if( most_ndf == plane1 ) best_plane = 1; 
	  if( most_ndf == plane2 ) best_plane = 2; 
	  //cout << "best_plane is " << best_plane << endl;
	  
	  dau_track_KE_RecoAlg[i][n_recoDauTracks_RecoAlg[i]][0] = trk_cal_RecoAlg[0]->KineticEnergy();
	  dau_track_KE_RecoAlg[i][n_recoDauTracks_RecoAlg[i]][1] = trk_cal_RecoAlg[1]->KineticEnergy();
	  dau_track_KE_RecoAlg[i][n_recoDauTracks_RecoAlg[i]][2] = trk_cal_RecoAlg[2]->KineticEnergy();
	  
	  std::vector<double> PIDAval;
	  std::vector<double> chi2;
	  PIDAcal( trk_cal_RecoAlg, PIDAval);
	  int idx=0;
	  for(auto const& val : PIDAval){
	    dau_track_PIDA_RecoAlg[i][n_recoDauTracks_RecoAlg[i]][idx] = val;
	    idx++;
	  }
	  dau_track_bestplane_RecoAlg[i][n_recoDauTracks_RecoAlg[i]] = best_plane;
	  
	  //save dE/dx & dQ/dx
	  for (unsigned iplane = 0; iplane < 3; iplane++ ) {
	    int npts = trk_cal_RecoAlg[iplane]->dEdx().size();
	    if (npts > MAX_CALO_PTS)
	      npts = MAX_CALO_PTS;
	    if (best_plane >= 0 && iplane == (unsigned)best_plane) {
	      n_dau_cal_points_RecoAlg[i][n_recoDauTracks_RecoAlg[i]] = npts;
	    }

	    n_dau_cal_points_byplane_RecoAlg[i][n_recoDauTracks_RecoAlg[i]][iplane] = npts;
	    
	    const vector<float> &dqdx = trk_cal_RecoAlg[iplane]->dQdx();
	    const vector<float> &dedx = trk_cal_RecoAlg[iplane]->dEdx();
	    const vector<float> &resr = trk_cal_RecoAlg[iplane]->ResidualRange();
	    const vector<float> &pitch = trk_cal_RecoAlg[iplane]->TrkPitchVec();
	    //auto &xyz = trk_cal_RecoAlg[iplane]->XYZ();
	    
	    size_t n_dqdx  = std::min(dqdx.size(),  static_cast<size_t>(MAX_CALO_PTS));
	    size_t n_dedx  = std::min(dedx.size(),  static_cast<size_t>(MAX_CALO_PTS));
	    size_t n_resr  = std::min(resr.size(),  static_cast<size_t>(MAX_CALO_PTS));
	    size_t n_pitch = std::min(pitch.size(), static_cast<size_t>(MAX_CALO_PTS));
	    
	    if (best_plane >= 0 && iplane == (unsigned)best_plane) {
	      /*
	      std::copy(dqdx.begin(), dqdx.end(), dau_track_dQ_dx_RecoAlg[i][n_recoDauTracks_RecoAlg[i]]);
	      std::copy(dedx.begin(), dedx.end(), dau_track_dE_dx_RecoAlg[i][n_recoDauTracks_RecoAlg[i]]);
	      std::copy(resr.begin(), resr.end(), dau_track_range_RecoAlg[i][n_recoDauTracks_RecoAlg[i]]);
	      std::copy(pitch.begin(), pitch.end(), dau_track_pitch_RecoAlg[i][n_recoDauTracks_RecoAlg[i]]);
	      */

	      std::copy(dqdx.begin(),  dqdx.begin()  + n_dqdx,  dau_track_dQ_dx_RecoAlg[i][n_recoDauTracks_RecoAlg[i]]);
	      std::copy(dedx.begin(),  dedx.begin()  + n_dedx,  dau_track_dE_dx_RecoAlg[i][n_recoDauTracks_RecoAlg[i]]);
	      std::copy(resr.begin(),  resr.begin()  + n_resr,  dau_track_range_RecoAlg[i][n_recoDauTracks_RecoAlg[i]]);
	      std::copy(pitch.begin(), pitch.begin() + n_pitch, dau_track_pitch_RecoAlg[i][n_recoDauTracks_RecoAlg[i]]);
	    }
	    /*
	    std::copy(dqdx.begin(), dqdx.end(), dau_track_dQ_dx_byplane_RecoAlg[i][n_recoDauTracks_RecoAlg[i]][iplane]);
	    std::copy(dedx.begin(), dedx.end(), dau_track_dE_dx_byplane_RecoAlg[i][n_recoDauTracks_RecoAlg[i]][iplane]);
	    std::copy(resr.begin(), resr.end(), dau_track_range_byplane_RecoAlg[i][n_recoDauTracks_RecoAlg[i]][iplane]);
	    std::copy(pitch.begin(), pitch.end(), dau_track_pitch_byplane_RecoAlg[i][n_recoDauTracks_RecoAlg[i]][iplane]);
	    */
	    std::copy(dqdx.begin(),  dqdx.begin()  + n_dqdx,  dau_track_dQ_dx_byplane_RecoAlg[i][n_recoDauTracks_RecoAlg[i]][iplane]);
	    std::copy(dedx.begin(),  dedx.begin()  + n_dedx,  dau_track_dE_dx_byplane_RecoAlg[i][n_recoDauTracks_RecoAlg[i]][iplane]);
	    std::copy(resr.begin(),  resr.begin()  + n_resr,  dau_track_range_byplane_RecoAlg[i][n_recoDauTracks_RecoAlg[i]][iplane]);
	    std::copy(pitch.begin(), pitch.begin() + n_pitch, dau_track_pitch_byplane_RecoAlg[i][n_recoDauTracks_RecoAlg[i]][iplane]);
	    
	    // save calo point's XYZ coords
	    /*
	    double* coords = (double*)track_calo_xyz_byplane_RecoAlg[n_recoDauTracks_RecoAlg[i]][iplane];
	    for(int k = 0; k < npts; j++) {
	      (coords+k*3)[0] = xyz[k].X();
	      (coords+k*3)[1] = xyz[k].Y();
	      (coords+k*3)[2] = xyz[k].Z();
	    }
	    */
	  }

	  //truth matcher
      
	  double dau_tmpEfrac_RecoAlg = 0;
	  const simb::MCParticle *dau_particle_RecoAlg;
	  double dau_tmpComplet_RecoAlg = 0;
	  std::vector<art::Ptr<recob::Hit>> dau_all_trackHits_RecoAlg = track_hits_RecoAlg.at(i);
	  truthMatcher( all_hits,  dau_all_trackHits_RecoAlg, dau_particle_RecoAlg, dau_tmpEfrac_RecoAlg, dau_tmpComplet_RecoAlg );
	  if(!dau_particle_RecoAlg) continue;
	  dau_track_mcID_RecoAlg[i][n_recoDauTracks_RecoAlg[i]] = dau_particle_RecoAlg->TrackId();
	  dau_track_mcPDG_RecoAlg[i][n_recoDauTracks_RecoAlg[i]] = dau_particle_RecoAlg->PdgCode();
	  dau_track_Efrac_RecoAlg[i][n_recoDauTracks_RecoAlg[i]] = dau_tmpEfrac_RecoAlg;
	  dau_track_complet_RecoAlg[i][n_recoDauTracks_RecoAlg[i]] = dau_tmpComplet_RecoAlg;
	  
	  //dau_track_mcPDG_RecoAlg[i][n_recoDauTracks_RecoAlg[i]] = dau_track_mcPDG_RecoAlg[j];

	  //n_recoDauTracks_RecoAlg[i]
	  n_recoDauTracks_RecoAlg[i]++;
	}
	
      }       
      
      // add track points
      n_track_points[i] = track->NPoints();
      for (int ipt=0; ipt < n_track_points[i]; ipt++) {
	// FIXME: get the track point coordinates!
	recob::tracking::Point_t pos_point = track->TrajectoryPoint(ipt).position;
	TVector3 pos(pos_point.x(),pos_point.y(),pos_point.z());
	//pos.GetXYZ(track_point_xyz[i][ipt]);
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
	  	  
	  for(auto pScore: pidScore){	    

	    double chi2value = pScore.fValue;
	    
	    // PIDA is always the last one and ndf there is -9999  
	    if(pScore.fAssumedPdg != 0){
	      track_pidndf[i][planenum] = pScore.fNdf; // This value is the same for each particle type, but different in each plane
	    }
	    switch(pScore.fAssumedPdg){
	    case 2212:
	      track_chi_pr[i][planenum] = chi2value;
	      break;
	    case 321:
	      track_chi_ka[i][planenum] = chi2value;
	      break;
	    case 211:
	      track_chi_pi[i][planenum] = chi2value;
	      break;
	    case 13:
	      track_chi_mu[i][planenum] = chi2value;
	      break;
	    case 0:
	      track_pida[i][planenum] = chi2value; 
	      break;      
	    }
	  }	  
	  
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
	//auto &xyz = trk_cal[iplane]->XYZ();
	
	size_t n_dqdx  = std::min(dqdx.size(),  static_cast<size_t>(MAX_CALO_PTS));
	size_t n_dedx  = std::min(dedx.size(),  static_cast<size_t>(MAX_CALO_PTS));
	size_t n_resr  = std::min(resr.size(),  static_cast<size_t>(MAX_CALO_PTS));
	size_t n_pitch = std::min(pitch.size(), static_cast<size_t>(MAX_CALO_PTS));

	if (best_plane >= 0 && iplane == (unsigned)best_plane) {
	  /*
	  std::copy(dqdx.begin(), dqdx.end(), track_dQ_dx[i]);
	  std::copy(dedx.begin(), dedx.end(), track_dE_dx[i]);
	  std::copy(resr.begin(), resr.end(), track_range[i]);
	  std::copy(pitch.begin(), pitch.end(), track_pitch[i]);
	  */

	  std::copy(dqdx.begin(),  dqdx.begin()  + n_dqdx,  track_dQ_dx[i]);
	  std::copy(dedx.begin(),  dedx.begin()  + n_dedx,  track_dE_dx[i]);
	  std::copy(resr.begin(),  resr.begin()  + n_resr,  track_range[i]);
	  std::copy(pitch.begin(), pitch.begin() + n_pitch, track_pitch[i]);
	}
	/*
	std::copy(dqdx.begin(), dqdx.end(), track_dQ_dx_byplane[i][iplane]);
	std::copy(dedx.begin(), dedx.end(), track_dE_dx_byplane[i][iplane]);
	std::copy(resr.begin(), resr.end(), track_range_byplane[i][iplane]);
	std::copy(pitch.begin(), pitch.end(), track_pitch_byplane[i][iplane]);
	*/
	std::copy(dqdx.begin(),  dqdx.begin()  + n_dqdx,  track_dQ_dx_byplane[i][iplane]);
	std::copy(dedx.begin(),  dedx.begin()  + n_dedx,  track_dE_dx_byplane[i][iplane]);
	std::copy(resr.begin(),  resr.begin()  + n_resr,  track_range_byplane[i][iplane]);
	std::copy(pitch.begin(), pitch.begin() + n_pitch, track_pitch_byplane[i][iplane]);

	// save calo point's XYZ coords
	/*
	double* coords = (double*)track_calo_xyz_byplane[i][iplane];
	for(int j = 0; j < npts; j++) {
	  (coords+j*3)[0] = xyz[j].X();
	  (coords+j*3)[1] = xyz[j].Y();
	  (coords+j*3)[2] = xyz[j].Z();
	}
	*/

      }     
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
  void RunningAnalyzer::PIDAcal( std::vector<const anab::Calorimetry*> cal, std::vector<double> &PIDA){
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
  void RunningAnalyzer::truthHitMatcher( art::Ptr<recob::Hit> hit, const simb::MCParticle *&MCparticle){

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
    std::unordered_map<int,double> TrackIDE;
    // amount of energy deposited by the particle that deposited more energy... tomato potato... blabla
    //!if the collection of hits have more than one particle associate save the particle w/ the highest energy deposition
    //!since we are looking for muons/pions/protons this should be enough
    if( !trkID_E.size() ) {
      MCparticle = 0;
      return; //Ghost track???
    }
    for(std::map<int,double>::iterator ii = trkID_E.begin(); ii!=trkID_E.end(); ++ii){
      total_E += ii->second;
      TrackIDE[ii->first] += ii->second;

      if( TrackIDE[ii->first] > max_E ){
	max_E = TrackIDE[ii->first];
	TrackID = ii->first;
	if( TrackID < 0 ) E_em += ii->second;
      }
      /*
      if((ii->second)>max_E){
	max_E = ii->second;
	TrackID = ii->first;
	if( TrackID < 0 ) E_em += ii->second;
      }
      */

    }
    // consider the most energetic hit to estimate PDG etc
    MCparticle = part_inv->TrackIdToParticle_P(TrackID);
  }


  //========================================================================
  void RunningAnalyzer::truthMatcher( std::vector<art::Ptr<recob::Hit>> all_hits, std::vector<art::Ptr<recob::Hit>> track_hits, const simb::MCParticle *&MCparticle, double &Efrac, double &Ecomplet){

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

  double RunningAnalyzer::truthLength( const simb::MCParticle *MCparticle ){
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
  bool RunningAnalyzer::insideFV( double vertex[4]){

    double x = vertex[0];
    double y = vertex[1];
    double z = vertex[2];

    if (x>fFidVolXmin && x<fFidVolXmax&&
	y>fFidVolYmin && y<fFidVolYmax&&
	z>fFidVolZmin && z<fFidVolZmax)
      return true;
    else
      return false;
  }

  void RunningAnalyzer::reset() {
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
    std::memset( track_isContained, 0, sizeof(track_isContained) );
    std::memset( track_ID, 0, sizeof(track_ID) );
    std::memset( vtxID_trk, 0, sizeof(vtxID_trk) );
    std::memset( track_vtxDir, 0, sizeof(track_vtxDir) );
    std::memset( track_vtx, 0, sizeof(track_vtx) );
    std::memset( track_end, 0, sizeof(track_end) );
    std::memset( track_length, 0, sizeof(track_length) );
    std::memset( dau_track_length, 0, sizeof(dau_track_length) );
    std::memset( dau_track_distance, 0, sizeof(dau_track_distance) );
    std::memset( dau_track_pdg, 0, sizeof(dau_track_pdg) );

    n_recoTracks_RecoAlg = 0;
    std::memset( n_recoDauTracks_RecoAlg, 0, sizeof(n_recoDauTracks_RecoAlg) );
    std::memset( dau_track_length_RecoAlg, 0, sizeof(dau_track_length_RecoAlg) );
    std::memset( dau_track_distance_RecoAlg, 0, sizeof(dau_track_distance_RecoAlg) );
    std::memset( dau_track_pdg_RecoAlg, 0, sizeof(dau_track_pdg_RecoAlg) );
    std::memset( dau_track_mcPDG_RecoAlg, 0, sizeof(dau_track_mcPDG_RecoAlg) );

    std::memset( n_recoRebDauTracks, 0, sizeof(n_recoRebDauTracks) );
    std::memset( rebdautracktrue_length, 0, sizeof(rebdautracktrue_length) );
    std::memset( rebdautracktruedir_length, 0, sizeof(rebdautracktruedir_length) );
    std::memset( rebdautrack_length, 0, sizeof(rebdautrack_length) );
    std::memset( rebdautrack_distance, 0, sizeof(rebdautrack_distance) );
    std::memset( rebdautrack_pdg, 0, sizeof(rebdautrack_pdg) );

    std::memset( best_peak_x, 0, sizeof(best_peak_x) );
    std::memset( best_peak_y, 0, sizeof(best_peak_y) );
    std::memset( best_peak_z, 0, sizeof(best_peak_z) );

    std::memset( best_peak_x_true, 0, sizeof(best_peak_x_true) );
    std::memset( best_peak_y_true, 0, sizeof(best_peak_y_true) );
    std::memset( best_peak_z_true, 0, sizeof(best_peak_z_true) );

    std::memset( dau_track_mcID, 0, sizeof(dau_track_mcID) );  
    std::memset( dau_track_mcPDG, 0, sizeof(dau_track_mcPDG) );
    std::memset( dau_track_Efrac, 0, sizeof(dau_track_Efrac) ); 
    std::memset( dau_track_complet, 0, sizeof(dau_track_complet) );
    std::memset( dau_track_PIDA, 0, sizeof(dau_track_PIDA) );
    std::memset( dau_track_PID_pdg, 0, sizeof(dau_track_PID_pdg) );
    std::memset( dau_track_KE, 0, sizeof(dau_track_KE) );
    std::memset( dau_track_bestplane, 0, sizeof(dau_track_bestplane) ); 
    std::memset( n_dau_track_points, 0, sizeof(n_dau_track_points) ); 
    std::memset( n_dau_cal_points, 0, sizeof(n_dau_cal_points) );
    std::memset( n_dau_cal_points_byplane, 0, sizeof(n_dau_cal_points_byplane) );
    std::memset( dau_track_dQ_dx, 0, sizeof(dau_track_dQ_dx) );
    std::memset( dau_track_dE_dx, 0, sizeof(dau_track_dE_dx) );
    std::memset( dau_track_range, 0, sizeof(dau_track_range) );
    std::memset( dau_track_pitch, 0, sizeof(dau_track_pitch) );
    std::memset( dau_track_chi_pr, 0, sizeof(dau_track_chi_pr) );
    std::memset( dau_track_chi_ka, 0, sizeof(dau_track_chi_ka) );
    std::memset( dau_track_chi_pi, 0, sizeof(dau_track_chi_pi) );
    std::memset( dau_track_chi_mu, 0, sizeof(dau_track_chi_mu) );
    std::memset( dau_track_pida, 0, sizeof(dau_track_pida) );
    std::memset( dau_track_pidndf, 0, sizeof(dau_track_pidndf) ); 
    std::memset( dau_track_dQ_dx_byplane, 0, sizeof(dau_track_dQ_dx_byplane) );
    std::memset( dau_track_dE_dx_byplane, 0, sizeof(dau_track_dE_dx_byplane) );
    std::memset( dau_track_range_byplane, 0, sizeof(dau_track_range_byplane) );
    std::memset( dau_track_pitch_byplane, 0, sizeof(dau_track_pitch_byplane) );

    std::memset( dau_track_mcID_RecoAlg, 0, sizeof(dau_track_mcID_RecoAlg) );  
    std::memset( dau_track_mcPDG_RecoAlg, 0, sizeof(dau_track_mcPDG_RecoAlg) );
    std::memset( dau_track_Efrac_RecoAlg, 0, sizeof(dau_track_Efrac_RecoAlg) ); 
    std::memset( dau_track_complet_RecoAlg, 0, sizeof(dau_track_complet_RecoAlg) );
    std::memset( dau_track_PIDA_RecoAlg, 0, sizeof(dau_track_PIDA_RecoAlg) );
    std::memset( dau_track_PID_pdg_RecoAlg, 0, sizeof(dau_track_PID_pdg_RecoAlg) );
    std::memset( dau_track_KE_RecoAlg, 0, sizeof(dau_track_KE_RecoAlg) );
    std::memset( dau_track_bestplane_RecoAlg, 0, sizeof(dau_track_bestplane_RecoAlg) ); 
    std::memset( n_dau_track_points_RecoAlg, 0, sizeof(n_dau_track_points_RecoAlg) ); 
    std::memset( n_dau_cal_points_RecoAlg, 0, sizeof(n_dau_cal_points_RecoAlg) );
    std::memset( n_dau_cal_points_byplane_RecoAlg, 0, sizeof(n_dau_cal_points_byplane_RecoAlg) );
    std::memset( dau_track_dQ_dx_RecoAlg, 0, sizeof(dau_track_dQ_dx_RecoAlg) );
    std::memset( dau_track_dE_dx_RecoAlg, 0, sizeof(dau_track_dE_dx_RecoAlg) );
    std::memset( dau_track_range_RecoAlg, 0, sizeof(dau_track_range_RecoAlg) );
    std::memset( dau_track_pitch_RecoAlg, 0, sizeof(dau_track_pitch_RecoAlg) );
    std::memset( dau_track_chi_pr_RecoAlg, 0, sizeof(dau_track_chi_pr_RecoAlg) );
    std::memset( dau_track_chi_ka_RecoAlg, 0, sizeof(dau_track_chi_ka_RecoAlg) );
    std::memset( dau_track_chi_pi_RecoAlg, 0, sizeof(dau_track_chi_pi_RecoAlg) );
    std::memset( dau_track_chi_mu_RecoAlg, 0, sizeof(dau_track_chi_mu_RecoAlg) );
    std::memset( dau_track_pida_RecoAlg, 0, sizeof(dau_track_pida_RecoAlg) );
    std::memset( dau_track_pidndf_RecoAlg, 0, sizeof(dau_track_pidndf) );
    std::memset( dau_track_dQ_dx_byplane_RecoAlg, 0, sizeof(dau_track_dQ_dx_byplane_RecoAlg) );
    std::memset( dau_track_dE_dx_byplane_RecoAlg, 0, sizeof(dau_track_dE_dx_byplane_RecoAlg) );
    std::memset( dau_track_range_byplane_RecoAlg, 0, sizeof(dau_track_range_byplane_RecoAlg) );
    std::memset( dau_track_pitch_byplane_RecoAlg, 0, sizeof(dau_track_pitch_byplane_RecoAlg) );

    //std::memset( track_dir_vtx, 0, sizeof(track_dir_vtx) );
    std::memset( track_PIDA, 0, sizeof(track_PIDA) );
    std::memset( track_PID_pdg, 0, sizeof(track_PID_pdg) );
    std::memset( track_KE, 0, sizeof(track_KE) );
    std::memset( track_bestplane, 0, sizeof(track_bestplane) );
    std::memset( track_Prange, 0, sizeof(track_Prange) );
    std::memset( track_Efrac, 0, sizeof(track_Efrac) );
    std::memset( track_complet, 0, sizeof(track_complet) );
    std::memset( track_mcID, 0, sizeof(track_mcID) );
    std::memset( track_mcPDG, 0, sizeof(track_mcPDG) );
    std::memset( dau_track_mcPDG, 0, sizeof(dau_track_mcPDG) );
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
    std::memset( track_chi_pr, 0, sizeof(track_chi_pr) );
    std::memset( track_chi_ka, 0, sizeof(track_chi_ka) );
    std::memset( track_chi_pi, 0, sizeof(track_chi_pi) );
    std::memset( track_chi_mu, 0, sizeof(track_chi_mu) );
    std::memset( track_pida, 0, sizeof(track_pida) );
    std::memset( track_pidndf, 0, sizeof(track_pidndf) ); 


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
