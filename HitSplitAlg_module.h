#ifndef HitSplitAlg_Module
#define HitSplitAlg_Module 1


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

#include "larevt/SpaceCharge/SpaceCharge.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"



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
  double 	fPIDA_endPoint;
  double	fMaxPIDAValue;
  double	fMinPIDAValue;
  bool	fSaveMCTree;

    std::vector<art::Ptr<recob::Hit>> fHits;
    std::vector<art::Ptr<recob::Track>> fTracks;
    std::vector<art::Ptr<recob::Shower>> fShowers;
    std::vector<art::Ptr<recob::SpacePoint>> fSpacePoints;
    std::map<art::Ptr<recob::Track>, std::vector<art::Ptr<recob::Hit>>> fTracksToHits; 
    std::map<art::Ptr<recob::Track>, std::vector<art::Ptr<recob::SpacePoint>>> fTracksToSpacePoints;
    std::map<art::Ptr<recob::Shower>, std::vector<art::Ptr<recob::Hit>>> fShowersToHits;
    std::map<art::Ptr<recob::Shower>, std::vector<art::Ptr<recob::SpacePoint>>> fShowersToSpacePoints;
    std::map<art::Ptr<recob::Hit>, art::Ptr<recob::SpacePoint>> fHitsToSpacePoints;
    //std::map<art::Ptr<recob::Hit>, art::Ptr<recob::SpacePoint>> fHitsToSpacePoints_recotrack;
    std::map<art::Ptr<recob::SpacePoint>, art::Ptr<recob::Hit>> fSpacePointsToHits;



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
    double best_peak_x[MAX_TRACKS][10];
    double best_peak_y[MAX_TRACKS][10];
    double best_peak_z[MAX_TRACKS][10];
    double best_peak_x_true[MAX_TRACKS];
    double best_peak_y_true[MAX_TRACKS];
    double best_peak_z_true[MAX_TRACKS];

    double track_dir_vtx[MAX_TRACKS][4];
    double track_PIDA[MAX_TRACKS][3];
    int	   track_PID_pdg[MAX_TRACKS][3];
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
