#include "CCKaonAnalyzerRebuild_module.h"

//#include "ubana/SinglePhotonAnalysis/SinglePhoton_module.h"                                                                                                                                                                                                                                                                 
#include "headers/analyze_Slice.h"
#include "headers/analyze_Tracks.h"
#include "headers/analyze_Showers.h"
#include "headers/analyze_MCTruth.h"
#include "headers/analyze_OpFlashes.h"

//#include "headers/particle_split_basetool_14Aug.h"                                                                                                                                                                                                                                                                          
#include "headers/particle_split_basetool.h"
#include "headers/track_production.h"


//#include "gallery/Event.h"                                                                                                                                                                                                                                                                                                  
//#include "gallery/ValidHandle.h"                                                                                                                                                                                                                                                                                            
//#include "headers/analyze_Template.h"                                                                                                                                                                                                                                                                                       
//#include "headers/second_shower_search.h"                                                                                                                                                                                                                                                                                   
//#include "headers/analyze_EventWeight.h"                                                                                                                                                                                                                                                                                    

/*                                                                                                                                                                                                                                                                                                                            
//#include "ubana/SinglePhotonAnalysis/fiducial_volume.h"                                                                                                                                                                                                                                                                     
 //#include "ubana/SinglePhotonAnalysis/isolation.h"                                                                                                                                                                                                                                                                           
*/

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"

#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Persistency/Common/TriggerResults.h"

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/ParameterSetRegistry.h"

#include "lardata/Utilities/AssociationUtil.h"

#include "larsim/EventWeight/Base/MCEventWeight.h"

#include "larcore/Geometry/Geometry.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include "larcoreobj/SummaryData/POTSummary.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"

#include "ubana/ParticleID/Algorithms/uB_PlaneIDBitsetHelperFunctions.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "PID/LLR_PID.h"
#include "PID/LLRPID_proton_muon_lookup.h"

#include "PID_K/LLR_PID_K.h"
#include "PID_K/LLRPID_kaon_proton_lookup.h"

//#include "larreco/RecoAlg/TrackMomentumCalculator.h"                                                                                                                                                                                                                                                                        
#include "larreco/RecoAlg/TrajectoryMCSFitter.h"
#include "ubana/ParticleID/Algorithms/uB_PlaneIDBitsetHelperFunctions.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"

#include "TTree.h"
#include "TMath.h"

#include <array>
#include <vector>
#include <map>
//#include "LinkDef.h"                                                                                                                                                                                                                                                                                                        


#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/View.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Utilities/make_tool.h"


#include "larcore/Geometry/Geometry.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/MCBase/MCStep.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/BeamInfo.h"
#include "lardataobj/RawData/TriggerData.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcoreobj/SummaryData/POTSummary.h"

#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larreco/Deprecated/BezierTrack.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "larsim/EventWeight/Base/MCEventWeight.h"
#include "ubobj/Trigger/ubdaqSoftwareTriggerData.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/FlashMatch.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "ubobj/Optical/UbooneOpticalFilter.h"
#include "ubana/AnalysisTree/MCTruth/IMCTruthMatching.h"

#include "lardataobj/RecoBase/PFParticleMetadata.h"


#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "ubana/ParticleID/Algorithms/uB_PlaneIDBitsetHelperFunctions.h"

#include "canvas/Persistency/Common/TriggerResults.h"
#include "fhiclcpp/ParameterSetRegistry.h"

#include <cstddef> // std::ptrdiff_t                                                                                                                                                                                                                                                                                          
#include <cstring> // std::memcpy()                                                                                                                                                                                                                                                                                           
#include <vector>
#include <map>
#include <iterator> // std::begin(), std::end()                                                                                                                                                                                                                                                                               
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <functional> // std::mem_fun_ref                                                                                                                                                                                                                                                                                     
#include <typeinfo>
#include <memory> // std::unique_ptr<>                                                                                                                                                                                                                                                                                        

#include "TTree.h"
#include "TTimeStamp.h"

//#include "ubana/SinglePhotonAnalysis/SinglePhoton_module.h"                                                                                                                                                                                                                                                                 

//#ifdef __MAKECINT__                                                                                                                                                                                                                                                                                                         
#ifdef __CLING__
#pragma link C++ class std::vector < std::vector<Float_t> >+;
#pragma link C++ class std::vector < std::vector< std::vector<Float_t> > >+;
#endif

//const int kMaxTracks=20;                                                                                                                                                                                                                                                                                                    

using namespace std;
namespace Kaon_Analyzer{

// ======================== Local Function Definition to get the reco origin ======================                                                                                                                                                                                                                           
  art::Ptr<simb::MCTruth>TrackIDToMCTruth(art::Event const & evt, std::string _geant_producer, int geant_track_id)
  {
    lar_pandora::MCTruthToMCParticles truthToParticles;
    lar_pandora::MCParticlesToMCTruth particlesToTruth;

    lar_pandora::LArPandoraHelper::CollectMCParticles(evt, _geant_producer, truthToParticles, particlesToTruth);

    for (auto iter : particlesToTruth) {
      if (iter.first->TrackId() == geant_track_id) {
	return iter.second;
      }
    }

    art::Ptr<simb::MCTruth> null_ptr;
    return null_ptr;
  }

  //========================================================================                                                                                                                                                                                                                                                    
  CCKaonAnalyzerRebuild::CCKaonAnalyzerRebuild(fhicl::ParameterSet const& pset) :
    EDAnalyzer(pset),
    fHitsModuleLabel          (pset.get< std::string >("HitsModuleLabel","gaushit")),
    fLArG4ModuleLabel         (pset.get< std::string >("LArG4ModuleLabel","largeant")),
    fGenieGenModuleLabel      (pset.get< std::string >("GenieGenModuleLabel","generator")),
    fTrackModuleLabel         (pset.get< std::string >("TrackModuleLabel","pandora")),
    fShowerModuleLabel        (pset.get< std::string >("ShowerModuleLabel","pandora")),
    //fCalorimetryModuleLabel   (pset.get< std::string >("CalorimetryModuleLabel","pandoracaliSCE")),                                                                                                                                                                                                                           
    fCalorimetryModuleLabel   (pset.get< std::string >("CalorimetryModuleLabel","pandoraKalmanShowercali")),
    fPIDLabel                 (pset.get< std::string >("PIDLabel","pandorapid")),
    fHitTruthAssns            (pset.get< std::string >("HitTruthAssn","gaushitTruthMatch")),
    fHitTrackAssns            (pset.get< std::string >("HitTrackAssn","pandora")),
    fHitShowerAssns            (pset.get< std::string >("HitShowerAssn","pandora")),
    m_pfp_producer            (pset.get< std::string >("pfp_producer","pandora")),
    fPFParticleLabel          (pset.get< std::string >("PFParticleLabel", "pandora")),
    fSpacePointproducer       (pset.get< std::string >("SpacePointproducer", "pandora")),
  //fSpacePointproducer  = p.get< art::InputTag >("SpacePointproducer");                                                                                                                                                                                                                                                      
  //  m_pandoraLabel            (pset.get< std::string >("PandoraLabel")),                                                                                                                                                                                                                                                      
  //  m_is_verbose              (pset.get<bool>("Verbose",false)),                                                                                                                                                                                                                                                              
    isMC                      (pset.get< bool >("IsMC",true))

    //reco_track_dEdx(nullptr),                                                                                                                                                                                                                                                                                                 
    //reco_track_ResRan(nullptr)                                                                                                                                                                                                                                                                                                
  {
    //  fm_piPFParticleLabel = pset.get<std::string>("PFParticleLabel");                                                                                                                                                                                                                                                        
    //  Kaon_Analyzer::CCKaonAnalyzerRebuild tmp;                                                                                                                                                                                                                                                                               
    //tmp.m_pandoraLabel = pset.get<std::string>("PandoraLabel");                                                                                                                                                                                                                                                               
    //tmp.m_is_verbose = pset.get<bool>("Verbose",false);                                                                                                                                                                                                                                                                       

    theDetector = lar::providerFrom<detinfo::DetectorPropertiesService>();
    detClocks   = lar::providerFrom<detinfo::DetectorClocksService>();
    SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
    geom = lar::providerFrom<geo::Geometry>();

    m_is_verbose = pset.get<bool>("Verbose",false);
    m_use_delaunay = pset.get<bool>("useDelaunay",false);
    m_is_data = pset.get<bool>("isData",false);
    m_is_overlayed = pset.get<bool>("isOverlayed",false);
    m_fill_trees = pset.get<bool>("FillTrees",true);
    m_run_pi0_filter = pset.get<bool>("RunPi0Filter",false);
    if(m_run_pi0_filter) m_is_data = true;// If running in filter mode, treat all as data                                                                                                                                                                                                                                 

    m_pandoraLabel = pset.get<std::string>("PandoraLabel");
    m_trackLabel = pset.get<std::string>("TrackLabel");
    //        m_trackLabel_old = pset.get<std::string>("TrackLabelOld");                                                                                                                                                                                                                                                  
    m_sliceLabel = pset.get<std::string>("SliceLabel","pandora");
    m_showerLabel = pset.get<std::string>("ShowerLabel");
    m_caloLabel = pset.get<std::string>("CaloLabel");
    m_flashLabel = pset.get<std::string>("FlashLabel");

    m_hitfinderLabel = pset.get<std::string>("HitFinderModule", "gaushit");
    m_badChannelLabel = pset.get<std::string>("BadChannelLabel","badmasks");
    m_badChannelProducer = pset.get<std::string>("BadChannelProducer","simnfspl1");

    m_generatorLabel = pset.get<std::string>("GeneratorLabel","generator");
    m_mcTrackLabel = pset.get<std::string>("MCTrackLabel","mcreco");
    m_mcShowerLabel = pset.get<std::string>("MCShowerLabel","mcreco");
    m_geantModuleLabel = pset.get<std::string>("GeantModule","largeant");
    m_backtrackerLabel = pset.get<std::string>("BackTrackerModule","gaushitTruthMatch");
    m_hitMCParticleAssnsLabel = pset.get<std::string>("HitMCParticleAssnLabel","gaushitTruthMatch");

    m_CRTTzeroLabel = pset.get<std::string>("CRTTzeroLabel","crttzero");
    m_runCRT = pset.get<bool>("runCRT",false);
    m_CRTHitProducer = pset.get<std::string>("CRTHitProducer", "crthitcorr");

    m_gain_mc =pset.get<std::vector<double>>("gain_mc");
    m_wire_spacing = pset.get<double>("wire_spacing");
    m_width_dqdx_box = pset.get<double>("width_box");
    m_length_dqdx_box = pset.get<double>("length_box");
    m_truthmatching_signaldef = pset.get<std::string>("truthmatching_signaldef");
    m_pidLabel = pset.get<std::string>("ParticleIDLabel","particleid");
    m_shower3dLabel = pset.get<std::string>("Shower3DLabel","shrreco3d");

    m_run_all_pfps = pset.get<bool>("runAllPFPs",false);
    rangen = new TRandom3(22);
    bool_make_sss_plots = true;

    // set dedx pdf parameters                                                                                                                                                                                                                                                                                                  
    llr_pid_calculator.set_dedx_binning(0, protonmuon_parameters.dedx_edges_pl_0);
    llr_pid_calculator.set_par_binning(0, protonmuon_parameters.parameters_edges_pl_0);
    llr_pid_calculator.set_lookup_tables(0, protonmuon_parameters.dedx_pdf_pl_0);

    llr_pid_calculator.set_dedx_binning(1, protonmuon_parameters.dedx_edges_pl_1);
    llr_pid_calculator.set_par_binning(1, protonmuon_parameters.parameters_edges_pl_1);
    llr_pid_calculator.set_lookup_tables(1, protonmuon_parameters.dedx_pdf_pl_1);

    llr_pid_calculator.set_dedx_binning(2, protonmuon_parameters.dedx_edges_pl_2);
    llr_pid_calculator.set_par_binning(2, protonmuon_parameters.parameters_edges_pl_2);
    llr_pid_calculator.set_lookup_tables(2, protonmuon_parameters.dedx_pdf_pl_2);


    llr_pid_calculator_k.set_dedx_binning(0, protonmuon_parameters_k.dedx_edges_pl_0);
    llr_pid_calculator_k.set_par_binning(0, protonmuon_parameters_k.parameters_edges_pl_0);
    llr_pid_calculator_k.set_lookup_tables(0, protonmuon_parameters_k.dedx_pdf_pl_0);

    llr_pid_calculator_k.set_dedx_binning(1, protonmuon_parameters_k.dedx_edges_pl_1);
    llr_pid_calculator_k.set_par_binning(1, protonmuon_parameters_k.parameters_edges_pl_1);
    llr_pid_calculator_k.set_lookup_tables(1, protonmuon_parameters_k.dedx_pdf_pl_1);

    llr_pid_calculator_k.set_dedx_binning(2, protonmuon_parameters_k.dedx_edges_pl_2);
    llr_pid_calculator_k.set_par_binning(2, protonmuon_parameters_k.parameters_edges_pl_2);
    llr_pid_calculator_k.set_lookup_tables(2, protonmuon_parameters_k.dedx_pdf_pl_2);
  }

  //========================================================================                                                                                                                                                                                                                                                    
  CCKaonAnalyzerRebuild::~CCKaonAnalyzerRebuild()
  {
    //destructor                                                                                                                                                                                                                                                                                                                
  }

  //========================================================================                                                                                                                                                                                                                                                    
  void CCKaonAnalyzerRebuild::endSubRun(const art::SubRun &subrun)
  {
    //if (!m_isData)                                                                                                                                                                                                                                                                                                            
    //{                                                                                                                                                                                                                                                                                                                         
    art::Handle<sumdata::POTSummary> potSummaryHandle;
    m_pot = subrun.getByLabel("generator", potSummaryHandle) ? static_cast<float>(potSummaryHandle->totpot) : 0.f;
    // -- std::cout << "[CCKaonAnalyzerRebuild::endSubRun] Storing POT info!" << std::endl;                                                                                                                                                                                                                                   
    //}                                                                                                                                                                                                                                                                                                                         

    m_run = subrun.run();
    m_subrun = subrun.subRun();
    fSubrunTree->Fill();
  }

  //========================================================================                                                                                                                                                                                                                                                    
  void CCKaonAnalyzerRebuild::filter(const art::Event &evt)
  {

    /*                                                                                                                                                                                                                                                                                                                          
  auto const TPC = (*geom).begin_TPC();                                                                                                                                                                                                                                                                                       
  auto ID = TPC.ID();                                                                                                                                                                                                                                                                                                         
  m_Cryostat = ID.Cryostat;                                                                                                                                                                                                                                                                                                   
  m_TPC = ID.TPC;                                                                                                                                                                                                                                                                                                             
  */

    _time2cm = theDetector->SamplingRate() / 1000.0 * theDetector->DriftVelocity( theDetector->Efield(), theDetector->Temperature() );//found in ProtoShowerPandora_tool.cc                                                                                                                                                     

    this->ClearVertex();


    //Collect the PFParticles from the event. This is the core!                                                                                                                                                                                                                                                                \
                                                                                                                                                                                                                                                                                                                              
    // Collect all the hits. We will need these. Lets grab both the handle as well as a vector of art::Ptr as I like both.                                                                                                                                                                                                
    art::ValidHandle<std::vector<recob::Hit>> const & hitHandle = evt.getValidHandle<std::vector<recob::Hit>>(m_hitfinderLabel);
    std::vector<art::Ptr<recob::Hit>> hitVector;
    art::fill_ptr_vector(hitVector,hitHandle);

    //Lets do "THE EXACT SAME STUFF" for Optical Flashes                                                                                                                                                                                                                                                                  
    art::ValidHandle<std::vector<recob::OpFlash>> const & flashHandle  = evt.getValidHandle<std::vector<recob::OpFlash>>(m_flashLabel);
    std::vector<art::Ptr<recob::OpFlash>> flashVector;
    art::fill_ptr_vector(flashVector,flashHandle);

    //tracks                                                                                                                                                                                                                                                                                                              
    art::ValidHandle<std::vector<recob::Track>> const & trackHandle  = evt.getValidHandle<std::vector<recob::Track>>(m_trackLabel);
    std::vector<art::Ptr<recob::Track>> trackVector;
    art::fill_ptr_vector(trackVector,trackHandle);

    //BadChannels                                                                                                                                                                                                                                                                                                         
    art::Handle<std::vector<int> > badChannelHandle;
    std::vector<int> badChannelVector;
    if(evt.getByLabel(m_badChannelProducer, m_badChannelLabel, badChannelHandle)){
      badChannelVector            = *(badChannelHandle);
    }


    //Collect the PFParticles from the event. This is the core!                                                                                                                                                                                                                                                           

    art::ValidHandle<std::vector<recob::PFParticle>> const & pfParticleHandle = evt.getValidHandle<std::vector<recob::PFParticle>>(m_pandoraLabel);
    std::vector<art::Ptr<recob::PFParticle>> pfParticleVector;
    art::fill_ptr_vector(pfParticleVector,pfParticleHandle);

    //get the cluster handle for the dQ/dx calc                                                                                                                                                                                                                                                                                 
    art::ValidHandle<std::vector<recob::Cluster>> const & clusterHandle = evt.getValidHandle<std::vector<recob::Cluster>>(m_pandoraLabel);
    std::vector< art::Ptr<recob::Cluster> > clusterVector;
    art::fill_ptr_vector(clusterVector,clusterHandle);

    PFParticleIdMap pfParticleMap;
    this->GetPFParticleIdMap(pfParticleHandle, pfParticleMap);

    //Slices                                                                                                                                                                                                                                                                                                                 \
                                                                                                                                                                                                                                                                                                                              
    art::ValidHandle<std::vector<recob::Slice>> const & sliceHandle  = evt.getValidHandle<std::vector<recob::Slice>>(m_pandoraLabel);
    std::vector<art::Ptr<recob::Slice>> sliceVector;
    art::fill_ptr_vector(sliceVector,sliceHandle);
    //And some associations                                                                                                                                                                                                                                                                                                  \
                                                                                                                                                                                                                                                                                                                              
    art::FindManyP<recob::PFParticle> pfparticles_per_slice(sliceHandle, evt, m_pandoraLabel);
    art::FindManyP<recob::Hit> hits_per_slice(sliceHandle, evt, m_pandoraLabel);

    std::map< art::Ptr<recob::Slice>, std::vector<art::Ptr<recob::PFParticle>> > sliceToPFParticlesMap;
    std::map<int, std::vector<art::Ptr<recob::PFParticle>> > sliceIDToPFParticlesMap;
    for(size_t i=0; i< sliceVector.size(); ++i){
      auto slice = sliceVector[i];
      sliceToPFParticlesMap[slice] =pfparticles_per_slice.at(slice.key());
      sliceIDToPFParticlesMap[slice->ID()] = pfparticles_per_slice.at(slice.key());
    }

    std::map< art::Ptr<recob::Slice>, std::vector<art::Ptr<recob::Hit>> > sliceToHitsMap;
    std::map<int, std::vector<art::Ptr<recob::Hit>> > sliceIDToHitsMap;
    for(size_t i=0; i< sliceVector.size(); ++i){
      auto slice = sliceVector[i];
      sliceToHitsMap[slice] =hits_per_slice.at(slice.key());
      sliceIDToHitsMap[slice->ID()] = hits_per_slice.at(slice.key());
    }

    //And some verticies.                                                                                                                                                                                                                                                                                                 
    art::ValidHandle<std::vector<recob::Vertex>> const & vertexHandle = evt.getValidHandle<std::vector<recob::Vertex>>(m_pandoraLabel);
    std::vector<art::Ptr<recob::Vertex>> vertexVector;
    art::fill_ptr_vector(vertexVector,vertexHandle);
    art::FindManyP<recob::Vertex> vertices_per_pfparticle(pfParticleHandle, evt, m_pandoraLabel);
    std::map< art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::Vertex>> > pfParticlesToVerticesMap;
    for(size_t i=0; i< pfParticleVector.size(); ++i){
      auto pfp = pfParticleVector[i];
      pfParticlesToVerticesMap[pfp] =vertices_per_pfparticle.at(pfp.key());
    }

    //------- 3D showers                                                                                                                                                                                                                                                                                                  

    art::FindOneP<recob::Shower> showerreco3D_per_pfparticle(pfParticleHandle, evt, m_shower3dLabel);
    std::map<art::Ptr<recob::PFParticle>, art::Ptr<recob::Shower>> pfParticlesToShowerReco3DMap;
    for(size_t i=0; i< pfParticleVector.size(); ++i){
      auto pfp = pfParticleVector[i];
      if(!showerreco3D_per_pfparticle.at(pfp.key()).isNull()){
	pfParticlesToShowerReco3DMap[pfp] = showerreco3D_per_pfparticle.at(pfp.key());
      }

    }


    // Once we have actual verticies, lets concentrate on JUST the neutrino PFParticles for now:                                                                                                                                                                                                                          
    //--------------------------------                                                                                                                                                                                                                                                                                    
    // Produce two PFParticle vectors containing final-state particles:                                                                                                                                                                                                                                                   
    // 1. Particles identified as cosmic-rays - recontructed under cosmic-hypothesis                                                                                                                                                                                                                                      
    // 2. Daughters of the neutrino PFParticle - reconstructed under the neutrino hypothesis                                                                                                                                                                                                                              
    std::vector< art::Ptr<recob::PFParticle> > crParticles;
    std::vector< art::Ptr<recob::PFParticle> > nuParticles;
    this->GetFinalStatePFParticleVectors(pfParticleMap, pfParticlesToVerticesMap, crParticles, nuParticles);


    if(m_is_verbose) std::cout<<"SinglePhoton::analyze() \t||\t Get Spacepoints"<<std::endl;
    //Look, here is a map that I just forced myself rather than build using helpers. Not that different is it. But for somereason I only use PFParticles.. huh,                                                                                                                                                           
    //Spacepoint associaitions                                                                                                                                                                                                                                                                                            
    art::FindManyP<recob::SpacePoint> spacePoints_per_pfparticle(pfParticleHandle, evt, m_pandoraLabel);
    std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::SpacePoint>> > pfParticleToSpacePointsMap;
    for(size_t i=0; i< nuParticles.size(); ++i){
      const art::Ptr<recob::PFParticle> pfp = nuParticles[i];
      pfParticleToSpacePointsMap[pfp] = spacePoints_per_pfparticle.at(pfp.key());
    }


    if(m_is_verbose) std::cout<<"SinglePhoton::analyze() \t||\t Get PandoraMetadata"<<std::endl;
    //add the associaton between PFP and metadata, this is important to look at the slices and scores                                                                                                                                                                                                                         

    art::FindManyP< larpandoraobj::PFParticleMetadata > pfPartToMetadataAssoc(pfParticleHandle, evt, m_pandoraLabel);
    std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> > pfParticleToMetadataMap;
    for(size_t i=0; i< pfParticleVector.size(); ++i){
      const art::Ptr<recob::PFParticle> pfp = pfParticleVector[i];
      pfParticleToMetadataMap[pfp] =  pfPartToMetadataAssoc.at(pfp.key());
    }


    if(m_is_verbose) std::cout<<"SinglePhoton::analyze() \t||\t Get Clusters"<<std::endl;

    //Get a map between the PFP's and the clusters. Although Mark isn't a fan of clusters, they're imporant for the shower dQ/dx                                                                                                                                                                                          
    //Also need a map between clusters and hits                                                                                                                                                                                                                                                                           
    art::FindManyP<recob::Cluster> clusters_per_pfparticle(pfParticleHandle, evt, m_pandoraLabel);
    art::FindManyP<recob::Hit> hits_per_cluster(clusterHandle, evt, m_pandoraLabel);
    std::map<art::Ptr<recob::PFParticle>,  std::vector<art::Ptr<recob::Cluster>> > pfParticleToClustersMap;
    std::map<art::Ptr<recob::Cluster>,  std::vector<art::Ptr<recob::Hit>> > clusterToHitsMap;
    //fill map PFP to Clusters                                                                                                                                                                                                                                                                                            
    for(size_t i=0; i< nuParticles.size(); ++i){
      auto pfp = nuParticles[i];
      pfParticleToClustersMap[pfp] = clusters_per_pfparticle.at(pfp.key());
      // pfParticleToSpacePointsMap[pfp] = spacePoints_per_pfparticle.at(pfp.key());                                                                                                                                                                                                                                    
      // pfParticleToMetadataMap[pfp] =  pfPartToMetadataAssoc.at(pfp.key());                                                                                                                                                                                                                                           
    }
    //fill map Cluster to Hits                                                                                                                                                                                                                                                                                            
    for(size_t i=0; i< clusterVector.size(); ++i){
      auto cluster = clusterVector[i];
      clusterToHitsMap[cluster] = hits_per_cluster.at(cluster.key());
    }

    if(m_is_verbose) std::cout<<"SinglePhoton::analyze() \t||\t Build hits to PFP Maps"<<std::endl;

    //taking out the Larpandora helper functions here because they don't match to non-neutrino slice hits for some reason                                                                                                                                                                                                 

    //OK Here we build two IMPORTANT maps for the analysis, (a) given a PFParticle get a vector of hits..                                                                                                                                                                                                                 
    //and (b) given a single hit, get the PFParticle it is in (MARK: is it only one? always? RE-MARK: Yes)                                                                                                                                                                                                                
    std::map<art::Ptr<recob::PFParticle>,  std::vector<art::Ptr<recob::Hit>> > pfParticleToHitsMap;
    //        std::map<art::Ptr<recob::Hit>, art::Ptr<recob::PFParticle>>                hitToPFParticleMap;                                                                                                                                                                                                              
    //Using a pandora helper here, but to be honest we should probably just build using normal associations so keep independant if pssoble                                                                                                                                                                                
    // lar_pandora::LArPandoraHelper::BuildPFParticleHitMaps(evt, m_pandoraLabel, pfParticleToHitsMap, hitToPFParticleMap, lar_pandora::LArPandoraHelper::kAddDaughters);                                                                                                                                                 

    //use pfp->cluster and cluster->hit to build pfp->hit map                                                                                                                                                                                                                                                             
    //for each PFP                                                                                                                                                                                                                                                                                                        
    for(size_t i=0; i<nuParticles.size(); ++i){
      auto pfp = nuParticles[i];

      // std::cout<<"starting to match to hits for pfp "<<pfp->Self()<<std::endl;                                                                                                                                                                                                                                       
      //get the associated clusters                                                                                                                                                                                                                                                                                     
      std::vector<art::Ptr<recob::Cluster>> clusters_vec  = pfParticleToClustersMap[pfp] ;

      //make empty vector to store hits                                                                                                                                                                                                                                                                                 
      std::vector<art::Ptr<recob::Hit>> hits_for_pfp = {};

      // std::cout<<"-- there are "<<clusters_vec.size()<<" associated clusters"<<std::endl;                                                                                                                                                                                                                            

      //for each cluster, get the associated hits                                                                                                                                                                                                                                                                       
      for (art::Ptr<recob::Cluster> cluster: clusters_vec){
	std::vector<art::Ptr<recob::Hit>> hits_vec =  clusterToHitsMap[cluster];

	//   std::cout<<"looking at cluster in pfp "<<pfp->Self()<<" with "<<hits_vec.size() <<" hits"<<std::endl;                                                                                                                                                                                                    
	//insert hits into vector                                                                                                                                                                                                                                                                                     
	hits_for_pfp.insert( hits_for_pfp.end(), hits_vec.begin(), hits_vec.end() );
      }

      //fill the map                                                                                                                                                                                                                                                                                                    
      pfParticleToHitsMap[pfp] = hits_for_pfp;
      //std::cout<<"saving a total of "<<hits_for_pfp.size()<<" hits for pfp "<<pfp->Self()<<std::endl;                                                                                                                                                                                                                 

    }//for each pfp                                                                                                                                                                                                                                                                                                       



    //these are all filled in analyze slice                                                                                                                                                                                                                                                                                  \
                                                                                                                                                                                                                                                                                                                              
    std::vector< art::Ptr<recob::Track> > tracks;
    std::vector< art::Ptr<recob::Shower> > showers;
    std::map< art::Ptr<recob::Track> , art::Ptr<recob::PFParticle >> trackToNuPFParticleMap;
    std::map< art::Ptr<recob::Shower> , art::Ptr<recob::PFParticle>> showerToNuPFParticleMap;


    if(m_is_verbose) std::cout<<"SinglePhoton::analyze() \t||\t Get Tracks and Showers"<<std::endl;
    this->CollectTracksAndShowers(nuParticles, pfParticleMap,  pfParticleHandle, evt, tracks, showers, trackToNuPFParticleMap, showerToNuPFParticleMap);

    //**********************************************************************************************/                                                                                                                                                                                                                     
    //**********************************************************************************************/                                                                                                                                                                                                                     
    //---------------------------------- MC TRUTH Data Only---------------------------                                                                                                                                                                                                                                    
    //**********************************************************************************************/                                                                                                                                                                                                                     
    //**********************************************************************************************/                                                                                                                                                                                                                     

    //Get the MCtruth handles and vectors                                                                                                                                                                                                                                                                                 
    std::vector<art::Ptr<simb::MCTruth>> mcTruthVector;
    std::vector<art::Ptr<simb::MCParticle>> mcParticleVector;

    //Then build a map from MCparticles to Hits and vice versa                                                                                                                                                                                                                                                            
    std::map< art::Ptr<simb::MCParticle>,  std::vector<art::Ptr<recob::Hit> >  >  mcParticleToHitsMap;
    std::map< art::Ptr<recob::Hit>, art::Ptr<simb::MCParticle> >                  hitToMCParticleMap;


    //Apparrently a MCParticle doesn't know its origin (thanks Andy!)                                                                                                                                                                                                                                                     
    //I would also like a map from MCparticle to MCtruth and then I will be done.  and Vice Versa                                                                                                                                                                                                                         
    //Note which map is which!       //First  is one-to-many.         //Second is one-to-one                                                                                                                                                                                                                              
    std::map< art::Ptr<simb::MCTruth>,    std::vector<art::Ptr<simb::MCParticle>>>  MCTruthToMCParticlesMap;
    std::map< art::Ptr<simb::MCParticle>, art::Ptr<simb::MCTruth>>                  MCParticleToMCTruthMap;
    std::map<int, art::Ptr<simb::MCParticle> >                                     MCParticleToTrackIdMap;

    std::vector<art::Ptr<sim::MCTrack>> mcTrackVector;
    std::vector<art::Ptr<sim::MCShower>> mcShowerVector;

    std::vector<art::Ptr<simb::MCParticle>> matchedMCParticleVector;
    std::map<art::Ptr<recob::Track>, art::Ptr<simb::MCParticle> > trackToMCParticleMap;
    std::map<art::Ptr<recob::Shower>, art::Ptr<simb::MCParticle> > showerToMCParticleMap;

    //Given a simb::MCParticle we would like a map to either a sim::MCTrack or sim::MCShower                                                                                                                                                                                                                              
    std::map< art::Ptr<simb::MCParticle>, art::Ptr<sim::MCTrack> > MCParticleToMCTrackMap;
    std::map< art::Ptr<simb::MCParticle>, art::Ptr<sim::MCShower> > MCParticleToMCShowerMap;


    //**********************************************************************************************/                                                                                                                                                                                                                     
    //**********************************************************************************************/                                                                                                                                                                                                                     
    //Some event based properties                                                                                                                                                                                                                                                                                         



    std::vector<std::pair<art::Ptr<recob::PFParticle>,int>> allPFPSliceIdVec; //stores a pair of all PFP's in the event and the slice ind                                                                                                                                                                                     
    std::vector<std::pair<art::Ptr<recob::PFParticle>,int>> primaryPFPSliceIdVec; //stores a pair of only the primary PFP's in the event and the slice ind                                                                                                                                                                   \
                                                                                                                                                                                                                                                                                                                              
    std::map<int, double> sliceIdToNuScoreMap; //map between a slice Id and neutrino score                                                                                                                                                                                                                                    
    std::map<art::Ptr<recob::PFParticle>, bool> PFPToClearCosmicMap; //returns true for clear cosmic, false otherwise                                                                                                                                                                                                         
    std::map<art::Ptr<recob::PFParticle>, int> PFPToSliceIdMap; //returns the slice id for all PFP's                                                                                                                                                                                                                          
    std::map<art::Ptr<recob::PFParticle>,bool> PFPToNuSliceMap;
    std::map<art::Ptr<recob::PFParticle>,double> PFPToTrackScoreMap;
    std::map<int, int> sliceIdToNumPFPsMap;
    std::cout<<"SinglePhoton::analyze::AnalyzeSlice()\t||\t Starting"<<std::endl;


    this->AnalyzeSlices(pfParticleToMetadataMap, pfParticleMap,  primaryPFPSliceIdVec, sliceIdToNuScoreMap, PFPToClearCosmicMap, PFPToSliceIdMap, PFPToNuSliceMap, PFPToTrackScoreMap);
    //std::cout<<"There are "<< allPFPSliceIdVec.size()<<" pfp-slice id matches stored in the vector"<<std::endl;                                                                                                                                                                                                             
    std::cout<<"SinglePhoton::analyze\t||\tthe number of PPF's with stored clear cosmic info is "<<PFPToClearCosmicMap.size()<<std::endl;
    std::cout<<"SinglePhoton::analyze\t||\tthe number of PFP's stored in the PFPToSliceIdMap is "<<PFPToSliceIdMap.size()<<std::endl;

    if (PFPToSliceIdMap.size() < 1){
      std::cout<<"ERROR, not storing PFP's in PFPToSliceIdMap"<<std::endl;
    }


    for (auto pair:sliceIDToPFParticlesMap){
      std::vector<art::Ptr<recob::PFParticle>> pfp_vec = pair.second;
      int slice_id = pair.first;
      //if (slice_vec[0]->Slice() != PFPToSliceIdMap[pfp] )                                                                                                                                                                                                                                                                  \
                                                                                                                                                                                                                                                                                                                              
      for(auto pfp: pfp_vec){
        if (slice_id != PFPToSliceIdMap[pfp] && PFPToSliceIdMap[pfp]>=0){
	  std::cout<<"sliceIDToPFParticlesMap[slice->ID()] for pfp "<<pfp->Self()<<" is slice "<< slice_id<< "but PFPToSliceIdMap[pfp] = "<<PFPToSliceIdMap[pfp]<<std::endl;
        }
      }

    }


    //if CRT info, get CRT hits                                                                                                                                                                                                                                                                                               
    art::Handle<std::vector<crt::CRTHit>> crthit_h; //only filled when there are hits, otherwise empty                                                                                                                                                                                                                        
    art::Handle<raw::DAQHeaderTimeUBooNE> rawHandle_DAQHeader;
    double evt_timeGPS_nsec = -999;
    if(m_runCRT){
      evt.getByLabel(m_DAQHeaderProducer, rawHandle_DAQHeader);
      evt.getByLabel(m_CRTHitProducer, crthit_h);
      raw::DAQHeaderTimeUBooNE const& my_DAQHeader(*rawHandle_DAQHeader);
      art::Timestamp evtTimeGPS = my_DAQHeader.gps_time();
      evt_timeGPS_nsec = evtTimeGPS.timeLow();

      std::cout<<"SinglePhoton::analyze \t||\t Got CRT hits"<<std::endl;
    }

    this->AnalyzeFlashes(flashVector, crthit_h, evt_timeGPS_nsec);

    this->AnalyzeTracks(tracks, trackToNuPFParticleMap, pfParticleToSpacePointsMap,  MCParticleToTrackIdMap, sliceIdToNuScoreMap, PFPToClearCosmicMap,  PFPToSliceIdMap,  PFPToTrackScoreMap, PFPToNuSliceMap,pfParticleMap);

    this->AnalyzeShowers(showers,showerToNuPFParticleMap, pfParticleToHitsMap, pfParticleToClustersMap, clusterToHitsMap,sliceIdToNuScoreMap, PFPToClearCosmicMap,  PFPToSliceIdMap, PFPToNuSliceMap, PFPToTrackScoreMap,pfParticleMap,pfParticlesToShowerReco3DMap);


    if(!m_is_data){

      art::ValidHandle<std::vector<simb::MCTruth>> const & mcTruthHandle= evt.getValidHandle<std::vector<simb::MCTruth>>(m_generatorLabel);
      art::fill_ptr_vector(mcTruthVector,mcTruthHandle);
      art::ValidHandle<std::vector<simb::MCParticle>> const & mcParticleHandle= evt.getValidHandle<std::vector<simb::MCParticle>>(m_geantModuleLabel);
      art::fill_ptr_vector(mcParticleVector,mcParticleHandle);

      this->CollectMCParticles(evt, m_geantModuleLabel, MCTruthToMCParticlesMap, MCParticleToMCTruthMap, MCParticleToTrackIdMap);

      art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData> mcparticles_per_hit(hitHandle, evt, m_hitMCParticleAssnsLabel);

      //mcc9 march miniretreat fix                                                                                                                                                                                                                                                                                      
      std::vector<art::Ptr<simb::MCParticle>> particle_vec; //vector of all MCParticles associated with a given hit in the reco PFP                                                                                                                                                                                     
      std::vector<anab::BackTrackerHitMatchingData const *> match_vec; //vector of some backtracker thing                                                                                                                                                                                                               

      m_test_matched_hits = 0;

      for(size_t j=0; j<hitVector.size();j++){
	const art::Ptr<recob::Hit> hit = hitVector[j];

	particle_vec.clear(); match_vec.clear(); //only store per hit                                                                                                                                                                                                                                                 

	mcparticles_per_hit.get(hit.key(), particle_vec, match_vec);

	if(particle_vec.size() > 0){
	  m_test_matched_hits++;
	}

      }

      this->BuildMCParticleHitMaps(evt, m_geantModuleLabel, hitVector,  mcParticleToHitsMap, hitToMCParticleMap, lar_pandora::LArPandoraHelper::kAddDaughters,  MCParticleToTrackIdMap);

      //std::cout<<"SinglePhoton\t||\t Starting backtracker on recob::track"<<std::endl;                                                                                                                                                                                                                                
      //std::vector<double> trk_overlay_vec = recoMCmatching<art::Ptr<recob::Track>>( tracks, trackToMCParticleMap, trackToNuPFParticleMap, pfParticleToHitsMap, mcparticles_per_hit, matchedMCParticleVector);                                                                                                         

      std::cout<<"SinglePhoton\t||\t Starting backtracker on recob::shower"<<std::endl;
      this->showerRecoMCmatching(showers, showerToMCParticleMap, showerToNuPFParticleMap, pfParticleToHitsMap, mcparticles_per_hit, matchedMCParticleVector, pfParticleMap,  MCParticleToTrackIdMap, sliceIdToNuScoreMap, PFPToClearCosmicMap,  PFPToSliceIdMap, PFPToNuSliceMap);



      std::cout<<"filling info in ncdelta slice tree"<<std::endl;
      this->AnalyzeRecoMCSlices( m_truthmatching_signaldef, MCParticleToTrackIdMap, showerToNuPFParticleMap , allPFPSliceIdVec, showerToMCParticleMap, trackToNuPFParticleMap, trackToMCParticleMap,  PFPToSliceIdMap);
    }

    for (auto pair:PFPToNuSliceMap){
      auto pfp = pair.first;
      auto is_nuslice = pair.second;
      if (is_nuslice){
	std::cout<<"pfp in nuslice "<<pfp->Self()<<std::endl;
      }

    }


    //---------------------- END OF LOOP, fill vertex ---------------------                                                                                                                                                                                                                                                   
    ncdelta_slice_tree->Fill();
    vertex_tree->Fill();
    eventweight_tree->Fill();


  }


  //========================================================================                                                                                                                                                                                                                                                    
  void CCKaonAnalyzerRebuild::beginJob()
  {

    //initialiseCanvas();                                                                                                                                                                                                                                                                                                       

    art::ServiceHandle<art::TFileService> tfs;

    fEventTree = tfs->make<TTree>("Event", "Event Tree from Reco");

    //gInterpreter->GenerateDictionary("vector<vector<vector<Float_t>> >", "vector");                                                                                                                                                                                                                                           
    fEventTree->Branch("event", &event, "event/I");
    fEventTree->Branch("run", &run, "run/I");
    fEventTree->Branch("subrun", &subrun, "surbrun/I");

    fEventTree->Branch("IsKaon", &IsKaon, "IsKaon/I");
    fEventTree->Branch("IsSingleKaon", &IsSingleKaon, "IsSingleKaon/I");
    fEventTree->Branch("IsAssociatedKaon", &IsAssociatedKaon, "IsAssociatedKaon/I");
    fEventTree->Branch("IsMuBR", &IsMuBR, "IsMuBR/I");
    fEventTree->Branch("IsPiBR", &IsPiBR, "IsPiBR/I");
    fEventTree->Branch("prip", &prip, "prip/I");
    fEventTree->Branch("prip_k_dau", &prip_k_dau, "prip_k_dau/I");


    fEventTree->Branch("process",&Process);
    fEventTree->Branch("endprocess",&EndProcess);

    fEventTree->Branch("event_weight" ,&event_weight);

    fEventTree->Branch("evtwgt_funcname" ,&evtwgt_funcname);
    fEventTree->Branch("evtwgt_weight" ,&evtwgt_weight);
    fEventTree->Branch("evtwgt_nweight" ,&evtwgt_nweight);

    fEventTree->Branch("true_nu_energy", &true_nu_energy, "true_nu_energy/F");
    fEventTree->Branch("true_nu_pdg", &true_nu_pdg, "true_nu_pdg/I");
    fEventTree->Branch("true_nu_mode", &true_nu_mode, "true_nu_mode/I");
    fEventTree->Branch("true_nu_ccnc", &true_nu_ccnc, "true_nu_ccnc/I");
    fEventTree->Branch("true_nu_vtx_x", &true_nu_vtx_x, "true_nu_vtx_x/F");
    fEventTree->Branch("true_nu_vtx_y", &true_nu_vtx_y, "true_nu_vtx_y/F");
    fEventTree->Branch("true_nu_vtx_z", &true_nu_vtx_z, "true_nu_vtx_z/F");
    fEventTree->Branch("true_nu_vtx_inTPC", &true_nu_vtx_inTPC, "true_nu_vtx_inTPC/O");
    fEventTree->Branch("true_nu_vtx_in5cmTPC", &true_nu_vtx_in5cmTPC, "true_nu_vtx_in5cmTPC/O");
    fEventTree->Branch("true_nu_vtx_inCCInclusiveTPC", &true_nu_vtx_inCCInclusiveTPC, "true_nu_vtx_inCCInclusiveTPC/O");

    fEventTree->Branch("true_lepton_pdg", &true_lepton_pdg, "true_lepton_pdg/I");
    fEventTree->Branch("true_lepton_p", &true_lepton_p, "true_lepton_p/F");
    fEventTree->Branch("true_lepton_ke", &true_lepton_ke, "true_lepton_ke/F");
    fEventTree->Branch("true_lepton_theta", &true_lepton_theta, "true_lepton_theta/F");
    fEventTree->Branch("true_lepton_costheta", &true_lepton_costheta, "true_lepton_costheta/F");
    fEventTree->Branch("true_lepton_phi", &true_lepton_phi, "true_lepton_phi/F");

    fEventTree->Branch("true_nkaons", &true_nkaons, "true_nkaons/I");

    fEventTree->Branch("true_kaon_length", &true_kaon_length, "true_kaon_length/F");
    fEventTree->Branch("true_kaon_p", &true_kaon_p, "true_kaon_p/F");
    fEventTree->Branch("true_kaon_ke", &true_kaon_ke, "true_kaon_ke/F");
    fEventTree->Branch("true_kaon_theta", &true_kaon_theta, "true_kaon_theta/F");
    fEventTree->Branch("true_kaon_costheta", &true_kaon_costheta, "true_kaon_costheta/F");
    fEventTree->Branch("true_kaon_phi", &true_kaon_phi, "true_kaon_phi/F");
    fEventTree->Branch("true_kaon_ccmuon_angle", &true_kaon_ccmuon_angle, "true_kaon_ccmuon_angle/F");
    fEventTree->Branch("true_kaon_ccmuon_cosangle", &true_kaon_ccmuon_cosangle, "true_kaon_ccmuon_cosangle/F");
    fEventTree->Branch("true_kaon_end_process", &true_kaon_end_process, "true_kaon_end_process/I");
    fEventTree->Branch("true_kaon_end_ke", &true_kaon_end_ke, "true_kaon_end_ke/F");
    fEventTree->Branch("true_kaon_start_x", &true_kaon_start_x, "true_kaon_start_x/F");
    fEventTree->Branch("true_kaon_start_y", &true_kaon_start_y, "true_kaon_start_y/F");
    fEventTree->Branch("true_kaon_start_z", &true_kaon_start_z, "true_kaon_start_z/F");
    fEventTree->Branch("true_kaon_end_x", &true_kaon_end_x, "true_kaon_end_x/F");
    fEventTree->Branch("true_kaon_end_y", &true_kaon_end_y, "true_kaon_end_y/F");
    fEventTree->Branch("true_kaon_end_z", &true_kaon_end_z, "true_kaon_end_z/F");
    fEventTree->Branch("true_kaon_end_inTPC", &true_kaon_end_inTPC, "true_kaon_end_inTPC/O");
    fEventTree->Branch("true_kaon_end_in5cmTPC", &true_kaon_end_in5cmTPC, "true_kaon_end_in5cmTPC/O");
    fEventTree->Branch("true_kaon_end_inCCInclusiveTPC", &true_kaon_end_inCCInclusiveTPC, "true_kaon_end_inCCInclusiveTPC/O");

    fEventTree->Branch("true_dau_muon_length", &true_dau_muon_length, "true_dau_muon_length/F");
    fEventTree->Branch("true_dau_muon_p", &true_dau_muon_p, "true_dau_muon_p/F");
    fEventTree->Branch("true_dau_muon_ke", &true_dau_muon_ke, "true_dau_muon_ke/F");
    fEventTree->Branch("true_dau_muon_theta", &true_dau_muon_theta, "true_dau_muon_theta/F");
    fEventTree->Branch("true_dau_muon_costheta", &true_dau_muon_costheta, "true_dau_muon_costheta/F");
    fEventTree->Branch("true_dau_muon_phi", &true_dau_muon_phi, "true_dau_muon_phi/F");
    fEventTree->Branch("true_dau_muon_ccmuon_angle", &true_dau_muon_ccmuon_angle, "true_dau_muon_ccmuon_angle/F");
    fEventTree->Branch("true_dau_muon_ccmuon_cosangle", &true_dau_muon_ccmuon_cosangle, "true_dau_muon_ccmuon_cosangle/F");
    fEventTree->Branch("true_dau_muon_end_process", &true_dau_muon_end_process, "true_dau_muon_end_process/I");
    fEventTree->Branch("true_dau_muon_end_ke", &true_dau_muon_end_ke, "true_dau_muon_end_ke/F");
    fEventTree->Branch("true_dau_muon_start_x", &true_dau_muon_start_x, "true_dau_muon_start_x/F");
    fEventTree->Branch("true_dau_muon_start_y", &true_dau_muon_start_y, "true_dau_muon_start_y/F");
    fEventTree->Branch("true_dau_muon_start_z", &true_dau_muon_start_z, "true_dau_muon_start_z/F");
    fEventTree->Branch("true_dau_muon_end_x", &true_dau_muon_end_x, "true_dau_muon_end_x/F");
    fEventTree->Branch("true_dau_muon_end_y", &true_dau_muon_end_y, "true_dau_muon_end_y/F");
    fEventTree->Branch("true_dau_muon_end_z", &true_dau_muon_end_z, "true_dau_muon_end_z/F");

    fEventTree->Branch("true_dau_pip_length", &true_dau_pip_length, "true_dau_pip_length/F");
    fEventTree->Branch("true_dau_pip_p", &true_dau_pip_p, "true_dau_pip_p/F");
    fEventTree->Branch("true_dau_pip_ke", &true_dau_pip_ke, "true_dau_pip_ke/F");
    fEventTree->Branch("true_dau_pip_theta", &true_dau_pip_theta, "true_dau_pip_theta/F");
    fEventTree->Branch("true_dau_pip_costheta", &true_dau_pip_costheta, "true_dau_pip_costheta/F");
    fEventTree->Branch("true_dau_pip_phi", &true_dau_pip_phi, "true_dau_pip_phi/F");
    fEventTree->Branch("true_dau_pip_ccmuon_angle", &true_dau_pip_ccmuon_angle, "true_dau_pip_ccmuon_angle/F");
    fEventTree->Branch("true_dau_pip_ccmuon_cosangle", &true_dau_pip_ccmuon_cosangle, "true_dau_pip_ccmuon_cosangle/F");
    fEventTree->Branch("true_dau_pip_end_process", &true_dau_pip_end_process, "true_dau_pip_end_process/I");
    fEventTree->Branch("true_dau_pip_end_ke", &true_dau_pip_end_ke, "true_dau_pip_end_ke/F");
    fEventTree->Branch("true_dau_pip_start_x", &true_dau_pip_start_x, "true_dau_pip_start_x/F");
    fEventTree->Branch("true_dau_pip_start_y", &true_dau_pip_start_y, "true_dau_pip_start_y/F");
    fEventTree->Branch("true_dau_pip_start_z", &true_dau_pip_start_z, "true_dau_pip_start_z/F");
    fEventTree->Branch("true_dau_pip_end_x", &true_dau_pip_end_x, "true_dau_pip_end_x/F");
    fEventTree->Branch("true_dau_pip_end_y", &true_dau_pip_end_y, "true_dau_pip_end_y/F");
    fEventTree->Branch("true_dau_pip_end_z", &true_dau_pip_end_z, "true_dau_pip_end_z/F");


    fEventTree->Branch("true_dau_pin_length", &true_dau_pin_length, "true_dau_pin_length/F");
    fEventTree->Branch("true_dau_pin_p", &true_dau_pin_p, "true_dau_pin_p/F");
    fEventTree->Branch("true_dau_pin_ke", &true_dau_pin_ke, "true_dau_pin_ke/F");
    fEventTree->Branch("true_dau_pin_theta", &true_dau_pin_theta, "true_dau_pin_theta/F");
    fEventTree->Branch("true_dau_pin_costheta", &true_dau_pin_costheta, "true_dau_pin_costheta/F");
    fEventTree->Branch("true_dau_pin_phi", &true_dau_pin_phi, "true_dau_pin_phi/F");
    fEventTree->Branch("true_dau_pin_ccmuon_angle", &true_dau_pin_ccmuon_angle, "true_dau_pin_ccmuon_angle/F");
    fEventTree->Branch("true_dau_pin_ccmuon_cosangle", &true_dau_pin_ccmuon_cosangle, "true_dau_pin_ccmuon_cosangle/F");
    fEventTree->Branch("true_dau_pin_end_process", &true_dau_pin_end_process, "true_dau_pin_end_process/I");
    fEventTree->Branch("true_dau_pin_end_ke", &true_dau_pin_end_ke, "true_dau_pin_end_ke/F");
    fEventTree->Branch("true_dau_pin_start_x", &true_dau_pin_start_x, "true_dau_pin_start_x/F");
    fEventTree->Branch("true_dau_pin_start_y", &true_dau_pin_start_y, "true_dau_pin_start_y/F");
    fEventTree->Branch("true_dau_pin_start_z", &true_dau_pin_start_z, "true_dau_pin_start_z/F");
    fEventTree->Branch("true_dau_pin_end_x", &true_dau_pin_end_x, "true_dau_pin_end_x/F");
    fEventTree->Branch("true_dau_pin_end_y", &true_dau_pin_end_y, "true_dau_pin_end_y/F");
    fEventTree->Branch("true_dau_pin_end_z", &true_dau_pin_end_z, "true_dau_pin_end_z/F");

    fEventTree->Branch("cheat_num_hits", &cheat_num_hits, "cheat_num_hits[20][10]/I");
    fEventTree->Branch("cheat_ini_hits", &cheat_ini_hits, "cheat_ini_hits[20][10]/I");
    fEventTree->Branch("cheat_closest_distance", &cheat_closest_distance, "cheat_closest_distance[20][10]/F");

    fEventTree->Branch("cheat_peak_pdg", &cheat_peak_pdg, "cheat_peak_pdg[20][10]/F");
    fEventTree->Branch("cheat_peak_theta", &cheat_peak_theta, "cheat_peak_theta[20][10]/F");
    fEventTree->Branch("cheat_pip_trkln", &cheat_pip_trkln, "cheat_pip_trkln/F");
    fEventTree->Branch("cheat_mup_trkln", &cheat_mup_trkln, "cheat_mup_trkln/F");

    fEventTree->Branch("cheat_peak_phi", &cheat_peak_phi, "cheat_peak_phi[20][10]/F");
    fEventTree->Branch("best_peak_theta", &best_peak_theta, "best_peak_theta[20][10]/F");
    fEventTree->Branch("best_peak_phi", &best_peak_phi, "best_peak_phi[20][10]/F");
    fEventTree->Branch("best_peak_trkln", &best_peak_trkln, "best_peak_trkln[20][10]/F");

    fEventTree->Branch("true_length", &true_length, "true_length[10]/F");
    fEventTree->Branch("true_p", &true_p, "true_p[10]/F");
    fEventTree->Branch("true_ke", &true_ke, "true_ke[10]/F");
    fEventTree->Branch("true_theta", &true_theta, "true_theta[10]/F");
    fEventTree->Branch("true_costheta", &true_costheta, "true_costheta[10]/F");
    fEventTree->Branch("true_phi", &true_phi, "true_phi[10]/F");
    fEventTree->Branch("true_ccmuon_angle", &true_ccmuon_angle, "true_ccmuon_angle[10]/F");
    fEventTree->Branch("true_ccmuon_cosangle", &true_ccmuon_cosangle, "true_ccmuon_cosangle[10]/F");
    fEventTree->Branch("true_end_process", &true_end_process, "true_end_process[10]/I");
    fEventTree->Branch("true_end_ke", &true_end_ke, "true_end_ke[10]/F");
    fEventTree->Branch("true_end_x", &true_end_x, "true_end_x[10]/F");
    fEventTree->Branch("true_end_y", &true_end_y, "true_end_y[10]/F");
    fEventTree->Branch("true_end_z", &true_end_z, "true_end_z[10]/F");
    fEventTree->Branch("true_end_inTPC", &true_end_inTPC, "true_end_inTPC[10]/O");
    fEventTree->Branch("true_end_in5cmTPC", &true_end_in5cmTPC, "true_end_in5cmTPC[10]/O");
    fEventTree->Branch("true_end_inCCInclusiveTPC", &true_end_inCCInclusiveTPC, "true_end_inCCInclusiveTPC[10]/O");

    fEventTree->Branch("true_kaon_ndaughters", &true_kaon_ndaughters, "true_kaon_ndaughters/I");
    fEventTree->Branch("true_kaon_ndaughters_decay", &true_kaon_ndaughters_decay, "true_kaon_ndaughters_decay/I");
    fEventTree->Branch("true_kaon_ndaughters_inelastic", &true_kaon_ndaughters_inelastic, "true_kaon_ndaughters_inelastic/I");
  fEventTree->Branch("true_kaon_ndecmup", &true_kaon_ndecmup, "true_kaon_ndecmup/I");
  fEventTree->Branch("true_kaon_ndecpip", &true_kaon_ndecpip, "true_kaon_ndecpip/I");
  fEventTree->Branch("true_kaon_ninekap", &true_kaon_ninekap, "true_kaon_ninekap/I");
  fEventTree->Branch("true_kaon_ninepip", &true_kaon_ninepip, "true_kaon_ninepip/I");
  fEventTree->Branch("true_kaon_ninepro", &true_kaon_ninepro, "true_kaon_ninepro/I");
  fEventTree->Branch("true_kaon_daughter_length", &true_kaon_daughter_length, "true_kaon_daughter_length/F");
  fEventTree->Branch("true_kaon_daughter_p", &true_kaon_daughter_p, "true_kaon_daughter_p/F");
  fEventTree->Branch("true_kaon_daughter_ke", &true_kaon_daughter_ke, "true_kaon_daughter_ke/F");
  fEventTree->Branch("true_kaon_daughter_theta", &true_kaon_daughter_theta, "true_kaon_daughter_theta/F");
  fEventTree->Branch("true_kaon_daughter_costheta", &true_kaon_daughter_costheta, "true_kaon_daughter_costheta/F");
  fEventTree->Branch("true_kaon_daughter_angle", &true_kaon_daughter_angle, "true_kaon_daughter_angle/F");
  fEventTree->Branch("true_kaon_daughter_cosangle", &true_kaon_daughter_cosangle, "true_kaon_daughter_cosangle/F");
  fEventTree->Branch("true_kaon_daughter_pdg", &true_kaon_daughter_pdg, "true_kaon_daughter_pdg/I");
  fEventTree->Branch("true_kaon_daughter_start_x", &true_kaon_daughter_start_x, "true_kaon_daughter_start_x/F");
  fEventTree->Branch("true_kaon_daughter_start_y", &true_kaon_daughter_start_y, "true_kaon_daughter_start_y/F");
  fEventTree->Branch("true_kaon_daughter_start_z", &true_kaon_daughter_start_z, "true_kaon_daughter_start_z/F");
  fEventTree->Branch("true_kaon_daughter_end_x", &true_kaon_daughter_end_x, "true_kaon_daughter_end_x/F");
  fEventTree->Branch("true_kaon_daughter_end_y", &true_kaon_daughter_end_y, "true_kaon_daughter_end_y/F");
  fEventTree->Branch("true_kaon_daughter_end_z", &true_kaon_daughter_end_z, "true_kaon_daughter_end_z/F");
  fEventTree->Branch("true_kaon_daughter_end_inTPC", &true_kaon_daughter_end_inTPC, "true_kaon_daughter_endInTPC/O");
  fEventTree->Branch("true_kaon_daughter_end_in5cmTPC", &true_kaon_daughter_end_in5cmTPC, "true_kaon_daughter_end_in5cmTPC/O");
  fEventTree->Branch("true_kaon_daughter_end_inCCInclusiveTPC", &true_kaon_daughter_end_inCCInclusiveTPC, "true_kaon_daughter_end_inCCInclusiveTPC/O");

  fEventTree->Branch("true_nhyperons", &true_nhyperons, "true_nhyperons/I");

  fEventTree->Branch("true_nkaons_anti", &true_nkaons_anti, "true_nkaons_anti/I");
  fEventTree->Branch("true_kaon_length_anti", &true_kaon_length_anti, "true_kaon_length_anti/F");
  fEventTree->Branch("true_kaon_p_anti", &true_kaon_p_anti, "true_kaon_p_anti/F");
  fEventTree->Branch("true_kaon_ke_anti", &true_kaon_ke_anti, "true_kaon_ke_anti/F");
  fEventTree->Branch("true_kaon_theta_anti", &true_kaon_theta_anti, "true_kaon_theta_anti/F");
  fEventTree->Branch("true_kaon_costheta_anti", &true_kaon_costheta_anti, "true_kaon_costheta_anti/F");
  fEventTree->Branch("true_kaon_phi_anti", &true_kaon_phi_anti, "true_kaon_phi_anti/F");
  fEventTree->Branch("true_kaon_ccmuon_angle_anti", &true_kaon_ccmuon_angle_anti, "true_kaon_ccmuon_angle_anti/F");
  fEventTree->Branch("true_kaon_ccmuon_cosangle_anti", &true_kaon_ccmuon_cosangle_anti, "true_kaon_ccmuon_cosangle_anti/F");
  fEventTree->Branch("true_kaon_end_process_anti", &true_kaon_end_process_anti, "true_kaon_end_process_anti/I");
  fEventTree->Branch("true_kaon_end_ke_anti", &true_kaon_end_ke_anti, "true_kaon_end_ke_anti/F");
  fEventTree->Branch("true_kaon_end_x_anti", &true_kaon_end_x_anti, "true_kaon_end_x_anti/F");
  fEventTree->Branch("true_kaon_end_y_anti", &true_kaon_end_y_anti, "true_kaon_end_y_anti/F");
  fEventTree->Branch("true_kaon_end_z_anti", &true_kaon_end_z_anti, "true_kaon_end_z_anti/F");
  fEventTree->Branch("true_kaon_end_inTPC_anti", &true_kaon_end_inTPC_anti, "true_kaon_end_inTPC_anti/O");
  fEventTree->Branch("true_kaon_end_in5cmTPC_anti", &true_kaon_end_in5cmTPC_anti, "true_kaon_end_in5cmTPC_anti/O");
  fEventTree->Branch("true_kaon_end_inCCInclusiveTPC_anti", &true_kaon_end_inCCInclusiveTPC_anti, "true_kaon_end_inCCInclusiveTPC_anti/O");

  fEventTree->Branch("true_kaon_ndaughters_anti", &true_kaon_ndaughters_anti, "true_kaon_ndaughters_anti/I");
  fEventTree->Branch("true_kaon_ndaughters_decay_anti", &true_kaon_ndaughters_decay_anti, "true_kaon_ndaughters_decay_anti/I");
  fEventTree->Branch("true_kaon_ndaughters_inelastic_anti", &true_kaon_ndaughters_inelastic_anti, "true_kaon_ndaughters_inelastic_anti/I");
  fEventTree->Branch("true_kaon_ndecmup_anti", &true_kaon_ndecmup_anti, "true_kaon_ndecmup_anti/I");
  fEventTree->Branch("true_kaon_ndecpip_anti", &true_kaon_ndecpip_anti, "true_kaon_ndecpip_anti/I");
  fEventTree->Branch("true_kaon_ninekap_anti", &true_kaon_ninekap_anti, "true_kaon_ninekap_anti/I");
  fEventTree->Branch("true_kaon_ninepip_anti", &true_kaon_ninepip_anti, "true_kaon_ninepip_anti/I");
  fEventTree->Branch("true_kaon_ninepro_anti", &true_kaon_ninepro_anti, "true_kaon_ninepro_anti/I");
  fEventTree->Branch("true_kaon_daughter_length_anti", &true_kaon_daughter_length_anti, "true_kaon_daughter_length_anti/F");
  fEventTree->Branch("true_kaon_daughter_p_anti", &true_kaon_daughter_p_anti, "true_kaon_daughter_p_anti/F");
  fEventTree->Branch("true_kaon_daughter_ke_anti", &true_kaon_daughter_ke_anti, "true_kaon_daughter_ke_anti/F");
  fEventTree->Branch("true_kaon_daughter_theta_anti", &true_kaon_daughter_theta_anti, "true_kaon_daughter_theta_anti/F");
  fEventTree->Branch("true_kaon_daughter_costheta_anti", &true_kaon_daughter_costheta_anti, "true_kaon_daughter_costheta_anti/F");
  fEventTree->Branch("true_kaon_daughter_angle_anti", &true_kaon_daughter_angle_anti, "true_kaon_daughter_angle_anti/F");
  fEventTree->Branch("true_kaon_daughter_cosangle_anti", &true_kaon_daughter_cosangle_anti, "true_kaon_daughter_cosangle_anti/F");
  fEventTree->Branch("true_kaon_daughter_pdg_anti", &true_kaon_daughter_pdg_anti, "true_kaon_daughter_pdg_anti/I");
  fEventTree->Branch("true_kaon_daughter_end_x_anti", &true_kaon_daughter_end_x_anti, "true_kaon_daughter_end_x_anti/F");
  fEventTree->Branch("true_kaon_daughter_end_y_anti", &true_kaon_daughter_end_y_anti, "true_kaon_daughter_end_y_anti/F");
  fEventTree->Branch("true_kaon_daughter_end_z_anti", &true_kaon_daughter_end_z_anti, "true_kaon_daughter_end_z_anti/F");
  fEventTree->Branch("true_kaon_daughter_end_inTPC_anti", &true_kaon_daughter_end_inTPC_anti, "true_kaon_daughter_endInTPC_anti/O");
  fEventTree->Branch("true_kaon_daughter_end_in5cmTPC_anti", &true_kaon_daughter_end_in5cmTPC_anti, "true_kaon_daughter_end_in5cmTPC_anti/O");
  fEventTree->Branch("true_kaon_daughter_end_inCCInclusiveTPC_anti", &true_kaon_daughter_end_inCCInclusiveTPC_anti, "true_kaon_daughter_end_inCCInclusiveTPC_anti/O");


  fEventTree->Branch("reco_nu_cc_filter", &reco_nu_cc_filter, "reco_nu_cc_filter/O");

  fEventTree->Branch("reco_nu_vtx_x", &reco_nu_vtx_x, "reco_nu_vtx_x/F");
  fEventTree->Branch("reco_nu_vtx_y", &reco_nu_vtx_y, "reco_nu_vtx_y/F");
  fEventTree->Branch("reco_nu_vtx_z", &reco_nu_vtx_z, "reco_nu_vtx_z/F");
  fEventTree->Branch("reco_nu_vtx_inTPC", &reco_nu_vtx_inTPC, "reco_nu_vtx_inTPC/O");
  fEventTree->Branch("reco_nu_vtx_in5cmTPC", &reco_nu_vtx_in5cmTPC, "reco_nu_vtx_in5cmTPC/O");
  fEventTree->Branch("reco_nu_vtx_inCCInclusiveTPC", &reco_nu_vtx_inCCInclusiveTPC, "reco_nu_vtx_inCCInclusiveTPC/O");
  fEventTree->Branch("reco_nu_ndaughters", &reco_nu_ndaughters, "reco_nu_ndaughters/I"          );
  fEventTree->Branch("reco_nu_cc_nmue", &reco_nu_cc_nmue, "reco_nu_cc_nmue/I"          );

  fEventTree->Branch("reco_ccmu_vtx_x", &reco_ccmu_vtx_x, "reco_ccmu_vtx_x/F");
  fEventTree->Branch("reco_ccmu_vtx_y", &reco_ccmu_vtx_y, "reco_ccmu_vtx_y/F");
  fEventTree->Branch("reco_ccmu_vtx_z", &reco_ccmu_vtx_z, "reco_ccmu_vtx_z/F");
  fEventTree->Branch("reco_ccmu_vtx_inTPC",&reco_ccmu_vtx_inTPC, "reco_ccmu_vtx_inTPC/O");
  fEventTree->Branch("reco_ccmu_vtx_in5cmTPC", &reco_ccmu_vtx_in5cmTPC, "reco_ccmu_vtx_in5cmTPC/O");
  fEventTree->Branch("reco_ccmu_vtx_inCCInclusiveTPC", &reco_ccmu_vtx_inCCInclusiveTPC, "reco_ccmu_vtx_inCCInclusiveTPC/O");
  fEventTree->Branch("reco_ccmu_true_pdg", &reco_ccmu_true_pdg, "reco_ccmu_true_pdg/I");
  fEventTree->Branch("reco_ccmu_true_origin", &reco_ccmu_true_origin, "reco_ccmu_true_origin/I");
  fEventTree->Branch("reco_ccmu_true_primary", &reco_ccmu_true_primary, "reco_ccmu_true_primary/O");
  fEventTree->Branch("reco_ccmu_true_end_inTPC", &reco_ccmu_true_end_inTPC, "reco_ccmu_true_end_inTPC/O");
  fEventTree->Branch("reco_ccmu_true_end_in5cmTPC", &reco_ccmu_true_end_in5cmTPC, "reco_ccmu_true_end_in5cmTPC/O");
  fEventTree->Branch("reco_ccmu_true_end_inCCInclusiveTPC", &reco_ccmu_true_end_inCCInclusiveTPC, "reco_ccmu_true_end_inCCInclusiveTPC/O");
  fEventTree->Branch("reco_ccmu_true_length", &reco_ccmu_true_length, "reco_ccmu_true_length/F");

  fEventTree->Branch("reco_ntracks", &reco_ntracks, "reco_ntracks/I");
  fEventTree->Branch("reco_nshowers", &reco_ntracks, "reco_nshowers/I");
  fEventTree->Branch("nTracks", &nTracks, "nTracks/I");
  fEventTree->Branch("nShowers", &nShowers, "nShowers/I");

  fEventTree->Branch("reco_track_start_x", &reco_track_start_x, "reco_track_start_x[20]/F");
  fEventTree->Branch("reco_track_start_y", &reco_track_start_y, "reco_track_start_y[20]/F");
  fEventTree->Branch("reco_track_start_z", &reco_track_start_z, "reco_track_start_z[20]/F");
  fEventTree->Branch("reco_track_end_x", &reco_track_end_x, "reco_track_end_x[20]/F");
  fEventTree->Branch("reco_track_end_y", &reco_track_end_y, "reco_track_end_y[20]/F");
  fEventTree->Branch("reco_track_end_z", &reco_track_end_z, "reco_track_end_z[20]/F");


  fEventTree->Branch("reco_track_distance", &reco_track_distance, "reco_track_distance[20]/F");
  fEventTree->Branch("reco_track_nhits0", &reco_track_nhits0, "reco_track_nhits0[20]/I");
  fEventTree->Branch("reco_track_nhits1", &reco_track_nhits1, "reco_track_nhits1[20]/I");
  fEventTree->Branch("reco_track_nhits2", &reco_track_nhits2, "reco_track_nhits2[20]/I");

  fEventTree->Branch("reco_track_kin0", &reco_track_kin0, "reco_track_kin0[20]/F");
  fEventTree->Branch("reco_track_kin1", &reco_track_kin1, "reco_track_kin1[20]/F");
  fEventTree->Branch("reco_track_kin2", &reco_track_kin2, "reco_track_kin2[20]/F");

  fEventTree->Branch("reco_track_dEdx_pl0", &reco_track_dEdx_pl0,"reco_track_dEdx_pl0[20][2000]");
  fEventTree->Branch("reco_track_ResRan_pl0", &reco_track_ResRan_pl0,"reco_track_ResRan_pl0[20][2000]");

  fEventTree->Branch("reco_track_dEdx_pl1", &reco_track_dEdx_pl1,"reco_track_dEdx_pl1[20][2000]");
  fEventTree->Branch("reco_track_ResRan_pl1", &reco_track_ResRan_pl1,"reco_track_ResRan_pl1[20][2000]");

  fEventTree->Branch("reco_track_dEdx_pl2", &reco_track_dEdx_pl2,"reco_track_dEdx_pl2[20][2000]");
  fEventTree->Branch("reco_track_ResRan_pl2", &reco_track_ResRan_pl2,"reco_track_ResRan_pl2[20][2000]");

  fEventTree->Branch("reco_track_length", &reco_track_length, "reco_track_length[20]/F");
  fEventTree->Branch("reco_track_theta", &reco_track_theta, "reco_track_theta[20]/F");
  fEventTree->Branch("reco_track_phi", &reco_track_phi, "reco_track_phi[20]/F");
  fEventTree->Branch("reco_track_dir", &reco_track_dir, "reco_track_dir[20]/O");

  fEventTree->Branch("reco_track_P_vtx", &reco_track_P_vtx, "reco_track_P_vtx[20]/F");
  fEventTree->Branch("reco_track_P_str", &reco_track_P_str, "reco_track_P_str[20]/F");
  fEventTree->Branch("reco_track_P_end", &reco_track_P_end, "reco_track_P_end[20]/F");

  fEventTree->Branch("reco_track_chi2ka_pl0", &reco_track_chi2ka_pl0, "reco_track_chi2ka_pl0[20]/F");
  fEventTree->Branch("reco_track_chi2pr_pl0", &reco_track_chi2pr_pl0, "reco_track_chi2pr_pl0[20]/F");
  fEventTree->Branch("reco_track_chi2pi_pl0", &reco_track_chi2pi_pl0, "reco_track_chi2pi_pl0[20]/F");
  fEventTree->Branch("reco_track_chi2mu_pl0", &reco_track_chi2mu_pl0, "reco_track_chi2mu_pl0[20]/F");
  fEventTree->Branch("reco_track_chi2ka_pl1", &reco_track_chi2ka_pl1, "reco_track_chi2ka_pl1[20]/F");
  fEventTree->Branch("reco_track_chi2pr_pl1", &reco_track_chi2pr_pl1, "reco_track_chi2pr_pl1[20]/F");
  fEventTree->Branch("reco_track_chi2pi_pl1", &reco_track_chi2pi_pl1, "reco_track_chi2pi_pl1[20]/F");
  fEventTree->Branch("reco_track_chi2mu_pl1", &reco_track_chi2mu_pl1, "reco_track_chi2mu_pl1[20]/F");
  fEventTree->Branch("reco_track_chi2ka_pl2", &reco_track_chi2ka_pl2, "reco_track_chi2ka_pl2[20]/F");
  fEventTree->Branch("reco_track_chi2pr_pl2", &reco_track_chi2pr_pl2, "reco_track_chi2pr_pl2[20]/F");
  fEventTree->Branch("reco_track_chi2pi_pl2", &reco_track_chi2pi_pl2, "reco_track_chi2pi_pl2[20]/F");
  fEventTree->Branch("reco_track_chi2mu_pl2", &reco_track_chi2mu_pl2, "reco_track_chi2mu_pl2[20]/F");

  fEventTree->Branch("reco_track_Bragg_fwd_ka_pl0", &reco_track_Bragg_fwd_ka_pl0, "reco_track_Bragg_fwd_ka_pl0[20]/F");
  fEventTree->Branch("reco_track_Bragg_fwd_pr_pl0", &reco_track_Bragg_fwd_pr_pl0, "reco_track_Bragg_fwd_pr_pl0[20]/F");
  fEventTree->Branch("reco_track_Bragg_fwd_pi_pl0", &reco_track_Bragg_fwd_pi_pl0, "reco_track_Bragg_fwd_pi_pl0[20]/F");
  fEventTree->Branch("reco_track_Bragg_fwd_mu_pl0", &reco_track_Bragg_fwd_mu_pl0, "reco_track_Bragg_fwd_mu_pl0[20]/F");
  fEventTree->Branch("reco_track_Bragg_fwd_ka_pl1", &reco_track_Bragg_fwd_ka_pl1, "reco_track_Bragg_fwd_ka_pl1[20]/F");
  fEventTree->Branch("reco_track_Bragg_fwd_pr_pl1", &reco_track_Bragg_fwd_pr_pl1, "reco_track_Bragg_fwd_pr_pl1[20]/F");
  fEventTree->Branch("reco_track_Bragg_fwd_pi_pl1", &reco_track_Bragg_fwd_pi_pl1, "reco_track_Bragg_fwd_pi_pl1[20]/F");
  fEventTree->Branch("reco_track_Bragg_fwd_mu_pl1", &reco_track_Bragg_fwd_mu_pl1, "reco_track_Bragg_fwd_mu_pl1[20]/F");
  fEventTree->Branch("reco_track_Bragg_fwd_ka_pl2", &reco_track_Bragg_fwd_ka_pl2, "reco_track_Bragg_fwd_ka_pl2[20]/F");
  fEventTree->Branch("reco_track_Bragg_fwd_pr_pl2", &reco_track_Bragg_fwd_pr_pl2, "reco_track_Bragg_fwd_pr_pl2[20]/F");
  fEventTree->Branch("reco_track_Bragg_fwd_pi_pl2", &reco_track_Bragg_fwd_pi_pl2, "reco_track_Bragg_fwd_pi_pl2[20]/F");
  fEventTree->Branch("reco_track_Bragg_fwd_mu_pl2", &reco_track_Bragg_fwd_mu_pl2, "reco_track_Bragg_fwd_mu_pl2[20]/F");


  fEventTree->Branch("reco_track_MIP_pl0", &reco_track_MIP_pl0, "reco_track_MIP_pl0[20]/F");
  fEventTree->Branch("reco_track_MIP_pl1", &reco_track_MIP_pl1, "reco_track_MIP_pl1[20]/F");
  fEventTree->Branch("reco_track_MIP_pl2", &reco_track_MIP_pl2, "reco_track_MIP_pl2[20]/F");


  fEventTree->Branch("reco_track_chi2ka_3pl", &reco_track_chi2ka_3pl, "reco_track_chi2ka_3pl[20]/F");
  fEventTree->Branch("reco_track_chi2pr_3pl", &reco_track_chi2pr_3pl, "reco_track_chi2pr_3pl[20]/F");
  fEventTree->Branch("reco_track_chi2pi_3pl", &reco_track_chi2pi_3pl, "reco_track_chi2pi_3pl[20]/F");
  fEventTree->Branch("reco_track_chi2mu_3pl", &reco_track_chi2mu_3pl, "reco_track_chi2mu_3pl[20]/F");
  fEventTree->Branch("reco_track_likepr_3pl", &reco_track_likepr_3pl, "reco_track_likepr_3pl[20]/F");
  fEventTree->Branch("reco_track_llrpid_3pl", &reco_track_llrpid_3pl, "reco_track_llrpid_3pl[20]/F");
  fEventTree->Branch("reco_track_total_llrpid_3pl", &reco_track_total_llrpid_3pl, "reco_track_total_llrpid_3pl[20]/F");
  fEventTree->Branch("reco_track_llrpid_k_3pl", &reco_track_llrpid_k_3pl, "reco_track_llrpid_k_3pl[20]/F");
  fEventTree->Branch("reco_track_vtx_inTPC", &reco_track_vtx_inTPC, "reco_track_vtx_inTPC[20]/O");
  fEventTree->Branch("reco_track_vtx_in5cmTPC", &reco_track_vtx_in5cmTPC, "reco_track_vtx_in5cmTPC[20]/O");
  fEventTree->Branch("reco_track_vtx_inCCInclusiveTPC", &reco_track_vtx_inCCInclusiveTPC, "reco_track_vtx_inCCInclusiveTPC[20]/O");
  fEventTree->Branch("reco_track_end_inTPC", &reco_track_end_inTPC, "reco_track_end_inTPC[20]/O");
  fEventTree->Branch("reco_track_end_in5cmTPC", &reco_track_end_in5cmTPC, "reco_track_end_in5cmTPC[20]/O");
  fEventTree->Branch("reco_track_end_inCCInclusiveTPC", &reco_track_end_inCCInclusiveTPC, "reco_track_end_inCCInclusiveTPC[20]/O");
  fEventTree->Branch("reco_track_true_pdg", &reco_track_true_pdg, "reco_track_true_pdg[20]/I");
  fEventTree->Branch("reco_track_true_origin", &reco_track_true_origin, "reco_track_true_origin[20]/I");
  fEventTree->Branch("reco_track_true_primary", &reco_track_true_primary, "reco_track_true_primary[20]/O");
  fEventTree->Branch("reco_track_true_end_inTPC", &reco_track_true_end_inTPC, "reco_track_true_end_inTPC[20]/O");
  fEventTree->Branch("reco_track_true_end_in5cmTPC", &reco_track_true_end_in5cmTPC, "reco_track_true_end_in5cmTPC[20]/O");
  fEventTree->Branch("reco_track_true_end_inCCInclusiveTPC", &reco_track_true_end_inCCInclusiveTPC, "reco_track_true_end_inCCInclusiveTPC[20]/O");
  fEventTree->Branch("reco_track_true_length", &reco_track_true_length, "reco_track_true_length[20]/F");


  fEventTree->Branch("reco_track_match_e", &reco_track_match_e, "reco_track_match_e[20][10]/F");
  fEventTree->Branch("reco_track_match_hit", &reco_track_match_hit, "reco_track_match_hit[20][10]/F");
  fEventTree->Branch("reco_track_match_epdg", &reco_track_match_epdg, "reco_track_match_epdg[20][10]/I");
  fEventTree->Branch("reco_track_match_hitpdg", &reco_track_match_hitpdg, "reco_track_match_hitpdg[20][10]/I");

  fEventTree->Branch("reco_track_daughter_match_e", &reco_track_daughter_match_e, "reco_track_daughter_match_e[20][20][10]/F");
  fEventTree->Branch("reco_track_daughter_match_hit", &reco_track_daughter_match_hit, "reco_track_daughter_match_hit[20][20][10]/F");
  fEventTree->Branch("reco_track_daughter_match_epdg", &reco_track_daughter_match_epdg, "reco_track_daughter_match_epdg[20][20][10]/I");
  fEventTree->Branch("reco_track_daughter_match_hitpdg", &reco_track_daughter_match_hitpdg, "reco_track_daughter_match_hitpdg[20][20][10]/I");

  fEventTree->Branch("reco_track_daughter_shower_match_e", &reco_track_daughter_shower_match_e, "reco_track_daughter_shower_match_e[20][20][10]/F");
  fEventTree->Branch("reco_track_daughter_shower_match_hit", &reco_track_daughter_shower_match_hit, "reco_track_daughter_shower_match_hit[20][20][10]/F");
  fEventTree->Branch("reco_track_daughter_shower_match_epdg", &reco_track_daughter_shower_match_epdg, "reco_track_daughter_shower_match_epdg[20][20][10]/I");
  fEventTree->Branch("reco_track_daughter_shower_match_hitpdg", &reco_track_daughter_shower_match_hitpdg, "reco_track_daughter_shower_match_hitpdg[20][20][10]/I");

  fEventTree->Branch("reco_track_true_pdg_sh", &reco_track_true_pdg_sh, "reco_track_true_pdg_sh[20]/I");
  fEventTree->Branch("reco_track_true_origin_sh", &reco_track_true_origin_sh, "reco_track_true_origin_sh[20]/I");
  fEventTree->Branch("reco_track_true_primary_sh", &reco_track_true_primary_sh, "reco_track_true_primary_sh[20]/O");
  fEventTree->Branch("reco_track_true_end_inTPC_sh", &reco_track_true_end_inTPC_sh, "reco_track_true_end_inTPC_sh[20]/O");
  fEventTree->Branch("reco_track_true_end_in5cmTPC_sh", &reco_track_true_end_in5cmTPC_sh, "reco_track_true_end_in5cmTPC_sh[20]/O");
  fEventTree->Branch("reco_track_true_end_inCCInclusiveTPC_sh", &reco_track_true_end_inCCInclusiveTPC_sh, "reco_track_true_end_inCCInclusiveTPC[20]_sh/O");
  fEventTree->Branch("reco_track_true_length_sh", &reco_track_true_length_sh, "reco_track_true_length_sh[20]/F");

  fEventTree->Branch("reco_track_daughter_true_pdg_sh", &reco_track_daughter_true_pdg_sh, "reco_track_daughter_true_pdg_sh[20][20]/I");
  fEventTree->Branch("reco_track_daughter_true_origin_sh", &reco_track_daughter_true_origin_sh, "reco_track_daughter_true_origin_sh[20][20]/I");
  fEventTree->Branch("reco_track_daughter_true_primary_sh", &reco_track_daughter_true_primary_sh, "reco_track_daughter_true_primary[20][20]_sh/O");
  fEventTree->Branch("reco_track_daughter_true_end_inTPC_sh", &reco_track_daughter_true_end_inTPC_sh, "reco_track_daughter_true_end_inTPC[20][20]_sh/O");
  fEventTree->Branch("reco_track_daughter_true_end_in5cmTPC_sh", &reco_track_daughter_true_end_in5cmTPC_sh, "reco_track_daughter_true_end_in5cmTPC[20][20]_sh/O");
  fEventTree->Branch("reco_track_daughter_true_end_inCCInclusiveTPC_sh", &reco_track_daughter_true_end_inCCInclusiveTPC_sh, "reco_track_daughter_true_end_inCCInclusiveTPC_sh[20][20]/O");
  fEventTree->Branch("reco_track_daughter_true_length_sh", &reco_track_daughter_true_length_sh, "reco_track_daughter_true_length_sh[20][20]/F");
  fEventTree->Branch("reco_track_daughter_true_mother_sh", &reco_track_daughter_true_mother_sh, "reco_track_daughter_true_mother_sh[20][20]/I");

  fEventTree->Branch("reco_track_daughter_start_x", &reco_track_daughter_start_x, "reco_track_daughter_start_x[20][20]/F");
  fEventTree->Branch("reco_track_daughter_start_y", &reco_track_daughter_start_y, "reco_track_daughter_start_y[20][20]/F");
  fEventTree->Branch("reco_track_daughter_start_z", &reco_track_daughter_start_z, "reco_track_daughter_start_z[20][20]/F");
  fEventTree->Branch("reco_track_daughter_end_x", &reco_track_daughter_end_x, "reco_track_daughter_end_x[20][20]/F");
  fEventTree->Branch("reco_track_daughter_end_y", &reco_track_daughter_end_y, "reco_track_daughter_end_y[20][20]/F");
  fEventTree->Branch("reco_track_daughter_end_z", &reco_track_daughter_end_z, "reco_track_daughter_end_z[20][20]/F");

  fEventTree->Branch("reco_shower_ndaughters", &reco_shower_ndaughters, "reco_shower_ndaughters[20]/I");
  fEventTree->Branch("reco_track_daughter_distance_sh", &reco_track_daughter_distance_sh, "reco_track_daughter_distance_sh[20][20]/F");
  fEventTree->Branch("reco_track_daughter_vtx_distance_sh", &reco_track_daughter_vtx_distance_sh, "reco_track_daughter_vtx_distance_sh[20][20]/F");
  fEventTree->Branch("reco_angle_track_daughter_sh", &reco_angle_track_daughter_sh, "reco_angle_track_daughter_sh[20][20]/F");
  fEventTree->Branch("reco_angle_daughter_track_daughter_sh", &reco_angle_daughter_track_daughter_sh, "reco_angle_daughter_track_daughter_sh[20][20]/F");
  fEventTree->Branch("reco_track_daughter_length_sh", &reco_track_daughter_length_sh, "reco_track_daughter_length_sh[20][20]/F");
  fEventTree->Branch("reco_track_daughter_theta_sh", &reco_track_daughter_theta_sh, "reco_track_daughter_theta_sh[20][20]/F");
  fEventTree->Branch("reco_track_daughter_phi_sh", &reco_track_daughter_phi_sh, "reco_track_daughter_phi_sh[20][20]/F");
  fEventTree->Branch("reco_track_daughter_open_angle_sh", &reco_track_daughter_open_angle_sh, "reco_track_daughter_open_angle_sh[20][20]/F");
  fEventTree->Branch("reco_track_daughter_dedx_pl0_sh", &reco_track_daughter_dedx_pl0_sh, "reco_track_daughter_dedx_pl0_sh[20][20]/F");
  fEventTree->Branch("reco_track_daughter_dedx_pl1_sh", &reco_track_daughter_dedx_pl1_sh, "reco_track_daughter_dedx_pl1_sh[20][20]/F");
  fEventTree->Branch("reco_track_daughter_dedx_pl2_sh", &reco_track_daughter_dedx_pl2_sh, "reco_track_daughter_dedx_pl2_sh[20][20]/F");
  fEventTree->Branch("reco_track_daughter_energy_pl0_sh", &reco_track_daughter_energy_pl0_sh, "reco_track_daughter_energy_pl0_sh[20][20]/F");
  fEventTree->Branch("reco_track_daughter_energy_pl1_sh", &reco_track_daughter_energy_pl1_sh, "reco_track_daughter_energy_pl1_sh[20][20]/F");
  fEventTree->Branch("reco_track_daughter_energy_pl2_sh", &reco_track_daughter_energy_pl2_sh, "reco_track_daughter_energy_pl2_sh[20][20]/F");

  fEventTree->Branch("reco_track_daughter_vtx_inTPC_sh", &reco_track_daughter_vtx_inTPC_sh, "reco_track_daughter_vtx_inTPC_sh[20][20]/O");
  fEventTree->Branch("reco_track_daughter_vtx_in5cmTPC_sh", &reco_track_daughter_vtx_in5cmTPC_sh, "reco_track_daughter_vtx_in5cmTPC_sh[20][20]/O");
  fEventTree->Branch("reco_track_daughter_vtx_inCCInclusiveTPC_sh", &reco_track_daughter_vtx_inCCInclusiveTPC_sh, "reco_track_daughter_vtx_inCCInclusiveTPC_sh[20][20]/O");
  fEventTree->Branch("reco_track_daughter_end_inTPC_sh", &reco_track_daughter_end_inTPC_sh, "reco_track_daughter_end_inTPC_sh[20][20]/O");
  fEventTree->Branch("reco_track_daughter_end_in5cmTPC_sh", &reco_track_daughter_end_in5cmTPC_sh, "reco_track_daughter_end_in5cmTPC_sh[20][20]/O");
  fEventTree->Branch("reco_track_daughter_end_inCCInclusiveTPC_sh", &reco_track_daughter_end_inCCInclusiveTPC_sh, "reco_track_daughter_end_inCCInclusiveTPC_sh[20][20]/O");

  fEventTree->Branch("reco_track_ndaughters", &reco_track_ndaughters, "reco_track_ndaughters[20]/I");
  fEventTree->Branch("reco_track_daughter_distance", &reco_track_daughter_distance, "reco_track_daughter_distance[20][20]/F");
  fEventTree->Branch("reco_track_daughter_vtx_distance", &reco_track_daughter_vtx_distance, "reco_track_daughter_vtx_distance[20][20]/F");
  fEventTree->Branch("reco_angle_track_daughter", &reco_angle_track_daughter, "reco_angle_track_daughter[20][20]/F");
  fEventTree->Branch("reco_track_daughter_nhits0", &reco_track_daughter_nhits0, "reco_track_daughter_nhits0[20][20]/I");
  fEventTree->Branch("reco_track_daughter_nhits1", &reco_track_daughter_nhits1, "reco_track_daughter_nhits1[20][20]/I");
  fEventTree->Branch("reco_track_daughter_nhits2", &reco_track_daughter_nhits2, "reco_track_daughter_nhits2[20][20]/I");

  //fEventTree->Branch("reco_track_daughter_dEdx", &reco_track_daughter_dEdx);                                                                                                                                                                                                                                                
  //fEventTree->Branch("reco_track_daughter_ResRan", &reco_track_daughter_ResRan);                                                                                                                                                                                                                                            
  //                                                                                                                                                                                                                                                                                                                          
  fEventTree->Branch("reco_track_daughter_dEdx_pl0", &reco_track_daughter_dEdx_pl0, "reco_track_daughter_dEdx_pl0[20][20][2000]/F");
  fEventTree->Branch("reco_track_daughter_ResRan_pl0", &reco_track_daughter_ResRan_pl0, "reco_track_daughter_ResRan_pl0[20][20][2000]/F");

  fEventTree->Branch("reco_track_daughter_dEdx_pl1", &reco_track_daughter_dEdx_pl1, "reco_track_daughter_dEdx_pl1[20][20][2000]/F");
  fEventTree->Branch("reco_track_daughter_ResRan_pl1", &reco_track_daughter_ResRan_pl1, "reco_track_daughter_ResRan_pl1[20][20][2000]/F");

  fEventTree->Branch("reco_track_daughter_dEdx_pl2", &reco_track_daughter_dEdx_pl2, "reco_track_daughter_dEdx_pl2[20][20][2000]/F");
  fEventTree->Branch("reco_track_daughter_ResRan_pl2", &reco_track_daughter_ResRan_pl2, "reco_track_daughter_ResRan_pl2[20][20][2000]/F");

  fEventTree->Branch("reco_track_daughter_length", &reco_track_daughter_length, "reco_track_daughter_length[20][20]/F");
  fEventTree->Branch("reco_track_daughter_theta", &reco_track_daughter_theta, "reco_track_daughter_theta[20][20]/F");
  fEventTree->Branch("reco_track_daughter_phi", &reco_track_daughter_phi, "reco_track_daughter_phi[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2ka_pl0", &reco_track_daughter_chi2ka_pl0, "reco_track_daughter_chi2ka_pl0[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2pr_pl0", &reco_track_daughter_chi2pr_pl0, "reco_track_daughter_chi2pr_pl0[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2pi_pl0", &reco_track_daughter_chi2pi_pl0, "reco_track_daughter_chi2pi_pl0[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2mu_pl0", &reco_track_daughter_chi2mu_pl0, "reco_track_daughter_chi2mu_pl0[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2ka_pl1", &reco_track_daughter_chi2ka_pl1, "reco_track_daughter_chi2ka_pl1[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2pr_pl1", &reco_track_daughter_chi2pr_pl1, "reco_track_daughter_chi2pr_pl1[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2pi_pl1", &reco_track_daughter_chi2pi_pl1, "reco_track_daughter_chi2pi_pl1[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2mu_pl1", &reco_track_daughter_chi2mu_pl1, "reco_track_daughter_chi2mu_pl1[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2ka_pl2", &reco_track_daughter_chi2ka_pl2, "reco_track_daughter_chi2ka_pl2[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2pr_pl2", &reco_track_daughter_chi2pr_pl2, "reco_track_daughter_chi2pr_pl2[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2pi_pl2", &reco_track_daughter_chi2pi_pl2, "reco_track_daughter_chi2pi_pl2[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2mu_pl2", &reco_track_daughter_chi2mu_pl2, "reco_track_daughter_chi2mu_pl2[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2ka_3pl", &reco_track_daughter_chi2ka_3pl, "reco_track_daughter_chi2ka_3pl[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2pr_3pl", &reco_track_daughter_chi2pr_3pl, "reco_track_daughter_chi2pr_3pl[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2pi_3pl", &reco_track_daughter_chi2pi_3pl, "reco_track_daughter_chi2pi_3pl[20][20]/F");
  fEventTree->Branch("reco_track_daughter_chi2mu_3pl", &reco_track_daughter_chi2mu_3pl, "reco_track_daughter_chi2mu_3pl[20][20]/F");
  fEventTree->Branch("reco_track_daughter_likepr_3pl", &reco_track_daughter_likepr_3pl, "reco_track_daughter_likepr_3pl[20][20]/F");
  fEventTree->Branch("reco_track_daughter_llrpid_3pl", &reco_track_daughter_llrpid_3pl, "reco_track_daughter_llrpid_3pl[20][20]/F");
  fEventTree->Branch("reco_track_daughter_llrpid_k_3pl", &reco_track_daughter_llrpid_k_3pl, "reco_track_daughter_llrpid_k_3pl[20][20]/F");
  fEventTree->Branch("reco_track_daughter_vtx_inTPC", &reco_track_daughter_vtx_inTPC, "reco_track_daughter_vtx_inTPC[20][20]/O");
  fEventTree->Branch("reco_track_daughter_vtx_in5cmTPC", &reco_track_daughter_vtx_in5cmTPC, "reco_track_daughter_vtx_in5cmTPC[20][20]/O");
  fEventTree->Branch("reco_track_daughter_vtx_inCCInclusiveTPC", &reco_track_daughter_vtx_inCCInclusiveTPC, "reco_track_daughter_vtx_inCCInclusiveTPC[20][20]/O");
  fEventTree->Branch("reco_track_daughter_end_inTPC", &reco_track_daughter_end_inTPC, "reco_track_daughter_end_inTPC[20][20]/O");
  fEventTree->Branch("reco_track_daughter_end_in5cmTPC", &reco_track_daughter_end_in5cmTPC, "reco_track_daughter_end_in5cmTPC[20][20]/O");
  fEventTree->Branch("reco_track_daughter_end_inCCInclusiveTPC", &reco_track_daughter_end_inCCInclusiveTPC, "reco_track_daughter_end_inCCInclusiveTPC[20][20]/O");
  fEventTree->Branch("reco_track_daughter_true_pdg", &reco_track_daughter_true_pdg, "reco_track_daughter_true_pdg[20][20]/I");
  fEventTree->Branch("reco_track_daughter_true_origin", &reco_track_daughter_true_origin, "reco_track_daughter_true_origin[20][20]/I");
  fEventTree->Branch("reco_track_daughter_true_primary", &reco_track_daughter_true_primary, "reco_track_daughter_true_primary[20][20]/O");
  fEventTree->Branch("reco_track_daughter_true_end_inTPC", &reco_track_daughter_true_end_inTPC, "reco_track_daughter_true_end_inTPC[20][20]/O");
  fEventTree->Branch("reco_track_daughter_true_end_in5cmTPC", &reco_track_daughter_true_end_in5cmTPC, "reco_track_daughter_true_end_in5cmTPC[20][20]/O");
  fEventTree->Branch("reco_track_daughter_true_end_inCCInclusiveTPC", &reco_track_daughter_true_end_inCCInclusiveTPC, "reco_track_daughter_true_end_inCCInclusiveTPC[20][20]/O");
  fEventTree->Branch("reco_track_daughter_true_length", &reco_track_daughter_true_length, "reco_track_daughter_true_length[20][20]/F");
  fEventTree->Branch("reco_track_daughter_true_mother", &reco_track_daughter_true_mother, "reco_track_daughter_true_mother[20][20]/I");


  fEventTree->Branch("reco_track_daughter_Bragg_fwd_ka_pl0", &reco_track_daughter_Bragg_fwd_ka_pl0, "reco_track_daughter_Bragg_fwd_ka_pl0[20][20]/F");
  fEventTree->Branch("reco_track_daughter_Bragg_fwd_pr_pl0", &reco_track_daughter_Bragg_fwd_pr_pl0, "reco_track_daughter_Bragg_fwd_pr_pl0[20][20]/F");
  fEventTree->Branch("reco_track_daughter_Bragg_fwd_pi_pl0", &reco_track_daughter_Bragg_fwd_pi_pl0, "reco_track_daughter_Bragg_fwd_pi_pl0[20][20]/F");
  fEventTree->Branch("reco_track_daughter_Bragg_fwd_mu_pl0", &reco_track_daughter_Bragg_fwd_mu_pl0, "reco_track_daughter_Bragg_fwd_mu_pl0[20][20]/F");
  fEventTree->Branch("reco_track_daughter_Bragg_fwd_ka_pl1", &reco_track_daughter_Bragg_fwd_ka_pl1, "reco_track_daughter_Bragg_fwd_ka_pl1[20][20]/F");
  fEventTree->Branch("reco_track_daughter_Bragg_fwd_pr_pl1", &reco_track_daughter_Bragg_fwd_pr_pl1, "reco_track_daughter_Bragg_fwd_pr_pl1[20][20]/F");
  fEventTree->Branch("reco_track_daughter_Bragg_fwd_pi_pl1", &reco_track_daughter_Bragg_fwd_pi_pl1, "reco_track_daughter_Bragg_fwd_pi_pl1[20][20]/F");
  fEventTree->Branch("reco_track_daughter_Bragg_fwd_mu_pl1", &reco_track_daughter_Bragg_fwd_mu_pl1, "reco_track_daughter_Bragg_fwd_mu_pl1[20][20]/F");
  fEventTree->Branch("reco_track_daughter_Bragg_fwd_ka_pl2", &reco_track_daughter_Bragg_fwd_ka_pl2, "reco_track_daughter_Bragg_fwd_ka_pl2[20][20]/F");
  fEventTree->Branch("reco_track_daughter_Bragg_fwd_pr_pl2", &reco_track_daughter_Bragg_fwd_pr_pl2, "reco_track_daughter_Bragg_fwd_pr_pl2[20][20]/F");
  fEventTree->Branch("reco_track_daughter_Bragg_fwd_pi_pl2", &reco_track_daughter_Bragg_fwd_pi_pl2, "reco_track_daughter_Bragg_fwd_pi_pl2[20][20]/F");
  fEventTree->Branch("reco_track_daughter_Bragg_fwd_mu_pl2", &reco_track_daughter_Bragg_fwd_mu_pl2, "reco_track_daughter_Bragg_fwd_mu_pl2[20][20]/F");


  fEventTree->Branch("reco_track_daughter_MIP_pl0", &reco_track_daughter_MIP_pl0, "reco_track_daughter_MIP_pl0[20][20]/F");
  fEventTree->Branch("reco_track_daughter_MIP_pl1", &reco_track_daughter_MIP_pl1, "reco_track_daughter_MIP_pl1[20][20]/F");
  fEventTree->Branch("reco_track_daughter_MIP_pl2", &reco_track_daughter_MIP_pl2, "reco_track_daughter_MIP_pl2[20][20]/F");



  fEventTree->Branch("k_can_trkid", &k_can_trkid,"k_can_trkid/I");
  fEventTree->Branch("mu_can_trkid", &mu_can_trkid,"mu_can_trkid/I");
  fEventTree->Branch("k_mu_can_dis", &k_mu_can_dis,"k_mu_can_dis/F");
  fEventTree->Branch("k_mu_open_angle", &k_mu_open_angle,"k_mu_open_angle/F");
  fEventTree->Branch("k_vtx_dis", &k_vtx_dis,"k_vtx_dis/F");

  fEventTree->Branch("cut_1", &cut_1,"cut_1/I");
  fEventTree->Branch("cut_2", &cut_2,"cut_2/I");
  fEventTree->Branch("cut_3", &cut_3,"cut_3/I");
  fEventTree->Branch("cut_4", &cut_4,"cut_4/I");
  fEventTree->Branch("cut_5", &cut_5,"cut_5/I");
  fEventTree->Branch("cut_6", &cut_6,"cut_6/I");
  fEventTree->Branch("cut_7", &cut_7,"cut_7/I");
  fEventTree->Branch("cut_8", &cut_8,"cut_8/I");
  fEventTree->Branch("cut_9", &cut_9,"cut_9/I");
  fEventTree->Branch("cut_10", &cut_10,"cut_10/I");
  fEventTree->Branch("cut_11", &cut_11,"cut_11/I");
  fEventTree->Branch("cut_12", &cut_12,"cut_12/I");
  fEventTree->Branch("cut_13", &cut_13,"cut_13/I");
  fEventTree->Branch("PFP_have_nuslice",& PFP_have_nuslice);
  fSubrunTree = tfs->make<TTree>("subruns", "SubRun Tree");
  fSubrunTree->Branch("run", &m_run, "run/i");
  fSubrunTree->Branch("subRun", &m_subrun, "subRun/i");
  //if (!m_isData)                                                                                                                                                                                                                                                                                                            
  fSubrunTree->Branch("pot", &m_pot, "pot/F");

  vertex_tree = tfs->make<TTree>("vertex_tree", "vertex_tree");
  pot_tree = tfs->make<TTree>("pot_tree", "pot_tree");
  eventweight_tree = tfs->make<TTree>("eventweight_tree", "eventweight_tree");
  ncdelta_slice_tree = tfs->make<TTree>("ncdelta_slice_tree", "ncdelta_slice_tree");

  this->CreateTrackBranches();
  this->CreateShowerBranches();
  this->CreateMCTruthBranches();
  this->CreateSliceBranches();

}
