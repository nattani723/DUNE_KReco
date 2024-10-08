include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
#link_libraries( -L$ENV{BOOST_LIB} -lboost_system -lboost_filesystem)
art_make( 
	  PandoraPFA::PandoraSDK
	  PandoraPFA::PandoraMonitoring
	  LIB_LIBRARIES 
	  larcore_Geometry_Geometry_service
	  larsim_Simulation lardataobj_Simulation 
	  lardataalg_DetectorInfo
	  lardataobj_RawData
	  lardataobj_RecoBase
	  lardataobj_AnalysisBase
	  lardata_Utilities
	  nusimdata::SimulationBase
	  
	  MODULE_LIBRARIES 
	  larcorealg::Geometry
	  larcore::Geometry_Geometry_service
	  larsim::Simulation
	  nug4::ParticleNavigation 
	  lardataobj::Simulation
	  larsim::MCCheater_BackTrackerService_service
	  larsim::MCCheater_ParticleInventoryService_service
          larsim::Utils
	  lardata::Utilities
	  larevt::Filters
	  lardataobj::RawData
	  lardataobj::RecoBase
	  lardataobj::AnalysisBase
	  lardata::RecoObjects
          lardata::ArtDataHelper
          lardataalg::DetectorInfo
	  larreco::RecoAlg
	  dunereco::HitFinderDUNE
	  larpandora::LArPandoraInterface
	  dunecore::SignalShapingServiceDUNE10kt_service	
	  dunecore::SignalShapingServiceDUNE35t_service
          nusimdata::SimulationBase
          ART_FRAMEWORK_CORE
	  ART_FRAMEWORK_PRINCIPAL
	  ART_FRAMEWORK_SERVICES_REGISTRY
	  ART_ROOT_IO_TFILE_SUPPORT
	  ROOT::Core
	  ART_ROOT_IO_TFILESERVICE_SERVICE
	  art::Persistency_Common canvas
	  art::Persistency_Provenance canvas
	  art::Utilities canvas
	  messagefacility::MF_MessageLogger
	  fhiclcpp::fhiclcpp
	  cetlib::cetlib cetlib_except
          ROOT_GEOM
          ROOT_XMLIO
          ROOT_GDML
	  ROOT_BASIC_LIB_LIST
	  dunecore::Geometry
	  ROOT::Geom
        )

install_headers()
install_fhicl()
install_source()
install_scripts()
