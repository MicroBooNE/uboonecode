
# for classes_def.xml!!
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
include_directories( $ENV{NUTOOLS_INC} )
cet_find_library( SQLITE3 NAMES sqlite3 PATHS ENV NO_DEFAULT_PATH )
art_make( BASENAME_ONLY
	    LIBRARY_NAME  uboonecode_uboone_MuCS
      MODULE_LIBRARIES larcore_Geometry
			larcore_Geometry_Geometry_service
			larsim_Simulation
			larsim_MCCheater_BackTracker_service
      lardata_Utilities
			larevt_Filters
			lardata_RawData
			lardata_RecoBase
			lardata_RecoObjects  
			larreco_RecoAlg
			uboonecode_uboone_MuCS
			SignalShapingServiceMicroBooNE_service
			lardata_AnalysisBase
      ${SIMULATIONBASE}
      ${LARRECO_LIB}
      ${LARDATA_LIB}
      ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			${ART_PERSISTENCY_COMMON}
			${ART_PERSISTENCY_PROVENANCE}
			${ART_UTILITIES}
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
      ${FHICLCPP}
			${CETLIB}
      ${ROOT_GEOM}
      ${ROOT_XMLIO}
      ${ROOT_GDML}
			${ROOT_BASIC_LIB_LIST}
			${SQLITE3}
        )

install_fhicl()
install_source()
install_headers()
