
# for classes_def.xml!!
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
include_directories( $ENV{NUTOOLS_INC} )
cet_find_library( SQLITE3 NAMES sqlite3 PATHS ENV NO_DEFAULT_PATH )
art_make( BASENAME_ONLY
	    LIBRARY_NAME  uboonecode_uboone_MuCS
      MODULE_LIBRARIES Geometry
			Geometry_service
			Simulation
			BackTracker_service
      Utilities
			Filters
			RawData
			RecoBase
			RecoObjects  
			RecoAlg
			uboonecode_uboone_MuCS
			SignalShapingServiceMicroBooNE_service
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

# install_headers()
install_fhicl()
install_source()
