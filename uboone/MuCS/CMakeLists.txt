
# for classes_def.xml!!
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
include_directories( $ENV{NUTOOLS_INC} )
cet_find_library( SQLITE3 NAMES sqlite3 PATHS ENV NO_DEFAULT_PATH )
cet_find_library( CRY NAMES CRY PATHS ENV CRY_LIB NO_DEFAULT_PATH )
cet_find_library( IFDH NAMES ifdh PATHS ENV IFDHC_LIB NO_DEFAULT_PATH )
include_directories ( $ENV{CRYHOME}/src )
include_directories( $ENV{IFDHC_FQ_DIR}/inc )

art_make( BASENAME_ONLY
	    LIBRARY_NAME  uboonecode_uboone_MuCS
      MODULE_LIBRARIES larcore_Geometry
			larcore_Geometry_Geometry_service
			larcoreobj_SummaryData
			larsim_Simulation larsimobj_Simulation
			larsim_MCCheater_BackTracker_service
      lardata_Utilities
			larevt_Filters
			lardataobj_RawData
			lardataobj_RecoBase
			lardata_RecoObjects  
			larreco_RecoAlg
			uboonecode_uboone_MuCS
			SignalShapingServiceMicroBooNE_service
			lardataobj_AnalysisBase
      nusimdata_SimulationBase
      ${LARRECO_LIB}
      ${LARDATA_LIB}
      ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			art_Persistency_Common canvas_Persistency_Common
			art_Persistency_Provenance canvas_Persistency_Provenance
			art_Utilities canvas_Utilities
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
      ${FHICLCPP}
			${CETLIB}
			${CLHEP}
			${ROOT_PHYSICS}
      ${ROOT_GEOM}
      ${ROOT_XMLIO}
      ${ROOT_GDML}
			${ROOT_BASIC_LIB_LIST}
			${ROOT_CORE}
			${ROOT_EG}
			${SQLITE3}
			${CRY}
			${IFDH}
)

install_fhicl()
install_source()
install_headers()
