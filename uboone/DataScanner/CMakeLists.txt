# define ROOTLIB
set( ROOTLIB -L$ENV{ROOTSYS}/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl)

set( LARLIGHT_DIR $ENV{MAKE_TOP_DIR}/core )

set( LARLIGHT_LIB -L$ENV{LARLIGHT_LIBDIR} -lBase -lDataFormat -lLArUtil)

set( UBOONECODE_LIB -L$ENV{UBOONECODE_LIB} )

include_directories ( ${LARLIGHT_DIR} )
include_directories ( ${LARLIGHT_DIR}/Base )
include_directories ( ${LARLIGHT_DIR}/DataFormat )

link_libraries( ${LIB_NAME} ${ROOTLIB} ${LARLIGHT_LIB}) # All execs need this.

add_subdirectory(job)

art_make( BASENAME_ONLY
          MODULE_LIBRARIES uboonecode_uboone_MCShowerReco
			   larcore_Geometry
			   larcore_Geometry_Geometry_service
			   larsim_Simulation
			   SimulationBase
			   larsim_MCCheater_BackTracker_service
			   lardata_Utilities
			   larevt_Filters
			   lardata_RawData
			   lardata_RecoBase
			   lardata_RecoObjects
			   lardata_Utilities_LArProperties_service
			   lardata_Utilities_DetectorProperties_service
			   ${UBOONECODE_LIB}
			   ${SIMULATIONBASE}
			   ${LARRECO_LIB}
			   ${LARDATA_LIB}
			   ${ART_FRAMEWORK_CORE}
			   ${ART_FRAMEWORK_PRINCIPAL}
                           ${ART_FRAMEWORK_BASIC}
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
			)

# install_headers()
install_fhicl()
install_source()
