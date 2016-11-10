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
			   larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
			   nusimdata_SimulationBase
			   larsim_MCCheater_BackTracker_service
			   lardata_Utilities
			   larevt_Filters
			   lardataobj_RawData
			   lardataobj_RecoBase
			   lardata_RecoObjects larreco_Deprecated
			   ${UBOONECODE_LIB}
			   nusimdata_SimulationBase
			   ${LARRECO_LIB}
			   ${LARDATA_LIB}
			   ${ART_FRAMEWORK_CORE}
			   ${ART_FRAMEWORK_PRINCIPAL}
                           ${ART_FRAMEWORK_BASIC}
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
			   ${ROOT_GEOM}
			   ${ROOT_XMLIO}
			   ${ROOT_GDML}
			   ${ROOT_BASIC_LIB_LIST}
			)

# install_headers()
install_fhicl()
install_source()
