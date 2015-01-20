# define ROOTLIB
set( ROOTLIB -L$ENV{ROOTSYS}/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl)

set( LARLITE_LIB -L$ENV{LARLITE_LIBDIR} -lLArLite_Base -lLArLite_DataFormat -lLArLite_LArUtil)

set( UBOONECODE_LIB -L$ENV{UBOONECODE_LIB} )

include_directories ( $ENV{LARLITE_COREDIR} )

link_libraries( ${LIB_NAME} ${ROOTLIB} ${LARLITE_LIB}) # All execs need this.

add_subdirectory(job)

art_make( BASENAME_ONLY
	  LIBRARY_NAME     uboonecode_uboone_LiteMaker
	  LIB_LIBRARIES
	  	   Geometry
		   Geometry_service
	   	   Simulation
		   SimulationBase
		   Utilities
		   RawData
		   RecoBase
		   MCBase
		   SummaryData
		   RecoObjects
		   LArProperties_service
		   DetectorProperties_service
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
          MODULE_LIBRARIES 	  		   
		   uboonecode_uboone_LiteMaker
		   Geometry
		   Geometry_service
		   Simulation
		   SimulationBase
		   Utilities
		   MCBase
		   SummaryData
		   RawData
		   RecoBase
		   RecoObjects
		   LArProperties_service
		   DetectorProperties_service
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
