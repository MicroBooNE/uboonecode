# define ROOTLIB
#set( ROOTLIB -L$ENV{ROOTSYS}/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -pthread -lm -ldl)

set( LARLITE_LIB -L$ENV{LARLITE_LIBDIR} -lLArLite_Base -lLArLite_DataFormat -lLArLite_LArUtil)

set( UBOONECODE_LIB -L$ENV{UBOONECODE_LIB} )

include_directories ( $ENV{LARLITE_COREDIR} )

link_libraries( ${LIB_NAME} -L$ENV{BOOST_LIB} -lboost_system ${ROOTLIB} ${LARLITE_LIB}) # All execs need this.

add_subdirectory(job)

art_make( EXCLUDE LLMetaMaker_service.cc
	  BASENAME_ONLY
	  LIBRARY_NAME     uboonecode_uboone_LiteMaker
	  LIB_LIBRARIES
	  	   larcore_Geometry
		   uboone_Geometry
		   uboonecode_uboone_MuCS
		   larcore_Geometry_Geometry_service
	   	   larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
		   lardataobj_RecoBase
		   nusimdata_SimulationBase
		   lardataobj_MCBase
		   larcoreobj_SummaryData
		   lardataobj_RawData
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

          MODULE_LIBRARIES 	  		   
		   uboonecode_uboone_LiteMaker	   
	      	   TFileMetadataMicroBooNE_service
	      	   FileCatalogMetadataMicroBooNE_service
		   LLMetaMaker_service
		   uboonecode_uboone_MuCS
		   larcore_Geometry
		   uboone_Geometry
		   larcore_Geometry_Geometry_service
		   larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
		   nusimdata_SimulationBase
		   lardata_Utilities
		   lardataobj_MCBase
		   larcoreobj_SummaryData
		   lardataobj_RawData
		   lardataobj_RecoBase
		   lardata_RecoObjects larreco_Deprecated
                   larevt_Filters
                   larevt_CalibrationDBI_IOVData
                   larevt_CalibrationDBI_Providers
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

simple_plugin(LLMetaMaker "service"
	      TFileMetadataMicroBooNE_service
	      FileCatalogMetadataMicroBooNE_service
	      ${ART_FRAMEWORK_SERVICES_SYSTEM_FILECATALOGMETADATA_SERVICE}	      
	      art_Utilities canvas_Utilities
	      ${MF_MESSAGELOGGER}
	      ${MF_UTILITIES}
	      ${FHICLCPP}
	      ${CETLIB}
	      ${ROOT_BASIC_LIB_LIST}
              BASENAME_ONLY
)

# install_headers()
install_fhicl()
install_source()
