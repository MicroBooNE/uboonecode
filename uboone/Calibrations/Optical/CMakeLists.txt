# Either set these here or set environment variables

include_directories ( $ENV{SWTRIGGER_INCDIR} )
find_library( SWTRIG_LIBBASE NAMES SWTriggerBase PATHS $ENV{SWTRIGGER_LIBDIR} )
find_library( SWTRIG_LIBFEMU NAMES FEMBeamTrigger PATHS $ENV{SWTRIGGER_LIBDIR} )

add_subdirectory(fcl)

art_make( LIBRARY_NAME SPEcalibrationCode
	  EXCLUDE SPEcalibration_module.cc FEMemulator_module.cc SWTrigger_module.cc CalibrationFit.C
	  lardataobj_RawData
	  UBOpReadoutMap_service
	  LIB_LIBRARIES ${ART_FRAMEWORK_CORE} 
	       ${ART_FRAMEWORK_IO_SOURCES}
               ${ART_FRAMEWORK_PRINCIPAL}
               ${ART_FRAMEWORK_SERVICES_REGISTRY}
               ${ART_FRAMEWORK_SERVICES_OPTIONAL}
               ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
               art_Persistency_Common canvas_Persistency_Common
               art_Persistency_Provenance canvas_Persistency_Provenance
               art_Utilities canvas_Utilities
	       ${ROOT_BASIC_LIB_LIST}	     
	       BASENAME_ONLY  
	)
		

simple_plugin( SPEcalibration "module"
	       SPEcalibrationCode
	       OpticalSubEvents
	       UBOpReadoutMap_service
	       lardataobj_RawData
	       ${ART_FRAMEWORK_CORE}
	       ${ART_FRAMEWORK_IO_SOURCES}
               ${ART_FRAMEWORK_PRINCIPAL}
               ${ART_FRAMEWORK_SERVICES_REGISTRY}
               ${ART_FRAMEWORK_SERVICES_OPTIONAL}
               ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
               art_Persistency_Common canvas_Persistency_Common
               art_Persistency_Provenance canvas_Persistency_Provenance
               art_Utilities canvas_Utilities
	       ${ROOT_BASIC_LIB_LIST}	       
               BASENAME_ONLY
             )

simple_plugin( FEMemulator "module"
	       SPEcalibrationCode
	       OpticalSubEvents
	       UBOpReadoutMap_service
	       lardataobj_RawData
	       ${ART_FRAMEWORK_CORE}
	       ${ART_FRAMEWORK_IO_SOURCES}
               ${ART_FRAMEWORK_PRINCIPAL}
               ${ART_FRAMEWORK_SERVICES_REGISTRY}
               ${ART_FRAMEWORK_SERVICES_OPTIONAL}
               ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
               art_Persistency_Common canvas_Persistency_Common
               art_Persistency_Provenance canvas_Persistency_Provenance
               art_Utilities canvas_Utilities
	       ${ROOT_BASIC_LIB_LIST}	       
               BASENAME_ONLY
             )

simple_plugin( SWTrigger "module"
	       UBOpReadoutMap_service
	       uboone_RawData
	       lardata_RawData
	       lardata_RecoBase
	       ${ART_FRAMEWORK_CORE}
	       ${ART_FRAMEWORK_IO_SOURCES}
               ${ART_FRAMEWORK_PRINCIPAL}
               ${ART_FRAMEWORK_SERVICES_REGISTRY}
               ${ART_FRAMEWORK_SERVICES_OPTIONAL}
               ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
               ${ART_PERSISTENCY_COMMON}
               ${ART_PERSISTENCY_PROVENANCE}
               ${ART_UTILITIES}
	       ${SWTRIG_LIBFEMU}
	       ${SWTRIG_LIBBASE}
	       ${ROOT_BASIC_LIB_LIST}	       
               BASENAME_ONLY
             )

install_fhicl()
install_headers()
install_source()

