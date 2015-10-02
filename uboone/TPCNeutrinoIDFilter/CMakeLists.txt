# basic source code CMakeLists.txt


art_make(BASENAME_ONLY
	LIBRARY_NAME TPCNeutrinoIDFilter
	LIB_LIBRARIES 
		     RecoBase
		     AnalysisBase
		     Geometry
		     Geometry_service
		     LArProperties_service
		     Utilities
		     ${ART_PERSISTENCY_COMMON}
		     ${ART_PERSISTENCY_PROVENANCE}
		     ${ART_FRAMEWORK_CORE}
		     ${ART_FRAMEWORK_PRINCIPAL}
		     ${ART_FRAMEWORK_SERVICES_REGISTRY}
		     ${ART_FRAMEWORK_SERVICES_OPTIONAL}
		     ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
		     ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILEDIRECTORY_SERVICE}
		     ${ART_UTILITIES}
		     ${MF_MESSAGELOGGER}
		     ${MF_UTILITIES}
		     ${FHICLCPP}
		     ${CETLIB}
		     ${ROOT_BASIC_LIB_LIST}
		     ${Boost_SYSTEM_LIBRARY}
	MODULE_LIBRARIES
		     TPCNeutrinoIDFilter
		     RecoBase
		     AnalysisBase
		     Geometry
		     Geometry_service
		     LArProperties_service
		     DetectorProperties_service
		     Utilities
		     ${ART_PERSISTENCY_COMMON}
		     ${ART_PERSISTENCY_PROVENANCE}
		     ${ART_FRAMEWORK_CORE}
		     ${ART_FRAMEWORK_PRINCIPAL}
		     ${ART_FRAMEWORK_SERVICES_REGISTRY}
		     ${ART_FRAMEWORK_SERVICES_OPTIONAL}
		     ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
		     ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILEDIRECTORY_SERVICE}
		     ${ART_UTILITIES}
		     ${MF_MESSAGELOGGER}
		     ${MF_UTILITIES}
		     ${FHICLCPP}
		     ${CETLIB}
		     ${ROOT_BASIC_LIB_LIST}
		     ${Boost_SYSTEM_LIBRARY}
         )

install_headers()
install_source()
install_fhicl()
