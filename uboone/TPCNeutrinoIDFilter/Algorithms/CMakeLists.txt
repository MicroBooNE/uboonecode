# basic source code CMakeLists.txt


art_make(BASENAME_ONLY
	LIBRARY_NAME TPCNeutrinoIDAlgs
	LIB_LIBRARIES 
		     lardata_RecoBase
		     lardata_AnalysisBase
		     larcore_Geometry
		     larcore_Geometry_Geometry_service
		     lardata_Utilities
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
		     lardata_RecoBase
		     lardata_AnalysisBase
		     larcore_Geometry
		     larcore_Geometry_Geometry_service
		     lardata_Utilities
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
