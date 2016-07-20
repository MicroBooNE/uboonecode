# basic source code CMakeLists.txt


art_make(BASENAME_ONLY
	LIBRARY_NAME TPCNeutrinoIDFilter
	LIB_LIBRARIES 
		     lardataobj_RecoBase
		     lardataobj_AnalysisBase
		     larcore_Geometry
		     larcore_Geometry_Geometry_service
		     lardata_Utilities
		     art_Persistency_Common canvas_Persistency_Common
		     art_Persistency_Provenance canvas_Persistency_Provenance
		     ${ART_FRAMEWORK_CORE}
		     ${ART_FRAMEWORK_PRINCIPAL}
		     ${ART_FRAMEWORK_SERVICES_REGISTRY}
		     ${ART_FRAMEWORK_SERVICES_OPTIONAL}
		     ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
		     ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILEDIRECTORY_SERVICE}
		     art_Utilities canvas_Utilities
		     ${MF_MESSAGELOGGER}
		     ${MF_UTILITIES}
		     ${FHICLCPP}
		     ${CETLIB}
		     ${ROOT_BASIC_LIB_LIST}
		     ${Boost_SYSTEM_LIBRARY}
	MODULE_LIBRARIES
		     TPCNeutrinoIDFilter
		     lardataobj_RecoBase
		     lardataobj_AnalysisBase
		     larcore_Geometry
		     larcore_Geometry_Geometry_service
		     lardata_Utilities
		     art_Persistency_Common canvas_Persistency_Common
		     art_Persistency_Provenance canvas_Persistency_Provenance
		     ${ART_FRAMEWORK_CORE}
		     ${ART_FRAMEWORK_PRINCIPAL}
		     ${ART_FRAMEWORK_SERVICES_REGISTRY}
		     ${ART_FRAMEWORK_SERVICES_OPTIONAL}
		     ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
		     ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILEDIRECTORY_SERVICE}
		     art_Utilities canvas_Utilities
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
