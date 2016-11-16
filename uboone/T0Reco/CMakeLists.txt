
art_make(	
	  LIB_LIBRARIES
                        lardata_Utilities_LArFFT_service
			${ART_FRAMEWORK_CORE}
			art_Utilities canvas_Utilities
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
			${CETLIB}
                        ${Boost_SYSTEM_LIBRARY}
			${ROOT_BASIC_LIB_LIST}
          MODULE_LIBRARIES 
			larcore_Geometry
			larcore_Geometry_Geometry_service
			larsim_MCCheater_BackTracker_service
			larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
			lardataobj_RecoBase
			lardataobj_AnalysisBase	
			lardata_RecoObjects larreco_Deprecated  
			larreco_RecoAlg
			lardata_Utilities
			larevt_Filters
			nusimdata_SimulationBase
			${LARRECO_LIB}
			${LARDATA_LIB}
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			art_Persistency_Common canvas_Persistency_Common
			art_Persistency_Provenance canvas_Persistency_Provenance
			art_Utilities canvas_Utilities
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
	                ${FHICLCPP}
			${CETLIB}
			${CLHEP}
                        ${ROOT_GEOM}
                        ${ROOT_XMLIO}
                        ${ROOT_GDML}
			${ROOT_BASIC_LIB_LIST}
        )

install_headers()
install_fhicl()
install_source()

