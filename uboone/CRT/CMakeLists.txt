include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

art_make( BASENAME_ONLY
	    LIBRARY_NAME  uboonecode_crt
      MODULE_LIBRARIES larcore_Geometry
			larcore_Geometry_Geometry_service
			larcoreobj_SummaryData
			larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
      lardata_Utilities
			larevt_Filters
			lardataobj_RawData
			lardataobj_RecoBase
			lardata_RecoObjects larreco_Deprecated
			larreco_RecoAlg
			uboonecode_uboone_MuCS
			lardataobj_AnalysisBase
      nusimdata_SimulationBase
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
			${IFDH}
)

install_fhicl()
install_source()
install_headers()
