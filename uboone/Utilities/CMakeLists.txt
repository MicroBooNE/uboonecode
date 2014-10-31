cet_enable_asserts()

set( uboone_util_lib_list LArFFT_service
                        LArProperties_service
			Geometry
			HitFinder
			RecoBase
			Geometry_service
			Utilities
  			${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_PERSISTENCY_COMMON}
			${ART_PERSISTENCY_PROVENANCE}
			${ART_FRAMEWORK_SERVICES_BASIC}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			${ART_UTILITIES}
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
			${CETLIB}
			${ROOT_BASIC_LIB_LIST}
    )


simple_plugin( SignalShapingServiceMicroBooNE  "service"
               ${uboone_util_lib_list}
	       BASENAME_ONLY
        )
simple_plugin( TinyHitFinderAna "module"
               ${uboone_util_lib_list}
	       BASENAME_ONLY
        )
             
 
install_headers()
# do not install duplicate standard_reco_uboone.fcl file
# hack until install_fhicl( EXCLUDES standard_reco_uboone.fcl ) does the right thing
install_fhicl( LIST databaseutil_microboone.fcl
		    detectorproperties_microboone.fcl
		    evdservices_microboone.fcl
		    tinyhitfinderana_microboone.fcl
		    hitfinder_ana_uboone.fcl
		    larfft_microboone.fcl
		    larproperties_microboone.fcl
		    optical_reco_uboone.fcl
		    pmtconfig_microboone.fcl
		    sam_microboone.fcl
		    services_microboone_1window.fcl
		    services_microboone_basic.fcl
		    services_microboone.fcl
		    services_microboone_simulation.fcl
		    signalservices_microboone.fcl
		    timeservice_microboone.fcl )
install_source()

add_subdirectory(test)

