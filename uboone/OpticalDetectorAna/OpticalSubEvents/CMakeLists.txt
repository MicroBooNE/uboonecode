# find_ups_product is meant to be called from the top CMakeLists.txt file

add_subdirectory( cfdiscriminator_algo )
add_subdirectory( subevent_algo )

simple_plugin( SubEventBuilder "module"
	       uboone_cfdiscriminator
	       uboone_subeventalgo
               ${ART_FRAMEWORK_CORE}
	       ${ART_FRAMEWORK_IO_SOURCES}
               ${ART_FRAMEWORK_PRINCIPAL}
               ${ART_FRAMEWORK_SERVICES_REGISTRY}
               ${ART_FRAMEWORK_SERVICES_OPTIONAL}
               ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
               ${ART_PERSISTENCY_COMMON}
               ${ART_PERSISTENCY_PROVENANCE}
               ${ART_UTILITIES}
	       ${ROOT_BASIC_LIB_LIST}	       
               BASENAME_ONLY
             )


install_headers()
install_source()

