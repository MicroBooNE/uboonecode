set( ROOTLIB -L$ENV{ROOTSYS}/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl)
link_libraries( ${LIB_NAME} -L$ENV{BOOST_LIB} ${ROOTLIB} )
art_make( BASENAME_ONLY
	  LIBRARY_NAME     uboonecode_uboone_BasicTool_GeoAlgo
	  LIB_LIBRARIES
	  ${ROOT_BASIC_LIB_LIST}
	)
install_headers()
install_source()
