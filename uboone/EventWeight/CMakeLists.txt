cet_find_library( IFDH NAMES ifdh PATHS ENV IFDHC_LIB NO_DEFAULT_PATH )
include_directories( $ENV{IFDHC_FQ_DIR}/inc )
file(GLOB WGHCALCFILES *WeightCalc* Factory* )
# include_directories ( $ENV{GENIE_INC}/GENIE )

art_make_library(LIBRARY_NAME WeightCalculators 
                 SOURCE	      ${WGHCALCFILES}
	         LIBRARIES
		       lardata_Utilities
#		       nutools_NuReweight_art
#                       nutools_NuReweight
		       larsim_Simulation
                       nutools_RandomUtils_NuRandomService_service
		       nusimdata_SimulationBase
              	      ${ART_FRAMEWORK_CORE}
		      ${ART_FRAMEWORK_PRINCIPAL}
		      ${ART_FRAMEWORK_SERVICES_REGISTRY}
		      ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
		      art_Persistency_Common canvas_Persistency_Common
                      art_Persistency_Provenance canvas_Persistency_Provenance
                      art_Utilities canvas_Utilities
                      ${CETLIB}
		      ${FHICLCPP}
                      ${MF_MESSAGELOGGER}
		      ${CLHEP}
                      ${ROOT_BASIC_LIB_LIST}
                      ${Boost_SYSTEM_LIBRARY}
                      ${IFDH}
	        ) 

art_make(LIBRARY_NAME   EventWeight
         EXCLUDE        ${WGHCALCFILES}
	 LIB_LIBRARIES
                      BNBPrimaryHadronWeightCalculators
		      WeightCalculators
#              	      ${ART_FRAMEWORK_CORE}
#		      ${ART_FRAMEWORK_IO_SOURCES}
#                      ${ART_FRAMEWORK_PRINCIPAL}
#                      ${ART_FRAMEWORK_SERVICES_REGISTRY}
#                      ${ART_FRAMEWORK_SERVICES_OPTIONAL}
#                      ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
#                      ${ART_PERSISTENCY_COMMON}
#                      ${ART_PERSISTENCY_PROVENANCE}
#                      ${ART_UTILITIES}
#		      nutools_RandomUtils_NuRandomService_service
                      ${MF_MESSAGELOGGER}
                      ${MF_UTILITIES}
                       ${IFDH}
#                      ${FHICLCPP}
#                      ${CETLIB}
#                      ${ROOT_GEOM}
#                      ${ROOT_XMLIO}
#                      ${ROOT_GDML}
#                      ${ROOT_BASIC_LIB_LIST}
                       ${Boost_SYSTEM_LIBRARY}
        MODULE_LIBRARIES EventWeight
			 WeightCalculators
                         BNBPrimaryHadronWeightCalculators

	)


install_headers()
install_fhicl()
install_source()
add_subdirectory(BNBPrimaryHadron)
