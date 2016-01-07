# project name
project(beamdaq)

# ----------------------------------------------------------------------
# Build this project's executables:

add_executable(bdaq_get bdaq_get.cpp)
add_executable(bdaq_readraw bdaq_readraw.cpp)

include_directories ( $ENV{UBOONEDAQ_DATATYPES_INC} )
set(CURLLIB "-lcurl")

find_library(UBDATATYPESLIB NAMES ubdata_types ubdata_share PATHS $ENV{UBOONEDAQ_DATATYPES_LIB})

cet_make(LIBRARY_NAME BeamDAQ LIBRARIES ${UBDATATYPESLIB}  
                                        ${MF_MESSAGELOGGER} ${MF_UTILITIES} ${CURLLIB} 
					${Boost_PROGRAM_OPTIONS_LIBRARY} ${Boost_FILESYSTEM_LIBRARY} ${Boost_DATE_TIME_LIBRARY} ${Boost_SERIALIZATION_LIBRARY} 
					${FHICLCPP} ${CETLIB} 
					EXCLUDE bdaq_get.cpp bdaq_get.h bdaq_readraw.cpp)

target_link_libraries( bdaq_get ${UBDATATYPESLIB}  
                                ${MF_MESSAGELOGGER} ${MF_UTILITIES} ${CURLLIB} 
				${Boost_PROGRAM_OPTIONS_LIBRARY} ${Boost_FILESYSTEM_LIBRARY} ${Boost_DATE_TIME_LIBRARY} ${Boost_SERIALIZATION_LIBRARY} 
				${FHICLCPP} ${CETLIB} BeamDAQ)
target_link_libraries( bdaq_readraw ${Boost_PROGRAM_OPTIONS_LIBRARY} ${Boost_FILESYSTEM_LIBRARY} ${Boost_DATE_TIME_LIBRARY} ${Boost_SERIALIZATION_LIBRARY} 
                                    ${UBDATATYPESLIB} BeamDAQ )


# ----------------------------------------------------------------------
# Install this project's headers, sources, library, and executables:

install_fhicl()

install_scripts()

install( TARGETS bdaq_get
         DESTINATION ${flavorqual_dir}/bin )

install( TARGETS bdaq_readraw
         DESTINATION ${flavorqual_dir}/bin )

