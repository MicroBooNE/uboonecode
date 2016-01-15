
cet_find_library( BOOST_SERIALIZATION NAMES boost_serialization PATHS ENV BOOST_LIB NO_DEFAULT_PATH )
cet_find_library( BOOST_DATE_TIME     NAMES boost_date_time     PATHS ENV BOOST_LIB NO_DEFAULT_PATH )


set( UBOONEDAQ_DATATYPES_DIR $ENV{UBOONEDAQ_DATATYPES_DIR} )

add_subdirectory(Geometry)
add_subdirectory(RawData)
if( UBOONEDAQ_DATATYPES_DIR )
	add_subdirectory(BeamData)
	add_subdirectory(BeamDAQ)
endif()
add_subdirectory(BeamAna)
add_subdirectory(Utilities)
add_subdirectory(Database)
add_subdirectory(Simulation)
add_subdirectory(LArG4)
add_subdirectory(PhotonPropagation)
add_subdirectory(DetSim)
add_subdirectory(EventGenerator)
add_subdirectory(CalData)
add_subdirectory(TriggerSim)
add_subdirectory(AnalysisTree)
add_subdirectory(Calibrations)
add_subdirectory(OpticalDetectorSim)
add_subdirectory(OpticalDetectorAna)
add_subdirectory(SpaceCharge)
add_subdirectory(TPCNeutrinoIDFilter)
#add_subdirectory(DataScanner)
#add_subdirectory(LiteMaker)

# Respect the order specified in lar_build, which is:
#
#     "SoftRelTools",
#     "SRT_LAR",
#     "SimpleTypesAndConstants",
#     "larcore_Geometry",
#     "larcore_SummaryData",
#     "lardata_RawData",
#     "lardata_RecoBase",
#     "lardata_OpticalDetectorData",
#     "lardata_AnalysisBase",
#     "lardata_Utilities",
#     "larsim_Simulation",
#     "larsim_MCCheater",
#     "lardata_AnalysisAlg",
#     "lardata_RecoObjects",
#     "larreco_RecoAlg",
#     "larsim_LArG4",
#     "larsim_PhotonPropagation",
#     "larsim_DetSim",
#     "larsim_EventGenerator",
#     "larevt_Filters",
#     "CalData",
#     "larreco_ClusterFinder",
#     "larreco_Genfit",
#     "larreco_HitFinder",
#     "larreco_TrackFinder",
#     "larreco_VertexFinder",
#     "larreco_EventFinder",
#     "larsim_TriggerAlgo",
#     "larana_OpticalDetector",
#     "larana_Calorimetry", 
#     "larreco_ShowerFinder",
#     "lareventdisplay_EventDisplay",
#     "larana_ParticleIdentification",
#     "larexamples_AnalysisExample",
#     "LArPandoraAlgorithms",
#     "larpandora_LArPandoraInterface"
