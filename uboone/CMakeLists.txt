
cet_find_library( BOOST_SERIALIZATION NAMES boost_serialization PATHS ENV BOOST_LIB NO_DEFAULT_PATH )
cet_find_library( BOOST_DATE_TIME     NAMES boost_date_time     PATHS ENV BOOST_LIB NO_DEFAULT_PATH )



add_subdirectory(Geometry)
add_subdirectory(RawData)
add_subdirectory(Utilities)
add_subdirectory(Simulation)
add_subdirectory(LArG4)
add_subdirectory(PhotonPropagation)
add_subdirectory(DetSim)
add_subdirectory(EventGenerator)
add_subdirectory(CalData)
add_subdirectory(TriggerAlgo)
add_subdirectory(AnalysisTree)

# Respect the order specified in lar_build, which is:
#
#     "SoftRelTools",
#     "SRT_LAR",
#     "SimpleTypesAndConstants",
#     "Geometry",
#     "SummaryData",
#     "RawData",
#     "RecoBase",
#     "OpticalDetectorData",
#     "AnalysisBase",
#     "Utilities",
#     "Simulation",
#     "MCCheater",
#     "AnalysisAlg",
#     "RecoObjects",
#     "RecoAlg",
#     "LArG4",
#     "PhotonPropagation",
#     "DetSim",
#     "EventGenerator",
#     "Filters",
#     "CalData",
#     "ClusterFinder",
#     "Genfit",
#     "HitFinder",
#     "TrackFinder",
#     "VertexFinder",
#     "EventFinder",
#     "TriggerAlgo",
#     "OpticalDetector",
#     "Calorimetry", 
#     "ShowerFinder",
#     "EventDisplay",
#     "ParticleIdentification",
#     "AnalysisExample",
#     "LArPandoraAlgorithms",
#     "LArPandoraInterface"