# various integration tests

include(CetTest)

cet_enable_asserts()



add_subdirectory(EventGenerator)
add_subdirectory(Geometry)
add_subdirectory(RecoObjects)
add_subdirectory(test_fcl)
add_subdirectory(Production)

add_subdirectory(ci)
