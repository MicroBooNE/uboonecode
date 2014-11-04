# various integration tests

include(CetTest)

cet_enable_asserts()



add_subdirectory(EventGenerator)
add_subdirectory(Geometry)
add_subdirectory(RecoObjects)

add_subdirectory(ci)
