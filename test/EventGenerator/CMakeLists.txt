
# simple examples using prodsingle

# use the general prodsingle_uboone.fcl file
# because OPTIONAL_GROUPS is defined, this test will not be run by default
# use mrb t --test-all to run all the tests
cet_test(prodsingle_uboone HANDBUILT
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config prodsingle_uboone.fcl
  OPTIONAL_GROUPS RELEASE
)

# use a local fhicl file
cet_test(prodsingle_uboone_max2 HANDBUILT
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config prodsingle_uboone_max2.fcl
  DATAFILES prodsingle_uboone_max2.fcl
)
