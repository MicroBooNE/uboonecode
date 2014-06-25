simple_plugin(KalmanFilterTest "module"
              RecoObjects
	      ${MF_MESSAGELOGGER}
	      ${MF_UTILITIES}
	      ${FHICLCPP}
	      ${CETLIB}
	      ${PQ}
              BASENAME_ONLY
             )

simple_plugin(PropTest "module"
              RecoObjects
	      ${MF_MESSAGELOGGER}
	      ${MF_UTILITIES}
	      ${FHICLCPP}
	      ${CETLIB}
	      ${PQ}
              BASENAME_ONLY
             )

cet_test( PropTest HANDBUILT 
  	  TEST_EXEC lar
	  TEST_ARGS -c PropTest.fcl
	  DATAFILES PropTest.fcl )
cet_test( KalmanFilterTest HANDBUILT 
  	  TEST_EXEC lar
	  TEST_ARGS -c KalmanFilterTest.fcl
	  DATAFILES KalmanFilterTest.fcl )

install_headers()
install_fhicl()
install_source()
